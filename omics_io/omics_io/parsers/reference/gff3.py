from __future__ import annotations
from dataclasses import asdict
from typing import Dict, Any, List, Tuple, Iterable, Optional
import pandas as pd
from ...errors import ParseError
from ...parse_obj import OmicData, FeatureRec, CrossRef
from ...identifiers import IDS
from ...utils import to_primitive

KEEP_TYPES = {"cds", "gene", "trna", "rrna", "ncrna", "misc_rna"}

def _open_text(path: str):
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path, "rt")

def _parse_attrs(field: str) -> Dict[str, str]:
    if not field or field == ".":
        return {}
    out: Dict[str, str] = {}
    for kv in field.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out

def _iter_dbxref(val: str):
    if not val:
        return
    for item in val.split(","):
        if ":" not in item:
            continue
        ns, target = item.split(":", 1)
        ns, target = ns.strip(), target.strip()
        if ns and target:
            yield ns, target

def _is_missing(x) -> bool:
    return x is None or x == "" or (pd.isna(x) if not isinstance(x, str) else False)

def _empty_matrix(index: List[str] | pd.Index) -> pd.DataFrame:
    idx = pd.Index(index, dtype="string")
    m = pd.DataFrame(index=idx, dtype="float64")
    m.columns.name = None
    return m

def _build_feature_meta(rows: List[Dict[str, Any]]) -> List[FeatureRec]:
    if not rows:
        return []
    df = (
        pd.DataFrame(rows)
        .drop_duplicates(subset=["locus_tag"])
        .set_index("locus_tag")
        .sort_index()
    )
    df.index = df.index.astype("string")
    for c in ("product", "type", "contig"):
        if c in df:
            df[c] = df[c].astype("string")
    for c in ("start", "end"):
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    if "strand" in df:
        df["strand"] = pd.to_numeric(df["strand"], errors="coerce").astype("Int8")

    keep = [c for c in ("product", "type", "contig", "start", "end", "strand") if c in df.columns]
    recs: List[FeatureRec] = []

    valid_entities = set(asdict(IDS.type).values())

    for fid, *vals in df[keep].itertuples(name=None):
        attrs = {k: to_primitive(v) for k, v in zip(keep, vals)}
        etype = str(attrs.pop("type", "")).lower()
        entity = etype if etype in valid_entities else IDS.type.dna
        recs.append(
            FeatureRec(
                id=str(fid),
                entity=entity,
                namespace=IDS.ns.RefSeq_Locus,
                attrs=attrs,
            )
        )
    return recs

def _build_cross_ref(
    cr_rows: List[Tuple[str, str, str, str]], valid_features: Iterable[str]
) -> Optional[List[CrossRef]]:
    if not cr_rows:
        return None
    valid = set(map(str, valid_features))
    out: List[CrossRef] = []
    for axis, src, ns, tgt in cr_rows:
        if axis not in ("feature", "column"):
            continue
        if axis == "feature" and str(src) not in valid:
            continue
        out.append(CrossRef(src=str(src), namespace=str(ns), target=str(tgt)))
    return out

def _safe_alias(k: str, v: str, qa: Dict[str, str], prefix: str = "cds:") -> str:
    """
    Normalize alias-like fields so they do not duplicate each other.
    Currently prevents Name == protein_id or gene == protein_id.
    Returns possibly prefixed value.
    """
    pid = qa.get("protein_id")
    if not v or not pid:
        return v
    # Only normalize textual identifiers, not numeric or complex.
    if v == pid and k != "protein_id":
        lt = qa.get("locus_tag")
        return f"{prefix}{pid}" if not lt else f"{prefix}{pid}|{lt}"
    return v


def parse_gff3(path: str) -> OmicData:
    """
    GFF3 → OmicData with:
      - matrix: empty float64, one row per locus_tag
      - feature_meta: List[FeatureRec] with attrs: product, contig, start, end, strand
      - cross_ref: List[CrossRef] from ID, Parent, locus_tag, gene, protein_id, Name, Dbxref[*]
                   plus product_label only if product exists
    """
    merged: Dict[str, Dict[str, Any]] = {}
    ivals: Dict[str, Dict[str, Any]] = {}
    cr_rows: List[Tuple[str, str, str, str]] = []

    def _update_interval(locus: str, contig: str, start_i: int, end_i: int, strand_i: int, ftype_l: str) -> None:
        cur = ivals.get(locus)
        if cur is None:
            ivals[locus] = {"contig": contig, "start": start_i, "end": end_i, "strand": strand_i, "type": ftype_l}
            return
        if cur["contig"] == contig:
            cur["start"] = min(cur["start"], start_i)
            cur["end"] = max(cur["end"], end_i)
            cur["strand"] = strand_i if cur["strand"] == strand_i else 0

    with _open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            ftype_l = str(ftype).lower()
            if ftype_l not in KEEP_TYPES:
                continue

            qa = _parse_attrs(attrs)
            locus = qa.get("locus_tag", "")
            if not locus:
                continue

            contig = seqid.split(".", 1)[0]
            start_i = int(start)
            end_i = int(end)
            strand_i = 1 if strand == "+" else (-1 if strand == "-" else 0)

            row = {
                "locus_tag": locus,
                "product": qa.get("product", ""),
                "type": ftype_l,
                "contig": contig,
                "start": start_i,
                "end": end_i,
                "strand": strand_i,
            }
            if locus not in merged:
                merged[locus] = row
            else:
                if _is_missing(merged[locus].get("product")) and not _is_missing(row.get("product")):
                    merged[locus]["product"] = row.get("product", "")
                if _is_missing(merged[locus].get("type")) and not _is_missing(row.get("type")):
                    merged[locus]["type"] = row.get("type")

            _update_interval(locus, contig, start_i, end_i, strand_i, ftype_l)

            for ns_key in ("ID", "Parent", "locus_tag", "gene", "protein_id", "Name"):
                val = qa.get(ns_key, "")
                val = _safe_alias(ns_key, val, qa, prefix="cds:")
                if val:
                    cr_rows.append(("feature", locus, ns_key, val))

            for ns, target in _iter_dbxref(qa.get("Dbxref", "")):
                cr_rows.append(("feature", locus, ns, target))

            prod = qa.get("product", "")
            if prod:
                cr_rows.append(("feature", locus, "product_label", prod.strip()))

    rows: List[Dict[str, Any]] = []
    if merged:
        for locus, base in merged.items():
            iv = ivals.get(locus, {})
            rows.append(
                {
                    "locus_tag": locus,
                    "product": base.get("product", ""),
                    "type": base.get("type", ""),
                    "contig": iv.get("contig", ""),
                    "start": iv.get("start", None),
                    "end": iv.get("end", None),
                    "strand": iv.get("strand", None),
                }
            )

    features = _build_feature_meta(rows)
    feature_ids = [r.id for r in features]
    matrix = _empty_matrix(feature_ids)
    cross_ref = _build_cross_ref(cr_rows, feature_ids)

    return OmicData(
        name="reference_gff3",
        omics_type=IDS.omics_type.genomics,
        matrix=matrix,
        feature_meta=features,
        column_meta=[],
        cross_ref=cross_ref,
    )
