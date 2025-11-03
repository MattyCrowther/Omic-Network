from __future__ import annotations
from dataclasses import asdict
from typing import Dict, Any, List, Tuple, Iterable, Optional
import pandas as pd
from ...errors import ParseError
from ...parse_obj import OmicData, FeatureRec, CrossRef
from ...identifiers import IDS
from ...utils import to_primitive
from Bio import SeqIO


def _q(q: dict, key: str) -> str | None:
    v = q.get(key)
    if not v:
        return None
    return v[0] if isinstance(v, list) else str(v)


def _is_missing(x) -> bool:
    return x is None or x == "" or pd.isna(x)


def _open_text(path: str):
    if str(path).endswith(".gz"):
        import gzip

        return gzip.open(path, "rt")
    return open(path, "rt")


def _iter_genbank(path: str):
    with _open_text(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            yield rec


def _parse_record(
    rec, keep_types: set[str], include_pseudogenes: bool
) -> Tuple[List[Dict[str, Any]], List[Tuple[str, str, str, str]]]:
    """Return merged feature rows and cross-ref rows for one record."""
    contig = (getattr(rec, "id", "") or getattr(rec, "name", "")).split(".", 1)[0]
    merged: Dict[str, Dict[str, Any]] = {}
    cr_rows: List[Tuple[str, str, str, str]] = []
    cr_seen: set[Tuple[str, str, str]] = set()

    for feat in rec.features:
        ftype = str(getattr(feat, "type", "")).lower()
        if ftype not in keep_types:
            continue

        q = feat.qualifiers or {}
        if not include_pseudogenes and ("pseudo" in q or _q(q, "pseudogene")):
            continue

        locus = _q(q, "locus_tag")
        if not locus:
            continue

        start = int(feat.location.start) + 1
        end = int(feat.location.end)
        strand = int(getattr(feat.location, "strand", 0) or 0)

        row = {
            "locus_tag": locus,
            "product": _q(q, "product"),
            "type": ftype,
            "contig": contig,
            "start": start,
            "end": end,
            "strand": strand,
        }
        if locus not in merged:
            merged[locus] = row
        else:
            for k, v in row.items():
                cur = merged[locus].get(k, None)
                if _is_missing(cur) and not _is_missing(v):
                    merged[locus][k] = v

        for ns_key in ("gene", "protein_id"):
            val = _q(q, ns_key)
            if val:
                trip = (locus, ns_key.lower(), val)
                if trip not in cr_seen:
                    cr_seen.add(trip)
                    cr_rows.append(("feature", locus, ns_key.lower(), val))

        for s in q.get("db_xref", []) or []:
            if ":" in s:
                ns, val = (t.strip() for t in s.split(":", 1))
                if ns and val:
                    trip = (locus, ns, val)
                    if trip not in cr_seen:
                        cr_seen.add(trip)
                        cr_rows.append(("feature", locus, ns, val))

        prod = _q(q, "product")
        if prod:
            cr_rows.append(("feature", locus, "product_label", prod.strip()))

    return list(merged.values()), cr_rows


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

    keep = [
        c
        for c in ("product", "type", "contig", "start", "end", "strand")
        if c in df.columns
    ]
    recs: List[FeatureRec] = []
    for fid, *vals in df[keep].itertuples(name=None):
        attrs = {k: to_primitive(v) for k, v in zip(keep, vals)}

        valid_entities = set(asdict(IDS.type).values())
        etype = str(attrs.pop("type", "")).lower()

        if etype in valid_entities:
            e_type = etype
        else:
            e_type = IDS.type.dna
        recs.append(
            FeatureRec(
                id=str(fid), entity=e_type, namespace=IDS.ns.RefSeq_Locus, attrs=attrs
            )
        )
    return recs


def _empty_matrix(index: List[str] | pd.Index) -> pd.DataFrame:
    idx = pd.Index(index, dtype="string")
    m = pd.DataFrame(index=idx, dtype="float64")
    m.columns.name = None
    return m


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
        out.append(CrossRef(str(src), str(ns), str(tgt)))
    return out


def parse_genbank(path: str, *, include_pseudogenes: bool = False) -> OmicData:
    """Parse GenBank into OmicData with per-feature metadata and cross-refs."""
    keep_types = {"cds", "gene", "trna", "rrna", "ncrna", "misc_rna"}
    all_rows: List[Dict[str, Any]] = []
    all_cr: List[Tuple[str, str, str, str]] = []

    try:
        for rec in _iter_genbank(path):
            rows, cr_rows = _parse_record(rec, keep_types, include_pseudogenes)
            all_rows.extend(rows)
            all_cr.extend(cr_rows)
    except Exception as e:
        raise ParseError(f"GenBank parsing failed: {e}", path=str(path), cause=e)

    features = _build_feature_meta(all_rows)
    feature_ids = [r.id for r in features]
    matrix = _empty_matrix(feature_ids)
    cross_ref = _build_cross_ref(all_cr, feature_ids)

    return OmicData(
        name="genbank_reference",
        omics_type=IDS.omics_type.genomics,
        matrix=matrix,
        feature_meta=features,
        column_meta=[],
        cross_ref=cross_ref,
    )
