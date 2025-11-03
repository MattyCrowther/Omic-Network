from __future__ import annotations
import re, gzip, io
from typing import Dict, Any, List, Optional, Tuple
import pandas as pd
from pyteomics import mztab
from ...parse_obj import OmicData, FeatureRec, ColumnRec, CrossRef
from ...errors import ParseError
from ...identifiers import IDS

_RUN_COL_PATS = (
    r"(best_search_engine_score\[\d+\])_ms_run\[(\d+)\]",
    r"(num_psms)_ms_run\[(\d+)\]",
    r"(num_peptides_(?:distinct|unique))_ms_run\[(\d+)\]",
)
_ABUND_PATS = (
    r"abundance_assay\[(\d+)\]",
    r"sub\[\d+\]_abundance_assay\[(\d+)\]",
)
_GN_RE = re.compile(r"\bGN=([A-Za-z0-9_.-]+)\b")
_ASSAY_NAME_RE = re.compile(r"assay\[(\d+)\]-name$")


def _open_mztab(path: str) -> mztab.MzTab:
    """Open .mzTab or .mzTab.gz safely and return a parsed MzTab object."""
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                text = fh.read()
        else:
            with open(path, "rt", encoding="utf-8") as fh:
                text = fh.read()
        return mztab.MzTab(io.StringIO(text))
    except Exception as e:
        raise ParseError(f"mzTab open failed: {e}", path=str(path), cause=e)
    

def _is_decoy(acc: str) -> bool:
    return bool(re.match(r"^(DECOY_|XXX_|REV__|CON__)", str(acc)))


def _rename_run_cols(cols: List[str]) -> Dict[str, str]:
    out = {}
    for c in cols:
        for pat in _RUN_COL_PATS:
            m = re.fullmatch(pat, c)
            if m:
                out[c] = f"{m.group(1)}_run{m.group(2)}"
                break
    return out


def _list_assay_cols(cols: List[str]) -> List[Tuple[str, str]]:
    out = []
    for c in cols:
        for pat in _ABUND_PATS:
            m = re.fullmatch(pat, c)
            if m:
                out.append((c, m.group(1)))
                break
    return out


def _gene_from_desc(x: Any) -> Optional[str]:
    if isinstance(x, str):
        m = _GN_RE.search(x)
        if m:
            return m.group(1)
    return None


def _normalize_protein_table(
    ptab: pd.DataFrame,
    *,
    filter_decoys: bool = True,
    min_psms: int = 0,
    min_unique: int = 0,
) -> pd.DataFrame:
    if ptab is None or ptab.empty:
        return pd.DataFrame(index=pd.Index([], name="accession"))
    ptab = ptab.copy().replace({"null": pd.NA, "NULL": pd.NA, "NaN": pd.NA})
    ptab = ptab.rename(columns=_rename_run_cols(list(ptab.columns)))
    if "accession" not in ptab.columns:
        raise ParseError("Protein table lacks 'accession'")
    ptab["accession"] = ptab["accession"].astype(str)
    ptab = ptab.set_index("accession", drop=True)
    if filter_decoys:
        ptab = ptab[~ptab.index.map(_is_decoy)]

    def _ok(row: pd.Series) -> bool:
        npsm = int(row.get("num_psms") or 0)
        nuni = int(row.get("num_peptides_unique") or 0)
        return (npsm >= min_psms) and (nuni >= min_unique)

    if min_psms or min_unique:
        ptab = ptab[ptab.apply(_ok, axis=1)]
    return ptab


def _select_matrix(
    ptab: pd.DataFrame, *, mztab_type: str, prefer_quant: bool = True
) -> Tuple[pd.DataFrame, Optional[str]]:
    assay_cols = _list_assay_cols(list(ptab.columns)) if not ptab.empty else []
    if prefer_quant and mztab_type.lower().startswith("quant") and assay_cols:
        cols, assays = zip(*assay_cols)
        mat = ptab.loc[:, list(cols)].copy()
        mat.columns = [f"assay_{a}" for a in assays]
        return mat.apply(pd.to_numeric, errors="coerce")

    prefer_cols = [
        c
        for c in ptab.columns
        if c.startswith(
            (
                "best_search_engine_score",
                "num_psms",
                "num_peptides_",
                "protein_coverage",
            )
        )
    ]
    fallback = [
        c
        for c in ["description", "search_engine", "database", "protein_coverage"]
        if c in ptab.columns
    ]
    keep = prefer_cols or fallback
    matrix = ptab.loc[:, keep].copy() if keep else pd.DataFrame(index=ptab.index)
    for c in matrix.columns:
        matrix[c] = pd.to_numeric(matrix[c], errors="ignore")
    return matrix


def _build_feature_meta_and_xrefs(
    ptab: pd.DataFrame,
) -> Tuple[List[FeatureRec], List[CrossRef]]:
    desc = ptab.get("description", pd.Series(index=ptab.index, dtype=object))
    amb = ptab.get("ambiguity_members", pd.Series(index=ptab.index, dtype=object))
    mods = ptab.get("modifications", pd.Series(index=ptab.index, dtype=object))
    taxs = ptab.get("taxid", pd.Series(index=ptab.index, dtype=object))

    feats: List[FeatureRec] = []
    xrefs: List[CrossRef] = []
    for acc in ptab.index:
        feats.append(
            FeatureRec(id=acc, entity=IDS.type.protein, namespace=IDS.ns.UniProtKB)
        )
        xrefs.append(CrossRef(acc, IDS.predicates.alias, acc))
        g = _gene_from_desc(desc.get(acc))
        if g:
            xrefs.append(CrossRef(acc, IDS.predicates.produced_by, g))
        amb_v = amb.get(acc)
        if isinstance(amb_v, str) and amb_v.strip() and str(amb_v).lower() != "null":
            for member in str(amb_v).split(","):
                if member.strip():
                    xrefs.append(
                        CrossRef(acc, IDS.predicates.contains, member.strip()
                        )
                    )
        mod_v = mods.get(acc)
        if isinstance(mod_v, str) and mod_v.strip() and str(mod_v).lower() != "null":
            for tok in str(mod_v).split(","):
                m = re.search(r"UNIMOD:\d+", tok)
                if m:
                    xrefs.append(CrossRef(acc, IDS.predicates.modification, m.group(0)))
        tax_v = taxs.get(acc)
        if pd.notna(tax_v):
            xrefs.append(CrossRef(acc, IDS.predicates.part_of, str(tax_v)))
    return feats, xrefs


def _parse_md_assays(md: Dict[str, Any]) -> Dict[str, str]:
    names = {}
    for k, v in (md or {}).items():
        m = _ASSAY_NAME_RE.search(str(k))
        if m:
            names[m.group(1)] = str(v)
    return names


def _build_column_meta(
    md: Dict[str, Any], assay_cols: List[Tuple[str, str]]
) -> Tuple[List[ColumnRec], List[CrossRef]]:
    assays = [a for _, a in assay_cols]
    assay_names = _parse_md_assays(md)

    cols: List[ColumnRec] = []
    for a in assays:
        sid = f"assay_{a}"
        attrs = {
            "assay_index": int(a) if str(a).isdigit() else a,
            "assay_name": assay_names.get(a),
        }
        cols.append(
            ColumnRec(
                id=sid,
                entity=IDS.type.sample,
                namespace=IDS.ns.SampleID,
                role="assay",
                attrs=attrs,
            )
        )
    return cols


def parse_mztab(
    path: str,
    *,
    prefer_quant: bool = True,
    filter_decoys: bool = True,
    min_psms: int = 0,
    min_unique: int = 0,
) -> OmicData:
    t = _open_mztab(path)
    md = t.metadata or {}

    ptab = _normalize_protein_table(
        t.protein_table,
        filter_decoys=filter_decoys,
        min_psms=min_psms,
        min_unique=min_unique,
    )
    matrix = _select_matrix(
        ptab, mztab_type=str(md.get("mzTab-type", "")), prefer_quant=prefer_quant
    )
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    feats, feat_xrefs = _build_feature_meta_and_xrefs(ptab)
    assay_cols = _list_assay_cols(list(ptab.columns))
    col_meta = _build_column_meta(md, assay_cols)

    col_meta = [
        ColumnRec(
            id=str(col),
            entity=("sample" if str(col).startswith("assay_") else "field"),
            namespace=("assay" if str(col).startswith("assay_") else "n/a"),
            role=("assay" if str(col).startswith("assay_") else "field"),
            attrs={}
        )
        for col in matrix.columns
    ]

    cross_ref: list[CrossRef] = feat_xrefs

    return OmicData(
        name="proteome_mztab",
        omics_type=IDS.omics_type.proteomics,
        matrix=matrix,
        feature_meta=feats,
        column_meta=col_meta,
        cross_ref=cross_ref,
    )
