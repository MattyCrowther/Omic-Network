from __future__ import annotations
import io, gzip
from typing import Optional, Dict, Any, List, Tuple
import pandas as pd
from ...errors import ParseError
from ...parse_obj import OmicData, FeatureRec, ColumnRec, CrossRef
from ...identifiers import IDS


COL_METAB = "metabolite_name"
COL_REFMET = "refmet_name"
COL_WB = "workbench_id"
COL_KEGG = "kegg_id"

NS_REFMET = IDS.ns.REFMET
NS_METAB = IDS.ns.METAB
NS_WB = IDS.ns.WB
NS_KEGG = IDS.ns.KEGG_Compound

FACTORS_ROW_TOKEN = "factors"
DEFAULT_SEP = "\t"

HEADER_ALIASES: Dict[str, str] = {
    "metabolite_name": COL_METAB,
    "metabolitename": COL_METAB,
    "metabolite": COL_METAB,
    "refmet_name": COL_REFMET,
    "refmetname": COL_REFMET,
    "refmet": COL_REFMET,
    "workbench_id": COL_WB,
    "workbenchid": COL_WB,
    "workbench": COL_WB,
    "kegg_id": COL_KEGG,
    "keggid": COL_KEGG,
    "kegg": COL_KEGG,
}
NON_SAMPLE_COLS = {COL_METAB, COL_REFMET}



def _read_table(path: str, sep: str) -> pd.DataFrame:
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                return pd.read_csv(fh, sep=sep, dtype=str)
        with open(path, "rt", encoding="utf-8") as fh:
            return pd.read_csv(fh, sep=sep, dtype=str)
    except Exception as e:
        raise ParseError(f"Failed to read table: {e}", path=path, cause=e)



def _normalize_headers(df: pd.DataFrame) -> pd.DataFrame:
    cols = []
    for c in df.columns:
        key = str(c).strip()
        lc = key.lower().replace(" ", "").replace("-", "_")
        cols.append(HEADER_ALIASES.get(lc, key))
    out = df.copy()
    out.columns = cols
    return out


def _pick_feature_name(df: pd.DataFrame) -> Tuple[pd.Series, str]:
    has_refmet = (COL_REFMET in df.columns) and df[COL_REFMET].notna().any()
    if has_refmet:
        s = df[COL_REFMET].fillna(
            df.get(COL_METAB, pd.Series(index=df.index, dtype=str))
        )
        ns = NS_REFMET
    else:
        if COL_METAB not in df.columns:
            raise ParseError(f"Missing required column '{COL_METAB}'")
        s, ns = df[COL_METAB], NS_METAB
    s = s.astype(str).str.strip()
    if s.duplicated().any():
        dups = s[s.duplicated()].unique()[:5]
        raise ParseError(f"Duplicate feature IDs after naming: {list(dups)}")
    s.name = "feature_id"
    return s, ns


def _extract_sample_meta(
    raw: pd.DataFrame,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame], List[str]]:
    first_col = raw.columns[0]
    factor_mask = (
        raw[first_col].astype(str).str.strip().str.lower().eq(FACTORS_ROW_TOKEN)
    )
    sample_cols = [c for c in raw.columns if c not in NON_SAMPLE_COLS]
    sample_meta = None
    df = raw
    if factor_mask.any():
        frow = raw.loc[factor_mask].iloc[0]
        records: Dict[str, Dict[str, str]] = {str(c): {} for c in sample_cols}
        for c in sample_cols:
            val = str(frow.get(c, "")).strip()
            if ":" in val:
                k, v = val.split(":", 1)
                k, v = k.strip(), v.strip()
                if k and v:
                    records[str(c)][k] = v
            elif val:
                records[str(c)]["factor"] = val
        sample_meta = pd.DataFrame.from_dict(records, orient="index")
        df = raw.loc[~factor_mask].reset_index(drop=True)
    return df, sample_meta, sample_cols


def _build_matrix(
    df: pd.DataFrame, sample_cols: List[str], chosen_ids: pd.Series
) -> pd.DataFrame:
    df = df.copy()
    df.index = chosen_ids
    drop_cols = [c for c in (COL_METAB, COL_REFMET) if c in df.columns]
    df = df.drop(columns=drop_cols)
    keep = [c for c in sample_cols if c in df.columns]
    mat = df[keep].apply(pd.to_numeric, errors="coerce")
    mat.index = mat.index.astype(str)
    mat.columns = pd.Index([str(c) for c in mat.columns], name="sample")
    return mat



def _load_mapping(mapping_path: Optional[str], sep: str) -> Optional[pd.DataFrame]:
    if not mapping_path:
        return None
    m = _read_table(mapping_path, sep=sep)
    if m.empty:
        return None
    m = _normalize_headers(m)
    for col in (COL_METAB, COL_REFMET, COL_WB, COL_KEGG):
        if col in m.columns:
            m[col] = m[col].astype(str).str.strip()
    return m


def _build_id_lookup(mapping: pd.DataFrame) -> Dict[Tuple[str, str], Dict[str, str]]:
    def _norm(s: pd.Series) -> pd.Series:
        return s.astype(str).str.strip().str.lower()

    lookups: Dict[Tuple[str, str], Dict[str, str]] = {}
    for key_col in (COL_REFMET, COL_METAB):
        if key_col in mapping.columns:
            k = _norm(mapping[key_col])
            if COL_WB in mapping.columns:
                lookups[(key_col, COL_WB)] = dict(zip(k, mapping[COL_WB].astype(str)))
            if COL_KEGG in mapping.columns:
                lookups[(key_col, COL_KEGG)] = dict(
                    zip(k, mapping[COL_KEGG].astype(str))
                )
    return lookups


def _emit_external_refs(
    raw_keys: pd.DataFrame, chosen_ids: pd.Series, mapping: Optional[pd.DataFrame]
) -> List[CrossRef]:
    if mapping is None or mapping.empty:
        return []

    def _norm(s: pd.Series) -> pd.Series:
        return s.astype(str).str.strip().str.lower()

    keys_cols = [c for c in (COL_METAB, COL_REFMET) if c in raw_keys.columns]
    keys = raw_keys[keys_cols].copy()
    keys.index = chosen_ids.astype(str)
    keys = keys.loc[chosen_ids.astype(str)]
    lookups = _build_id_lookup(mapping)
    key_ref = (
        _norm(keys[COL_REFMET])
        if COL_REFMET in keys.columns
        else pd.Series(index=keys.index, dtype=str)
    )
    key_met = (
        _norm(keys[COL_METAB])
        if COL_METAB in keys.columns
        else pd.Series(index=keys.index, dtype=str)
    )
    ids = pd.DataFrame(index=keys.index)
    for col in (COL_WB, COL_KEGG):
        s = pd.Series(index=keys.index, dtype=object)
        if (COL_REFMET, col) in lookups:
            s = key_ref.map(lookups[(COL_REFMET, col)])
        if (COL_METAB, col) in lookups:
            s = s.combine_first(key_met.map(lookups[(COL_METAB, col)]))
        ids[col] = s
    xrefs: List[CrossRef] = []
    for acc, row in ids.fillna("").iterrows():
        wb = str(row.get(COL_WB, "")).strip()
        kg = str(row.get(COL_KEGG, "")).strip()
        if wb:
            xrefs.append(CrossRef(str(acc), NS_WB, wb))
        if kg:
            xrefs.append(CrossRef(str(acc), NS_KEGG, kg))
    return xrefs



def _build_feature_meta(mat: pd.DataFrame, chosen_ns: str) -> List[FeatureRec]:
    ns = NS_REFMET if chosen_ns == NS_REFMET else NS_METAB
    return [
        FeatureRec(id=str(fid), entity=IDS.type.metabolite, namespace=ns, attrs={})
        for fid in mat.index
    ]


def _build_column_meta(
    mat: pd.DataFrame, sample_meta: Optional[pd.DataFrame]
) -> List[ColumnRec]:
    col_meta: List[ColumnRec] = []
    if sample_meta is not None:
        sample_meta = sample_meta.reindex(mat.columns)
    for sid in mat.columns:
        attrs = {}
        if sample_meta is not None and sid in sample_meta.index:
            attrs = {str(k): str(v) for k, v in sample_meta.loc[sid].dropna().items()}
        col_meta.append(
            ColumnRec(
                id=str(sid),
                entity=IDS.type.sample,
                namespace=IDS.ns.SampleID,
                role="sample",
                attrs=attrs,
            )
        )
    return col_meta


def _build_cross_refs(
    mat: pd.DataFrame,
    chosen_ns: str,
    external: List[CrossRef],
    sample_meta: Optional[pd.DataFrame] = None,
) -> List[CrossRef]:
    x: List[CrossRef] = [CrossRef(str(fid), chosen_ns, str(fid)) for fid in mat.index]
    x.extend(external)
    if sample_meta is not None and not sample_meta.empty:
        sample_meta = sample_meta.reindex(mat.columns)
        for sid in mat.columns:
            if sid in sample_meta.index:
                for k, v in sample_meta.loc[sid].dropna().items():
                    x.append(CrossRef(str(sid), str(k), str(v)))
    return x



def parse_metabolomics(
    intensity_path: str,
    *,
    mapping_path: Optional[str] = None,
    sep: str = DEFAULT_SEP,
) -> OmicData:
    raw = _read_table(intensity_path, sep=sep)
    if raw.empty:
        raise ParseError("Intensity table is empty", path=intensity_path)
    raw = _normalize_headers(raw)
    if COL_METAB not in raw.columns and COL_REFMET not in raw.columns:
        raise ParseError(
            f"Missing required column '{COL_METAB}' (or '{COL_REFMET}')",
            path=intensity_path,
        )

    df, sample_meta, sample_cols = _extract_sample_meta(raw)
    chosen_ids, chosen_ns = _pick_feature_name(df)
    mat = _build_matrix(df, sample_cols, chosen_ids)
    if mat.index.has_duplicates or mat.columns.has_duplicates:
        raise ValueError("Duplicate feature or sample IDs in matrix")

    mapping = _load_mapping(mapping_path, sep) if mapping_path else None
    external_refs = _emit_external_refs(df, chosen_ids, mapping)

    feature_meta = _build_feature_meta(mat, chosen_ns)
    column_meta  = _build_column_meta(mat, sample_meta)
    cross_ref    = _build_cross_refs(mat, chosen_ns, external_refs, sample_meta)

    return OmicData(
        name="metabolomics_matrix",
        omics_type=IDS.omics_type.metabolomics,
        matrix=mat,
        feature_meta=feature_meta,
        column_meta=column_meta,
        cross_ref=cross_ref or None,
    )
