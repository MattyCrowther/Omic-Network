from __future__ import annotations
import pandas as pd
from typing import Optional, List
from ...errors import ParseError
from ...parse_obj import OmicData, FeatureRec, ColumnRec, CrossRef
from ...identifiers import IDS
from ...utils import to_primitive

DROP_COLS = {
    "contact", "creator", "run_date", "alignment",
    "R1", "R2", "BAM", "Notes",
    "passed_fastqc", "passed_pct_reads_mapped",
    "passed_reads_mapped_to_CDS", "passed_global_correlation",
    "passed_similar_replicates", "passed_number_replicates",
    "n_replicates", "Published", "full_name",
    "LibraryLayout", "Platform",
}
REF_COLS = ["GEO", "SRX", "SRR", "Run", "BioSample", "BioProject"]


def _load_matrix(path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, index_col=0)
    except Exception as e:
        raise ParseError(f"Failed to load matrix: {e}", path=path)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.index.has_duplicates or df.columns.has_duplicates:
        raise ValueError("Duplicate feature or sample IDs in matrix")
    return df

def _build_feature_meta(expr: pd.DataFrame) -> List[FeatureRec]:
    """Build one FeatureRec per locus_tag."""
    return [
        FeatureRec(
            id=str(fid),
            entity=IDS.type.rna,
            namespace=IDS.ns.RefSeq_Locus,
            attrs={},
        )
        for fid in expr.index
    ]

def _load_sample_metadata(path: str, sample_ids: List[str]) -> pd.DataFrame:
    try:
        smeta = pd.read_csv(path, index_col=0, dtype=str)
    except Exception as e:
        raise ParseError(f"Failed to load sample metadata: {e}", path=path)
    smeta.index = smeta.index.astype(str)
    smeta = smeta.reindex(index=sample_ids)
    smeta.index.name = "sample"
    keep_cols = [c for c in smeta.columns if c not in DROP_COLS]
    return smeta[keep_cols]

def _build_column_meta_and_xrefs(smeta: pd.DataFrame) -> tuple[List[ColumnRec], List[CrossRef]]:
    """Create ColumnRec list and sample CrossRef entries."""
    column_meta: List[ColumnRec] = []
    xrefs: List[CrossRef] = []

    xref_fields = [c for c in smeta.columns if c in REF_COLS]
    meta_fields = [c for c in smeta.columns if c not in REF_COLS]

    for sid, row in smeta.iterrows():
        attrs = {k: to_primitive(v) for k, v in row.items() if pd.notna(v) and k in meta_fields}
        column_meta.append(
            ColumnRec(
                id=str(sid),
                entity=IDS.type.sample,
                namespace=IDS.ns.SampleID,
                attrs=attrs,
            )
        )
        for ns in xref_fields:
            val = row.get(ns)
            if pd.notna(val) and str(val).strip():
                xrefs.append(CrossRef(str(sid), ns, str(val).strip()))
    return column_meta, xrefs

def _feature_selfrefs(dataframe: pd.DataFrame) -> List[CrossRef]:
    """Feature self-mapping cross-refs."""
    refs = []

    '''
    for gene in dataframe.index:
        for col in dataframe.columns:
            if dataframe.at[gene, col] != 0:
                refs.append(CrossRef(gene, IDS.predicates.has_feature, col))
    '''
    return refs

def parse_precise2(
    expression_path: str,
    *,
    presence_path: Optional[str] = None,
    metadata_path: Optional[str] = None,
) -> OmicData:
    """PRECISE2 → OmicData (modular version)."""
    expr = _load_matrix(expression_path)
    features = _build_feature_meta(expr)
    cross_ref = _feature_selfrefs(expr)

    column_meta: List[ColumnRec] = []
    if metadata_path:
        smeta = _load_sample_metadata(metadata_path, list(expr.columns))
        col_meta, sample_refs = _build_column_meta_and_xrefs(smeta)
        column_meta.extend(col_meta)
        cross_ref.extend(sample_refs)
    else:
        column_meta = [
            ColumnRec(
                id=str(cid),
                entity=IDS.type.sample,
                namespace=IDS.ns.SampleID,
            )
            for cid in expr.columns
        ]

    return OmicData(
        name="precise2_transcriptome",
        omics_type=IDS.omics_type.transcriptomics,
        matrix=expr,
        feature_meta=features,
        column_meta=column_meta,
        cross_ref=cross_ref,
    )
