from __future__ import annotations
from dataclasses import dataclass
from typing import Final


@dataclass(frozen=True, slots=True)
class _Entity:
    dna: Final[str] = "dna"
    gene: Final[str] = "gene"
    cds: Final[str] = "cds"
    rna: Final[str] = "rna"
    protein: Final[str] = "protein"
    metabolite: Final[str] = "metabolite"
    sample: Final[str] = "sample"
    other: Final[str] = "other"
    unknown: Final[str] = "unknown"


@dataclass(frozen=True, slots=True)
class _Role:
    measurement: Final[str] = "measurement"
    covariate: Final[str] = "covariate"
    qc: Final[str] = "qc"
    derived: Final[str] = "derived"


@dataclass(frozen=True, slots=True)
class _Predicate:
    produced_by: Final[str] = "produced_by"
    contains: Final[str] = "contains"
    modification: Final[str] = "modification"
    part_of: Final[str] = "part_of"
    has_feature: Final[str] = "has_feature"
    confidence: Final[str] = "confidence"
    alias: Final[str] = "aliases"
    product: Final[str] = "product"

@dataclass(frozen=True, slots=True)
class _OmicsType:
    genomics: Final[str] = "genomics"
    transcriptomics: Final[str] = "transcriptomics"
    proteomics: Final[str] = "proteomics"
    metabolomics: Final[str] = "metabolomics"
    epigenomics: Final[str] = "epigenomics"

@dataclass(frozen=True, slots=True)
class _Namespace:

    Ensembl_Gene: Final[str] = "Ensembl.Gene"
    Ensembl_Transcript: Final[str] = "Ensembl.Transcript"
    RefSeq_Locus: Final[str] = "RefSeq.locus_tag"
    SGD_Yeast: Final[str] = "SGD.Orf"

    UniProtKB: Final[str] = "UniProtKB"

    HMDB: Final[str] = "HMDB"
    KEGG_Compound: Final[str] = "KEGG.Compound"

    SampleID: Final[str] = "sample_id"

    NA: Final[str] = "n/a"
    REFMET: Final[str] = "RefMet"
    METAB: Final[str] = "Metabolite"
    WB: Final[str] = "WorkBench"


@dataclass(frozen=True, slots=True)
class Identifiers:
    """Dot-access constants. No runtime cost, no typos."""
    type: _Entity = _Entity()
    role: _Role = _Role()
    ns: _Namespace = _Namespace()
    predicates: _Predicate = _Predicate()
    omics_type: _OmicsType = _OmicsType()


IDS = Identifiers()
