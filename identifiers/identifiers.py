from __future__ import annotations
from dataclasses import dataclass
from typing import Final

@dataclass(frozen=True, slots=True)
class _Entity:
    dna: Final[str] = "dna"
    gene: Final[str] = "gene"
    transcript: Final[str] = "transcript"
    protein: Final[str] = "protein"
    metabolite: Final[str] = "metabolite"
    sample: Final[str] = "sample"
    other: Final[str] = "other"

@dataclass(frozen=True, slots=True)
class _Role:
    measurement: Final[str] = "measurement"
    covariate: Final[str] = "covariate"
    qc: Final[str] = "qc"
    derived: Final[str] = "derived"

@dataclass(frozen=True, slots=True)
class _OmicsType:
    genomics: Final[str] = "genomics"
    transcriptomics: Final[str] = "transcriptomics"
    proteomics: Final[str] = "proteomics"
    metabolomics: Final[str] = "metabolomics"
    epigenomics: Final[str] = "epigenomics"

@dataclass(frozen=True, slots=True)
class _Namespace:
    # genes / transcripts
    Ensembl_Gene: Final[str] = "Ensembl.Gene"
    Ensembl_Transcript: Final[str] = "Ensembl.Transcript"
    RefSeq_Locus: Final[str] = "RefSeq.locus_tag"
    SGD_Yeast: Final[str] = "SGD.Orf"
    # proteins
    UniProtKB: Final[str] = "UniProtKB"
    # metabolites
    HMDB: Final[str] = "HMDB"
    KEGG_Compound: Final[str] = "KEGG.Compound"
    # samples
    SampleID: Final[str] = "sample_id"
    # generic
    NA: Final[str] = "n/a"

@dataclass(frozen=True, slots=True)
class Identifiers:
    """Dot-access constants. No runtime cost, no typos."""
    type: _Entity = _Entity()
    role: _Role = _Role()
    ns: _Namespace = _Namespace()
    omics_type: _OmicsType = _OmicsType()

IDS = Identifiers()