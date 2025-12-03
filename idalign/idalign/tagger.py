from omics_io.omics_io.identifiers import IDS

TYPE_MAP = {
    # Aliases (unique IDs / synonyms per entity)
    "asap": "alias",
    "ecocyc": "alias",
    "geneid": "alias",
    "id": "alias",
    "kegg.compound": "alias",
    "name": "alias",
    "refmet": "alias",
    "workbench": "alias",
    "alias": "alias",
    "gene": "alias",
    "locus_tag": "alias",
    "run": "alias",            # SRR accession (unique)
    "srx": "alias",            # SRX accession (unique)

    # Relations (shared target nodes → structure, not equivalence)
    "geo": "relation",         # GSE series ID (shared parent)
    "genotype": "relation",    # shared genotype label/entity
    IDS.predicates.modification: "relation",
    IDS.predicates.produced_by: "relation",
    IDS.predicates.contains: "relation",
    IDS.predicates.has_feature: "relation",
    IDS.predicates.product: "relation",
    "product_label": "relation",
    "protein_id": "relation",
    "uniprotkb/swiss-prot": "relation",
    "parent": "relation",
    "genbank": "relation",
}
TYPE_MAP = {
    "asap": "alias",
    "ecocyc": "alias",
    "geneid": "alias",
    "id": "alias",
    "kegg.compound": "alias",
    "name": "alias",
    "refmet": "alias",
    "workbench": "alias",
    "alias": "alias",
    "gene": "alias",
    "locus_tag": "alias",
    "run": "alias",
    "srx": "alias", 
    "geo": "relation",
    "genotype": "relation",
    IDS.predicates.modification: "relation",
    IDS.predicates.produced_by: "relation",
    IDS.predicates.contains: "relation",
    IDS.predicates.has_feature: "relation",
    IDS.predicates.product: "relation",
    "product_label": "relation",
    "protein_id": "relation",
    "uniprotkb/swiss-prot": "relation",
    "parent": "relation",
    "genbank": "relation",
}


_REL_TYPE_MAP = {
    "protein_id": IDS.type.protein,
    "uniprotkb/swiss-prot": IDS.type.protein,
    "kegg.compound": IDS.type.metabolite,
    "modification": lambda a: a,
    "parent": IDS.type.dna,
    "contains": lambda a: a,
    "produced_by": lambda a: IDS.type.protein if a == IDS.type.gene else "unknown",
    "product_label": IDS.type.protein,
    "product": IDS.type.protein,
    "geo" : IDS.type.study
}

REGISTRY_NS = {
    "geneid","uniprot","uniprotkb","ensembl","refseq","ecocyc","asap",
    "chebi","kegg","hmdb","pubchem","genbank_reference","refseq.locus_tag","kegg.compound",
}
LOCAL_NS = {
    "gene","gene_name","gene_symbol","symbol","synonym","gene_synonym",
    "locus_tag","name","id","workbench",
}

REL_NS_TO_REGISTRY = {
    "uniprotkb/swiss-prot": "uniprot",
    "protein_id": "protein",
    "genbank": "dna",
    "parent": "dna",
    "contains": "dna",
    "modification": "metabolite",
    "produced_by": "protein",
    "product_label": "protein",
}

REL_NS_TO_PREDICATE = {
    "uniprotkb/swiss-prot": "product",
    "protein_id": "product",
    "produced_by": "produced_by",
    "product_label": "product",
    "parent": "part_of",
    "contains": "contains",
    "genbank": "contains",
    "modification": "modification",
}

def tag_record(namespace):
    if namespace.lower() not in TYPE_MAP:
        return "unclassified"
    return TYPE_MAP[namespace.lower()]

def derive_type(namespace: str, source_type: str) -> str:
    rule = _REL_TYPE_MAP.get(namespace.lower())
    if rule is None:
        return "unknown"
    return rule(source_type) if callable(rule) else rule

def normalize_feature_scope(od, ns: str | None) -> str:
    """Feature scope: registries are global; 'local' labels map to dataset scope; else use given ns."""
    ns = (ns or "").strip().lower()
    if not ns:
        return od.name.strip().lower()
    if ns in REGISTRY_NS or "." in ns:
        return ns
    if ns in LOCAL_NS:
        return od.name.strip().lower()
    return ns

def alias_target_scope(rel_ns: str, src_scope: str, fallback: str) -> str:
    """Alias targets: registry/specific ns → that registry; local labels stay 
    in the source’s dataset scope."""
    ns = (rel_ns or "").lower()
    if ns.startswith("db_xref:"):
        return ns.split(":", 1)[1] or fallback
    if ns in REGISTRY_NS:
        return ns
    if ns in LOCAL_NS:
        return src_scope
    return fallback

def relation_target_scope(rel_ns: str, fallback: str) -> str:
    ns = (rel_ns or "").lower()
    if ns in REL_NS_TO_REGISTRY:
        return REL_NS_TO_REGISTRY[ns]
    if ns in REGISTRY_NS:
        return ns
    if ns.startswith("db_xref:"):
        return ns.split(":", 1)[1] or fallback
    return fallback