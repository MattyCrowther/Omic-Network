import unittest
import sys
import tempfile
from pathlib import Path
import pandas as pd
from collections import Counter

sys.path.insert(0, "..")
sys.path.insert(0, "../..")


from idalign.aligner import match_references
from idalign.mapping_data import Node, AlignmentResult
from omics_io.omics_io.parsers.reference.genbank import parse_genbank
from omics_io.omics_io.parsers.reference.gff3 import parse_gff3
from omics_io.omics_io.parsers.transcriptomics.precise2 import parse_precise2
from omics_io.omics_io.parsers.proteomics.mztab import parse_mztab
from omics_io.omics_io.parsers.metabolomics.table import parse_metabolomics
from idalign.aligner import match_references


from collections import Counter, defaultdict
from dataclasses import asdict
from typing import Dict, Set, Tuple, List, Iterable, Optional
from pathlib import Path

# assumes: parse_genbank, match_references, AlignmentResult, Node imported

def summarize_alignment(res) -> Dict[str, float]:
    n_groups = len(res.groups)
    sizes = [len(m) for m in res.groups.values()]
    n_nodes = sum(sizes)
    singletons = sum(1 for s in sizes if s == 1)
    max_size = max(sizes) if sizes else 0
    unknown_nodes = sum(1 for n in res.node_to_gid if getattr(n, "type", "unknown") == "unknown")
    return {
        "n_groups": n_groups,
        "n_nodes": n_nodes,
        "mean_group_size": (n_nodes / n_groups) if n_groups else 0.0,
        "singletons": singletons,
        "max_group_size": max_size,
        "unknown_type_nodes": unknown_nodes,
        "unknown_type_frac": (unknown_nodes / n_nodes) if n_nodes else 0.0,
        "n_rel_triplets": len(res.group_relations),
        "n_rel_edges_total": sum(res.group_relations.values()),
    }

def group_size_hist(res) -> Dict[int, int]:
    return dict(Counter(len(m) for m in res.groups.values()))

def namespaces_per_group(res) -> Dict[int, Set[str]]:
    out = {}
    for gid, members in res.groups.items():
        out[gid] = {n.namespace for n in members}
    return out

def mixed_namespace_groups(res) -> List[int]:
    return [gid for gid, nss in namespaces_per_group(res).items() if len(nss) > 1]

def groups_with_unknown_type(res) -> List[int]:
    out = []
    for gid, members in res.groups.items():
        types = {getattr(n, "type", "unknown") for n in members}
        if len(types) != 1 or "unknown" in types:
            out.append(gid)
    return out

def top_groups_by_size(res, k: int = 10) -> List[Tuple[int, int]]:
    sizes = [(gid, len(members)) for gid, members in res.groups.items()]
    return sorted(sizes, key=lambda x: x[1], reverse=True)[:k]

def print_group(res, gid: int, limit: Optional[int] = 20) -> None:
    members = sorted(res.groups[gid], key=lambda n: (n.namespace, n.identifier))
    print(f"G{gid} size={len(members)} type={_group_type(res, gid)}")
    for i, n in enumerate(members):
        if limit is not None and i >= limit:
            print("  ...")
            break
        print(f"  {n.namespace}\t{n.identifier}\t{getattr(n,'type','unknown')}")

def _group_type(res, gid: int) -> str:
    types = {getattr(n, "type", "unknown") for n in res.groups[gid]}
    return next(iter(types)) if len(types) == 1 else "mixed/unknown"

def relation_stats(res) -> Dict[str, int]:
    by_rel = defaultdict(int)
    for (_, rel, _), c in res.group_relations.items():
        by_rel[rel] += c
    return dict(sorted(by_rel.items(), key=lambda x: x[1], reverse=True))

def degree_by_group(res) -> Dict[int, Dict[str, int]]:
    out = {gid: {"out": 0, "in": 0} for gid in res.groups}
    for (ga, _, gb), c in res.group_relations.items():
        out[ga]["out"] += c
        out[gb]["in"] += c
    return out

def verify_cross_dataset_alias(res, namespace: str, ident: str) -> List[int]:
    """Return all gids that contain (namespace, ident). Expect length 1 if merged."""
    gids = []
    for gid, members in res.groups.items():
        if any((n.namespace == namespace and n.identifier == ident) for n in members):
            gids.append(gid)
    return gids
def assert_unique_scoped_ids(res):
    seen = {}
    for gid, members in res.groups.items():
        for n in members:
            k = (n.namespace, n.identifier)
            if k in seen and seen[k] != gid:
                raise AssertionError(f"Scoped ID {k} in multiple groups: {seen[k]} and {gid}")
            seen[k] = gid
    return len(seen)

def same_group_relations(res):
    return [ (ga, rel, gb, c)
             for (ga, rel, gb), c in res.group_relations.items()
             if ga == gb ]

def groups_missing_registry(res, registry="geneid"):
    return [ (gid, sorted((n.namespace, n.identifier) for n in members))
             for gid, members in res.groups.items()
             if any(n.namespace in {"gene","locus_tag"} for n in members)
             and all(n.namespace != registry for n in members) ]

def coverage_by_registry(res):
    from collections import Counter
    cnt = Counter()
    for gid, members in res.groups.items():
        regs = {n.namespace for n in members}
        for r in regs:
            cnt[r] += 1
    return dict(cnt)



# --- quick smoke test ---
if __name__ == "__main__":
    seq = parse_genbank(Path("data/sequence.gb"))
    res = match_references(seq)

    print("Summary:", summarize_alignment(res))
    print("Group size hist:", group_size_hist(res))
    print("Top groups:", top_groups_by_size(res, k=5))
    print("Relation stats:", relation_stats(res))

    # Inspect the largest group
    tg = top_groups_by_size(res, 1)
    if tg:
        print_group(res, tg[0][0], limit=30)

    # Example cross-merge check for GeneID
    gids = verify_cross_dataset_alias(res, "geneid", "947440")
    print("GeneID 947440 GIDs:", gids)

    # Example use
    print("unique scoped IDs:", assert_unique_scoped_ids(res))
    print("same-group edges:", same_group_relations(res))
    print("groups missing GeneID (first 5):", groups_missing_registry(res)[:5])
    print("registry coverage (top 10):", sorted(coverage_by_registry(res).items(), key=lambda x: x[1], reverse=True)[:10])