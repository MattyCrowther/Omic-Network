
from __future__ import annotations

from typing import Dict, Iterable, List, Set, Tuple
from collections import defaultdict

from .mapping_data import AlignmentResult, Node
from .tagger import (
    tag_record,
    derive_type,
    alias_target_scope,
    normalize_feature_scope,
    relation_target_scope,
    REL_NS_TO_PREDICATE,
)
from .union_find import UnionFind
from omics_io.omics_io.parse_obj import OmicData

# ---------------- Types ----------------
ScopedID = Tuple[str, str]
AliasKey = Tuple[ScopedID, ScopedID]
RelTriple = Tuple[ScopedID, str, ScopedID]


# ---------------- Small helpers ----------------
def key(scope: str, ident: str) -> ScopedID:
    return (str(scope), str(ident).strip())


def alias_ok(a, b, id_type, inferred):
    ta = id_type.get(a) or inferred.get(a)
    tb = id_type.get(b) or inferred.get(b)
    
    if ta is None and tb is None:
        return False
    return ta is None or tb is None or ta == tb


# ---------------- Stage 1: feature meta ----------------
def collect_feature_scopes_and_types(
    datasets: Iterable[OmicData],
) -> Tuple[Dict[Tuple[str, str], str], Dict[ScopedID, str]]:
    feat_scope: Dict[Tuple[str, str], str] = {}
    id_type: Dict[ScopedID, str] = {}
    for od in datasets:
        for f in od.feature_meta:
            scope = normalize_feature_scope(od, f.namespace)
            feat_scope[(od.name, f.id)] = scope
            if f.entity:
                id_type[key(scope, f.id)] = f.entity

        for f in od.column_meta:
            scope = normalize_feature_scope(od, f.namespace)
            feat_scope[(od.name, f.id)] = scope
            if f.entity:
                id_type[key(scope, f.id)] = f.entity

    return feat_scope, id_type


# ---------------- Stage 2: cross-ref ingest ----------------
def ingest_cross_refs(
    datasets: Iterable[OmicData],
    feat_scope: Dict[Tuple[str, str], str],
    id_type: Dict[ScopedID, str],
) -> Tuple[Set[AliasKey], List[RelTriple], List[RelTriple], Dict[ScopedID, str]]:
    alias_edges: Set[AliasKey] = set()
    rel_edges: List[RelTriple] = []
    unknown_edges: List[RelTriple] = []
    inferred_type: Dict[ScopedID, str] = {}

    for od in datasets:
        ds_default_scope = normalize_feature_scope(od, None)
        for cr in od.cross_ref:
            src_raw = str(cr.src).strip()
            tgt_raw = str(cr.target).strip()
            rel_ns = str(cr.namespace).lower().strip()
            tag = tag_record(rel_ns)

            s_src = feat_scope.get((od.name, src_raw), ds_default_scope)

            if tag == "alias":
                s_tgt_fb = feat_scope.get((od.name, tgt_raw), ds_default_scope)
                s_tgt = alias_target_scope(rel_ns, s_src, s_tgt_fb)

                k_src, k_tgt = key(s_src, src_raw), key(s_tgt, tgt_raw)

                # propagate explicit types across alias
                ts, tt = id_type.get(k_src), id_type.get(k_tgt)
                if ts and k_tgt not in id_type:
                    inferred_type[k_tgt] = ts
                if tt and k_src not in id_type:
                    inferred_type[k_src] = tt

                a, b = sorted((k_src, k_tgt))
                alias_edges.add((a, b))
                continue

            # relation or unknown
            s_tgt_fb = feat_scope.get((od.name, tgt_raw), ds_default_scope)
            s_tgt = relation_target_scope(rel_ns, s_tgt_fb)
            k_src, k_tgt = key(s_src, src_raw), key(s_tgt, tgt_raw)

            if tag == "relation":
                predicate = REL_NS_TO_PREDICATE.get(rel_ns, rel_ns)
                src_t = id_type.get(k_src, "unknown")
                tgt_t = derive_type(rel_ns, src_t)
                if k_tgt not in id_type and tgt_t != "unknown":
                    inferred_type[k_tgt] = tgt_t
                rel_edges.append((k_src, predicate, k_tgt))
            else:
                print(f'Unknown Edge: {k_src,rel_ns,k_tgt}')
                unknown_edges.append((k_src, rel_ns, k_tgt))

    return alias_edges, rel_edges, unknown_edges, inferred_type


# ---------------- Stage 3: union-find build ----------------
def build_union_find(
    alias_edges: Set[AliasKey],
    rel_edges: List[RelTriple],
    id_type: Dict[ScopedID, str],
    inferred_type: Dict[ScopedID, str],
) -> Tuple[UnionFind, Set[ScopedID]]:
    uf = UnionFind()
    all_keys: Set[ScopedID] = set()

    def touch(k: ScopedID) -> None:
        all_keys.add(k)
        uf.find(k)

    for a, b in alias_edges:
        touch(a)
        touch(b)
        if alias_ok(a, b, id_type, inferred_type):
            uf.union(a, b)
    for a, _, b in rel_edges:
        touch(a)
        touch(b)

    return uf, all_keys


# ---------------- Stage 4: resolve types ----------------
def resolve_component_types(
    uf: UnionFind,
    all_keys: Set[ScopedID],
    id_type: Dict[ScopedID, str],
    inferred_type: Dict[ScopedID, str],
) -> Dict[ScopedID, str]:
    root_to_members: Dict[ScopedID, List[ScopedID]] = defaultdict(list)
    for ksid in all_keys:
        root_to_members[uf.find(ksid)].append(ksid)

    final_type: Dict[ScopedID, str] = {}
    for _, members in root_to_members.items():
        known = {id_type.get(m) or inferred_type.get(m) for m in members}
        known.discard(None)
        ty = next(iter(known)) if len(known) == 1 else "unknown"
        for m in members:
            final_type[m] = ty
    return final_type


# ---------------- Stage 5: materialize ----------------
def materialize_alignment(
    uf: UnionFind,
    all_keys: Set[ScopedID],
    final_type: Dict[ScopedID, str],
    rel_edges: List[RelTriple],
    unknown_edges,
) -> AlignmentResult:
    roots: Dict[ScopedID, int] = {}
    next_gid = 0

    def gid_for(k: ScopedID) -> int:
        nonlocal next_gid
        r = uf.find(k)
        if r not in roots:
            roots[r] = next_gid
            next_gid += 1
        return roots[r]

    node_to_gid: Dict[Node, int] = {}
    groups: Dict[int, Set[Node]] = defaultdict(set)

    for ksid in all_keys:
        ns, ident = ksid
        n = Node(
            identifier=str(ident),
            namespace=str(ns),
            type=final_type.get(ksid, "unknown"),
        )
        gid = gid_for(ksid)
        node_to_gid[n] = gid
        groups[gid].add(n)

    group_relations: Dict[Tuple[int, str, int], int] = defaultdict(int)
    for a, predicate, b in rel_edges:
        ga, gb = gid_for(a), gid_for(b)
        group_relations[(ga, predicate, gb)] += 1
    
    return AlignmentResult(
        node_to_gid=node_to_gid,
        groups=dict(groups),
        group_relations=dict(group_relations),
        unknown_edges=unknown_edges,
    )


def match_references(*omic_data_sets: OmicData) -> AlignmentResult:
    feat_scope, id_type = collect_feature_scopes_and_types(omic_data_sets)
    alias_edges, rel_edges, unknown_edges, inferred_type = ingest_cross_refs(
        omic_data_sets, feat_scope, id_type
    )
    uf, all_keys = build_union_find(alias_edges, rel_edges, id_type, inferred_type)
    final_type = resolve_component_types(uf, all_keys, id_type, inferred_type)
    return materialize_alignment(uf, all_keys, final_type, rel_edges, unknown_edges)
