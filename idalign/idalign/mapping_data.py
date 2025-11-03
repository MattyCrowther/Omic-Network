from dataclasses import dataclass, field
from typing import Dict, Set, Tuple, List, Optional, Iterator
import pandas as pd
import os
import json
import datetime as dt


@dataclass(frozen=True)
class Node:
    identifier: str
    type: str 
    namespace: str

@dataclass
class AlignmentResult:
    """Unified alignment result connecting identifiers across omics datasets.

    Represents all alias groups and their inter-group relations as derived
    from one or more parsed cross-reference tables.

    Attributes:
        node_to_gid (dict[Node, int]): Mapping from node to its assigned group ID.
        group_relations (dict[tuple[int, str, int], int]): Relation edges between
            groups, with counts for how many times each was observed.
        groups (dict[int, set[Node]]): Grouped sets of alias-equivalent nodes.
        alias_meta (dict[tuple[Node, Node], AliasEvidence]): Optional alias evidence
            metadata between specific node pairs.
    """

    node_to_gid: Dict[str, int]                     
    groups: Dict[int, dict]                        
    group_relations: Dict[Tuple[int, str, int], int]
    unknown_edges: List[Tuple[Node, str, Node]] = None
    
    def gid_of(self, namespace: str, ident: str) -> Optional[int]:
        """Return the group ID containing the given identifier, if present.

        Args:
            namespace: Namespace of the identifier.
            ident: Identifier value.

        Returns:
            int | None: The numeric group ID, or None if not present.
        """
        return self.node_to_gid.get(Node(namespace, ident))

    def members(self, gid: int) -> Set[Node]:
        """Return all member nodes belonging to a given group ID.

        Args:
            gid: Group identifier.

        Returns:
            set[Node]: All nodes unified in that alias group.
        """
        return self.groups.get(gid, set())

    def relations_from(self, gid: int) -> List[Tuple[str, int, int]]:
        """Return all outgoing relations originating from a group.

        Args:
            gid: Source group ID.

        Returns:
            list[tuple[str, int, int]]: Triplets of (relation_label, target_gid, count).
        """
        return [
            (rel, g2, cnt)
            for (g1, rel, g2), cnt in self.group_relations.items()
            if g1 == gid
        ]

    def relations_to(self, gid: int) -> List[Tuple[str, int, int]]:
        """Return all incoming relations targeting a group.

        Args:
            gid: Target group ID.

        Returns:
            list[tuple[str, int, int]]: Triplets of (relation_label, source_gid, count).
        """
        return [
            (rel, g1, cnt)
            for (g1, rel, g2), cnt in self.group_relations.items()
            if g2 == gid
        ]

    def structural_orphans(self) -> Set[int]:
        """Return group IDs with no incoming or outgoing relations.

        These are groups that are completely disconnected at the relation level,
        even if they contain multiple aliases internally.
        """
        linked = {g1 for (g1, _, _) in self.group_relations} | {
            g2 for (_, _, g2) in self.group_relations
        }
        return set(self.groups.keys()) - linked

    def isolated_groups(self) -> Set[int]:
        """Return group IDs containing only a single node.

        These represent entities that have no alias relationships.
        """
        return {gid for gid, nodes in self.groups.items() if len(nodes) == 1}

    def connected_groups(self) -> Set[int]:
        """Return group IDs that participate in at least one relation."""
        linked = {g1 for (g1, _, _) in self.group_relations} | {
            g2 for (_, _, g2) in self.group_relations
        }
        return linked

    def to_frames(self):
        """Convert alias and relation data into DataFrame representations.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]:
                - alias_df: with columns [gid, namespace, ident]
                - rel_df: with columns [gid_src, relation, gid_tgt, count]
        """
        alias_rows = [
            {"gid": g, "namespace": n.namespace, "ident": n.ident}
            for g, nodes in self.groups.items()
            for n in nodes
        ]
        rel_rows = [
            {"gid_src": a, "relation": r, "gid_tgt": b, "count": c}
            for (a, r, b), c in self.group_relations.items()
        ]
        return pd.DataFrame(alias_rows), pd.DataFrame(rel_rows)

    @classmethod
    def from_frames(cls, alias_df, rel_df):
        """Reconstruct an AlignmentResult from alias and relation DataFrames.

        Args:
            alias_df (pd.DataFrame): DataFrame of alias mappings.
            rel_df (pd.DataFrame): DataFrame of group relations.

        Returns:
            AlignmentResult: A reconstructed instance.
        """
        node_to_gid: Dict[Node, int] = {}
        groups: Dict[int, Set[Node]] = {}

        for row in alias_df.itertuples(index=False):
            node = Node(str(row.namespace), str(row.ident))
            gid = int(row.gid)
            node_to_gid[node] = gid
            groups.setdefault(gid, set()).add(node)

        group_relations: Dict[Tuple[int, str, int], int] = {}
        for row in rel_df.itertuples(index=False):
            src = int(row.gid_src)
            rel = str(row.relation)
            tgt = int(row.gid_tgt)
            count = int(row.count)
            group_relations[(src, rel, tgt)] = count

        return cls(
            node_to_gid=node_to_gid, group_relations=group_relations, groups=groups
        )

    @classmethod
    def load(cls, path: str) -> "AlignmentResult":
        """Load an AlignmentResult from a saved directory.

        Args:
            path: Directory containing saved parquet and metadata files.

        Returns:
            AlignmentResult: Loaded alignment result.
        """
        alias_df = pd.read_parquet(os.path.join(path, "alias.parquet"))
        rel_df = pd.read_parquet(os.path.join(path, "relations.parquet"))
        return cls.from_frames(alias_df, rel_df)

    def save(self, path):
        """Save the current alignment result to a directory.

        Args:
            path (str or Path): Directory path for output files.

        Notes:
            Produces two parquet files ('alias.parquet', 'relations.parquet')
            and one JSON metadata file ('meta.json') containing version info,
            creation timestamp, and summary counts.
        """
        os.makedirs(path, exist_ok=True)
        alias_df, rel_df = self.to_frames()
        alias_df.to_parquet(os.path.join(path, "alias.parquet"))
        rel_df.to_parquet(os.path.join(path, "relations.parquet"))
        meta = {
            "version": "1",
            "created": dt.datetime.now(dt.timezone.utc).isoformat(),
            "n_groups": len(self.groups),
            "n_relations": len(self.group_relations),
        }
        with open(os.path.join(path, "meta.json"), "w") as f:
            json.dump(meta, f)


    def __iter__(self) -> Iterator[Tuple[int, Set[Node]]]:
        """Iterate over groups only: (gid, nodes)."""
        for gid, nodes in self.groups.items():
            yield gid, nodes

    def relationships(self) -> Iterator[Tuple[int, str, int, int]]:
        """Yield all edges as (gid_src, relation, gid_tgt, count)."""
        for (g1, rel, g2), cnt in self.group_relations.items():
            yield g1, rel, g2, cnt