class UnionFind:
    def __init__(self):
        """Initialize the union–find data structure.

        Creates two internal mappings:
        - parent: maps each element to its parent (or itself if it is a root)
        - rank: stores a heuristic depth used for efficient merging.
        """

        self.parent = {}
        self.rank = {}

    def find(self, x):
        # initialize if unseen
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        # path compression
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        """Merge the two disjoint sets containing a and b.

        Args:
            a (Hashable): First element.
            b (Hashable): Second element.

        Notes:
            Applies union by rank to maintain shallow trees.
            If ranks are equal, one root becomes parent and its rank increases by one.
        """
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        rank_a = self.rank.get(ra, 0)
        rank_b = self.rank.get(rb, 0)
        if rank_a < rank_b:
            self.parent[ra] = rb
        elif rank_a > rank_b:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] = rank_a + 1

    def groups(self):
        """Return all disjoint sets currently tracked.

        Returns:
            list[list[Hashable]]: A list of connected components, where each
            inner list contains elements that share a common root.

        Notes:
            The order of groups and members is arbitrary. Path compression
            is applied before grouping to ensure consistent roots.
        """
        out = {}
        for node in list(self.parent.keys()):
            root = self.find(node)
            out.setdefault(root, []).append(node)
        return list(out.values())

