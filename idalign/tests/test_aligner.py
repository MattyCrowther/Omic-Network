import unittest
import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from idalign.aligner import match_references

# ---- lightweight dummies matching matcher input shape ----
class XRef:
    __slots__ = ("src", "namespace", "target")
    def __init__(self, src, namespace, target):
        self.src = src
        self.namespace = namespace
        self.target = target

class Feat:
    __slots__ = ("id", "entity", "namespace")
    def __init__(self, id, entity=None, namespace=None):
        self.id = id
        self.entity = entity
        self.namespace = namespace

class DummyOD:
    __slots__ = ("name", "feature_meta", "cross_ref")
    def __init__(self, name, feature_meta=None, cross_ref=None):
        self.name = name
        self.feature_meta = feature_meta or []
        self.cross_ref = cross_ref or []

# ---- helpers against AlignmentResult ----
def find_gid(res, namespace, ident):
    for gid, members in res.groups.items():
        for n in members:
            if n.namespace == namespace and n.identifier == ident:
                return gid
    return None

def group_type(res, gid):
    types = {getattr(n, "type", "unknown") for n in res.groups[gid]}
    return next(iter(types)) if len(types) == 1 else "mixed/unknown"

def scoped_id_dupes(res):
    seen = {}
    dupes = []
    for gid, members in res.groups.items():
        for n in members:
            k = (n.namespace, n.identifier)
            if k in seen and seen[k] != gid:
                dupes.append((k, seen[k], gid))
            else:
                seen[k] = gid
    return dupes

def same_group_relations(res):
    return [ (ga, rel, gb, c) for (ga, rel, gb), c in res.group_relations.items() if ga == gb ]


class TestAlignment(unittest.TestCase):

    def test_registry_alias_merges_across_datasets(self):
        # A: local gene under dataset scope "rna" (matcher lowercases dataset names)
        A = DummyOD(
            name="rna",
            feature_meta=[Feat("sad", entity="gene", namespace="gene")],
            cross_ref=[
                XRef("sad", "geneid", "947440"),    # alias → scopes target to "geneid"
                XRef("sad", "locus_tag", "b1525"),  # alias → stays in source dataset scope ("rna")
            ],
        )
        # B: same GeneID appears as a feature in another dataset
        B = DummyOD(
            name="protein",
            feature_meta=[Feat("947440", entity="gene", namespace="geneid")],
            cross_ref=[],
        )

        res = match_references(A, B)

        gid_geneid = find_gid(res, "geneid", "947440")
        gid_sad    = find_gid(res, "rna", "sad")      # dataset-local scope (lowercased)
        gid_btag   = find_gid(res, "rna", "b1525")    # dataset-local scope

        self.assertIsNotNone(gid_geneid)
        self.assertEqual(gid_geneid, gid_sad)
        self.assertEqual(gid_geneid, gid_btag)
        self.assertEqual(group_type(res, gid_geneid), "gene")

    def test_synonyms_stay_local_without_registry_bridge(self):
        A = DummyOD(
            name="set1",
            feature_meta=[Feat("sad", entity="gene", namespace="gene")],
            cross_ref=[XRef("sad", "locus_tag", "b1525")],
        )
        B = DummyOD(
            name="set2",
            feature_meta=[Feat("sad", entity="gene", namespace="gene")],
            cross_ref=[XRef("sad", "locus_tag", "b1525")],
        )

        res = match_references(A, B)

        gid1 = find_gid(res, "set1", "sad")  # dataset-local
        gid2 = find_gid(res, "set2", "sad")  # dataset-local
        self.assertIsNotNone(gid1)
        self.assertIsNotNone(gid2)
        self.assertNotEqual(gid1, gid2)

    def test_unknown_types_do_not_block_alias_merge(self):
        A = DummyOD(
            name="setA",
            feature_meta=[Feat("X", entity="gene", namespace="gene")],
            cross_ref=[XRef("X", "geneid", "111")],
        )
        B = DummyOD(
            name="setB",
            feature_meta=[Feat("111", entity=None, namespace="geneid")],  # unknown type
            cross_ref=[],
        )

        res = match_references(A, B)
        gid = find_gid(res, "geneid", "111")
        self.assertIsNotNone(gid)
        members = {(n.namespace, n.identifier, n.type) for n in res.groups[gid]}
        # dataset scopes are lowercased by the matcher
        self.assertIn(("seta", "X", "gene"), members)
        self.assertIn(("geneid", "111", "gene"), members)

    def test_relations_aggregate_by_group(self):
        # Gene g1 produces protein P99999 (via produced_by); also alias g1 to a registry GeneID
        A = DummyOD(
            name="set",
            feature_meta=[
                Feat("g1", entity="gene", namespace="gene"),   # local → scope "set"
            ],
            cross_ref=[
                XRef("g1", "geneid", "100"),                  # alias to registry
                XRef("g1", "produced_by", "P99999"),          # relation gene -> protein
            ],
        )

        res = match_references(A)

        gid_g1 = find_gid(res, "set", "g1")            # dataset-local
        gid_p  = find_gid(res, "protein", "P99999")    # relation target normalised to "protein" scope

        self.assertIsNotNone(gid_g1)
        self.assertIsNotNone(gid_p)

        total_rel_edges = sum(
            c for (ga, rel, gb), c in res.group_relations.items()
            if ga == gid_g1 and gb == gid_p
        )
        self.assertGreaterEqual(total_rel_edges, 1)

    def test_no_self_loops_after_grouping(self):
        A = DummyOD(
            name="set",
            feature_meta=[Feat("g1", entity="gene", namespace="gene")],
            cross_ref=[XRef("g1", "geneid", "200")],
        )
        res = match_references(A)
        self.assertEqual(same_group_relations(res), [])

    def test_scoped_id_uniqueness(self):
        A = DummyOD(
            name="set",
            feature_meta=[Feat("g1", entity="gene", namespace="gene")],
            cross_ref=[XRef("g1", "geneid", "300")],
        )
        res = match_references(A)
        self.assertEqual(scoped_id_dupes(res), [])

    def test_alias_target_scoped_by_registry(self):
        A = DummyOD(
            name="A",
            feature_meta=[Feat("sad", entity="gene", namespace="gene")],
            cross_ref=[XRef("sad", "GeneID", "947440")],   # case-insensitive namespace
        )
        B = DummyOD(
            name="B",
            feature_meta=[Feat("947440", entity="gene", namespace="GeneID")],
            cross_ref=[],
        )

        res = match_references(A, B)
        gid_geneid = find_gid(res, "geneid", "947440")
        gid_sad    = find_gid(res, "a", "sad")  # dataset-local AND lowercased
        self.assertIsNotNone(gid_geneid)
        self.assertIsNotNone(gid_sad)
        self.assertEqual(gid_geneid, gid_sad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
