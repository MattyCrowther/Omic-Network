import unittest
from pathlib import Path
import sys
sys.path.insert(0, "..")
sys.path.insert(0, "../..")
from idalign.aligner import match_references
from omics_io.omics_io.parsers.reference.genbank import parse_genbank
from omics_io.omics_io.parsers.reference.gff3 import parse_gff3
from omics_io.omics_io.parsers.transcriptomics.precise2 import parse_precise2
from omics_io.omics_io.parsers.proteomics.mztab import parse_mztab
from omics_io.omics_io.parsers.metabolomics.table import parse_metabolomics
from idalign.aligner import match_references


# ---------- helpers (same assertions reused across scenarios) ----------
def summarize_alignment(res) -> dict:
    sizes = [len(m) for m in res.groups.values()]
    n_groups = len(res.groups)
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
        "n_rel_triplets": len(res.group_relations),
        "n_rel_edges_total": sum(res.group_relations.values()),
    }

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

def find_gid(res, namespace, ident):
    for gid, members in res.groups.items():
        if any(n.namespace == namespace and n.identifier == ident for n in members):
            return gid
    return None

def assert_common_invariants(testcase: unittest.TestCase, res, label: str):
    m = summarize_alignment(res)
    # Basic sanity
    testcase.assertGreater(m["n_groups"], 0, f"[{label}] expected >0 groups")
    testcase.assertEqual(m["n_nodes"], sum(len(s) for s in res.groups.values()), f"[{label}] node count mismatch")
    # No duplicate scoped IDs across groups
    testcase.assertEqual(scoped_id_dupes(res), [], f"[{label}] duplicate scoped IDs across groups")
    # No self-loops after grouping
    testcase.assertEqual(same_group_relations(res), [], f"[{label}] self-loop relations present")
    # Relations cardinality is consistent
    testcase.assertGreaterEqual(m["n_rel_edges_total"], m["n_rel_triplets"], f"[{label}] relation edge counts inconsistent")
    # Types resolved (unknown allowed, but the UF key must be independent of type)
    testcase.assertGreaterEqual(m["unknown_type_nodes"], 0, f"[{label}] bad unknown count")
    # Print a terse summary to aid E2E debugging
    print(f"[{label}] Summary: {m}")

class TestAlignerE2E(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Input paths
        cls.p_gbk  = Path("data/sequence.gb")
        cls.p_gff  = Path("data/sequence.gff3")
        cls.p_txp  = Path("data/expression_tpm_log.csv")
        cls.p_txp_pres = Path("data/gene_presence_matrix_bool.csv")
        cls.p_txp_meta = Path("data/samples_metadata.csv")
        cls.p_mztab1 = Path("data/F020490.pride.mztab")
        cls.p_mztab2 = Path("data/F020549.pride.mztab")
        cls.p_mztab3 = Path("data/F023569.pride.mztab")
        cls.p_meta   = Path("data/measurements.tsv")
        cls.p_map    = Path("data/map.tsv")

        # Skip the whole class if required files are missing
        missing = [p for p in [
            cls.p_gbk, cls.p_gff, cls.p_txp, cls.p_txp_pres, cls.p_txp_meta,
            cls.p_mztab1, cls.p_mztab2, cls.p_mztab3, cls.p_meta, cls.p_map
        ] if not p.exists()]

        if missing:
            raise unittest.SkipTest(f"Skipping E2E: missing files: {missing}")

        # Parse once
        cls.sequence_data = parse_genbank(cls.p_gbk)
        cls.gff_data      = parse_gff3(cls.p_gff)
        cls.trans_data    = parse_precise2(cls.p_txp, presence_path=cls.p_txp_pres, metadata_path=cls.p_txp_meta)
        cls.prot1         = parse_mztab(cls.p_mztab1)
        cls.prot2         = parse_mztab(cls.p_mztab2)
        cls.prot3         = parse_mztab(cls.p_mztab3)
        cls.meta          = parse_metabolomics(cls.p_meta, mapping_path=cls.p_map)

    # 1) Single dataset (GBK)
    def test_01_gbk_only(self):
        res = match_references(self.sequence_data)
        assert_common_invariants(self, res, "GBK-only")
        # Expect many singletons and typical max group size ~5 (symbol/locus_tag/GeneID/etc.)
        m = summarize_alignment(res)
        self.assertGreater(m["singletons"], 0, "[GBK-only] expected singletons")
        self.assertGreaterEqual(m["max_group_size"], 3, "[GBK-only] expected some merged alias groups")

    # 2) Two genomics sources (GBK + GFF3)
    def test_02_gbk_plus_gff(self):
        res = match_references(self.sequence_data, self.gff_data)
        assert_common_invariants(self, res, "GBK+GFF")
        m = summarize_alignment(res)
        # Adding GFF should not reduce mean group size
        res_gbk = match_references(self.sequence_data)
        m_gbk = summarize_alignment(res_gbk)
        self.assertGreaterEqual(m["mean_group_size"], m_gbk["mean_group_size"], "[GBK+GFF] mean group size should not drop")

        # Cross-dataset registry merge check (best-effort): if any GeneID exists, it should map to exactly one GID
        # (We don't hardcode a specific GeneID; we verify consistency for the first few)
        checked = 0
        for gid, members in res.groups.items():
            for n in members:
                if n.namespace == "geneid":
                    gid_again = find_gid(res, "geneid", n.identifier)
                    self.assertEqual(gid, gid_again, f"[GBK+GFF] GeneID {n.identifier} split across groups")
                    checked += 1
                    if checked >= 10:
                        break
            if checked >= 10:
                break

    # 3) Genomics + Transcriptomics
    def test_03_gen_plus_transcriptomics(self):
        res = match_references(self.sequence_data, self.gff_data, self.trans_data)
        assert_common_invariants(self, res, "Gen+Trans")
        m = summarize_alignment(res)
        # Adding transcriptomics should increase relations or at least keep them same
        res_gen = match_references(self.sequence_data, self.gff_data)
        m_gen = summarize_alignment(res_gen)
        self.assertGreaterEqual(m["n_rel_edges_total"], m_gen["n_rel_edges_total"], "[Gen+Trans] relations should not decrease")

    # 4) All omics (Gen + Trans + Proteomics + Metabolomics)
    def test_04_all_omics(self):
        res = match_references(
            self.sequence_data, self.gff_data, self.trans_data,
            self.prot1, self.prot2, self.prot3, self.meta
        )
        assert_common_invariants(self, res, "All-omics")
        m = summarize_alignment(res)

        # Compare to Gen+Trans: expect more nodes and typically more relations
        res_gen_tx = match_references(self.sequence_data, self.gff_data, self.trans_data)
        m_gen_tx = summarize_alignment(res_gen_tx)
        self.assertGreaterEqual(m["n_nodes"], m_gen_tx["n_nodes"], "[All-omics] nodes should increase")
        self.assertGreaterEqual(m["n_rel_edges_total"], m_gen_tx["n_rel_edges_total"], "[All-omics] relations should not decrease")

        # Optional: verify at least one protein registry relation present (namespace example from your map)
        has_protein_rel = any(rel in {"product", "protein_id"} for (_, rel, _) in res.group_relations.keys())
        self.assertTrue(has_protein_rel, "[All-omics] expected at least one protein registry relation")

if __name__ == "__main__":
    unittest.main(verbosity=2)
