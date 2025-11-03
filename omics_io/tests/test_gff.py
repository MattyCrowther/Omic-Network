import sys
import unittest
import tempfile
from pathlib import Path
import warnings
import pandas as pd

sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from omics_io.parsers.reference.gff3 import parse_gff3
from omics_io.parse_obj import FeatureRec, CrossRef


GFF_SINGLE = """
NC_000913.3\tRefSeq\tgene\t6\t20\t.\t+\t.\tID=gene-b0001;locus_tag=b0001;gene=gA;product=pA
NC_000913.3\tRefSeq\tCDS\t6\t20\t.\t+\t0\tID=cds-b0001;Parent=gene-b0001;locus_tag=b0001;gene=gA;product=pA;Dbxref=GeneID:946000,UniProtKB:P0AAA0
"""

GFF_MULTI = """
CHR1.1\tRefSeq\tCDS\t1\t10\t.\t+\t0\tID=b0002;locus_tag=b0002;gene=gB;product=PB
CHR1.1\tRefSeq\trepeat_region\t30\t40\t.\t+\t.\tID=rr1
CHR2.1\tRefSeq\ttRNA\t5\t25\t.\t+\t.\tID=t001;locus_tag=t001;product=tRNA-Lys
CHR2.1\tRefSeq\tCDS\t11\t20\t.\t+\t0\tID=b0003;locus_tag=b0003;gene=gC;product=PC
"""

GFF_EMPTY = """
NC_000913.3\tRefSeq\trepeat_region\t10\t20\t.\t+\t.\tID=r001
"""


class TestGFF3(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)

    def _write_tmp(self, text: str) -> str:
        td = tempfile.TemporaryDirectory()
        p = Path(td.name) / "tmp.gff"
        p.write_text(text)
        self.addCleanup(td.cleanup)
        return str(p)

    def _assert_invariants(self, od):
        self.assertEqual(od.omics_type, "genomics")
        self.assertIsInstance(od.matrix, pd.DataFrame)
        self.assertTrue(od.matrix.index.is_unique)
        self.assertTrue(all(pd.api.types.is_numeric_dtype(dt) for dt in od.matrix.dtypes))

        self.assertIsInstance(od.feature_meta, list)
        self.assertTrue(all(isinstance(r, FeatureRec) for r in od.feature_meta))
        self.assertListEqual([r.id for r in od.feature_meta], list(od.matrix.index))

        if od.cross_ref is not None:
            self.assertTrue(all(isinstance(r, CrossRef) for r in od.cross_ref))
            for r in od.cross_ref:
                self.assertIn(r.axis, {"feature", "column"})
                if r.axis == "feature":
                    self.assertIn(r.src, od.matrix.index)

    def test_single_record_merges_gene_and_cds(self):
        path = self._write_tmp(GFF_SINGLE)
        od = parse_gff3(path)
        self._assert_invariants(od)

        self.assertEqual(len(od.matrix), 1)
        fid = od.feature_meta[0]
        self.assertEqual(fid.id, "b0001")
        attrs = fid.attrs
        self.assertEqual(attrs["product"], "pA")
        self.assertEqual(attrs["contig"], "NC_000913")
        self.assertIn("start", attrs)
        self.assertIn("end", attrs)

        cr = od.cross_ref
        have = {(r.axis, r.src, r.namespace, r.target) for r in cr}
        need = {
            ("feature", "b0001", "locus_tag", "b0001"),
            ("feature", "b0001", "gene", "gA"),
            ("feature", "b0001", "GeneID", "946000"),
            ("feature", "b0001", "UniProtKB", "P0AAA0"),
        }
        for w in need:
            self.assertIn(w, have)

    def test_multiple_records_includes_all_ftypes(self):
        path = self._write_tmp(GFF_MULTI)
        od = parse_gff3(path)
        self._assert_invariants(od)

        fids = [r.id for r in od.feature_meta]
        self.assertSetEqual(set(fids), {"b0002", "t001", "b0003"})
        contigs = {r.attrs["contig"] for r in od.feature_meta}
        self.assertSetEqual(contigs, {"CHR1", "CHR2"})
        types = {r.entity for r in od.feature_meta}
        self.assertTrue(types.issubset({"cds", "trna", "dna"}))

    def test_empty_file_yields_empty(self):
        path = self._write_tmp(GFF_EMPTY)
        od = parse_gff3(path)
        self._assert_invariants(od)
        self.assertEqual(len(od.matrix), 0)
        self.assertEqual(len(od.feature_meta), 0)


    def test_empty_file_yields_empty(self):
        path = self._write_tmp(GFF_EMPTY)
        od = parse_gff3(path)
        self._assert_invariants(od)

        self.assertEqual(len(od.matrix), 0)

if __name__ == "__main__":
    unittest.main()
