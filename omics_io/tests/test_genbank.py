import sys
import unittest
import tempfile
from pathlib import Path
import warnings
import pandas as pd

sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from omics_io.parsers.reference.genbank import parse_genbank
from omics_io.errors import ParseError

SEQ60 = "atgcatgcat" * 6
SEQ80 = "atgcatgcat" * 8

GB_SINGLE = f"""\
LOCUS       CHR1                60 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  test record 1.
ACCESSION   CHR1
VERSION     CHR1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
FEATURES             Location/Qualifiers
     source          1..60
     gene            6..20
                     /gene="gA"
                     /locus_tag="b0001"
     CDS             6..20
                     /gene="gA"
                     /locus_tag="b0001"
                     /product="pA"
                     /protein_id="protA"
                     /db_xref="GeneID:946000"
                     /db_xref="UniProtKB:P0AAA0"
ORIGIN
        1 {SEQ60}
//
"""

GB_MULTI = f"""\
LOCUS       CHR1                80 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  test record chr1.
ACCESSION   CHR1
VERSION     CHR1.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
FEATURES             Location/Qualifiers
     source          1..80
     CDS             1..10
                     /locus_tag="b0002"
                     /gene="gB"
                     /product="PB"
     gene            1..10
                     /locus_tag="b0002"
                     /gene="gB"
                     /product="PB"
     repeat_region   30..40
ORIGIN
        1 {SEQ80}
//
LOCUS       CHR2                60 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  test record chr2.
ACCESSION   CHR2
VERSION     CHR2.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
FEATURES             Location/Qualifiers
     source          1..60
     tRNA            5..25
                     /locus_tag="t001"
                     /product="tRNA-Lys"
     CDS             11..20
                     /locus_tag="b0003"
                     /gene="gC"
                     /product="PC"
ORIGIN
        1 {SEQ60}
//
"""

GB_NO_GENE = f"""\
LOCUS       CHR0                60 bp    DNA     linear   BCT 01-JAN-2000
DEFINITION  no kept features.
ACCESSION   CHR0
VERSION     CHR0.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
FEATURES             Location/Qualifiers
     source          1..60
     repeat_region   5..10
ORIGIN
        1 {SEQ60}
//
"""

CR_FIELDS = ["axis", "src", "namespace", "target"]


class TestGenbank(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.filterwarnings("ignore", category=UserWarning)
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        try:
            from Bio import BiopythonParserWarning

            warnings.filterwarnings("ignore", category=BiopythonParserWarning)
        except Exception:
            pass

    def _write_tmp(self, text: str) -> str:
        td = tempfile.TemporaryDirectory()
        p = Path(td.name) / "tmp.gb"
        p.write_text(text)
        self.addCleanup(td.cleanup)
        return str(p)

    def _feature_df(self, od) -> pd.DataFrame:
        """Flatten list[FeatureRec] -> DataFrame indexed by id with attrs expanded."""
        recs = od.feature_meta
        if not recs:
            return pd.DataFrame(columns=["entity", "namespace"]).set_index(
                pd.Index([], name="id")
            )
        rows = []
        for r in recs:
            row = {"id": r.id, "entity": r.entity, "namespace": r.namespace}
            row.update(r.attrs or {})
            rows.append(row)
        return pd.DataFrame(rows).set_index("id").sort_index()

    def _assert_invariants(self, od):
        self.assertEqual(od.omics_type, "genomics")
        self.assertIsInstance(od.matrix, pd.DataFrame)

        feat_ids = [r.id for r in od.feature_meta]
        self.assertListEqual(list(od.matrix.index), feat_ids)
        self.assertEqual(od.matrix.shape[1], 0)

        if od.cross_ref is not None:
            cr = od.cross_ref
            self.assertIsInstance(cr, list)
            self.assertTrue(
                all(
                    hasattr(r, "axis")
                    and hasattr(r, "src")
                    and hasattr(r, "namespace")
                    and hasattr(r, "target")
                    for r in cr
                )
            )
            self.assertTrue({r.axis for r in cr}.issubset({"feature"}))
            self.assertTrue({r.src for r in cr}.issubset(set(feat_ids)))

    def test_single_record_merges_gene_and_cds(self):
        path = self._write_tmp(GB_SINGLE)
        od = parse_genbank(path)
        self._assert_invariants(od)

        self.assertEqual(len(od.feature_meta), 1)
        fdf = self._feature_df(od)
        self.assertIn("b0001", fdf.index)

        self.assertEqual(fdf.loc["b0001", "product"], "pA")

        self.assertEqual(fdf.loc["b0001", "contig"], "CHR1")

        cr = od.cross_ref
        want = {
            ("feature", "b0001", "gene", "gA"),
            ("feature", "b0001", "protein_id", "protA"),
            ("feature", "b0001", "GeneID", "946000"),
            ("feature", "b0001", "UniProtKB", "P0AAA0"),
        }
        got = {(r.axis, r.src, r.namespace, r.target) for r in cr}
        for w in want:
            self.assertIn(w, got)

    def test_multiple_records_includes_all_feature_types(self):
        path = self._write_tmp(GB_MULTI)
        od = parse_genbank(path)
        self._assert_invariants(od)

        fdf = self._feature_df(od)
        self.assertSetEqual(set(fdf.index), {"b0002", "b0003", "t001"})
        self.assertSetEqual(set(fdf["contig"]), {"CHR1", "CHR2"})

    def test_non_gene_records_yield_empty(self):
        path = self._write_tmp(GB_NO_GENE)
        try:
            od = parse_genbank(path)
        except ParseError:
            self.skipTest("Parser configured to raise on empty table.")
            return
        self._assert_invariants(od)
        self.assertEqual(len(od.feature_meta), 0)
        self.assertEqual(len(od.matrix), 0)

        if od.cross_ref is not None:
            self.assertEqual(len(od.cross_ref), 0)


if __name__ == "__main__":
    unittest.main()
