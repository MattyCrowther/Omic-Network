import sys
import unittest
import tempfile
from pathlib import Path
import pandas as pd

sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from omics_io.parsers.transcriptomics.precise2 import parse_precise2
from omics_io.parse_obj import FeatureRec, ColumnRec, CrossRef


class TestPrecise2(unittest.TestCase):
    def _write_csv(self, df: pd.DataFrame, name: str, tmpdir: Path) -> Path:
        p = tmpdir / name
        df.to_csv(p)
        return p

    def test_expression_only(self):
        with tempfile.TemporaryDirectory() as td:
            tmpdir = Path(td)

            expr = pd.DataFrame(
                [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]],
                index=["b0001", "b0002"],
                columns=["s1", "s2", "s3"],
            )
            exp_path = self._write_csv(expr, "expression.csv", tmpdir)

            od = parse_precise2(str(exp_path))

            self.assertEqual(od.omics_type, "transcriptomics")
            pd.testing.assert_frame_equal(od.matrix, expr)

            self.assertTrue(all(isinstance(f, FeatureRec) for f in od.feature_meta))
            self.assertEqual([f.id for f in od.feature_meta], ["b0001", "b0002"])
            self.assertTrue(all(f.entity == "rna" for f in od.feature_meta))

            self.assertTrue(all(isinstance(c, ColumnRec) for c in od.column_meta))
            self.assertEqual([c.id for c in od.column_meta], ["s1", "s2", "s3"])

            self.assertTrue(all(isinstance(x, CrossRef) for x in od.cross_ref))
            srcs = [x.src for x in od.cross_ref]
            self.assertSetEqual(set(srcs), {"b0001", "b0002"})
            ns = [x.namespace for x in od.cross_ref]
            self.assertIn("locus_tag", ns)

    def test_metadata_integration(self):
        with tempfile.TemporaryDirectory() as td:
            tmpdir = Path(td)

            expr = pd.DataFrame(
                [[1, 2], [3, 4]],
                index=["g1", "g2"],
                columns=["S1", "S2"],
            )
            meta = pd.DataFrame(
                {
                    "GEO": ["GSE1", "GSE2"],
                    "SRR": ["SRR1", "SRR2"],
                    "SRX": ["SRX1", None],
                    "project": ["P1", "P1"],
                    "contact": ["x", "y"],
                },
                index=["S1", "S2"],
            )

            exp_path = self._write_csv(expr, "expr.csv", tmpdir)
            md_path = self._write_csv(meta, "meta.csv", tmpdir)

            od = parse_precise2(str(exp_path), metadata_path=str(md_path))

            kept_cols = [c.attrs for c in od.column_meta]
            self.assertTrue(all("project" in k for k in kept_cols))
            self.assertTrue(all("contact" not in k for k in kept_cols))

            cr_samples = [x for x in od.cross_ref if x.axis == "column"]
            self.assertTrue(any(x.namespace == "GEO" for x in cr_samples))
            self.assertTrue(any(x.namespace == "SRR" for x in cr_samples))

    def test_duplicate_ids_raises(self):
        with tempfile.TemporaryDirectory() as td:
            tmpdir = Path(td)
            expr = pd.DataFrame(
                [[1, 2], [3, 4]],
                index=["g1", "g1"],
                columns=["s1", "s2"],
            )
            exp_path = self._write_csv(expr, "exp.csv", tmpdir)
            with self.assertRaises(ValueError):
                parse_precise2(str(exp_path))

    def test_returns_expected_types(self):
        with tempfile.TemporaryDirectory() as td:
            tmpdir = Path(td)
            expr = pd.DataFrame(
                [[1, 2, 3]],
                index=["b0001"],
                columns=["a", "b", "c"],
            )
            exp_path = self._write_csv(expr, "exp.csv", tmpdir)
            od = parse_precise2(str(exp_path))
            self.assertIsInstance(od.matrix, pd.DataFrame)
            self.assertIsInstance(od.feature_meta, list)
            self.assertIsInstance(od.column_meta, list)
            self.assertIsInstance(od.cross_ref, list)


if __name__ == "__main__":
    unittest.main()
