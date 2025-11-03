import gzip
import sys
import os
import unittest
import tempfile
from pathlib import Path
import pandas as pd

sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from omics_io.parsers.metabolomics.table import parse_metabolomics
from omics_io.errors import ParseError


def write_tsv(df: pd.DataFrame, path: Path, sep="\t", index=False):
    df.to_csv(path, sep=sep, index=index)


def write_tsv_gz(df: pd.DataFrame, path: Path, sep="\t", index=False):
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        df.to_csv(fh, sep=sep, index=index)


from dataclasses import asdict


def cr_df(od):
    return pd.DataFrame([asdict(r) for r in (od.cross_ref or [])])


class TestIntensityTable(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="metab_ut_"))

        self.samples = ["ArcA_1", "ArcA_2", "WT_1", "WT_2"]
        self.metabs = ["Alanine", "Arginine"]
        rows = [
            [
                "Factors",
                "-",
                "Genotype:Mutant",
                "Genotype:Mutant",
                "Genotype:Wildtype",
                "Genotype:Wildtype",
            ],
        ]
        for m in self.metabs:
            rows.append([m, m, 1.0, 2.0, 3.0, 4.0])
        self.intensity_df = pd.DataFrame(
            rows, columns=["Metabolite_name", "RefMet_name"] + self.samples, dtype=str
        )

        self.map_df = pd.DataFrame(
            {
                "Metabolite_Name": ["Alanine", "Arginine"],
                "RefMet_Name": ["Alanine", "Arginine"],
                "WorkBench_ID": ["ME0001", "ME0002"],
                "KEGG_ID": ["C00041", "C00062"],
            }
        )

        self.map_df_ci = self.map_df.rename(
            columns={
                "Metabolite_Name": "metabolite",
                "RefMet_Name": "refmet",
                "WorkBench_ID": "WORKBENCH",
                "KEGG_ID": "kEgG",
            }
        )

    def tearDown(self):
        for p in self.tmp.glob("*"):
            try:
                p.unlink()
            except Exception:
                pass
        try:
            os.rmdir(self.tmp)
        except Exception:
            pass

    def test_basic_parse_with_mapping(self):
        p_int = self.tmp / "intensity.tsv"
        p_map = self.tmp / "map.tsv"
        write_tsv(self.intensity_df, p_int, index=False)
        write_tsv(self.map_df, p_map, index=False)

        od = parse_metabolomics(str(p_int), mapping_path=str(p_map))
        self.assertEqual(od.omics_type, "metabolomics")

        self.assertEqual(od.feature_meta[0].namespace, "RefMet")
        self.assertEqual(od.matrix.shape, (2, 4))
        self.assertTrue(pd.api.types.is_numeric_dtype(od.matrix.dtypes.iloc[0]))

        cr = cr_df(od)
        counts = cr.groupby(["axis", "namespace"]).size().to_dict()

        self.assertEqual(counts.get(("feature", "RefMet"), 0), 2)
        self.assertEqual(counts.get(("feature", "WorkBench"), 0), 2)
        self.assertEqual(counts.get(("feature", "KEGG.Compound"), 0), 2)
        self.assertEqual(counts.get(("column", "Genotype"), 0), 4)

        feat_ok = set(cr.query("axis == 'feature'")["src"]) <= set(od.matrix.index)
        col_ok = set(cr.query("axis == 'column'")["src"]) <= set(od.matrix.columns)
        self.assertTrue(feat_ok and col_ok)

    def test_case_insensitive_mapping_headers(self):
        p_int = self.tmp / "intensity.tsv"
        p_map = self.tmp / "map_ci.tsv"
        write_tsv(self.intensity_df, p_int, index=False)
        write_tsv(self.map_df_ci, p_map, index=False)

        od = parse_metabolomics(str(p_int), mapping_path=str(p_map))
        cr = cr_df(od)
        alanine_rows = cr[cr["src"] == "Alanine"]
        self.assertTrue(
            (
                (alanine_rows["namespace"] == "WorkBench")
                & (alanine_rows["target"] == "ME0001")
            ).any()
        )
        self.assertTrue(
            (
                (alanine_rows["namespace"] == "KEGG.Compound")
                & (alanine_rows["target"] == "C00041")
            ).any()
        )

    def test_gzip_input(self):
        p_int = self.tmp / "intensity.tsv.gz"
        p_map = self.tmp / "map.tsv.gz"
        write_tsv_gz(self.intensity_df, p_int, index=False)
        write_tsv_gz(self.map_df, p_map, index=False)

        od = parse_metabolomics(str(p_int), mapping_path=str(p_map))
        cr = cr_df(od)
        self.assertEqual(od.matrix.shape, (2, 4))
        self.assertIn("KEGG.Compound", set(cr["namespace"]))

    def test_no_factors_row(self):
        df = self.intensity_df[self.intensity_df["Metabolite_name"] != "Factors"].copy()
        p_int = self.tmp / "nofactors.tsv"
        p_map = self.tmp / "map.tsv"
        write_tsv(df, p_int, index=False)
        write_tsv(self.map_df, p_map, index=False)

        od = parse_metabolomics(str(p_int), mapping_path=str(p_map))
        cr = cr_df(od)

        self.assertFalse((cr["axis"] == "column").any())

        self.assertTrue((cr["namespace"] == "KEGG.Compound").any())

    def test_refmet_absent_uses_metabolite_namespace(self):
        df = self.intensity_df.drop(columns=["RefMet_name"])
        p_int = self.tmp / "no_refmet.tsv"
        write_tsv(df, p_int, index=False)

        od = parse_metabolomics(str(p_int))
        self.assertEqual(od.feature_meta[0].namespace, "Metabolite")
        cr = cr_df(od)
        self.assertTrue(
            ((cr["namespace"] == "Metabolite") & (cr["src"] == "Alanine")).any()
        )


if __name__ == "__main__":
    unittest.main()
