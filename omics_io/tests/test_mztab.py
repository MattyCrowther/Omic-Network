import io
import sys
import gzip
import unittest
import tempfile
from pathlib import Path
import pandas as pd

sys.path.insert(0, "..")
sys.path.insert(0, "../..")

from omics_io.parsers.proteomics.mztab import parse_mztab
from omics_io.errors import ParseError
from omics_io.parse_obj import CrossRef

_MIN_ID_MZTAB = """MTD\tmzTab-version\t1.0
MTD\tmzTab-mode\tComplete
MTD\tmzTab-type\tIdentification
PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tbest_search_engine_score[1]\tnum_psms\tnum_peptides_unique\tambiguity_members\tmodifications\tprotein_coverage
PRT\tsp|P0A8V2|YBAX_ECOLI\tProtein YbaX OS=Escherichia coli (strain K12) GN=ybaX\t562\tEscherichia coli\tuniprot-ecoli.fasta\tver\t[MS, MS:1001207, Mascot, ]\t200.0\t5\t2\tsp|P0A8V2|YBAX_ECOLI,sp|P0A8W0|XXX\t152-UNIMOD:35\t12.3
"""

_MIN_QUANT_MZTAB = """MTD\tmzTab-version\t1.0
MTD\tmzTab-mode\tComplete
MTD\tmzTab-type\tQuantification
MTD\tassay[1]-name\tCondA_rep1
MTD\tassay[2]-name\tCondA_rep2
PRH\taccession\tdescription\tdatabase\tabundance_assay[1]\tabundance_assay[2]\tmodifications
PRT\tsp|P0A7V8|RPLA_ECOLI\t50S ribosomal protein L1 OS=Escherichia coli GN=rplA\tuniprot\t10.0\t15.0\t83-UNIMOD:35
PRT\tDECOY_REV__XYZ\tDecoy protein\tuniprot\tnull\tnull\tnull
"""

_QUANT_NO_ASSAY = """MTD\tmzTab-version\t1.0
MTD\tmzTab-mode\tComplete
MTD\tmzTab-type\tQuantification
PRH\taccession\tdescription\tdatabase
PRT\tsp|Q9XYZ1|TEST_ECOLI\tTest protein GN=tst\tuniprot
"""


def _write(tmpdir: Path, name: str, text: str) -> str:
    p = tmpdir / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def _write_gz(tmpdir: Path, name: str, text: str) -> str:
    p = tmpdir / name
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(text)
    return str(p)


class TestMzTabParser(unittest.TestCase):
    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp(prefix="mztab_ut_"))

    def tearDown(self):
        try:
            for f in self.tmpdir.glob("*"):
                f.unlink(missing_ok=True)
            self.tmpdir.rmdir()
        except Exception:
            pass

    def test_identification_basic_matrix_and_crossrefs(self):
        path = _write(self.tmpdir, "id.mztab", _MIN_ID_MZTAB)
        od = parse_mztab(path, prefer_quant=True)

        self.assertEqual(od.omics_type, "proteomics")

        acc = "sp|P0A8V2|YBAX_ECOLI"
        cr = od.cross_ref

        has_gene = any(r.namespace == "produced_by" and r.target == "ybaX" for r in cr)
        self.assertTrue(has_gene)

    def test_quant_abundance_matrix_and_decoy_filter(self):
        path = _write(self.tmpdir, "quant.mztab", _MIN_QUANT_MZTAB)
        od = parse_mztab(path, prefer_quant=True, filter_decoys=True)

        self.assertListEqual(sorted(od.matrix.columns.tolist()), ["assay_1", "assay_2"])
        self.assertEqual(od.matrix.shape[0], 1)

    def test_gzip_input(self):
        path = _write_gz(self.tmpdir, "quant.mztab.gz", _MIN_QUANT_MZTAB)
        od = parse_mztab(path, prefer_quant=True)
        self.assertIn("assay_2", od.matrix.columns)

    def test_numeric_matrix(self):
        path = _write(self.tmpdir, "quant_nulls.mztab", _MIN_QUANT_MZTAB)
        od = parse_mztab(path)

        self.assertTrue(
            all(pd.api.types.is_numeric_dtype(dt) for dt in od.matrix.dtypes)
        )

    def test_evidence_thresholds_filter(self):
        path = _write(self.tmpdir, "id_thresh.mztab", _MIN_ID_MZTAB)
        od = parse_mztab(path, min_psms=10, min_unique=10)
        self.assertEqual(od.matrix.shape[0], 0)

        self.assertTrue(od.cross_ref == [] or od.cross_ref is None)

    def test_sample(self):
        prot1 = parse_mztab(Path("data/F020490.pride.mztab"))
        for e in prot1.column_meta:
            print(e)


if __name__ == "__main__":
    unittest.main()
