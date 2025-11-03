import unittest
import sys
from pathlib import Path

sys.path.insert(0, "..")
sys.path.insert(0, "../..")
sys.path.insert(0, "../neo4j_interface")

from omics_io.omics_io.parsers.reference.genbank import parse_genbank
from omics_io.omics_io.parsers.reference.gff3 import parse_gff3
from omics_io.omics_io.parsers.transcriptomics.precise2 import parse_precise2
from omics_io.omics_io.parsers.proteomics.mztab import parse_mztab
from omics_io.omics_io.parsers.metabolomics.table import parse_metabolomics

from idalign.idalign.aligner import match_references

from network_builder.builder import OmicGraphBuilder
from network_builder.builder import username, url, password, _count_to_conf
from network_builder.neo4j_interface.storage import Neo4jStorage
from network_builder.utils.network_builder_convention import NetworkBuilderConvention
from omics_io.omics_io.identifiers import IDS


class TestNetworkBuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass
        """
        cls.sequence_data = parse_genbank(Path("data/sequence.gb"))

        cls.gff_data = parse_gff3(Path("data/sequence.gff3"))
        cls.trans_data = parse_precise2(
            Path("data/expression_tpm_log.csv"),
            presence_path=Path("data/gene_presence_matrix_bool.csv"),
            metadata_path=Path("data/samples_metadata.csv"),
        )
        cls.prot1 = parse_mztab(Path("data/F020490.pride.mztab"))
        cls.prot2 = parse_mztab(Path("data/F020549.pride.mztab"))
        cls.prot3 = parse_mztab(Path("data/F023569.pride.mztab"))
        cls.meta = parse_metabolomics(
            Path("data/measurements.tsv"),
            mapping_path=Path("data/map.tsv"),
        )
        """
        cls.storage = Neo4jStorage(url, username=username, 
                                   password=password,
                                   convention=NetworkBuilderConvention())

    def setUp(self):
        self.storage.drop()

    def assert_identifiers_present(self, dataset):
        # Every element from each dataset must appear either
        # as a concrete node ID or inside a stub’s aliases.
        ids = [ele for ele, _, _ in dataset]
        found = {n.id for n in self.storage.find_nodes(ids=ids)}
        stubs = self.storage.find_nodes(label=IDS.type.unknown)
        alias_in_stubs = set()
        for s in stubs:
            alias_in_stubs.update(s.properties.get(IDS.predicates.alias, []))
        missing = [i for i in ids if i not in found and i not in alias_in_stubs]
        self.assertFalse(missing, f"Missing identifiers: {missing}")

    def assert_groups_materialized(self, result):
        # Each alignment group must have a focal node id,
        # either merged real or a stub UNK:{gid},
        # and all group members are present as node
        # IDs or aliases on that focal.
        for gid, nodes in result:
            r_node_ids = list(set([e.identifier for e in nodes]))
            s_nodes = self.storage.find_nodes(ids=r_node_ids)
            if len(s_nodes) == 0:
                uk_node_ids = [f"UNK:{r}" for r in r_node_ids]
                s_nodes = self.storage.find_nodes(ids=uk_node_ids)
                self.assertTrue(
                    len(s_nodes) == 1, f"Given {r_node_ids}, cant find a node. {gid}"
                )
                self.assertCountEqual(
                    r_node_ids,
                    s_nodes[0].properties.get(IDS.predicates.alias)
                )
            else:
                self.assertEqual(len(s_nodes), 1, f"Given {r_node_ids}, cant find a node. {gid}")
                self.assertCountEqual(
                    r_node_ids,
                    s_nodes[0].properties.get(IDS.predicates.alias) + [s_nodes[0].id],
                    f"{r_node_ids},{s_nodes[0].properties.get(IDS.predicates.alias) + [s_nodes[0].id]}",
                )

    def assert_relations_projected(self, result):
        # Edges exist for every alignment edge, with confidence in $(0,1)$.
        for g1, rel, g2, cnt in result.relationships():
            n1 = self.storage.find_nodes(ids=[n.identifier for n in result.members(g1)])
            n2 = self.storage.find_nodes(ids=[n.identifier for n in result.members(g2)])
            if len(n2) == 0:
                uk_node_ids = [f"UNK:{r}" for r in [n.identifier for n in result.members(g2)]]
                n2 = self.storage.find_nodes(ids=uk_node_ids)
            self.assertEqual(len(n1) , 1,f'{n1}')
            self.assertEqual(len(n2) , 1,f'{n2}')
            res = self.storage.find_relationships(rel, n1[0].id, n2[0].id)
            self.assertEqual(len(res),1)

    def assert_idempotent(self, datasets, result,gb):
        n1 = len(self.storage.find_nodes())
        r1 = len(self.storage.find_relationships())
        for ds in datasets:
            gb.add_omic_set(ds)
        gb.add_alignment_data(result)
        n2 = len(self.storage.find_nodes())
        r2 = len(self.storage.find_relationships())
        self.assertEqual((n1, r1), (n2, r2))

    def test_build_single_datafile(self):
        gb = OmicGraphBuilder()
        sequence_data = parse_genbank(Path("data/sequence.gb"))
        result = match_references(sequence_data)

        gb.add_omic_set(sequence_data)
        gb.add_alignment_data(result)

        self.assert_identifiers_present(sequence_data)
        self.assert_groups_materialized(result)
        self.assert_relations_projected(result)
        self.assert_idempotent([sequence_data], result,gb)

    def test_build_single_omic(self):
        gb = OmicGraphBuilder()
        sequence_data = parse_genbank(Path("data/sequence.gb"))
        gff_data = parse_gff3(Path("data/sequence.gff3"))
        result = match_references(sequence_data,gff_data)

        gb.add_omic_set(sequence_data)
        gb.add_omic_set(gff_data)
        gb.add_alignment_data(result)

        self.assert_identifiers_present(sequence_data)
        self.assert_identifiers_present(gff_data)
        self.assert_groups_materialized(result)
        self.assert_relations_projected(result)
        self.assert_idempotent([sequence_data,gff_data], result,gb)

    def test_build_double_omic(self):
        gb = OmicGraphBuilder()
        sequence_data = parse_genbank(Path("data/sequence.gb"))
        gff_data = parse_gff3(Path("data/sequence.gff3"))
        trans_data = parse_precise2(
            Path("data/expression_tpm_log.csv"),
            presence_path=Path("data/gene_presence_matrix_bool.csv"),
            metadata_path=Path("data/samples_metadata.csv"),
        )
        result = match_references(sequence_data,gff_data,trans_data)

        gb.add_omic_set(sequence_data)
        gb.add_omic_set(gff_data)
        gb.add_omic_set(trans_data)
        gb.add_alignment_data(result)

        self.assert_identifiers_present(sequence_data)
        self.assert_identifiers_present(gff_data)
        self.assert_identifiers_present(trans_data)
        self.assert_groups_materialized(result)
        self.assert_relations_projected(result)
        self.assert_idempotent([sequence_data,gff_data,trans_data], result,gb)

    def test_build_all_omic(self):
        result = match_references(
            self.sequence_data,
            self.gff_data,
            self.trans_data,
            self.prot1,
            self.prot2,
            self.prot3,
            self.meta,
        )
        build_network(self.sequence_data, alignment_data=result)
