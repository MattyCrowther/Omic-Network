import unittest
import sys

sys.path.insert(0, "..")
sys.path.insert(0, "../..")
sys.path.insert(0, "../neo4j_interface")

from network_builder.neo4j_interface.storage import Neo4jStorage
from network_builder.builder import username, url, password
from omics_io.omics_io.identifiers import IDS
from network_builder.utils.network_builder_convention import NetworkBuilderConvention
from network_builder.integration_cleanup import clean_network

class TestNetworkBuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storage = Neo4jStorage(url, username=username, 
                                   password=password,
                                   convention=NetworkBuilderConvention())

    def test_run(self):
        clean_network(self.storage)


