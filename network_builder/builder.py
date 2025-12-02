import yaml
from pathlib import Path
import math
import os
from collections import defaultdict

from neo4j_interface.storage import Neo4jStorage
from neo4j_interface.utils.storage_objects import RelationshipRow, NodeRow
from network_builder.utils.network_builder_convention import NetworkBuilderConvention
from omics_io.omics_io.identifiers import IDS
from omics_io.omics_io.parse_obj import OmicData

config_path = Path(__file__).parent / "config.yaml"

with open(config_path, "r") as f:
    config = yaml.safe_load(f)

url = config["STORAGE"]["uri"]
username = config["STORAGE"]["username"]
password = config["STORAGE"]["password"]

matrix_store = "measure_store"
if not os.path.isdir(matrix_store):
    os.mkdir(matrix_store)

# Things to consider (No order)

# Currently, name isnt exactly a good metrix for an identifier.
# This is more a failing of omics_io. We could also request an ID,
# which could be the filename perhaps.

# Measurement data, What do we do with this? Its not very useful in a network.
# I reckon, youd end up file system storing it and linking in the MD.
# Like Matrix: "blabla_filename". Decide if attach to node above or below.

c50 = 5
K = math.log(2) / max(1e-9, c50)


def _count_to_conf(count):
    return 1 - math.exp(-K * count)


class OmicGraphBuilder:
    def __init__(self):
        self._storage = Neo4jStorage(
            url,
            username=username,
            password=password,
            convention=NetworkBuilderConvention(),
        )

    def add_omic_set(self, dataset : OmicData):
        nodes_by_label = defaultdict(dict)
        rel_rows = []
        mat_fn = Path(f"{matrix_store}/{dataset.name}")
        self._storage.upsert_node(dataset.omics_type, 
                                  dataset.name, 
                                  props={"location": str(mat_fn)})
        for ele, md in dataset.features():
            label = md.entity
            nodes_by_label[label][ele] = md.properties
            rel_rows.append(RelationshipRow(dataset.name, 
                                            ele, 
                                            None,
                                            dataset.omics_type,
                                            label))
        for ele,md in dataset.columns():
            nodes_by_label[md.entity][ele] = md.properties
            rel_rows.append(RelationshipRow(dataset.name, 
                                            ele, 
                                            None,
                                            dataset.omics_type,
                                            md.entity))

        for label, idmap in nodes_by_label.items():
            rows = [NodeRow(_id,p) for _id, p in idmap.items()]
            self._storage.upsert_nodes(label, rows)
        self._storage.upsert_relationships(IDS.predicates.has_feature, 
                                           rel_rows)

    def add_alignment_data(self, alignment_data):
        gid_to_nid = {}
        unknown_rows = []
        unknown_ids = []

        # A) batch lookup for all candidate ids
        all_ids = set()
        groups = []
        for gid, nodes in alignment_data:
            ids = list(set([n.identifier for n in nodes]))
            types = list(set([n.type for n in nodes]))
            groups.append((gid, ids,types))
            all_ids.update(ids)

        existing = {n.id: n for n in 
                    self._storage.find_nodes(ids=list(all_ids))}
        existing_set = set(existing.keys())

        # B) resolve groups; defer unknown creation and alias adds
        alias_adds = defaultdict(list)
        for gid, ids, types in groups:
            ids = list(set(ids))
            have = [i for i in ids if i in existing_set]
            type = [t for t in types if t != "unknown"]
            assert(len(type) <= 1)
            if len(type) == 0:
                type = None
            else:
                type = type[0]
            if not have:
                # stable unknown id from aliases
                f_id = "UNK:" + ids[0]
                gid_to_nid[gid] = [f_id,type]
                unknown_ids.append(f_id)
                unknown_rows.append(NodeRow(f_id, {IDS.predicates.alias: ids}))
            else:
                f_node = self._storage.merge_nodes(have,label=type)
                gid_to_nid[gid] = (f_node.id,f_node.label)
                missing = [i for i in ids if i not in have]
                if missing:
                    print(f_node.id)
                    alias_adds[f_node.id].extend(missing)

        if unknown_rows:
            self._storage.upsert_nodes(IDS.type.unknown, unknown_rows)

        # batch alias additions per node
        for nid, aliases in alias_adds.items():
            self._storage.add_property(nid, IDS.predicates.alias, 
                                       list(set(aliases)))

        # C) pre-aggregate relationships to dedupe
        rel_bins = defaultdict(list)
        agg = defaultdict(int)
        for g1, rel, g2, cnt in alignment_data.relationships():
            if g1 == g2:
                continue
            s_id,s_type = gid_to_nid.get(g1)
            t_id,t_type = gid_to_nid.get(g2)
            if s_id is None or t_id is None:
                raise ValueError(f"Cant find:{s_id},{t_id}")
            agg[(rel, s_id, t_id, s_type, t_type)] += cnt

        for (rel, s_id, t_id,s_type, t_type), tot in agg.items():
            element = RelationshipRow(s_id,
                                     t_id,
                                     {IDS.predicates.confidence: 
                                      _count_to_conf(tot)},
                                     _normalize_labels(s_type),
                                     _normalize_labels(t_type))

            rel_bins[rel].append(element)

        for rel_label, rows in rel_bins.items():
            try:
                self._storage.upsert_relationships(rel_label, rows)
            except Exception as ex:
                print(ex)
                print(rel_label)
                print(rows)
                exit()
        return unknown_ids

def _normalize_labels(x):
    if not x:
        return []
    if isinstance(x, str):
        candidates = [x]
    else:
        candidates = list(x)
    labels = [str(v).strip() for v in candidates if v is not None and str(v).strip() != ""]
    return labels