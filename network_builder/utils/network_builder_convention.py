from dataclasses import fields
from neo4j_interface.utils.conventions import DefaultConvention
from omics_io.omics_io.identifiers import IDS

def _collect_labels():
    labs = []
    for f in fields(IDS.type):
        labs.append(getattr(IDS.type, f.name))
    for f in fields(IDS.omics_type):
        labs.append(getattr(IDS.omics_type, f.name))
    return sorted(set(labs))

class NetworkBuilderConvention(DefaultConvention):
    def __init__(self):
        super().__init__()

    def get_constraints(self) -> dict:
        labels = _collect_labels()
        return {
            "constraints": [
                {
                    "name": "uniq_id",
                    "entity": "node",
                    "kind": "unique",
                    "properties": [self.id_prop],
                    "labels": labels,
                }
            ],
            "indexes": [
                {
                    "name": "idx_alias",
                    "entity": "node",
                    "kind": "btree",
                    "properties": [IDS.predicates.alias],
                    "labels": labels,
                }
            ],
        }
