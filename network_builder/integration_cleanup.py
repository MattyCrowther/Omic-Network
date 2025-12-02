from network_builder.neo4j_interface.storage import Neo4jStorage
from network_builder.neo4j_interface.utils.storage_objects import NodeObject
from omics_io.omics_io.identifiers import IDS

P_ALIAS = IDS.predicates.alias

def clean_network(storage:Neo4jStorage):
    '''
    Idea here is to apply simple rules for cleaning up artifacts or 
    inconsistencies inserted into the network when it was built. 
    These are likely things that would be unknownable until the 
    network was built such as missed matches due to type.
    '''
    _handle_unknown_nodes(storage)
    _handle_multiple_products(storage)



def _handle_no_node(storage:Neo4jStorage,node:NodeObject):
    # Note see we have with_relationships
    raise ValueError("WHY DOES A NODE HAVE NO EDGES????")
    if len(node.relationships) != 0:
        print(f"Found Unknown with no similar nodes {node} but with rels...")
        exit()
    
    rels = storage.find_relationships(right_id=node.id)
    if len(rels) != 0:
        print(f"Found Unknown with no similar nodes {node} but with incoming rels...")
        exit()

    lab_node = storage.find_nodes(node.id)
    if len(lab_node) != 0:
        print(f"Found Unknown {node} but with different label. ")
        exit()

    alias_nodes = storage.find_nodes(node.properties.get(P_ALIAS))
    if len(alias_nodes) != 0:
        print(f"Found Unknown {node} but with different label. ")
        exit()

    

def _handle_unknown_nodes(storage:Neo4jStorage):
    '''
    Case 1:
        UN doesnt have another node with the same ID.
    Case 2: 
        UN has a single node with the same ID.
        Merge with node irrespective of its layer.
    Case 3:
        UN has multiple nodes with the same ID.
        We try find one of the nodes that is in the same layer.
    '''
    direct_merges = []
    results = {"Direct Label Merge" : direct_merges}
    for i in storage.find_nodes(label=IDS.type.unknown,with_relationships=True):
        existing_nodes = list_diff(storage.find_nodes([i.id]), [i])
        if len(existing_nodes) == 0:
            _handle_no_node(storage,i)

        if len(existing_nodes) == 1:
            storage.merge_nodes([i.id],
                                canonical_label=existing_nodes[0].label)
            direct_merges.append(i.id)
        else:
            print(f"Found Unknown with multiple similar nodes {i} - {existing_nodes}")
            exit()
            i_layer = _find_layer(storage,i)
            for e in existing_nodes:
                e_layer = _find_layer(storage,e)

    return results


def _handle_multiple_products(storage:Neo4jStorage):
        # Assume gene has single product (Prokaryotes),
        for i in storage.find_nodes(label=IDS.type.gene,with_relationships=True):
            if IDS.predicates.product in i.relationships:
                products = i.relationships[IDS.predicates.product]
                if len(products) > 1:
                    storage.merge_nodes(products)

def _find_layer(storage:Neo4jStorage,node):
    o_t = None
    for r in storage.find_relationships(None,node.id):
        if r.end_id in IDS.omics_type:
            if o_t is not None:
                raise ValueError(f'{node} has multiple layers.')
            o_t = r.end_id.type

def list_diff(l1, l2):
    return [x for x in l1 if x not in l2]