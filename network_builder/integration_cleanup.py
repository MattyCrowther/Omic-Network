from network_builder.neo4j_interface.storage import Neo4jStorage
from network_builder.neo4j_interface.utils.storage_objects import NodeObject
from omics_io.omics_io.identifiers import IDS
from neo4j_interface.utils.storage_objects import RelationshipRow

P_ALIAS = IDS.predicates.alias

def clean_network(storage:Neo4jStorage):
    '''
    Idea here is to apply simple rules for cleaning up artifacts or 
    inconsistencies inserted into the network when it was built. 
    These are likely things that would be unknownable until the 
    network was built such as missed matches due to type.
    '''
    #_connect_layers(storage)
    #_handle_unknown_nodes(storage)
    _handle_multiple_products(storage)



def _connect_layers(storage:Neo4jStorage):
    def _get_entities(label):
        omics = storage.find_nodes(label=label, 
                                   with_relationships=True)
        eles = []
        for omic in omics:
            for rel_list in omic.relationships.values():
                eles.extend(rel_list)
        return set(eles)
    
    def _build_rel_row(set1, set2, set1_labels, set2_labels):
        for x in set1 & set2:
            rel_rows.append(RelationshipRow(x, x, None,
                                            set1_labels,
                                            set2_labels))

    g_elems = _get_entities(IDS.omics_type.genomics)
    t_elems = _get_entities(IDS.omics_type.transcriptomics)
    p_elems = _get_entities(IDS.omics_type.proteomics)
    m_elems = _get_entities(IDS.omics_type.metabolomics)
    
    rel_rows = []
    g_labs = [IDS.type.dna,IDS.type.gene,IDS.type.cds]
    t_labs = [IDS.type.rna]
    p_labs = [IDS.type.protein]
    m_labs = [IDS.type.metabolite]
    _build_rel_row(g_elems,t_elems,g_labs,t_labs)
    _build_rel_row(t_elems,p_elems,t_labs,p_labs)
    _build_rel_row(p_elems,m_elems,p_labs,m_labs)

    storage.upsert_relationships(IDS.predicates.maps,rel_rows)

def _handle_no_node(storage:Neo4jStorage,node:NodeObject):
    # Note see we have with_relationships
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
        prot_fams = []
        for i in storage.find_nodes(label=IDS.type.rna,with_relationships=True):
            if IDS.predicates.product in i.relationships:
                products = i.relationships[IDS.predicates.product]
                to_add = []
                to_merge = []
                if len(products) > 1:
                    for p in products:
                        p_rels = storage.find_relationships(right_id=p)
                        if len(p_rels) > 1:
                            to_add.append(p)
                            prot_fams.append(p)
                        else:
                            to_merge.append(p)
                    m_node = storage.merge_nodes(to_merge)
                    for a in to_add:
                        storage.add_property(m_node.id,
                                             IDS.predicates.alias,
                                             a,m_node.label)
        storage.remove_nodes(prot_fams,label=IDS.type.protein)

def _find_layer(storage:Neo4jStorage,node):
    o_t = None
    for r in storage.find_relationships(None,node.id):
        if r.end_id in IDS.omics_type:
            if o_t is not None:
                raise ValueError(f'{node} has multiple layers.')
            o_t = r.end_id.type

def list_diff(l1, l2):
    return [x for x in l1 if x not in l2]