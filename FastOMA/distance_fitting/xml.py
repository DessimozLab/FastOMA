from xml.etree import ElementTree as ET


def xml_to_newick(node, _toplevel=True):
    '''
    Converts HOG to Newick format and returns it. Based on code from hogprop.HOGParser but
    for non-orthoXML defined XML in "orthostub"-like format.
    '''
    if node.tag == 'geneRef':
        name = node.attrib['id']
        children = None
    elif node.tag == 'orthologGroup' or node.tag == 'paralogGroup':
        name = ''
        children = list(filter(lambda c: c is not None,
                               map(lambda n: xml_to_newick(n, _toplevel=False),
                                   node)))
    else:
        return

    newick = (name if not children else
              '({:s}){:s}'.format(', '.join(children), name))
    if _toplevel:
        newick += ';'
    return newick


def map_distances_to_xml(t_node, xml_node, tag_name='BranchLength'):
    '''
    Mapping distances from tree back to hog.
    '''
    valid_node = lambda n: (n.tag in {'geneRef', 'orthologGroup', 'paralogGroup'})

    if valid_node(xml_node):
        d = t_node.edge.length
        if d is not None:
            ET.SubElement(xml_node,
                          "score",
                          attrib={"id": tag_name, "value": f'{d:.03f}'})

        for c in zip(t_node.child_nodes(), filter(valid_node, xml_node)):
            map_distances_to_xml(*c, tag_name=tag_name)


def save_stats_to_xml(xml_node, tag_name='DistanceFitting', **stats):
    '''
    Saving statistics to the xml
    '''
    for (k, v) in stats.items():
        ET.SubElement(xml_node,
                      "score",
                      attrib={"id": f"{tag_name}_{k}", "value": f'{v}'})
