from dendropy import Tree
from .xml import xml_to_newick, map_distances_to_xml, save_stats_to_xml
from .topological_matrix import TopologicalMatrix


def fit_distances(hog_xml, dist_df, delta=1e-3, tag_name='BranchLength'):
    '''
    Fit distances to HOG structure using matrix method, using dendropy intermediary
    1. convert to newick
    2. load with dendropy
    3. convert to matrix
    '''
    nwk = xml_to_newick(hog_xml)
    t = Tree.get_from_string(nwk, schema='newick')

    pairs = dist_df[['geneRef1', 'geneRef2']].values
    d = dist_df['Distance'].values

    tm = TopologicalMatrix(t)
    (it, r, t_solve) = tm.solve(pairs, d, delta=delta)

    map_distances_to_xml(tm.t.seed_node, hog_xml, tag_name=tag_name)
    save_stats_to_xml(hog_xml, tag_name=tag_name, NrIter=it, Resid=r, SolveTime=t_solve)
    #return tm.t.as_string('newick')
