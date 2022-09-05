
from os import listdir
from ete3 import Tree
from orthoxml_to_newick import orthoxml_to_newick


swiss_trees_folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/swisstree_raw/"
list_swiss_trees = listdir(swiss_trees_folder)
swiss_trees = []
swiss_leaves_list = []
for file in list_swiss_trees:
    if file.endswith("_v5_sens.nhx"):
        swiss_tree = Tree(swiss_trees_folder+file)
        swiss_trees.append(swiss_tree)
        swiss_leaves = [leaf.name for leaf in swiss_tree.get_leaves()]
        swiss_leaves_list.append(swiss_leaves)

print("swiss tree read")



working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
gene_id_pickle_file = working_folder + "gene_id_30aug_s500.pickle"

import pickle

with open(gene_id_pickle_file, 'rb') as handle:
    gene_id_name = pickle.load(handle)

gene_dic_all = {}
for query_species_name, list_prots in gene_id_name.items():
    for (gene_idx_integer, query_prot_name) in list_prots:
        prot_short = query_prot_name.split("|")[1]
        gene_dic_all[prot_short] = prot_short+"_"+query_species_name[:-1]

print("gene_dic_all is ready")

# test_xml = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/small_output.xml"
test_xml = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/out_3sep_all_new_6.xml"

trees_dic = orthoxml_to_newick(test_xml)

print("trees are  prepared from orthoxml file")





hog_leaves_list =[]
hog_trees =[]
for sub_hog_name, tree_nhx in trees_dic.items():
    hog_tree = Tree(tree_nhx)
    hog_trees.append(hog_tree)
    hog_leaves = [leaf.name for leaf in hog_tree.get_leaves()]
    hog_leaves_list.append(hog_leaves)
print("hog leaves ready ")

print("first swiss, then hog")
swiss_tree_idx = 0
swiss_leaves = swiss_leaves_list[swiss_tree_idx]
for swiss_tree_idx, swiss_leaves in enumerate(swiss_leaves_list):
    print("\n"+"**"*4, swiss_tree_idx, "\n")

    max_intersection = 0
    max_idx = -1
    for hog_idx, hog_leaves in enumerate(hog_leaves_list):
        intersection_leaves = set(hog_leaves) & set(swiss_leaves)
        len_intersection = len(intersection_leaves)
        if len_intersection > max_intersection:
            max_idx = hog_idx
            max_intersection = len_intersection
            max_intersection_leaves = list(set(hog_leaves) & set(swiss_leaves))


    if max_idx > 0:
        print("max ", max_intersection)
        swiss_tree = swiss_trees[swiss_tree_idx]
        swiss_tree.prune(max_intersection_leaves)
        for node in swiss_tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node.name = gene_dic_all[node.name]

        print(swiss_tree.write(),"\n")
        hog_tree = hog_trees[max_idx]
        hog_tree.prune(max_intersection_leaves)
        for node in hog_tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node.name = gene_dic_all[node.name]

        print(hog_tree.write())

    else:
        print("no intersection")


        for node in swiss_tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node.name = gene_dic_all[node.name]


#
# for node in species_tree.traverse(strategy="postorder"):
# node_name = node.name
# node_children = node.children


# import  pickle
# with open(address_group_xml_ortho, 'wb') as handle:
#     # dill_pickle.dump(gene_id_name, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
#     pickle.dump((groups_xml, gene_id_name, orthoxml_file), handle, protocol=pickle.HIGHEST_PROTOCOL)


