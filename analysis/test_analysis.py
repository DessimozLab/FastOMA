
from os import listdir
from ete3 import Tree
from orthoxml_to_newick import orthoxml_to_newick
import pickle
#
# swiss_trees_folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/swisstree_raw/"
# list_swiss_trees = listdir(swiss_trees_folder)
# swiss_trees = []
# swiss_leaves_list = []
# for file in list_swiss_trees:
#     if file.endswith("_v5_sens.nhx"):
#         swiss_tree = Tree(swiss_trees_folder+file)
#         swiss_trees.append(swiss_tree)
#         swiss_leaves = [leaf.name for leaf in swiss_tree.get_leaves()]
#         swiss_leaves_list.append(swiss_leaves)
#
# print("swiss tree read")
#
#
#
# in_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# gene_id_pickle_file = in_folder + "gene_id_30aug_s500.pickle"
#
#
#
# with open(gene_id_pickle_file, 'rb') as handle:
#     gene_id_name = pickle.load(handle)
#
# gene_dic_all = {}
# for query_species_name, list_prots in gene_id_name.items():
#     for (gene_idx_integer, query_prot_name) in list_prots:
#         prot_short = query_prot_name.split("|")[1]
#         gene_dic_all[prot_short] = prot_short+"_"+query_species_name[:-1]
#
# print("gene_dic_all is ready")
#
# # test_xml = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/small_output.xml"
gethog33_xml_file = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/xml_output/out_3sep_all_new_6.xml"

trees_dic_gethog33 = orthoxml_to_newick(gethog33_xml_file)

print("trees are  prepared from orthoxml file")

hog_leaves_list =[]
hog_trees =[]
for sub_hog_name, tree_nhx in trees_dic_gethog33.items():
    hog_tree = Tree(tree_nhx)
    hog_trees.append(hog_tree)
    hog_leaves = [leaf.name for leaf in hog_tree.get_leaves()]
    hog_leaves_list.append(hog_leaves)
print("hog leaves ready ")

if 0:  # with_swiss
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


if 1:  # with gethog2

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gethog2_xml_file = working_folder+ "/archive/public_projects_2020_benchmarkservice_org/HierarchicalGroups_OMA_GETHOGs.orthoxml"

    trees_dic_gethog2 = orthoxml_to_newick(gethog2_xml_file)

    hog_leaves_list_gethog2 =[]
    hog_trees_gethog2 =[]
    for sub_hog_name, tree_nhx in trees_dic_gethog2.items():
        hog_tree = Tree(tree_nhx)
        hog_trees_gethog2.append(hog_tree)
        hog_leaves = [leaf.name for leaf in hog_tree.get_leaves()]
        hog_leaves_list_gethog2.append(hog_leaves)
    print("hog gethog2 leaves ready ", len(hog_leaves_list_gethog2))

    count_list =[ ]
    for hog_v2 in hog_leaves_list_gethog2[:5000]:
        count_list_ = []
        max_intersection = 0
        for hog_v3 in hog_leaves_list:
            len_intersection = len(set(hog_v2) & set(hog_v3))
            if len_intersection != 0:
                count_list_.append(len_intersection)
        count_list.append(count_list_)

for count1 in count_list:
    if sum(count1) > 20:
        sorted_c = sorted(count1, reverse=True)
        print(sorted_c[:20])

    print("trees are  prepared from orthoxml file")
#
#
#
# if len_intersection > max_intersection:
#     max_idx = hog_idx
#     max_intersection = len_intersection
#     max_intersection_leaves = list(set(hog_leaves) & set(swiss_leaves))
