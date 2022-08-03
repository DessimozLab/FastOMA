# # def traverse_tree_recursively(sub_species_tree):
# #     children_nodes = sub_species_tree.children
# #     for node in children_nodes:
# #         if not node.is_leaf():
# #             traverse_tree_recursively(node)
# #             infer_hog_a_level(node)
# #     if sub_species_tree.is_root():
# #         infer_hog_a_level(sub_species_tree)
# #     return 1
# #
# # def infer_hog_a_level(node):
# #     print(node.name)
# #     return 1
# # from ete3 import Tree
# #mytree = Tree('(((H,K)HK,(F,I)FI)FIHK,E)FIHKE;', format=1)
# #print(mytree)
# #traverse_tree_recursively(mytree)
#
#
#
# from ete3 import Tree
# mytree = Tree('(((H,K)HK,(F,I)FI)FIHK,E)FIHKE;', format=1)
#
#
# def infer_hogs_for_a_rhog(sub_species_tree, rhog_i, species_names_rhog, rhogid_num, gene_trees_folder):
#     if sub_species_tree.is_leaf():
#         children_nodes = []
#     else:
#         children_nodes =  sub_species_tree.children
#
#     for node_species_tree_child in children_nodes:
#        infer_hogs_for_a_rhog(node_species_tree_child, rhog_i, species_names_rhog, rhogid_num, gene_trees_folder)
#
#     HOG_this_level_list = infer_HOG_thisLevel(sub_species_tree)
#
#     return HOG_this_level_list
#
#
# def infer_HOG_thisLevel(sub_species_tree):
#     print(sub_species_tree.name)
#     return 10
#
#
# sub_species_tree = mytree
#
# infer_hogs_for_a_rhog(sub_species_tree, [], [], [ ], [])
#

#
# for i in range(6):
#
#     if i > 3:
#         continue
#     print(i)

