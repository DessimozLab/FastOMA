



from ete3 import Tree


sub_species_tree= Tree("((H:0.3,I:0.1):0.5, A:1, (B:0.4,(C:1,D:1):0.5):0.5);")

for node in sub_species_tree.traverse(strategy="postorder"):
  if not node.is_leaf():
    node_children = [i.name for i in node.children]
    node.name = "_".join(node_children)
#
# for node_species_tree in sub_species_tree.traverse(strategy="levelorder", is_leaf_fn=lambda n: hasattr(n,"processed") and n.processed == True):
#   species_leaves_names = [i.name for i in node_species_tree.get_leaves()]
#   print(node_species_tree.name)
#
#   if len(species_leaves_names) <= 3:
#     #print("here")
#     for node in node_species_tree.traverse():
#       node.processed = True
#



def infer_hogs_for_rhog_levels_recursively(sub_species_tree):

    if sub_species_tree.is_leaf():
      if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
        print("* ",sub_species_tree.name)
      return 1
    children_nodes = sub_species_tree.children
    hogs_children_level_list = []
    for node_species_tree_child in children_nodes:
        hogs_children_level_list_i = infer_hogs_for_rhog_levels_recursively(node_species_tree_child)
    if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
      print("* ", sub_species_tree.name)
    return 1


for node_species_tree in sub_species_tree.traverse(strategy="postorder"):
  species_leaves_names = [i.name for i in node_species_tree.get_leaves()]
  if len(species_leaves_names) > 1:
    infer_hogs_for_rhog_levels_recursively(node_species_tree)
    for node in node_species_tree.traverse():
      node.processed = True







    # hogs_a_rhog_1 = infer_hogs_for_rhog_levels_recursively(node_species_tree, recursive_4inputs)

    #secede()
    #hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively, node_species_tree, recursive_4inputs)
    #hogs_a_rhog = hogs_a_rhog_future.result()
    #rejoin()


  # improvment ?? i don't need to write all hogs of differetn taxonomic level as  picke files




# import xml.etree.ElementTree as ET
# # import dill as dill_pickle
# import pickle
# from os import listdir
# from xml.dom import minidom
# import os
#
# print("started ")
# in_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# pickle_folder = in_folder + "/pickle_folder_3sep10pm_allv3_new/"
#
#
# output_xml_name = "out_3sep_all_new_6.xml"
# gene_id_pickle_file = in_folder + "gene_id_28aug.pickle"
#
#
# if os.path.getsize(gene_id_pickle_file) > 0:
#     with open(gene_id_pickle_file, 'rb') as handle:
#         gene_id_name = pickle.load(handle)   # (groups_xml, gene_id_name, orthoxml_file)
#
# print("gene_id_name read ")
#
#
# orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
#                                                 "originVersion": "Nov 2021", "version": "0.3"})  #
#
#
# for query_species_name, list_prots in gene_id_name.items():
#     species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
#     database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
#     genes_xml = ET.SubElement(database_xml, "genes")
#
#     for (gene_idx_integer, query_prot_name) in list_prots:
#         query_prot_name_pure = query_prot_name.split("||")[0].strip().split("|")[1]
#         gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
# print("gene_xml created ")
#
#
# pickle_files_adress = listdir(pickle_folder)
# print("number of hog is", len(pickle_files_adress))
# hogs_a_rhog_xml_all = []
# print("Following files are with 0 byte (if any, so hog inferecnce should be repeated after deleting the files):")
# for pickle_file_adress in pickle_files_adress:
#     if os.path.getsize(pickle_folder + pickle_file_adress) > 0:
#         with open(pickle_folder + pickle_file_adress, 'rb') as handle:
#             hogs_a_rhog_xml_batch = pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
#             hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
#     else:
#         print(pickle_file_adress)
#
# print("after reading number of hogs is  ", len(hogs_a_rhog_xml_all))
#
#
#
# groups_xml = ET.SubElement(orthoxml_file, "groups")
#
# for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
#     groups_xml.append(hogs_a_rhog_xml)
#
# xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
#
#
# with open(in_folder +output_xml_name, "w") as file_xml:
#     file_xml.write(xml_str)
# file_xml.close()
#
# print("orthoxml is written in  "+ in_folder +output_xml_name)
#
#
#
# #
#
#
# # from xml.dom import minidom
# # import xml.etree.ElementTree as ET
# # import _utils
# # import _inferhog
# #
# # # from _utils import logger_hog
# # import FastOMA._utils_rhog as _utils_rhog
# #
# # # from distributed import get_client
# # # from dask.distributed import rejoin, secede
# #
# # if __name__ == '__main__':
# #     in_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# #     gene_trees_folder = "" # in_folder + "/gene_trees_/"
# #     # check gene_trees_folder exist otherwise mkdir this
# #
# #     address_rhogs_folder = in_folder + "/rhog_g10k/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
# #     file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)
# #
# #
# #     oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/FastOMA/archive/OmaServer.h5"
# #     # in_folder+"omamer_database/oma_path/OmaServer.h5"
# #     print("rHOG inferece has started. The oma database address is in ", oma_database_address)
# #     (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
# #     (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, in_folder)
# #     hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, in_folder)
# #
# #     (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
# #     prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
# #
# #     for prot_i, prot in enumerate(query_prot_names_species_mapped):
# #
# #         orthoxml_to_newick.py=1
# #
# #
# #
# #
# #     rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# #
# #
# #
# #
#
