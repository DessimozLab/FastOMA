


import _utils
import _inferhog
import sys



address_rhogs_folder =  "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf/rhogs/204/" # sys.argv[1]
inferhog_concurrent_on = True  # sys.argv[2]

print("input is", address_rhogs_folder)

list_rhog_fastas_files = _utils.list_rhog_fastas(address_rhogs_folder)
print("there are ",len(list_rhog_fastas_files),"rhogs in the input folder")
folder = address_rhogs_folder.split("/")[-2]
print("rhogs in the input folder",folder)

list_done_rhogid = []

list_rhog_fastas_files_rem = [i for i in list_rhog_fastas_files if i not in list_done_rhogid]

print("there are ",len(list_rhog_fastas_files_rem),"rhogs remained in the input folder")

hogs_rhog_xml_batch = _inferhog.read_infer_xml_rhogs_batch(list_rhog_fastas_files_rem,inferhog_concurrent_on, folder)


print("finsihed ",address_rhogs_folder)



# to do   FileNotFoundError: [Errno 2] No such file or directory: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3//gethog3_nf_//pickles_rhog//file_69696.pickle'





# import pyham # pyham package for ham analysis
# import pandas as pd # pandas for dataframes
#
#
#
# working_folder="/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/"
# nwk_path = working_folder+"concatanted_363.fasta.contree_edited.nwk"
#
# tree_str = pyham.utils.get_newick_string(nwk_path, type="nwk")
# tree_str[:20]
# # Then you select your favorite orthoXML file
# orthoxml_path =  working_folder+ "pickle_0.5_3000_6Nov_.xml"
# # removing from tree internal node    <property name="TaxRange" value="internal_247_rhg602692"/>
# #
# # with open(orthoxml_path + "_ed.xml", 'w') as orthoxml_file_edited:
# #     with open(orthoxml_path, 'r') as orthoxml_file:
# #         for val in orthoxml_file:
# #             val_edited = val
# #             if val.strip().startswith("<property name=\"TaxRange\" value="):
# #                 val_edited = "_".join(val.split("_")[:2]) + "\"/>\n"
# #
# #             orthoxml_file_edited.write(val_edited)
# #
# # from ete3 import Tree
# # species_tree = Tree(nwk_path)
# # counter_internal = 0
# # for node in species_tree.traverse(strategy="postorder"):
# #     node_name = node.name
# #     num_leaves_no_name = 0
# #     if len(node_name) < 1:
# #         if node.is_leaf():
# #             node.name = "leaf_" + str(num_leaves_no_name)
# #         else:
# #             node.name = "internal_" + str(counter_internal)
# #             counter_internal += 1
# # # print("Working on the following species tree.")
# # # print(species_tree)
# # species_tree.set_outgroup("SPHPU_")
# #
# # species_tree.write(format=0, outfile=nwk_path+"_edit.nwk")
#
#
#
# ham_analysis = pyham.Ham(nwk_path+"_edit.nwk", orthoxml_path+"_ed.xml", tree_format='newick') # , use_internal_name=False
#
#
# a=2
#
# print(a)
#
#
# # import xml.etree.ElementTree as ET
# # import dill as dill_pickle
# # from os import listdir
# # from xml.dom import minidom
# # import os
# #
# #
# #
# #
#
#
#
# #
# # print("started ")
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/ali_code_31aug/"
# # # gene_trees_folder = ""  # working_folder + "/gene_trees_/"
# # # check gene_trees_folder exist otherwise mkdir this
# #
# # #address_rhogs_folder = working_folder + "/rhog_g501_done/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
# # #species_tree_address = working_folder + "/archive/lineage_tree_qfo.phyloxml"
# # pickle_folder = working_folder + "/pickle_2sep5pm/"
# # # add warning when pickle folder is not empty
# # output_xml_name = "out_ali_2sep5pm_test.xml"
# # gene_id_pickle_file = working_folder + "group_xml_ortho_adjusted_family_40_2sep5pm.pickle"
# #
# # #
# # # orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
# # #                                                "originVersion": "Nov 2021", "version": "0.3"})  #
# #
# # if os.path.getsize(gene_id_pickle_file) > 0:
# #     with open(gene_id_pickle_file, 'rb') as handle:
# #         (groups_xml, gene_id_name, orthoxml_file) = dill_pickle.load(handle)
# #     # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
# # print("gene_id_name read ")
# #
# # #
# # #
# # # for query_species_name, list_prots in gene_id_name.items():
# # #
# # #     species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
# # #     database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
# # #     genes_xml = ET.SubElement(database_xml, "genes")
# # #
# # #     for (gene_idx_integer, query_prot_name) in list_prots:
# # #         query_prot_name_pure = query_prot_name.split("||")[0].strip().split("|")[1]
# # #         gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
# #
# # print("gene_xml created ")
# # pickle_files_adress = listdir(pickle_folder)
# # print("number of hog is", len(pickle_files_adress))
# # hogs_a_rhog_xml_all = []
# # list_len = []
# # for pickle_file_adress in pickle_files_adress:
# #     if os.path.getsize(pickle_folder + pickle_file_adress) > 0:
# #         with open(pickle_folder + pickle_file_adress, 'rb') as handle:
# #             hogs_a_rhog_xml_batch = dill_pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
# #             hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
# #             list_len.append(len(hogs_a_rhog_xml_batch))
# #     else:
# #         print("empty ", pickle_file_adress)
# #         # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.
# #
# # print("after reading number of hogs is  ", len(hogs_a_rhog_xml_all))
# # # orthoxml_file2 = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
# # #                                                "originVersion": "Nov 2021", "version": "0.3"})  #
# #
# # # groups_xml2 = ET.SubElement(orthoxml_file2, "groups")
# #
# # for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
# #     groups_xml.append(hogs_a_rhog_xml)
# #
# # #xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
# # # print(xml_str[:-1000])
# #
# # with open(working_folder +output_xml_name, "w") as file_xml:
# #     file_xml.write(xml_str)
# # file_xml.close()
# #
# # print("orthoxml is written in  "+ working_folder +output_xml_name)
# #
# #
# # a=232
# #
# #
#
# #
# # import xml.etree.ElementTree as ET
# # import dill as dill_pickle
# # from os import listdir
# # from xml.dom import minidom
# # import os
# # from Bio import SeqIO
# # #import dill as dill_pickle
# # import dill as pickle
# # #import pickle
# #
# #
# #
# # address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/ali_code_31aug/"
# #
# # address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/rhog_all_v3_g2_s500/"
# # address_group_xml_ortho = address_working_folder+"group_xml_ortho_adjusted_family_40_2sep5pm_dill.pickle"
# #
# #
# # rhog_files = listdir(address_rhogs_folder)[:]
# #
# # rhog_files = listdir(address_rhogs_folder)
# # rhogid_num_list = []
# # for rhog_file in rhog_files:
# #     if rhog_file.split(".")[-1] == "fa":
# #         rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
# #         rhogid_num_list.append(rhogid_num)
# #
# # rhogid_num_list_temp = rhogid_num_list
# #
# # species_prot_dic = {}
# # # all_prot_temp_list= []
# # for rhogid_num in rhogid_num_list_temp:
# #     prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
# #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# #     for prot_i in rhog_i:
# #         prot_i_name = prot_i.id  # .split("||")[0] # .split("|")[1]  # tr|E3JPS4|E3JPS4_PUCGT or new || ||
# #         species_i = prot_i.id.split("||")[1][:-1]  # prot_i.id.split("|")[-1].split("_")[-1]
# #         if species_i in species_prot_dic:
# #             species_prot_dic[species_i].append(prot_i_name)
# #         else:
# #             species_prot_dic[species_i] = [prot_i_name]
# #         # all_prot_temp_list.append(prot_i.id)
# #
# # print("there are species ", len(species_prot_dic))
# # orthoxml_file = ET.Element("orthoXML",
# #                            attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA", "originVersion": "Nov 2021",
# #                                    "version": "0.3"})  #
# #
# # gene_counter = 100000
# # gene_id_name = {}
# # query_species_names_rHOGs = list(species_prot_dic.keys())
# # for species_name in query_species_names_rHOGs:
# #     no_gene_species = True  # for code develop ment
# #     species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
# #     database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
# #     genes_xml = ET.SubElement(database_xml, "genes")
# #
# #     prot_list = species_prot_dic[species_name]
# #     for prot_itr in range(len(prot_list)):  # [12:15]
# #         prot_i_name = prot_list[prot_itr]
# #         gene_id_name[prot_i_name] = gene_counter
# #         prot_i_name_short = prot_i_name.split("||")[0].split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
# #         gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
# #         gene_counter += 1
# #
# # groups_xml = ET.SubElement(orthoxml_file, "groups")
# #
# #
# #
# # with open(address_group_xml_ortho, 'wb') as handle:
# #     # dill_pickle.dump(gene_id_name, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
# #     pickle.dump((groups_xml, gene_id_name, orthoxml_file), handle, protocol=pickle.HIGHEST_PROTOCOL)
# #
# # print("saved as ", address_group_xml_ortho)
# #
# #
# #
# #
# #
#
#
#
# # import _utils
# # from Bio import SeqIO
# # import os
# #
# #
# # print("started ")
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # address_rhogs_folder = working_folder + "rhog_all_v3_s500/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
# # pickl_out_folder = working_folder + "pickle_folder_g2_s500/"  # rhog_all_v3_g2_s500
# # pickle_folder = working_folder + "pickle_folder_30aug3pm/"
# #
# # not_c = 0
# # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# # for rhogid_num in rhogid_num_list:
# #     prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
# #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# #
# #     if len(rhog_i) > 2:
# #         os.popen('cp '+pickle_folder+'file_'+str(rhogid_num)+".pickle " + pickl_out_folder+'file_'+str(rhogid_num)+".pickle" )
# #
# #     else:
# #         not_c += 1
# # print("not part of ", not_c)
#
#     # import dill as dill_pickle
# # import os
# # from os import listdir
# # working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# # pickle_folder = working_folder + "/pickle_folder_29aug/"
# #
# #
# # orthoxml_to_newick.py=2
# # with open(pickle_folder + "file_9114.pickle", 'rb') as handle:
# #     hogs_a_rhog_xml_batch = dill_pickle.load(handle)
# #     orthoxml_to_newick.py=23
# #     orthoxml_to_newick.py=5
#
#
#
#
# #pickle_files_adress = listdir(pickle_folder)
# # print("number of hog is", len(pickle_files_adress))
# # hogs_a_rhog_xml_all = []
# # for pickle_file_adress in pickle_files_adress:
# #     if os.path.getsize(pickle_folder + pickle_file_adress) > 0:
# #         with open(pickle_folder + pickle_file_adress, 'rb') as handle:
# #             hogs_a_rhog_xml_batch = dill_pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
# #             hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
# #
#
#
#
# # from _utils import logger_hog
# # import _utils_rhog
# # import _utils
# # # from distributed import get_client
# # # from dask.distributed import rejoin, secede
# # from Bio import SeqIO
# #
# #
# # if __name__ == '__main__':
# #     working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
# #     gene_trees_folder = "" # working_folder + "/gene_trees_/"
# #     # check gene_trees_folder exist otherwise mkdir this
# #     pickle_folder = ""
# #     species_tree_address = ""
# #
# #
# #     address_rhogs_folder = working_folder + "/rhog_all_v2_filbig/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
# #     address_rhogs_folder_fil = working_folder + "/"
# #     file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)
# #
# #
# #     oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastoma/archive/OmaServer.h5"
# #     # working_folder+"omamer_database/oma_path/OmaServer.h5"
# #     print("rHOG inferece has started. The oma database address is in ", oma_database_address)
# #     (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
# #     (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
# #     hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)
# #
# #     (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_fscore_allspecies,
# #     prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
# #
# #     print("start finding filtered proteins")
# #     prot_name_filt = [ ]
# #     for species_i, query_prot_names in enumerate(query_prot_names_species_mapped):
# #         print(species_i)
# #         prots_hogmap_seqlen = prots_hogmap_seqlen_allspecies[species_i]
# #         prots_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen_allspecies[species_i]
# #         prots_hogmap_fscore = prots_hogmap_fscore_allspecies[species_i]
# #         for prot_i, prot_name in enumerate(query_prot_names):
# #             prot_hogmap_seqlen = prots_hogmap_seqlen[prot_i]
# #             prot_hogmap_subfmedseqlen = prots_hogmap_subfmedseqlen[prot_i]
# #             prot_hogmap_fscore = prots_hogmap_fscore[prot_i]
# #             if prot_hogmap_seqlen != "na" and prot_hogmap_subfmedseqlen != "na":
# #                 prot_hogmap_subfmedseqlen = int(prot_hogmap_subfmedseqlen)
# #                 prot_hogmap_seqlen = int(prot_hogmap_seqlen)
# #                 prot_hogmap_fscore= float(prot_hogmap_fscore)
# #                 if prot_hogmap_subfmedseqlen * 0.8 < prot_hogmap_seqlen < prot_hogmap_subfmedseqlen * 1.8 and prot_hogmap_fscore > 0.8:
# #                     prot_name_filt.append(prot_name)
# #     print("finished finding filtered proteins")
# #     rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
# #
# #     for rhogid_num in [811184]:# rhogid_num_list
# #         print("working on rhog", rhogid_num )
# #         prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
# #         rhog_i = list(SeqIO.parse(prot_address, "fasta"))
# #         rhog_i_filtered = []
# #         for rhog in rhog_i:
# #             if rhog.name.split("||")[0] in prot_name_filt:
# #                 rhog_i_filtered.append(rhog)
# #         print("writing rhog", rhogid_num)
# #         SeqIO.write(rhog_i_filtered, address_rhogs_folder_fil + "HOG_B" + str(rhogid_num).zfill(7) + ".fa", "fasta")
# #     print("finished writing filtered rhogs")
# #
# #
