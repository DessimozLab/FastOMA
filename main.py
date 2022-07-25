
"""
import os
from xml.dom import minidom
#import concurrent.futures

import dask
from dask.distributed import LocalCluster
from dask_jobqueue import SLURMCluster
import gc
"""

from dask.distributed import Client

import _utils
import _inferhog
from _utils import logger_hog


if __name__ == '__main__':

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"

    gene_trees_folder = working_folder + "/gene_trees_test/"
    address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # "/rhog_size_g2_s500/" sample_rootHOG
    species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"

    step = "hog"
    print("we are here ")
    if step == "roothog":
        """
        Structure of folders:
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address= the address to the oma databases

        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is a HOG.
        """
        # import pyoma.browser.db as db
        # oma_database_address = address_working_folder+"omamer_database/oma_path/OmaServer.h5"
        # print("program has started. The oma database address is in ",oma_database_address)
        # (oma_db, list_oma_speices) = parse_oma_db(oma_database_address)
        # (query_species_names, query_prot_records_species) = parse_proteome(list_oma_speices)
        # query_prot_records_species = add_species_name(query_prot_records_species,query_species_names)
        # hogmap_allspecies_elements = parse_hogmap_omamer(query_species_names)
        # (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
        # prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
        # query_prot_records_species_filtered =  filter_prot_mapped(query_species_names, query_prot_records_species, query_prot_names_species_mapped)
        # print(len(query_prot_records_species_filtered),len(query_prot_records_species_filtered[0]))
        # (rhogid_num_list, rhogids_prot_records_query) = group_prots_rootHOGs(prots_hogmap_hogid_allspecies, address_rhogs_folder)

    rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
    logger_hog.info("Number of root hog is "+str(len(rhogid_num_list))+".")
    print(rhogid_num_list[:2])
    # rhogid_num_list_temp = [836500]  # rhogid_num_list[23]  # [833732]

    from Bio import SeqIO
    # def traverse_tree_recursively(sub_species_tree):
    #     """
    #     test function:
    #     from ete3 import Tree
    #     mytree = Tree('(((H,K)HK,(F,I)FI)FIHK,E)FIHKE;', format=1)
    #     print(mytree)
    #     traverse_tree_recursively(mytree)
    #     """
    #     children_nodes = sub_species_tree.children
    #     for node in children_nodes:
    #         if not node.is_leaf():
    #             traverse_tree_recursively(node)
    #             infer_hog_a_level(node)
    #
    #     if sub_species_tree.is_root():
    #         infer_hog_a_level(sub_species_tree)
    #     return 1


    # HOG_thisLevel_list = []
    # len_HOG_thisLevel_list = []
    # HOG_thisLevel_xml_all = []
    #
    # for rhogid_num in rhogid_num_list[4:6]:
    #
    #     logger_hog.info("\n"+"="*50+"\n"+"Working on root hog: "+str(rhogid_num)+". \n")  # +", ",rhogid_num_i,"-th. \n"
    #     prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
    #     rhog_i = list(SeqIO.parse(prot_address, "fasta"))
    #     logger_hog.info("number of proteins in the rHOG is "+str(len(rhog_i))+".")
    #
    #     (species_tree) = utils.read_species_tree(species_tree_address)
    #     (species_tree, species_names_rhog, prot_names_rhog) = utils.prepare_species_tree(rhog_i, species_tree)
    #     # species_tree.write();  print(species_tree.write())
    #
    #     dic_sub_hogs = {}
    #
    #     # for node_species_tree in species_tree.traverse(strategy="postorder"):
    #     #     if node_species_tree.is_leaf():
    #     #         # each leaf itself is a subhog
    #     #         continue
    #     #
    #     #     dic_sub_msas = []
    #     #     print(node_species_tree.name)
    #     #     (dic_sub_hogs) = utils.infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num, gene_trees_folder)
    #     #     exit()
    #
    #
    #
    #     (dic_sub_hogs) = infer_hogs_for_a_rhog(species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    #                                                        rhogid_num, gene_trees_folder)
    #





    """
        to do :

            check whether delay work on infer_hogs_for_a_rhog

            create a function  for running   ( for rhogid_num in rhogid_num_list[4:6]: )
                input list of rhg num
                write xml

                think how pickle per level ?

                think how to distribute rhog into list


    """

    rhogid_num_list_input =  rhogid_num_list[6:9]
    (groups_xml, gene_id_name, orthoxml_file) = _utils.prepare_xml(rhogid_num_list_input, address_rhogs_folder)
    # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
    #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
    # len(gene_id_name)

    HOGs_rhogs_xml_all = []
    for rhogid_num in  rhogid_num_list_input:
        HOGs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder)

        HOGs_rhogs_xml_all.append(HOGs_a_rhog_xml_all)

    print("here")


    for HOGs_a_rhog_xml in HOGs_a_rhog_xml_all:
        groups_xml.append(HOGs_a_rhog_xml)

    from xml.dom import minidom
    import xml.etree.ElementTree as ET
    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    print(xml_str)


    #client = Client(processes=False)  # start local workers as processes
    #future_1 = client.submit(infer_hogs_for_a_rhog, species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num, gene_trees_folder)
    #(dic_sub_hogs)= future_1.result()

    # futures = client.map(inc, range(1000))
    # as completed


    #print(dic_sub_hogs)

    print("**")

