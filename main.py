
"""
import os
from xml.dom import minidom
#import concurrent.futures

from dask_jobqueue import SLURMCluster
import gc
"""

import dask
from dask.distributed import Client
from dask.distributed import LocalCluster

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
        # (oma_db, list_oma_species) = parse_oma_db(oma_database_address)
        # (query_species_names, query_prot_records_species) = parse_proteome(list_oma_species)
        # query_prot_records_species = add_species_name(query_prot_records_species,query_species_names)
        # hogmap_allspecies_elements = parse_hogmap_omamer(query_species_names)
        # (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_subfscore_allspecies,
        # prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements
        # query_prot_records_species_filtered =  filter_prot_mapped(query_species_names, query_prot_records_species,
        # query_prot_names_species_mapped)
        # print(len(query_prot_records_species_filtered),len(query_prot_records_species_filtered[0]))
        # (rhogid_num_list, rhogids_prot_records_query) = group_prots_rootHOGs(prots_hogmap_hogid_allspecies,
        # address_rhogs_folder)

    rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
    logger_hog.info("Number of root hog is "+str(len(rhogid_num_list))+".")
    print(rhogid_num_list[:2])

    rhogid_num_list_input = rhogid_num_list[6:9]
    (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_num_list_input, address_rhogs_folder)
    # # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
    # #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
    # # len(gene_id_name)


    ###  Running dask

    print("*** client **** ")
    cluster = LocalCluster()
    client = Client(cluster)
    print(cluster.dashboard_link)
    print(cluster.get_logs())

    len_tresh = 100

    rhogid_num = rhogid_num_list_input[0]
    for rhogid_num_i in range(len(rhogid_num_list_input)):
        rhogid_num = rhogid_num_list_input[rhogid_num_i]
        rhogid_len = rhogid_len_list[rhogid_num_i]

        if rhogid_len  < len_tresh:

            dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name,
                                     address_rhogs_folder, species_tree_address, gene_trees_folder)
        else:








    #     out = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address,
    #                                   gene_trees_folder)
    #     print("done", out)


    # dask_working.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')
    # print("visualized.")
    # dask_result = dask_working.compute()
    # print(dask_result)  # prints "55"

    # dask_a = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address,
    #  gene_trees_folder)
    # print("this is dask_a: \n ",dask_a)
    # dask_a.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')

    # dask.compute()

    # HOGs_rhogs_xml_all = []
    # for rhogid_num in  rhogid_num_list_input:
    #     HOGs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder,
    #     species_tree_address, gene_trees_folder)
    #     HOGs_rhogs_xml_all.append(HOGs_a_rhog_xml_all)
    # print("here")
    # for HOGs_a_rhog_xml in HOGs_a_rhog_xml_all:
    #     groups_xml.append(HOGs_a_rhog_xml)
    # from xml.dom import minidom
    # import xml.etree.ElementTree as ET
    # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str)

    #client = Client(processes=False)  # start local workers as processes
    # future_1 = client.submit(infer_hogs_for_a_rhog, species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    # rhogid_num, gene_trees_folder)
    # future = client.scatter(parameters)

    dask.visualize(gene_id_name)



    #out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder)

    # (dic_sub_hogs)= future_1.result()

    # futures = client.map(inc, range(1000))
    # as completed
    # future.cancel()

    # rhogid_num_list_temp = [836500]  # rhogid_num_list[23]  # [833732]

    #     (dic_sub_hogs) = infer_hogs_for_a_rhog(species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    #                                                        rhogid_num, gene_trees_folder)




    """
        to do :
                input list of rhg num
                think how pickle per level ?
                think how to distribute rhog into list        
                # Your functions should not change the inputs directly.
                https://docs.dask.org/en/stable/delayed-best-practices.html#
                dic hog ??
    """
    print("**")
