
"""
import os
#import concurrent.futures


import gc
"""

from _dask_env import client_dask

from xml.dom import minidom
import xml.etree.ElementTree as ET

import _utils
import _inferhog
from _utils import logger_hog

if __name__ == '__main__':

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gene_trees_folder = working_folder + "/gene_trees_test/"
    address_rhogs_folder = working_folder + "/rhog_size_g2_s500/"  # "  old3/rhog_all/  /rhog_size_g2_s500/" sample_rootHOG
    species_tree_address = working_folder + "lineage_tree_qfo.phyloxml"
    pickle_address = working_folder + "/pickle_folder/"



    # format_prot_name = 0 # bird dataset   TYTALB_R04643
    format_prot_name = 1  # qfo dataset   # 'tr|E3JPS4|E3JPS4_PUCGT

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
        # (rhogid_num_list, rhogids_prot_records_query) = group_prots_roothogs(prots_hogmap_hogid_allspecies,
        # address_rhogs_folder)




    # format_prot_name = 0  # bird dataset   TYTALB_R04643
    format_prot_name = 1  # qfo dataset   # 'tr|E3JPS4|E3JPS4_PUCGT

    rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
    logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list)) + ".")
    print(rhogid_num_list[:2])

    rhogid_num_list = rhogid_num_list[:10]
    number_roothog = len(rhogid_num_list)
    num_per_parralel = 2
    parralel_num = int(number_roothog / num_per_parralel)
    rhogid_batch_list = []
    for list_idx in range(parralel_num + 1):
        if list_idx == parralel_num:
            rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
        else:
            rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
        rhogid_batch_list.append(rhogid_num_list_portion)

    dask_future = False

    # if dask_future:


    dask_out_list = []
    for rhogid_batch_idx in range(len(rhogid_batch_list)):
        rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
        # rhogid_num_list_input = rhogid_batch
        logger_hog.info("Number of working root hog is " + str(len(rhogid_batch)) + ".")
        (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_batch,
                                                                                        address_rhogs_folder,
                                                                                        format_prot_name,
                                                                                        rhogid_batch_idx)

        # # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
        # #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
        print("length of gene_id_name ", len(gene_id_name))

        if dask_future:
            # len_tresh = 1000
            # for rhogid_num_i in range(len(rhogid_num_list_input)):
            #    rhogid_num = rhogid_num_list_input[rhogid_num_i]
            #    rhogid_len = rhogid_len_list[rhogid_num_i]
            #    if rhogid_len < len_tresh:

            dask_future_taxon = True
            vars_input = (
                gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address,
                dask_future, dask_future_taxon, format_prot_name)

            # vars_input_future = client_dask.scatter(vars_input)
            vars_input_future = vars_input
            # read_infer_xml_rhogs(rhogid_batch_list, vars_input)
            dask_out = client_dask.submit(_inferhog.read_infer_xml_rhogs, rhogid_batch, vars_input_future)
            dask_out_list.append(dask_out)

            # dask_out = client_dask.submit(_inferhog.read_infer_xml_rhog, rhogid_num, vars_input)
            # dask_out_list.append(dask_out)
            print("*a*" * 5)
        #             else:
        #                 print("*b*" * 100)
        #                 dask_future_taxon = True  # second level of parralelizion
        #                 print(rhogid_num_i, rhogid_num, rhogid_len)
        #                 # dask_out = client_dask.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name,
        #                 #                          address_rhogs_folder, species_tree_address, gene_trees_folder,
        #                 #                          pickle_address, dask_future, dask_future_taxon)
        #                 # dask_out_list.append(dask_out)
        #                 print("here")
        #

        #
        else:

            dask_future_taxon = True

            vars_input = (
                gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address,
                dask_future, dask_future_taxon, format_prot_name)

            out = _inferhog.read_infer_xml_rhogs(rhogid_batch, vars_input)
    print("working ")
    if dask_future:
        for dask_out in dask_out_list:
            hogs_a_rhog_xml_all = dask_out.result()
            print(hogs_a_rhog_xml_all)

    # rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
    # logger_hog.info("Number of root hogs is "+str(len(rhogid_num_list))+".")
    # print(rhogid_num_list[:2])
    #
    # rhogid_num_list_input = rhogid_num_list[:100]
    # logger_hog.info("Number of working root hog is " + str(len(rhogid_num_list_input)) + ".")
    # (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list) = _utils.prepare_xml(rhogid_num_list_input, address_rhogs_folder, format_prot_name )
    # # # with open(address_working_folder + "/group_xml_ortho.pickle", 'rb') as handle:
    # # #     (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)
    # len(gene_id_name)
    #
    # dask_future = False
    #
    # if dask_future:
    #     # print("*** client **** ")
    #
    #     # print(cluster.dashboard_link)
    #     # print(cluster.get_logs())
    #     ncore = 1 # Total number of cores per job
    #     njobs = 2  # Cut the job up into this many processes.
    #     # # By default, process ~= sqrt(cores) so that the number of processes = the number of threads per process
    #     nproc = ncore
    #
    #     cluster = LocalCluster()
    #     # cluster = SLURMCluster(cores=ncore, processes=nproc, memory="20GB", walltime="01:00:00")
    #     cluster.scale(njobs)  # # ask for one jobs
    #     client_dask = Client(cluster)
    #
    #     # futures = client.map(score, x_values)
    #     # results = client.gather(futures)
    #     # hogs_a_rhog_xml_all = results
    #
    #     len_tresh = 1000
    #     dask_out_list =[]
    #     for rhogid_num_i in range(len(rhogid_num_list_input)):
    #         rhogid_num = rhogid_num_list_input[rhogid_num_i]
    #         rhogid_len = rhogid_len_list[rhogid_num_i]
    #         if rhogid_len < len_tresh:
    #
    #             dask_future_taxon = False
    #             vars_input = (gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future, dask_future_taxon, format_prot_name)
    #             vars_input_future = client.scatter(vars_input)
    #             dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, vars_input_future)
    #             dask_out_list.append(dask_out)
    #
    #             # dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, vars_input)
    #             # dask_out_list.append(dask_out)
    #             print("*a*" * 100)
    #         else:
    #             print("*b*" * 100)
    #             dask_future_taxon = True  # second level of parralelizion
    #             print(rhogid_num_i, rhogid_num, rhogid_len)
    #             # dask_out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name,
    #             #                          address_rhogs_folder, species_tree_address, gene_trees_folder,
    #             #                          pickle_address, dask_future, dask_future_taxon)
    #             # dask_out_list.append(dask_out)
    #             print("here")
    #
    #     # for dask_out in dask_out_list :
    #     #     hogs_a_rhog_xml_all = dask_out.result()
    #     #     print(hogs_a_rhog_xml_all)
    #
    # else:
    #     print("*d*" * 100)
    #     dask_future_taxon = False
    #     for rhogid_num_i in range(len(rhogid_num_list_input)):
    #         rhogid_num = rhogid_num_list_input[rhogid_num_i]
    #         rhogid_len = rhogid_len_list[rhogid_num_i]
    #         if rhogid_len > 50:
    #             vars_input = (
    #             gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future,
    #             dask_future_taxon, format_prot_name)
    #             hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, vars_input)
    #             exit
    #
    #
    # # for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
    # #     groups_xml.append(hogs_a_rhog_xml)
    # # xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # # print(xml_str)
    #
    #
    # print("test")

    # dask_working.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')
    # print("visualized.")
    # dask_result = dask_working.compute()
    # print(dask_result)  # prints "55"

    # dask_a = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address,
    #  gene_trees_folder)
    # print("this is dask_a: \n ",dask_a)
    # dask_a.visualize(filename='/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/out_4.svg')

    # dask.compute()

    # hogs_rhogs_xml_all = []
    # for rhogid_num in  rhogid_num_list_input:
    #     hogs_a_rhog_xml_all = _inferhog.read_infer_xml_rhog(rhogid_num, gene_id_name, address_rhogs_folder,
    #     species_tree_address, gene_trees_folder)
    #     hogs_rhogs_xml_all.append(hogs_a_rhog_xml_all)
    # print("here")

    #client = Client(processes=False)  # start local workers as processes
    # future_1 = client.submit(infer_hogs_for_a_rhog, species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    # rhogid_num, gene_trees_folder)
    # future = client.scatter(parameters)

    #dask.visualize(gene_id_name)



    #out = client.submit(_inferhog.read_infer_xml_rhog, rhogid_num, gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder)

    # (dic_sub_hogs)= future_1.result()

    # futures = client.map(inc, range(1000))
    # as completed
    # future.cancel()

    # rhogid_num_list_temp = [836500]  # rhogid_num_list[23]  # [833732]

    #     (dic_sub_hogs) = infer_hogs_for_a_rhog(species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
    #                                                        rhogid_num, gene_trees_folder)



    print("**")

"""
to do :
    input list of rhg num
    think how pickle per level ?
    think how to distribute rhog into list        
    # Your functions should not change the inputs directly.
    https://docs.dask.org/en/stable/delayed-best-practices.html#
    dic hog ??
    
    dobule check function merge_subhogs 
            
"""

