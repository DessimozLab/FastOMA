# from xml.dom import minidom
# import xml.etree.ElementTree as ET
import _utils
import _inferhog

from _utils import logger_hog
import _utils_rhog
from os import listdir
import os
# from distributed import get_client
# from dask.distributed import rejoin, secede

if __name__ == '__main__':
    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/"
    gene_trees_folder = working_folder+"gene_tree_test2/" # working_folder + "gene_trees_/"
    # check gene_trees_folder exist otherwise mkdir this
    oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/omafast/archive/OmaServer.h5"

    address_rhogs_folder = working_folder + "rhog_all_v3_filbig/"  # old3/rhog_all/ /rhog_size_g2_s500/" sample_rootHOG
    species_tree_address = working_folder + "archive/lineage_tree_qfo.phyloxml"
    pickle_folder = working_folder + "pickle_folder_test2/"
    gene_id_pickle_file = working_folder + "gene_id_30aug_s500.pickle"
    # add warning when pickle folder is not empty
    output_xml_name = "out_t_aug2a_property_test.xml"

    if not os.path.exists(pickle_folder):
        os.mkdir(pickle_folder)
    if not os.path.exists(gene_trees_folder):
        os.mkdir(gene_trees_folder)

    # format_prot_name = 1  # 0:bird(TYTALB_R04643)  1:qfo(tr|E3JPS4|E3JPS4_PUCGT)
    file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)

    # step = "rhog"  # to infer roothogs when you have the proteome & hogmap.
    # step = "rhog"     # to infersubhogs when roothogs are ready.
    step = "hog"  # collect pickle file and write xml file

    print("we are here line25")
    if step == "rhog":
        """
        Structure of folders:
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address= the address to the oma databases
        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is orthoxml_to_newick.py HOG.
        """

        # working_folder+"omamer_database/oma_path/OmaServer.h5"

        print("rHOG inferece has started. The oma database address is in ", oma_database_address)
        (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
        (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
        query_prot_recs = _utils_rhog.add_species_name_gene_id(query_prot_recs,
                                                               query_species_names, gene_id_pickle_file)
        hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)

        (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_fscore_allspecies,
        prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

        query_prot_recs_filt = _utils_rhog.filter_prot_mapped(query_species_names,
                                                              query_prot_recs,
                                                              query_prot_names_species_mapped)

        print(len(query_prot_recs_filt), len(query_prot_recs_filt[0]))
        (rhogid_num_list, rhogids_prot_records_query) = _utils_rhog.group_prots_roothogs(prots_hogmap_hogid_allspecies,
                                                                                         address_rhogs_folder,
                                                                                         query_species_names,
                                                                                         query_prot_recs_filt)
        #step = "hog"

    if step == "hog":
        print("we are here line 60")
        rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder)
        logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list)) + ".")

        # rhogid_num_list =  [3339] #rhogid_num_list[:200]

        rhogid_num_list_raw = [ 834261 ] #  811161  # rhogid_num_list[:10] # [613860]  # , 618939, 615514, 834209 ]  #rhogid_num_list



        list_done_raw = listdir(pickle_folder)
        list_done = []
        for file in list_done_raw:
            numr = int(file.split(".")[0].split("_")[1])
            list_done.append(numr)

        rhogid_num_list = [i for i in rhogid_num_list_raw if i not in list_done]
        print("number of remained is ", len(rhogid_num_list))
        print(rhogid_num_list)
        a=2


        dask_level = 0  # 1:one level (rhog), 2:both levels (rhog+taxonomic)  3:only taxonomic level  0: no dask

        print(dask_level)
        if dask_level != 0:
            from _dask_env import client_dask

        print(rhogid_num_list[:7])
        number_roothog = len(rhogid_num_list)
        num_per_parralel = 1
        parralel_num = int(number_roothog/num_per_parralel)
        if number_roothog != parralel_num*num_per_parralel:
            parralel_num += 1
        rhogid_batch_list = []

        for list_idx in range(parralel_num):
            if list_idx == parralel_num:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:]
            else:
                rhogid_num_list_portion = rhogid_num_list[list_idx * num_per_parralel:(list_idx + 1) * num_per_parralel]
            rhogid_batch_list.append(rhogid_num_list_portion)

        if dask_level == 1 or dask_level == 2:
            dask_out_list = []
        dask_out_list = []
        hogs_rhogs_xml_all =[]
        for rhogid_batch_idx in range(len(rhogid_batch_list)):
            rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
            # logger_hog.info("\n *==* \nNumber of working root hog in the batchid:"+str(rhogid_batch_idx)+" is " +
            # str(len(rhogid_batch)) + ".")

            if dask_level == 1 or dask_level == 2:
                # vars_input_future = client_dask.scatter(vars_input)
                # client_dask = get_client()

                dask_out = client_dask.submit(_inferhog.read_infer_xml_rhogs_batch, rhogid_batch, file_folders, dask_level)
                dask_out_list.append(dask_out)

            else:
                hogs_rhog_xml_batch = _inferhog.read_infer_xml_rhogs_batch(rhogid_batch, file_folders, dask_level)
                # hogs_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
                hogs_rhogs_xml_all.extend(hogs_rhog_xml_batch)
                # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.
                print(hogs_rhogs_xml_all)

        if dask_level == 1 or dask_level == 2:
            for dask_out in dask_out_list:
                hogs_rhog_xml_batch = dask_out.result()
                hogs_rhogs_xml_all.extend(hogs_rhog_xml_batch)
            print(hogs_rhogs_xml_all)
            print("dask out gathered")

        step = "collect"

    if step == "collect":
        _inferhog.collect_write_xml(working_folder, pickle_folder, output_xml_name, gene_id_pickle_file)

    print("main py is finished s !.")

"""
to do : 
    - dobule check function merge_subhogs 
"""
