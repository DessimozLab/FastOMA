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


"""
Hard coded parameters

    max_num_prot = int(1e9)
    max_num_prot_per_sp = int(1e6) 


    _hog_class.py
    max_num_seq = 5
     
    _inferhog
    if (len(merged_msa) > 1000 and len(merged_msa[0]) > 3000) or (len(merged_msa) > 500 and len(merged_msa[0]) > 5000):
    tresh_ratio_gap_row = 0.1   # by 0.6 the whole row with few domains will
    tresh_ratio_gap_col = 0.2
    
    
    if len(rhog_i) > 20 and (dask_level == 2 or dask_level == 3): 

"""




if __name__ == '__main__':
    oma_database_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/omafast/archive/OmaServer.h5"

    working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/" #fastget/qfo2/"
    gene_id_pickle_file = working_folder + "gene_id_v2_bird.pickle"
    species_tree_address = working_folder + "concatanted_363.fasta.contree_edited.nwk"

    omamer_fscore_treshold_big_rhog = 0.5  # 0.2
    treshold_big_rhog_szie = 3000

    name = str(omamer_fscore_treshold_big_rhog)+"_"+str(treshold_big_rhog_szie)

    address_rhogs_folder_raw = working_folder + "rhogs_v2_raw/"
    address_rhogs_folder_filt = working_folder + "rhogs_v2_" + name + "/"
    pickle_folder = working_folder + "pickle_"+name+"/"
    gene_trees_folder = "no_write_tree_no"  #  working_folder+"genetree_"+name+"/"
    output_xml_name = "out_xml__"+name+"_.xml"


    # format_prot_name = 1  # 0:bird(TYTALB_R04643)  1:qfo(tr|E3JPS4|E3JPS4_PUCGT)
    file_folders = (address_rhogs_folder_filt, gene_trees_folder, pickle_folder, species_tree_address)

    # step = "find_rhog"  # to infer roothogs when you have the proteome & hogmap.
    # step = "find_subhog"     # to infer subhogs when roothogs are ready.

    step = "find_subhog"
    # find_subhog  find_rhog
    # collect pickle file and write xml file

    if step == "find_rhog":

        """
        Structure of folders:
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address= the address to the oma databases
        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is orthoxml_to_newick.py HOG.
        """

        # working_folder+"omamer_database/oma_path/OmaServer.h5"
        logger_hog.info("rHOG inferece has started. The oma database address is in "+oma_database_address)
        (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(oma_database_address)
        (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species, working_folder)
        query_prot_recs = _utils_rhog.add_species_name_gene_id(query_prot_recs,
                                                               query_species_names, gene_id_pickle_file)
        hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names, working_folder)

        # (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_fscore_allspecies,
        # prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

        (prots_hogmap_name_allspecies, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
        prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

        query_prot_names_species_mapped = prots_hogmap_name_allspecies # double check ?

        query_prot_recs_filt = _utils_rhog.filter_prot_mapped(query_species_names,
                                                              query_prot_recs,
                                                              query_prot_names_species_mapped)

        logger_hog.info("size of query_prot_recs_filt is "+str(len(query_prot_recs_filt))+" "+str(len(query_prot_recs_filt[0])))


        rhogids_list, rhogids_prot_records_query = _utils_rhog.group_prots_roothogs(prots_hogmap_hogid_allspecies, query_species_names, query_prot_recs_filt)
        rhogid_num_list_raw = _utils_rhog.write_rhog(rhogids_list, rhogids_prot_records_query, address_rhogs_folder_raw, 2)  # min_rhog_size=1, max_rhog_size=1e100

        rhogids_list_filt, rhogids_prot_records_query_filt = _utils_rhog.filter_rhog(rhogids_list, rhogids_prot_records_query, prots_hogmap_fscore_allspecies, query_species_names,  prots_hogmap_name_allspecies, omamer_fscore_treshold_big_rhog, treshold_big_rhog_szie)

        rhogid_num_list_filt = _utils_rhog.write_rhog(rhogids_list_filt, rhogids_prot_records_query_filt, address_rhogs_folder_filt, 2)  # min_rhog_size=1, max_rhog_size=1e100


       #step = "find_subhog"


    if step == "find_subhog":

        if not os.path.exists(gene_trees_folder):
            os.mkdir(gene_trees_folder)
        if not os.path.exists(pickle_folder):
            os.mkdir(pickle_folder)

        rhogid_num_list = _utils.list_rhog_fastas(address_rhogs_folder_filt)
        logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list)) + ".")

        rhogid_num_list_raw = rhogid_num_list[:10]  # 605945 # 560403 [570080] #
        # rhog_num_input = sys.argv[1]; rhogid_num_list_raw = [int(rhog_num_input)]

        # small size [614128, 599704,839732, 581211, 594354, 606190, 581722]
        # 613986 337 prots
        # 0589674 56 prots
        # [606409, 575384, 834730, 606033, 618436, 620754, 614327, 613986  ] #    # rhogid_num_list[:10] # [613860]  # , 618939, 615514, 834209 ]  #rhogid_num_list
        # itermediate size 834261 614102
        # very big 811161 811184
        print("***** pickle_folder ", pickle_folder)
        list_done_raw = listdir(pickle_folder)
        list_done = []
        for file in list_done_raw:
            numr = int(file.split(".")[0].split("_")[1])
            list_done.append(numr)

        #rhogid_num_list = [i for i in rhogid_num_list_raw if i not in list_done]
        rhogid_num_list = rhogid_num_list_raw[:10]
        logger_hog.info("number of remained is " + str(len(rhogid_num_list)))
        if not rhogid_num_list:
            exit


        dask_level = 2   # 1:one level (rhog), 2:both levels (rhog+taxonomic)  3:only taxonomic level  0: no dask

        logger_hog.info("Dask level is "+str(dask_level))
        if dask_level != 0:
            from _dask_env import client_dask
            import dask.distributed
            # export DASK_DISTRIBUTED__SCHEDULER__EVENTS_CLEANUP_DELAY=10h
            print(dask.config.get("distributed.scheduler"))
            dask.config.set({'distributed.scheduler.events-cleanup-delay': "10h"})
            print(dask.config.get("distributed.scheduler"))



        logger_hog.info("Few of rhog num ids are "+" ".join([str(i) for i in  rhogid_num_list[:7]]))
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




        if dask_level == 1 or dask_level == 2:
            logger_hog.info("start gathering dask")
            for dask_out in dask_out_list:
                hogs_rhog_xml_batch = dask_out.result()
                hogs_rhogs_xml_all.extend(hogs_rhog_xml_batch)



            logger_hog.info("dask out gathered")
            client_dask.close()
            client_dask.shutdown()
            logger_hog.info("client dask closed and shut down.")

        # step = "collect"

    if step == "collect":
        logger_hog.info("start writing xml")
        _inferhog.collect_write_xml(working_folder, pickle_folder, output_xml_name, gene_id_pickle_file)
        logger_hog.info("writing xml finished")

    logger_hog.info("main py is finished !.")

"""
to do : 
    -  
"""

