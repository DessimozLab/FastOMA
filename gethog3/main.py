# from xml.dom import minidom
# import xml.etree.ElementTree as ET
import _utils
import _inferhog

from _utils import logger_hog
import _utils_rhog
from os import listdir
import os

import _config
# from distributed import get_client
# from dask.distributed import rejoin, secede


"""

input species tree should have internal node name

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

    # step = "find_rhog"   #  to infer roothogs when you have the proteome & hogmap.
    # step = "find_subhog" # to infer subhogs when roothogs are ready.

    step = "find_subhog"  # find_subhog  find_rhog

    if step == "find_rhog":

        """
        Structure of folders in working_folder
        Put proteomes of species as fasta files in /omamer_search/proteome/
        Run omamer and put the output of omamer in /omamer_search/hogmap/
        oma_database_address = the address to the oma databases
        hog and HOG are used interchangeably here. 
        rHOG=rootHOG.  A subHOG itself is orthoxml_to_newick.py HOG.
        """

        if not os.path.exists(_config.working_folder):
            os.mkdir(_config.working_folder)

        # working_folder+"omamer_database/oma_path/OmaServer.h5"
        logger_hog.info("rHOG inferece has started. The oma database address is in "+_config.oma_database_address)
        (oma_db, list_oma_species) = _utils_rhog.parse_oma_db(_config.oma_database_address)
        (query_species_names, query_prot_recs) = _utils_rhog.parse_proteome(list_oma_species)
        query_prot_recs = _utils_rhog.add_species_name_gene_id(query_prot_recs, query_species_names)
        hogmap_allspecies_elements = _utils_rhog.parse_hogmap_omamer(query_species_names)

        (query_prot_names_species_mapped, prots_hogmap_hogid_allspecies, prots_hogmap_overlp_allspecies,
        prots_hogmap_fscore_allspecies, prots_hogmap_seqlen_allspecies, prots_hogmap_subfmedseqlen_allspecies) = hogmap_allspecies_elements

        query_prot_recs_filt = _utils_rhog.filter_prot_mapped(query_species_names, query_prot_recs,
                                                              query_prot_names_species_mapped)

        logger_hog.info("size of query_prot_recs_filt is "+str(len(query_prot_recs_filt))+" "+str(len(query_prot_recs_filt[0])))

        rhogids_list, rhogids_prot_records_query = _utils_rhog.group_prots_roothogs(prots_hogmap_hogid_allspecies, query_species_names, query_prot_recs_filt)
        # rhogid_num_list_raw=utils_rhog.write_rhog(rhogids_list,rhogids_prot_records_query, _config.working_folder+"rhogs/". "rhogs_raw",2)

        rhogids_list_filt, rhogids_prot_records_query_filt = _utils_rhog.filter_rhog(rhogids_list, rhogids_prot_records_query, prots_hogmap_fscore_allspecies, query_species_names,  query_prot_names_species_mapped)

        rhogid_num_list_filt1 = _utils_rhog.write_rhog(rhogids_list_filt, rhogids_prot_records_query_filt,
                                                      _config.working_folder+"rhogs/", 2)  # min_rhog_size, max_rhog_size


        #step = "find_subhog"

    if step == "find_subhog":

        pickle_folder = _config.working_folder + "pickles_rhog"
        if _config.gene_trees_write:
            gene_trees_folder = _config.working_folder + "genetrees/"
            if not os.path.exists(gene_trees_folder):
                os.mkdir(gene_trees_folder)
        if not os.path.exists(pickle_folder):
            os.mkdir(pickle_folder)

        address_rhogs_folder_filt = _config.working_folder + "rhogs/"
        rhogid_num_list_raw = _utils.list_rhog_fastas(address_rhogs_folder_filt)
        logger_hog.info("Number of root hogs is " + str(len(rhogid_num_list_raw)) + ".")

        print("***** pickle_folder ", pickle_folder)
        list_done_raw = listdir(pickle_folder)
        list_done_rhogid = []
        for file in list_done_raw:
            numr = int(file.split(".")[0].split("_")[1])
            list_done_rhogid.append(numr)

        # rhogid_num_list = rhogid_num_list_raw
        rhogid_num_list = [i for i in rhogid_num_list_raw if i not in list_done_rhogid]

        logger_hog.info("number of remained is " + str(len(rhogid_num_list)))

        rhogid_num_list = rhogid_num_list[:5]
        # print(rhogid_num_list[:4])
        logger_hog.info("working on a list with number of " + str(len(rhogid_num_list)))
        if not rhogid_num_list:
            exit()


        logger_hog.info("Dask level is "+str(_config.dask_level))
        if _config.dask_level != 0:
            from _dask_env import client_dask
            import dask.distributed
            dask.config.set({'distributed.scheduler.events-cleanup-delay': "10h"})
            dask.config.set({'distributed.logging.distributed': "debug"})
            dask.config.set({'distributed.logging.client': "debug"})
            dask.config.set({'distributed.logging.scheduler': "debug"})
            #print(dask.config.get("distributed"))

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

        if _config.dask_level == 1 or _config.dask_level == 2:
            dask_out_list = []
        dask_out_list = []
        hogs_rhogs_xml_all =[]
        for rhogid_batch_idx in range(len(rhogid_batch_list)):
            rhogid_batch = rhogid_batch_list[rhogid_batch_idx]
            if _config.dask_level == 1 or _config.dask_level == 2:
                # vars_input_future = client_dask.scatter(vars_input)
                # client_dask = get_client()
                dask_out = client_dask.submit(_inferhog.read_infer_xml_rhogs_batch, rhogid_batch)
                dask_out_list.append(dask_out)
            else:
                hogs_rhog_xml_batch = _inferhog.read_infer_xml_rhogs_batch(rhogid_batch)
                # hogs_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
                hogs_rhogs_xml_all.extend(hogs_rhog_xml_batch)

        if _config.dask_level == 1 or _config.dask_level == 2:
            logger_hog.info("start gathering dask")
            for dask_out in dask_out_list:
                hogs_rhog_xml_batch = dask_out.result()
                hogs_rhogs_xml_all.extend(hogs_rhog_xml_batch)

            logger_hog.info("Dask out gathered")
            client_dask.close()
            client_dask.shutdown()
            logger_hog.info("Client dask closed and shut down.")

        # step = "collect"

    if step == "collect":
        logger_hog.info("start writing xml")
        _inferhog.collect_write_xml(_config.working_folder)
        logger_hog.info("writing xml finished")

    logger_hog.info("main py is finished !.")

"""
to do : 


- as an argument differnt format of ete3 for reading species tree

run mani again
- double check 

    elif (dask_level == 2 or dask_level == 3): # 200


    if len(rhog_i) > 20000 and (dask_level == 2 or dask_level == 3): # 200
        # dask_future_taxon = True
        logger_hog.



to improve:


-  use  uniq internal node naming, and save as nwk, code below

    - dask=2: i don't need to write all hogs of different taxonomic level as  pickle
    - dask=0: or  rhog<200,  i don't need to write all hogs of different taxonomic level as  pickle
    
    - how decide  subtree 
    if len(species_leaves_names) <= 10:
    
    
    
from ete3 import Tree 
species_tree = Tree(nwk_path)
counter_internal = 0
for node in species_tree.traverse(strategy="postorder"):
    node_name = node.name
    num_leaves_no_name = 0
    if len(node_name) < 1:
        if node.is_leaf():
            node.name = "leaf_" + str(num_leaves_no_name)
        else:
            node.name = "internal_" + str(counter_internal)
            counter_internal += 1
# print("Working on the following species tree.")
# print(species_tree)
species_tree.write(format=1, outfile=nwk_path+"_edit.nwk")

    
"""

