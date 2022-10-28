
from ete3 import Tree
from Bio import SeqIO
# import dill as dill_pickle
import pickle
import gc
from distributed import get_client
import os

from os import listdir
import xml.etree.ElementTree as ET
from xml.dom import minidom
from Bio.Align import MultipleSeqAlignment

import _wrappers
import _utils
from _hog_class import HOG
from _utils import logger_hog
# from _dask_env import client_dask
# from dask.distributed import as_completed

from dask.distributed import rejoin, secede

import networkx as nx
import matplotlib.pyplot as plt


def read_infer_xml_rhogs_batch(rhogid_batch_list, file_folders, dask_level):
    # file_folders = (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address)

    hogs_rhog_xml_batch = []
    #print("There are "+str(len(rhogid_batch_list))+" rhogs in the batch.")
    for rhogid_num in rhogid_batch_list:
        hogs_rhogs_xml = read_infer_xml_rhog(rhogid_num, file_folders, dask_level)  # orthoxml_to_newick.py list of hog object
        hogs_rhog_xml_batch.extend(hogs_rhogs_xml)

    return hogs_rhog_xml_batch  # orthoxml_to_newick.py list of hog object


def read_infer_xml_rhog(rhogid_num, file_folders, dask_level):
    (address_rhogs_folder, gene_trees_folder, pickle_folder, species_tree_address) = file_folders
    hogs_children_level_pickle_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/gethog3_27oct/pickle_hog_children2/"
    hogs_children_level_pickle_folder_rhog = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num)

    if not os.path.exists(hogs_children_level_pickle_folder):
        os.mkdir(hogs_children_level_pickle_folder)
    if not os.path.exists(hogs_children_level_pickle_folder_rhog):
        os.mkdir(hogs_children_level_pickle_folder_rhog)


    logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
    prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
    rhog_i = list(SeqIO.parse(prot_address, "fasta"))
    logger_hog.debug("number of proteins in the rHOG is "+str(len(rhog_i))+".")
    (species_tree) = _utils.read_species_tree(species_tree_address)
    (species_tree, species_names_rhog, prot_names_rhog) = _utils.prepare_species_tree(rhog_i, species_tree, rhogid_num)
    species_names_rhog = list(set(species_names_rhog))
    logger_hog.debug("The number of unique species in the rHOG " + str(rhogid_num) + "is " + str(len(species_names_rhog)) + ".")
    # species_tree.write();  print(species_tree.write())

    # hogs_a_rhog = infer_hogs_for_rhog_levels_future(species_tree, recursive_input)

    recursive_4inputs = (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder)

    if len(rhog_i) > 1 and (dask_level == 2 or dask_level == 3): # 200
        # dask_future_taxon = True
        logger_hog.debug("Dask future taxon is on for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))

        client_dask_working = get_client()
        secede()
        # #recursive_input = (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder)
        hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_future_v2, species_tree,recursive_4inputs)

        #hogs_a_rhog_1 = infer_hogs_for_rhog_future_v2(species_tree, recursive_4inputs)
        # hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, species_tree, recursive_4inputs)
        #
        # #rhog_i_future = client_dask_working.scatter(rhog_i)
        # #recursive_input_future = (rhog_i_future, species_names_rhog, rhogid_num, gene_trees_folder)
        # #hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, species_tree, recursive_input_future)
        #
        hogs_a_rhog = hogs_a_rhog_future.result()
        # rejoin()




    else:
        # recursive_input = (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder)
        # ?? we can have recursive_input includign rhog_i  for  small rhogs
        # dask_future_taxon = False
        logger_hog.debug("Dask future taxon is off for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))

        hogs_a_rhog_1 = infer_hogs_for_rhog_levels_recursively(species_tree, recursive_4inputs)
        # hogs_a_rhog_1  is len

    root_node_name= species_tree.name
    hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_"+str(rhogid_num) + "/_"+ str(root_node_name)
    with open(hogs_children_level_pickle_file+".pickle", 'rb') as handle:
        hogs_a_rhog = pickle.load(handle)
    # ?? logger_hog.debug("subhogs in this level are "+' '.join(["[" + str(i) + "]" for i in ?? ])+".") # hogs_a_rhog
    keep_intermediate_files = True
    if not keep_intermediate_files:
     os.rmdir(hogs_children_level_pickle_folder_rhog)

    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        # print(hog_i)
        if len(hog_i._members) > 1:
            # could be improved # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_a_rhog_xml = hog_i.to_orthoxml()
            hogs_rhogs_xml.append(hogs_a_rhog_xml)
    logger_hog.debug("we are not reporting single tone hogs in the output xml.")
    # how to handle empty hogs !? why happening and if not save pickle, file not exist error from rhog id list
    with open(pickle_folder + '/file_' + str(rhogid_num) + '.pickle', 'wb') as handle:
        #dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("***** hogs are written as orthoxml_to_newick.py pickle " + pickle_folder + '/file_' + str(rhogid_num) + '.pickle')

    del hogs_a_rhog
    gc.collect()

    return hogs_rhogs_xml

#############

def infer_hogs_for_rhog_future_v2(sub_species_tree, recursive_4inputs):

    out_1 = 0
    species_leaves_names = [i.name for i in sub_species_tree.get_leaves()]
    if len(species_leaves_names) <= 10:
        out_1 = infer_hogs_for_rhog_levels_recursively_v2(sub_species_tree, recursive_4inputs)

    for node_species_tree in sub_species_tree.traverse(strategy="postorder"):
        species_leaves_names = [i.name for i in node_species_tree.get_leaves()]

        if len(species_leaves_names) > 10:
            # out_1 = infer_hogs_for_rhog_levels_recursively_v2(node_species_tree, recursive_4inputs)
            client_dask_working = get_client()
            secede()
            hogs_level_list_futures = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_v2, node_species_tree, recursive_4inputs)
            out_1 = client_dask_working.gather(hogs_level_list_futures)
            rejoin()

            for node in node_species_tree.traverse():
                node.processed = True


    return out_1



def infer_hogs_for_rhog_NO_future_v2(sub_species_tree, recursive_4inputs):
    out_1 = 0
    species_leaves_names = [i.name for i in sub_species_tree.get_leaves()]
    if len(species_leaves_names) <= 10:
        out_1 = infer_hogs_for_rhog_levels_recursively_v2(sub_species_tree, recursive_4inputs)

    for node_species_tree in sub_species_tree.traverse(strategy="postorder"):
        species_leaves_names = [i.name for i in node_species_tree.get_leaves()]

        if len(species_leaves_names) > 10:
            out_1 = infer_hogs_for_rhog_levels_recursively_v2(node_species_tree, recursive_4inputs)
            for node in node_species_tree.traverse():
                node.processed = True
    return out_1




# def infer_hogs_for_rhog_levels_recursively_future(sub_species_tree, recursive_4inputs):
#     #(rhog_i, species_names_rhog, rhogid_num, gene_trees_folder) = recursive_input
#     #logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
#
#     if sub_species_tree.is_leaf():
#
#         (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder) = recursive_4inputs
#         singletone_hog_out = singletone_hog_(sub_species_tree, species_names_rhog, rhogid_num, address_rhogs_folder)
#         return singletone_hog_out
#
#     children_nodes = sub_species_tree.children
#
#     client_dask_working = get_client()
#     secede()
#     hogs_children_level_list_futures = [client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, child, recursive_4inputs) for child in children_nodes ]
#
#     hogs_children_level_list_futures = client_dask_working.gather(hogs_children_level_list_futures)
#     rejoin()
#     # hogs_children_level_list = hogs_children_level_list_futures
#     # hogs_children_level_list = []
#     # for future in hogs_children_level_list_futures:
#     #    hogs_children_level_list.extend(future.result())
#
#
#     # if hogs_children_level_list_futures:
#     #     if isinstance(hogs_children_level_list_futures[0], list):
#     #         hogs_children_level_list_flatten = []
#     #         for hog_ in hogs_children_level_list_futures:
#     #             # for hog in hogs_list:
#     #             hogs_children_level_list_flatten.extend(hog_)
#
#     # hogs_children_level_list = hogs_children_level_list_flatten
#
#
#     infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, recursive_4inputs) # hogs_children_level_list
#
#     return infer_hogs_this_level_out

def infer_hogs_for_rhog_levels_recursively(sub_species_tree, recursive_4inputs):

    if sub_species_tree.is_leaf():
        (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder) = recursive_4inputs
        #hogs_this_level_list = singletone_hog(sub_species_tree, rhog_i, species_names_rhog, rhogid_num)
        singletone_hog_out = singletone_hog_(sub_species_tree, species_names_rhog, rhogid_num, address_rhogs_folder)
        # out 1 succesful
        return singletone_hog_out

    children_nodes = sub_species_tree.children

    hogs_children_level_list = []
    for node_species_tree_child in children_nodes:
        hogs_children_level_list_i = infer_hogs_for_rhog_levels_recursively(node_species_tree_child, recursive_4inputs)
        # hogs_children_level_list_i should be 1
        #hogs_children_level_list.extend(hogs_children_level_list_i)
    infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, recursive_4inputs) # ,hogs_children_level_list
    # hogs_this_level_list should be one

    return infer_hogs_this_level_out




def infer_hogs_for_rhog_levels_recursively_v2(sub_species_tree, recursive_4inputs):

    infer_hogs_this_level_out= 0
    if sub_species_tree.is_leaf():
        if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
            (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder) = recursive_4inputs
            infer_hogs_this_level_out = singletone_hog_(sub_species_tree, species_names_rhog, rhogid_num, address_rhogs_folder)
        return infer_hogs_this_level_out

    children_nodes = sub_species_tree.children
    for node_species_tree_child in children_nodes:
        hogs_children_level_list_i = infer_hogs_for_rhog_levels_recursively_v2(node_species_tree_child, recursive_4inputs)
    infer_hogs_this_level_out = 0
    if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
        infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, recursive_4inputs)  # output is an integer, which is the length

    return infer_hogs_this_level_out



# def singletone_hog(node_species_tree, rhog_i, species_names_rhog, rhogid_num):
#
#     node_species_name = node_species_tree.name  # there is only one species (for the one protein)
#     prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog)) if
#                                  species_names_rhog[idx] == node_species_name]
#     rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]
#
#     hogs_this_level_list = []
#     for prot in rhog_part:
#         hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
#         hogs_this_level_list.append(hog_leaf)
#     return hogs_this_level_list


def singletone_hog_(node_species_tree, species_names_rhog, rhogid_num, address_rhogs_folder):

    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    this_level_node_name = node_species_name
    logger_hog.debug("* reading prot address  " + str(this_level_node_name))

    prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
    # logger_hog.debug("* reading prot address  " + prot_address)
    rhog_i = list(SeqIO.parse(prot_address, "fasta"))

    prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog)) if
                                 species_names_rhog[idx] == node_species_name]
    rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]

    hogs_this_level_list = []
    for prot in rhog_part:
        hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
        hogs_this_level_list.append(hog_leaf)


    hogs_children_level_pickle_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/gethog3_27oct/pickle_hog_children2/"
    hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num) + "/_" + str(this_level_node_name)
    with open(hogs_children_level_pickle_file+".pickle", 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs was written as pickle file.")

    return len(hogs_this_level_list) # hogs_this_level_list _


def infer_hogs_this_level(sub_species_tree, recursive_4inputs):  # hogs_children_level_list
    (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder) = recursive_4inputs
    hogs_children_level_pickle_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/bird_hog/gethog3_27oct/pickle_hog_children2/"

    node_species_tree = sub_species_tree
    this_level_node_name = node_species_tree.name
    if node_species_tree.is_leaf():
        print(" issue 1235, *** * * * ** * * * ",rhogid_num )
        exit
        # assert hogs_children_level_list == []
        #hogs_this_level_list = singletone_hog_(node_species_tree, species_names_rhog, rhogid_num)
        # we shouldnt be here ???
        #
        # child_name = node_species_tree.name
        # hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num) + "/_" + str(child_name)
        # with open(hogs_children_level_pickle_file+".pickle", 'wb') as handle:
        #     pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return 1  # [child_name] # hogs_this_level_list

    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list =[]
    for child_name in children_name:
        hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_"+str(rhogid_num) + "/_"+ str(child_name)
        with open(hogs_children_level_pickle_file+".pickle", 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))
     # check file doenst exist ? how to handle, if not


    #if len(node_species_tree.name.split("_")) > 1:
    #if "internal" in node_species_tree.name:
    logger_hog.debug("Finding hogs for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name)+"\n"+str(node_species_tree.write())+"\n")



    if len(hogs_children_level_list) == 1:

        assert len(children_name) == 1

        hogs_this_level_list = hogs_children_level_list
        hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num) + "/_" + str(this_level_node_name)
        with open(hogs_children_level_pickle_file+".pickle", 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return len(hogs_children_level_list)

    # print("*7*", len(hogs_children_level_list), hogs_children_level_list)
    # if len(hogs_children_level_list)>0:
    #    print("**",hogs_children_level_list[0], len(hogs_children_level_list))

    # hogs_children_level_list
    # if isinstance(hogs_children_level_list[0], list):
    #     hogs_children_level_list_flatten = []
    #     for hogs_list in hogs_children_level_list:
    #         for hog in hogs_list:
    #         hogs_children_level_list_flatten.extend(hog)
    #
    # hogs_children_level_list = hogs_children_level_list_flatten


    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list]
    gene_tree_file_addr = gene_trees_folder + "/tree_" + str(rhogid_num) + "_" + str(
        node_species_tree.name) + ".nwk"

    if len(gene_tree_file_addr) > 245:
        import random
        rand_num = random.randint(1, 10000)
        gene_tree_file_addr = gene_tree_file_addr[:245] + str(rand_num)+".nwk"
    logger_hog.debug("Merging "+str(len(sub_msa_list_lowerLevel_ready))+" MSAs for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name))
    merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready, gene_tree_file_addr)

    if merged_msa:
        logger_hog.debug("All sub-hogs are merged, merged msa is with length of " + str(len(merged_msa)) + " " + str(
        len(merged_msa[0])) + " for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
                node_species_tree.name))

        # merged_msa_filt = merged_msa
        # 893*4839, 10 mins
        tresh_ratio_gap_row = 0.4
        tresh_ratio_gap_col = 0.2
        min_cols_msa_to_filter = 3000      # used for msa before gene tree inference and  saving msa in hog class

        if len(merged_msa[0]) >= min_cols_msa_to_filter:
            # (len(merged_msa) > 10000 and len(merged_msa[0]) > 3000) or (len(merged_msa) > 500 and len(merged_msa[0]) > 5000) or (len(merged_msa) > 200 and len(merged_msa[0]) > 9000):
            # for very big MSA, gene tree is slow. if it is full of gaps, let's trim the msa.
            logger_hog.debug("We are doing MSA trimming "+str(rhogid_num)+", for taxonomic level:"+str(node_species_tree.name))

            #print(len(merged_msa), len(merged_msa[0]))
            msa_filt_col = _utils.msa_filter_col(merged_msa, tresh_ratio_gap_col, gene_tree_file_addr)
            #print(len(msa_filt_col), len(msa_filt_col[0]))
            msa_filt_row_col = _utils.msa_filter_row(msa_filt_col, tresh_ratio_gap_row, gene_tree_file_addr)
            #print(len(msa_filt_row_col), len(msa_filt_row_col[0]))
            merged_msa_filt = msa_filt_row_col
        else:
            msa_filt_row_col = merged_msa
            msa_filt_col = merged_msa

            # the msa may be empty
        if len(msa_filt_row_col) < 2:
            msa_filt_row_col = msa_filt_col[:2]
    else:
        logger_hog.info("Issue 1455, merged_msa is empty " + str(rhogid_num) + ", for taxonomic level:" + str(node_species_tree.name))





    gene_tree_raw = _wrappers.infer_gene_tree(msa_filt_row_col, gene_tree_file_addr)
    gene_tree = Tree(gene_tree_raw + ";", format=0)
    logger_hog.debug("Gene tree is inferred with length of " + str(len(gene_tree)) + " for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name))
    R_outgroup = gene_tree.get_midpoint_outgroup()
    gene_tree.set_outgroup(R_outgroup)  # print("Midpoint rooting is done for gene tree.")
    gene_tree = _utils.lable_sd_internal_nodes(gene_tree)
    # print("Overlap speciation is done for internal nodes of gene tree, as following:")
    # print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")
    logger_hog.debug("Merging sub-hogs of children started  for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name))

    # merging two filtered msa

    # we may use this merged_msa here
    hogs_this_level_list = merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, msa_filt_col)


    logger_hog.debug("Hogs of this level is found for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(this_level_node_name))



    # check for conflicts in merging
    #     for i in range(subHOG_to_be_merged_set_other_Snodes):  if
    #         for i in range(subHOG_to_be_merged_set_other_Snodes):  print("*&*& ",node_species_tree.name)
    # # dvelopmnet mode  logger info
    # prot_list_sbuhog = [i._members for i in hogs_this_level_list]
    # prot_list_sbuhog_short = []
    # for prot_sub_list_sbuhog in prot_list_sbuhog:
    #     if format_prot_name == 0:  # bird dataset TYTALB_R04643
    #         prot_list_sbuhog_short = prot_sub_list_sbuhog
    #     elif format_prot_name == 1:  # qfo dataset  'tr|E3JPS4|E3JPS4_PUCGT
    #         prot_list_sbuhog_short.append([prot.split('|')[2] for prot in prot_sub_list_sbuhog])
    # logger_hog.debug(str(len(hogs_this_level_list))+" hogs are inferred at the level "+node_species_tree.name+": "+' '.join(
    #     [str(i) for i in prot_list_sbuhog_short]))

    hogs_children_level_pickle_file = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num) + "/_" + str(this_level_node_name)
    with open(hogs_children_level_pickle_file + ".pickle", 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return len(hogs_children_level_list)



def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa):

    """
    this function could be improved
    """

    # hogs_children_level_list_flatten = []
    # for hogs_list in hogs_children_level_list:
    #     hogs_children_level_list_flatten += hogs_list
    # hogs_children_level_list = hogs_children_level_list_flatten

    subhogs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
    hogs_this_level_list = []
    subHOG_to_be_merged_set_other_Snodes = []
    subHOG_to_be_merged_set_other_Snodes_flattned_temp = []

    # hoggraph_node_name = [i._hogid.split("_")[1][3:] for i in hogs_children_level_list]
    # hog_size_dic = {}
    # dic_hog = {}
    # for hog in hogs_children_level_list:
    #     hog_id_short = hog._hogid.split("_")[1][3:]
    #     for prot in hog._members:
    #         dic_hog[prot] = hog_id_short + "_len" + str(len(hog._members))
    #     hog_size_dic[hog_id_short] = len(hog._members)
    # for hog in hogs_children_level_list:
    #     hog_id_short = hog._hogid.split("_")[1][3:] + "_len" + str(len(hog._members))
    #     # print(hog_id_short, ":", hog._members)

    # hoggraph_node_name_len = [i+"_len"+str(hog_size_dic[i]) for i in hoggraph_node_name]
    # hoggraph = nx.Graph()
    # hoggraph.add_nodes_from(hoggraph_node_name_len)

    for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n,"processed") and n.processed==True):
        # print("Leaves assigned to hog are ", assigned_leaves_to_hog)   #print("Traversing gene tree. Now at node", node.name)
        if node.is_leaf():
            continue
        node_leaves_name = [i.name for i in node.get_leaves()]   # if node.name[0] == "D":  print(2)
        # s_gene_tree_leaves = [i.name for i in node.get_leaves()]
        #s_gene_tree_leaves_update = [i.split("_")[1][3:] for i in s_gene_tree_leaves]
        if node.name[0] == "S":  # this is orthoxml_to_newick.py sub-hog.
            # num_prot = len(s_gene_tree_leaves)
            # for i in range(num_prot):
            #     hog_i = dic_hog[s_gene_tree_leaves[i]]
            #     for j in range(i):
            #         hog_j = dic_hog[s_gene_tree_leaves[j]]
            #         if hoggraph.has_edge(hog_i, hog_j):
            #             hoggraph[hog_i][hog_j]['weight'] += 1
            #         else:
            #             hoggraph.add_edge(hog_i, hog_j, weight=1)

            # print(hoggraph.edges(data=True))

            subHOG_to_be_merged = []
            for node_leave_name in node_leaves_name:  # print(node_leave_name)
                for subHOG in hogs_children_level_list:
                    subHOG_members = subHOG._members
                    if node_leave_name in subHOG_members:  # could be improved
                        if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp:
                            subHOG_to_be_merged.append(subHOG)
                            subhogs_id_children_assigned.append(subHOG._hogid)
                        else:  # this hog is already decided to be merged  #print("issue 184", node.name, subHOG._hogid, node_leave_name)
                            if "processed" in node:
                                print("issue 1863", node.name, subHOG._hogid, node_leave_name) #print("processed", node.name) #else: #    print("processed not in ", node.name)  # print(node_leave_name,"is in ",subHOG._hogid)
            if subHOG_to_be_merged:
                subHOG_to_be_merged_set = set(subHOG_to_be_merged)
                taxnomic_range = node_species_tree.name
                HOG_this_node = HOG(subHOG_to_be_merged_set, taxnomic_range, rhogid_num, msa=merged_msa)
                hogs_this_level_list.append(HOG_this_node)
                subHOG_to_be_merged_set_other_Snodes.append([i._hogid for i in subHOG_to_be_merged_set])
                subHOG_to_be_merged_set_other_Snodes_flattned_temp = [item for items in
                                                                      subHOG_to_be_merged_set_other_Snodes for
                                                                      item in items]
                #  I don't need to traverse deeper in this clade
            node.processed = True  # print("?*?*  ", node.name)


        subHOG_to_be_merged_set_other_Snodes_flattned = [item for items in subHOG_to_be_merged_set_other_Snodes for
                                                         item in items]
        if [i._hogid for i in hogs_children_level_list] == subHOG_to_be_merged_set_other_Snodes_flattned:
            break
        # print("node name ", node.name)

    # fig = plt.figure(figsize=(300, 200), dpi=60)
    # pos = nx.spring_layout(hoggraph, k=0.25, iterations=30)  # For better example looking  # smaller k, biger space between
    # nx.draw(hoggraph, pos, with_labels=True, node_color='y', node_size=500, font_size=16) # , alpha=0.4
    # # nx.draw(G, pos,, edge_color="r", font_size=16, with_labels=True)
    # labels = {e: hoggraph.edges[e]['weight'] for e in hoggraph.edges}
    # nx.draw_networkx_edge_labels(hoggraph, pos, edge_labels=labels, font_size=16)
    # import random
    # num = random.randint(3, 1000000)
    # plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/hoggraph/" + hogs_children_level_list[0]._hogid[4:] + "file_rndm"+str(num)+".jpg")
    # # plt.show()


    for subHOG in hogs_children_level_list:  # for the single branch  ( D include orthoxml_to_newick.py  subhog and orthoxml_to_newick.py S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)

    # this could be improved
    # we expect to see orthoxml_to_newick.py list of list as ooutput
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    return hogs_this_level_list


def collect_write_xml(working_folder, pickle_folder, output_xml_name, gene_id_pickle_file):

    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        #gene_id_name = dill_pickle.load(handle)
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)

    for query_species_name, list_prots in gene_id_name.items():

        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")

        for (gene_idx_integer, query_prot_name) in list_prots:
            query_prot_name_pure1 = query_prot_name.split("||")[0].strip()
            if "|" in query_prot_name_pure1:
                query_prot_name_pure = query_prot_name_pure1.split("|")[1]
            else:
                query_prot_name_pure = query_prot_name
            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})

    pickle_files_adress = listdir(pickle_folder)

    hogs_a_rhog_xml_all = []
    for pickle_file_adress in pickle_files_adress:
        with open(pickle_folder + pickle_file_adress, 'rb') as handle:
            # hogs_a_rhog_xml_batch = dill_pickle.load(handle)
            hogs_a_rhog_xml_batch = pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
            # hogs_rhogs_xml_all is orthoxml_to_newick.py list of hog object.

    print("number of hogs in all batches is ", len(hogs_a_rhog_xml_all))

    groups_xml = ET.SubElement(orthoxml_file, "groups")

    for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
        groups_xml.append(hogs_a_rhog_xml)

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str[:-1000])

    with open(working_folder+output_xml_name, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    print("orthoxml is written in "+ working_folder+output_xml_name)
    return 1


