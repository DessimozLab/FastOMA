
from dask.distributed import rejoin, secede
from ete3 import Tree
from Bio import SeqIO
from distributed import get_client
import concurrent.futures

import time
import os
import shutil

import pickle
import gc

import _wrappers
import _utils
from _hog_class import HOG
from _utils import logger_hog
import _config
# import networkx as nx
# import matplotlib.pyplot as plt

from ete3 import PhyloTree



def read_infer_xml_rhogs_batch(rhogid_batch_list, inferhog_concurrent_on, folder=""):

    # for now each contain one rhog
    hogs_rhog_xml_batch = []
    for rhogid_num in rhogid_batch_list:
        hogs_rhogs_xml = read_infer_xml_rhog_v2(rhogid_num, inferhog_concurrent_on, folder)  # orthoxml_to_newick.py list of hog object
        hogs_rhog_xml_batch.extend(hogs_rhogs_xml)

    return hogs_rhog_xml_batch  # orthoxml_to_newick.py list of hog object


def read_infer_xml_rhog_v2(rhogid_num , inferhog_concurrent_on, folder=""):
    pickles_subhog_folder_all = _config.working_folder + "/pickles_subhog/"
    # if not os.path.exists(pickles_subhog_folder_all):
    #     os.makedirs(pickles_subhog_folder_all)

    pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
    if not os.path.exists(pickles_subhog_folder):
        os.makedirs(pickles_subhog_folder)

    logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
    rhog_i_prot_address = _config.working_folder + "rhogs/" + folder+"/HOG_B" + str(rhogid_num).zfill(7) + ".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    logger_hog.debug("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")
    (species_tree) = _utils.read_species_tree_add_internal(_config.species_tree_address)
    (species_tree, species_names_rhog, prot_names_rhog) = _utils.prepare_species_tree(rhog_i, species_tree, rhogid_num)
    species_names_rhog = list(set(species_names_rhog))
    logger_hog.debug(
        "The number of unique species in the rHOG " + str(rhogid_num) + "is " + str(len(species_names_rhog)) + ".")
    # species_tree.write();  print(species_tree.write())

    logger_hog.debug("Dask future taxon is off for hogid " + str(rhogid_num) + " with length " + str(len(rhog_i)))

    if inferhog_concurrent_on:  # _config.inferhog_concurrent_on:
        hogs_a_rhog_1 = infer_hogs_concurrent(species_tree, rhogid_num, folder )
    else:
        hogs_a_rhog_1 = infer_hogs_for_rhog_levels_recursively(species_tree, rhogid_num, folder)
    # hogs_a_rhog_1  is an integeer as the length

    root_node_name = species_tree.name
    pickle_subhog_file = pickles_subhog_folder + str(root_node_name)
    with open(pickle_subhog_file + ".pickle", 'rb') as handle:
        hogs_a_rhog = pickle.load(handle)
    # ?? logger_hog.debug("subhogs in this level are "+' '.join(["[" + str(i) + "]" for i in ?? ])+".") # hogs_a_rhog

    if not _config.keep_subhog_each_pickle:
        shutil.rmtree(pickles_subhog_folder)

    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        # print(hog_i)
        if len(hog_i._members) >= _config.inferhog_min_hog_size_xml:  # previously > 1
            # could be improved   # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_a_rhog_xml = hog_i.to_orthoxml()
            hogs_rhogs_xml.append(hogs_a_rhog_xml)
    logger_hog.debug("we are not reporting single tone hogs in the output xml.")
    # how to handle empty hogs !? why happening and if not save pickle, file not exist error from rhog id list
    pickles_rhog_folder = _config.working_folder + "/pickles_rhog/"
    pickles_rhog_file = pickles_rhog_folder + '/file_' + str(rhogid_num) + '.pickle'
    with open(pickles_rhog_file, 'wb') as handle:
        # dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("***** hogs are written as orthoxml_to_newick.py pickle " + pickles_rhog_file)

    del hogs_a_rhog
    gc.collect()

    return hogs_rhogs_xml


max_workers_num = 10  # config


def infer_hogs_concurrent(species_tree, rhogid_num, folder ="" ):
    pending_futures = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=_config.inferhog_max_workers_num) as executor:

        for node in species_tree.traverse(strategy="preorder"):
            node.dependencies_fulfilled = set()  # a set
            # node.infer_submitted = False

            if node.is_leaf():
                future_id = executor.submit(singletone_hog_, node, rhogid_num, folder)
                pending_futures[future_id] = node.name
                # node.infer_submitted = True

        while len(pending_futures) > 0:
            time.sleep(0.01)
            future_id_list = list(pending_futures.keys())
            for future_id in future_id_list:

                if future_id.done():
                    species_node_name = pending_futures[future_id]
                    del pending_futures[future_id]
                    species_node = species_tree.search_nodes(name=species_node_name)[0]

                    # print(future_id)
                    parent_node = species_node.up
                    if not parent_node:  # we reach the root
                        # assert len(pending_futures) == 0, str(species_node_name)+" "+str(rhogid_num)
                        assert species_node.name == species_tree.name
                        break

                    parent_node.dependencies_fulfilled.add(species_node_name)  # a set

                    childrend_parent_nodes = set(node.name for node in parent_node.get_children())
                    if parent_node.dependencies_fulfilled == childrend_parent_nodes:
                        # print("here", species_node_name)
                        #if not parent_node.infer_submitted:
                        future_id_parent = executor.submit(infer_hogs_this_level, parent_node, rhogid_num)
                            # parent_node.infer_submitted = True
                        # future_id_parent= parent_node.name+"aaa"
                        pending_futures[future_id_parent] = parent_node.name
                        #for future_id
                        #    del pending_futures[future_id]
                        # i need another dictionary the other way arround to removes this futures

            # print(" ", pending_futures)

    return len(pending_futures) + 1


#
# def read_infer_xml_rhog(rhogid_num):
#
#     pickles_subhog_folder_all = _config.working_folder + "/pickles_subhog/"
#
#     pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num)+"/"
#     if not os.path.exists(pickles_subhog_folder):
#         os.makedirs(pickles_subhog_folder)
#
#
#     logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
#     rhog_i_prot_address = _config.working_folder + "rhogs/" +"HOG_B"+str(rhogid_num).zfill(7)+".fa"
#     rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
#     logger_hog.debug("number of proteins in the rHOG is "+str(len(rhog_i))+".")
#     (species_tree) = _utils.read_species_tree(_config.species_tree_address)
#     (species_tree, species_names_rhog, prot_names_rhog) = _utils.prepare_species_tree(rhog_i, species_tree, rhogid_num)
#     species_names_rhog = list(set(species_names_rhog))
#     logger_hog.debug("The number of unique species in the rHOG " + str(rhogid_num) + "is " + str(len(species_names_rhog)) + ".")
#     # species_tree.write();  print(species_tree.write())
#
#     if len(rhog_i) > _config.inferhog_dask_2nd_rhogsize and (_config.dask_level == 2 or _config.dask_level == 3): # 200
#         logger_hog.debug("Dask future taxon is on for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))
#         client_dask_working = get_client()
#         secede()
#
#         # #recursive_input = (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder)
#         #hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_future_v2, species_tree, recursive_4inputs)
#         hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, species_tree, rhogid_num)
#         # #rhog_i_future = client_dask_working.scatter(rhog_i)
#         # #recursive_input_future = (rhog_i_future, species_names_rhog, rhogid_num, gene_trees_folder)
#         # #hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, species_tree, species_names_rhog, rhogid_num)
#         hogs_a_rhog = hogs_a_rhog_future.result()
#         rejoin()
#
#     elif (_config.dask_level == 2 or _config.dask_level == 3):
#         # logger_hog.debug("Dask future taxon is off for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))
#         client_dask_working = get_client()
#         secede()
#         hogs_a_rhog_future = client_dask_working.submit(infer_hogs_for_rhog_levels_recursively, species_tree, rhogid_num)
#         hogs_a_rhog = hogs_a_rhog_future.result()
#         rejoin()
#
#     else:
#         logger_hog.debug("Dask future taxon is off for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))
#         hogs_a_rhog_1 = infer_hogs_for_rhog_levels_recursively(species_tree, rhogid_num)
#         # hogs_a_rhog_1  is an integeer as the length
#
#     root_node_name= species_tree.name
#     pickle_subhog_file = pickles_subhog_folder + str(root_node_name)
#     with open(pickle_subhog_file+".pickle", 'rb') as handle:
#         hogs_a_rhog = pickle.load(handle)
#     # ?? logger_hog.debug("subhogs in this level are "+' '.join(["[" + str(i) + "]" for i in ?? ])+".") # hogs_a_rhog
#
#     if not _config.keep_subhog_each_pickle:
#         shutil.rmtree(pickles_subhog_folder)
#
#     hogs_rhogs_xml = []
#     for hog_i in hogs_a_rhog:
#         # print(hog_i)
#         if len(hog_i._members) >= _config.inferhog_min_hog_size_xml:  # previously > 1
#             # could be improved   # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
#             hogs_a_rhog_xml = hog_i.to_orthoxml()
#             hogs_rhogs_xml.append(hogs_a_rhog_xml)
#     logger_hog.debug("we are not reporting single tone hogs in the output xml.")
#     # how to handle empty hogs !? why happening and if not save pickle, file not exist error from rhog id list
#     pickles_rhog_folder = _config.working_folder + "/pickles_rhog/"
#     pickles_rhog_file = pickles_rhog_folder + '/file_' + str(rhogid_num) + '.pickle'
#     with open(pickles_rhog_file, 'wb') as handle:
#         #dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
#         pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     logger_hog.debug("***** hogs are written as orthoxml_to_newick.py pickle " + pickles_rhog_file)
#
#     del hogs_a_rhog
#     gc.collect()
#
#     return hogs_rhogs_xml

#############
#
# def infer_hogs_for_rhog_future_v2(sub_species_tree, species_names_rhog, rhogid_num):
#
#     out_1 = 0
#     species_leaves_names = [i.name for i in sub_species_tree.get_leaves()]
#     if len(species_leaves_names) <= 10:
#         out_1 = infer_hogs_for_rhog_subtree_v2(sub_species_tree, species_names_rhog, rhogid_num)
#     # hogs_level_list_futures_all = []
#     client_dask_working = get_client()
#     for node_species_tree in sub_species_tree.traverse(strategy="postorder"):
#         if len(species_leaves_names) > 10:
#             # secede()
#             hogs_level_list_futures = client_dask_working.submit(infer_hogs_for_rhog_subtree_v2, node_species_tree, species_names_rhog, rhogid_num)
#             # hogs_level_list_futures_all.append(hogs_level_list_futures)
#             # out_1 = client_dask_working.gather(hogs_level_list_futures)
#             # out_1 = hogs_level_list_futures.result()
#             # rejoin()
#             for node in node_species_tree.traverse():
#                 node.processed = True
#     # hogs_level_list_futures_all_result = client_dask_working.gather(hogs_level_list_futures_all)
#     # rejoin()
#
#     return out_1
#    # recursive_4inputs = (species_names_rhog, rhogid_num, gene_trees_folder, address_rhogs_folder)
# def infer_hogs_for_rhog_NOfuture_v2(sub_species_tree, recursive_4inputs):
#     out_1 = 0
#     species_leaves_names = [i.name for i in sub_species_tree.get_leaves()]
#     if len(species_leaves_names) <= 10:
#         out_1 = infer_hogs_for_rhog_subtree_v2(sub_species_tree, recursive_4inputs)
#
#     for node_species_tree in sub_species_tree.traverse(strategy="postorder"):
#         species_leaves_names = [i.name for i in node_species_tree.get_leaves()]
#
#         if len(species_leaves_names) > 10:
#             out_1 = infer_hogs_for_rhog_subtree_v2(node_species_tree, recursive_4inputs)
#             for node in node_species_tree.traverse():
#                 node.processed = True
#     return out_1

#
# def infer_hogs_for_rhog_levels_recursively_future(sub_species_tree, rhogid_num):
#     # (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder) = recursive_input
#     # logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
#
#     if sub_species_tree.is_leaf():
#
#         singletone_hog_out = singletone_hog_(sub_species_tree, rhogid_num)
#         return singletone_hog_out
#
#     children_nodes = sub_species_tree.children
#
#     client_dask_working = get_client()
#     secede()
#     hogs_children_level_list_futures = [client_dask_working.submit(infer_hogs_for_rhog_levels_recursively_future, child, rhogid_num) for child in children_nodes]
#     hogs_children_level_list_futures = client_dask_working.gather(hogs_children_level_list_futures)
#     rejoin()
#     # hogs_children_level_list = hogs_children_level_list_futures
#     # hogs_children_level_list = []
#     # for future in hogs_children_level_list_futures:
#     #    hogs_children_level_list.extend(future.result())
#     # if hogs_children_level_list_futures:
#     #     if isinstance(hogs_children_level_list_futures[0], list):
#     #         hogs_children_level_list_flatten = []
#     #         for hog_ in hogs_children_level_list_futures:
#     #             # for hog in hogs_list:
#     #             hogs_children_level_list_flatten.extend(hog_)
#     # hogs_children_level_list = hogs_children_level_list_flatten
#     infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, rhogid_num) # hogs_children_level_list
#
#     return infer_hogs_this_level_out

def infer_hogs_for_rhog_levels_recursively(sub_species_tree, rhogid_num, folder=""):

    if sub_species_tree.is_leaf():
        # hogs_this_level_list = singletone_hog(sub_species_tree, rhog_i, species_names_rhog, rhogid_num)
        singletone_hog_out = singletone_hog_(sub_species_tree, rhogid_num, folder)
        # out 1 succesful
        return singletone_hog_out
    children_nodes = sub_species_tree.children

    hogs_children_level_list = []
    for node_species_tree_child in children_nodes:
        hogs_children_level_list_i = infer_hogs_for_rhog_levels_recursively(node_species_tree_child, rhogid_num,folder)
        # hogs_children_level_list_i should be 1
        # hogs_children_level_list.extend(hogs_children_level_list_i)
    infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, rhogid_num)  # ,hogs_children_level_list
    # hogs_this_level_list should be one

    return infer_hogs_this_level_out
#
#
# def infer_hogs_for_rhog_subtree_v2(sub_species_tree, species_names_rhog, rhogid_num):
#
#     infer_hogs_this_level_out= 0
#     if sub_species_tree.is_leaf():
#         if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
#
#             infer_hogs_this_level_out = singletone_hog_(sub_species_tree, species_names_rhog, rhogid_num)
#         return infer_hogs_this_level_out
#
#     children_nodes = sub_species_tree.children
#     for node_species_tree_child in children_nodes:
#         hogs_children_level_list_i = infer_hogs_for_rhog_subtree_v2(node_species_tree_child, species_names_rhog, rhogid_num)
#     infer_hogs_this_level_out = 0
#     if not (hasattr(sub_species_tree, "processed") and sub_species_tree.processed == True):
#         infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, rhogid_num)  # output is an integer, which is the length
#
#     return infer_hogs_this_level_out
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


def singletone_hog_(node_species_tree, rhogid_num, folder=""):
    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    this_level_node_name = node_species_name

    if _config.inferhog_resume_subhog:
        pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)
        # open already calculated subhogs , but not completed till root in previous run
        if os.path.exists(pickle_subhog_file + ".pickle"):
            if os.path.getsize(pickle_subhog_file + ".pickle") > 3:  # 3 bytes
                with open(pickle_subhog_file + ".pickle", 'rb') as handle:
                    # i don't even need to open this even
                    hogs_children_level_list = pickle.load(handle)
                    if hogs_children_level_list:  # hogs_children_level_list is sth like [an object of class HOG of hogID=HOG:B0574027_sub10001, length=1, taxonomy= PSETE_]
                        return len(hogs_children_level_list)

    logger_hog.debug("* reading prot address  " + str(this_level_node_name))
    rhog_i_prot_address = _config.working_folder + "rhogs/"+folder+"/HOG_B"+str(rhogid_num).zfill(7)+".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))

    species_names_rhog_nonuniq = [seq.id.split("||")[1] for seq in rhog_i]
    prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog_nonuniq)) if
                                 species_names_rhog_nonuniq[idx] == node_species_name]
    rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]

    hogs_this_level_list = []
    for prot in rhog_part:
        hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
        hogs_this_level_list.append(hog_leaf)
    pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)
    with open(pickle_subhog_file+".pickle", 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs was written as pickle file.")

    return len(hogs_this_level_list)


def infer_hogs_this_level(sub_species_tree, rhogid_num):  # hogs_children_level_list

    node_species_tree = sub_species_tree

    this_level_node_name = node_species_tree.name
    if node_species_tree.is_leaf():
        print(" issue 1235, single tone hogs are treated in another function, the code shouldn't reach here.", rhogid_num )
        exit()
        # assert hogs_children_level_list == []
        #hogs_this_level_list = singletone_hog_(node_species_tree, species_names_rhog, rhogid_num)
        # we shouldnt be here ???
        #
        # child_name = node_species_tree.name
        # pickle_subhog_file = hogs_children_level_pickle_folder + "rhog_" + str(rhogid_num) + "/" + str(child_name)
        # with open(pickle_subhog_file+".pickle", 'wb') as handle:
        #     pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 1  # [child_name] # hogs_this_level_list

    pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)

    if _config.inferhog_resume_subhog:
        # open already calculated subhogs , but not completed till root in previous run
        if  os.path.exists(pickle_subhog_file + ".pickle"):
            if os.path.getsize(pickle_subhog_file + ".pickle") > 3: # bytes
                with open(pickle_subhog_file + ".pickle", 'rb') as handle:
                    # i don't even need to open this even
                    hogs_children_level_list = pickle.load(handle)
                    if hogs_children_level_list:
                        return len(hogs_children_level_list)

    # print(this_level_node_name, rhogid_num)
    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list = []
    pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
    for child_name in children_name:
        pickle_subhog_file = pickles_subhog_folder + str(child_name)
        with open(pickle_subhog_file+".pickle", 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))
      # to do, check file doenst exist ? how to handle, if not
    # if len(node_species_tree.name.split("_")) > 1:
    logger_hog.debug("Finding hogs for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name)+"\n"+str(node_species_tree.write())+"\n")

    if len(hogs_children_level_list) == 1:
        assert len(children_name) == 1
        hogs_this_level_list = hogs_children_level_list
        pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)
        with open(pickle_subhog_file+".pickle", 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return len(hogs_children_level_list)

    # hogs_children_level_list
    # if isinstance(hogs_children_level_list[0], list):
    #     hogs_children_level_list_flatten = []
    #     for hogs_list in hogs_children_level_list:
    #         for hog in hogs_list:
    #         hogs_children_level_list_flatten.extend(hog)
    #
    # hogs_children_level_list = hogs_children_level_list_flatten

    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list]
    gene_tree_file_addr = _config.working_folder + "/genetrees/tree_" + str(rhogid_num) + "_" + str(
        node_species_tree.name) + ".nwk"
    if len(gene_tree_file_addr) > 245:
        # there is a limitation on length of file name. I want to  keep it consistent ,msa and gene tree names.
        import random
        rand_num = random.randint(1, 10000)
        gene_tree_file_addr = gene_tree_file_addr[:245] + str(rand_num)+".nwk"
    logger_hog.debug("Merging "+str(len(sub_msa_list_lowerLevel_ready))+" MSAs for rhogid_num: "+
                     str(rhogid_num)+", for taxonomic level:"+str(node_species_tree.name))

    merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready, gene_tree_file_addr)

    if merged_msa:
        logger_hog.debug("All sub-hogs are merged, merged msa is with length of " + str(len(merged_msa)) + " " + str(
        len(merged_msa[0])) + " for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
                node_species_tree.name))
        # merged_msa_filt = merged_msa
        # 893*4839, 10 mins

        if len(merged_msa[0]) >= _config.inferhog_min_cols_msa_to_filter:
            # (len(merged_msa) > 10000 and len(merged_msa[0]) > 3000) or (len(merged_msa) > 500 and len(merged_msa[0]) > 5000) or (len(merged_msa) > 200 and len(merged_msa[0]) > 9000):
            # for very big MSA, gene tree is slow. if it is full of gaps, let's trim the msa.
            logger_hog.debug("We are doing MSA trimming "+str(rhogid_num)+", for taxonomic level:"+str(node_species_tree.name))
            # print(len(merged_msa), len(merged_msa[0]))
            msa_filt_col = _utils.msa_filter_col(merged_msa, _config.inferhog_tresh_ratio_gap_col, gene_tree_file_addr)
            # print(len(msa_filt_col), len(msa_filt_col[0]))
            if msa_filt_col and msa_filt_col[0] and len(msa_filt_col[0]):
                msa_filt_row_col = _utils.msa_filter_row(msa_filt_col, _config.inferhog_tresh_ratio_gap_row, gene_tree_file_addr)
            else:
                msa_filt_row_col = msa_filt_col
            # print(len(msa_filt_row_col), len(msa_filt_row_col[0]))
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

    #gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
    # outliers = find_outlier_leaves(gene_tree)
    # R = midpoint_rooting_outgroup(gene_tree, leaves_to_exclude=outliers)
    # gene_tree.set_outgroup(R)

    gene_tree = _utils.lable_sd_internal_nodes(gene_tree)
    # print("Overlap speciation is done for internal nodes of gene tree, as following:")
    # print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")
    logger_hog.debug("Merging sub-hogs of children started  for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name))

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

    pickles_subhog_folder = _config.working_folder + "/pickles_subhog/rhog_" + str(rhogid_num) + "/"
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)
    with open(pickle_subhog_file + ".pickle", 'wb') as handle:
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

