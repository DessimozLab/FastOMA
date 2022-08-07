
from ete3 import Tree
from Bio import SeqIO
import dill as dill_pickle
import gc
from distributed import get_client

import _wrappers
import _utils
from _hog_class import HOG
from _utils import logger_hog
# from _dask_env import client_dask

def read_infer_xml_rhogs(rhogid_batch_list, vars_input):
    # (gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future,
    # dask_future_taxon, format_prot_name) = vars_input

    hogs_rhogs_xml_all = []
    print("There are "+str(len(rhogid_batch_list))+" rhogs in the batch.")
    for rhogid_num in rhogid_batch_list:
        hogs_rhogs_xml = read_infer_xml_rhog(rhogid_num, vars_input)
        hogs_rhogs_xml_all += hogs_rhogs_xml

    return hogs_rhogs_xml_all

def read_infer_xml_rhog(rhogid_num, vars_input):
    (gene_id_name, address_rhogs_folder, species_tree_address, gene_trees_folder, pickle_address, dask_future, dask_future_taxon, format_prot_name) = vars_input

    logger_hog.info("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")

    prot_address = address_rhogs_folder+"HOG_B"+str(rhogid_num).zfill(7)+".fa"
    rhog_i = list(SeqIO.parse(prot_address, "fasta"))
    logger_hog.info("number of proteins in the rHOG is "+str(len(rhog_i))+".")
    (species_tree) = _utils.read_species_tree(species_tree_address)
    (species_tree, species_names_rhog, prot_names_rhog) = _utils.prepare_species_tree(rhog_i, species_tree, format_prot_name)
    # species_tree.write();  print(species_tree.write())

    input_vars2 = (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder, format_prot_name)
    if len(rhog_i) > 30000:
        # dask_future_taxon = True
        print("Dask future taxon is on for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))
        hogs_a_rhog = infer_hogs_for_rhog_levels_recursively_future(species_tree, input_vars2)

    else:
        # dask_future_taxon = False
        print("Dask future taxon is off for hogid "+str(rhogid_num)+" with length "+str(len(rhog_i)))
        hogs_a_rhog = infer_hogs_for_rhog_levels_recursively(species_tree, input_vars2)

    logger_hog.info("subhogs in this level are "+' '.join(["[" + str(i) + "]" for i in hogs_a_rhog])+".")
    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        print(hog_i)
        if len(hog_i._members) > 1:
            # could be improved
            hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_rhogs_xml.append(hogs_a_rhog_xml)

    # for hogs in hogs_a_rhog:
    #     if isinstance(hogs, list): # bwecause of dask future, sometimes it is  a list of list, sometimes list of hogs. need to be imporved
    #         for hog_i in hogs:
    #             print(hog_i)
    #             if len(hog_i._members) > 1:
    #                 # could be improved
    #                 hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
    #                 hogs_rhogs_xml.append(hogs_a_rhog_xml)
    #     # else:  # bwecause of dask future, sometimes it is  a list of list, sometimes list of hogs. need to be imporved
    #     #     if len(hogs._members) > 1:
    #     #         hogs_a_rhog_xml = hogs.to_orthoxml(**gene_id_name)
    #     #         hogs_rhogs_xml.append(hogs_a_rhog_xml)
    print(hogs_rhogs_xml)
    logger_hog.info("we are not reporting single tone hogs in the output xml.")
    # how to handle empty hogs !? why happening and if not save pickle, file not exist error from rhog id list

    with open(pickle_address + '/file_' + str(rhogid_num) + '.pickle', 'wb') as handle:
        dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
    logger_hog.info("***** hogs are written as a pickle " + pickle_address + '/file_' + str(rhogid_num) + '.pickle')

    del hogs_a_rhog
    gc.collect()


    return hogs_rhogs_xml


# only one level parralelization
def infer_hogs_for_rhog_levels_recursively_future(sub_species_tree, input_vars2):
    # (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder, format_prot_name) = input_vars2

    if sub_species_tree.is_leaf():
        children_nodes = []
    else:
        children_nodes = sub_species_tree.children

    hogs_children_level_list_future = []
    # for node_species_tree_child in children_nodes:
    #     hogs_this_level_list_future_i = infer_hogs_for_rhog_levels_recursively_future(node_species_tree_child, rhog_i, species_names_rhog, rhogid_num, gene_trees_folder, format_prot_name)
    #     hogs_children_level_list_future += hogs_this_level_list_future_i
    # hogs_this_level_list = infer_hogs_this_level(sub_species_tree, rhog_i, species_names_rhog, hogs_children_level_list, rhogid_num, gene_trees_folder, format_prot_name)

    client_dask_working = get_client()
    hogs_children_level_list_futures = client_dask_working.map(infer_hogs_for_rhog_levels_recursively_future, children_nodes, [input_vars2]* len(children_nodes) )
    hogs_children_level_list = client_dask_working.gather(hogs_children_level_list_futures)

    print(sub_species_tree.name)
    hogs_this_level_list = infer_hogs_this_level(sub_species_tree, input_vars2, hogs_children_level_list)

    hogs_children_level_list_flatten = []
    for hogs_list in hogs_children_level_list:
        hogs_children_level_list_flatten += hogs_list

    return hogs_children_level_list_flatten


# only one level parralelization
def infer_hogs_for_rhog_levels_recursively(sub_species_tree, input_vars2):
    # input_vars2 = (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder, format_prot_name) = input_vars2

    if sub_species_tree.is_leaf():
        children_nodes = []
    else:
        children_nodes = sub_species_tree.children

    hogs_children_level_list = []
    for node_species_tree_child in children_nodes:
        hogs_children_level_list_i = infer_hogs_for_rhog_levels_recursively(node_species_tree_child, input_vars2)
        hogs_children_level_list.extend(hogs_children_level_list_i)
    hogs_this_level_list = infer_hogs_this_level(sub_species_tree, input_vars2, hogs_children_level_list)

    return hogs_this_level_list


def singletone_hog(node_species_tree, rhog_i, species_names_rhog, rhogid_num):

    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog)) if
                                 species_names_rhog[idx] == node_species_name]
    rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]

    hogs_this_level_list = []
    for prot in rhog_part:
        hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
        hogs_this_level_list.append(hog_leaf)
    return hogs_this_level_list


def infer_hogs_this_level(sub_species_tree, input_vars2, hogs_children_level_list):

    (rhog_i, species_names_rhog, rhogid_num, gene_trees_folder, format_prot_name) = input_vars2
    node_species_tree = sub_species_tree
    if len(node_species_tree.name.split("_"))>1:
        logger_hog.info("Finding hogs for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(
            node_species_tree.name)+"\n"+str(node_species_tree.write())+"\n")

    if node_species_tree.is_leaf():
        assert hogs_children_level_list == []
        hogs_this_level_list = singletone_hog(node_species_tree, rhog_i, species_names_rhog, rhogid_num)

        return hogs_this_level_list

    if len(hogs_children_level_list) == 1:
        hogs_this_level_list = hogs_children_level_list
        return hogs_this_level_list

    # print("*7*", len(hogs_children_level_list), hogs_children_level_list)
    # if len(hogs_children_level_list)>0:
    #    print("**",hogs_children_level_list[0], len(hogs_children_level_list))

    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list]
    merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready)
    # logger_hog.info("All subhogs are merged, merged msa is with length of " + str(len(merged_msa)) + " " + str(
    #   len(merged_msa[0])) + ".")

    gene_tree_file_addr = gene_trees_folder + "/tree_" + str(rhogid_num) + "_" + str(
        node_species_tree.name) + ".nwk"
    gene_tree_raw = _wrappers.infer_gene_tree(merged_msa, gene_tree_file_addr)
    gene_tree = Tree(gene_tree_raw + ";", format=0)
    logger_hog.info("Gene tree is inferred with length of " + str(len(gene_tree)) + ".")
    R_outgroup = gene_tree.get_midpoint_outgroup()
    gene_tree.set_outgroup(R_outgroup)  # print("Midpoint rooting is done for gene tree.")
    gene_tree = _utils.lable_SD_internal_nodes(gene_tree)
    # print("Overlap speciation is done for internal nodes of gene tree, as following:")
    print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")
    hogs_this_level_list = merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa)

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
    # logger_hog.info(str(len(hogs_this_level_list))+" hogs are inferred at the level "+node_species_tree.name+": "+' '.join(
    #     [str(i) for i in prot_list_sbuhog_short]))

    return hogs_this_level_list




def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa):

    """
    this function could be improved
    """
    subhogs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
    hogs_this_level_list = []
    subHOG_to_be_merged_set_other_Snodes = []
    subHOG_to_be_merged_set_other_Snodes_flattned_temp = []
    for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n,
                                                                                     "processed") and n.processed == True):  # start from root
        # print("Leaves assigned to hog are ", assigned_leaves_to_hog)   #print("Traversing gene tree. Now at node", node.name)
        if node.is_leaf():
            continue
        node_leaves_name = [i.name for i in node.get_leaves()]
        if node.name[0] == "S":  # this is a sub-hog.
            subHOG_to_be_merged = []
            for node_leave_name in node_leaves_name:  # print(node_leave_name)
                for subHOG in hogs_children_level_list:
                    subHOG_members = subHOG._members
                    if node_leave_name in subHOG_members:  # could be improved
                        if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp:
                            subHOG_to_be_merged.append(subHOG)
                            subhogs_id_children_assigned.append(subHOG._hogid)
                        else:
                            print("issue 184", node.name, subHOG._hogid, node_leave_name)
                            if "processed" in node:
                                print(node.name)
                            else:
                                print("processed not in ", node.name)  # print(node_leave_name,"is in ",subHOG._hogid)

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
    for subHOG in hogs_children_level_list:  # for the single branch  ( D include a  subhog and a S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)

    # this could be improved
    # we expect to see a list of list as ooutput
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    return hogs_this_level_list

