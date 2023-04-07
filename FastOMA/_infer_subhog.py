
from ete3 import Tree
from ete3 import PhyloTree
from Bio import SeqIO
import concurrent.futures
import time
import os
import shutil
import pickle
import gc
import random
# import networkx as nx
# import matplotlib.pyplot as plt
from . import _wrappers
from . import _utils_subhog
from ._hog_class import HOG
from ._utils_subhog import logger_hog
from . import _config


def read_infer_xml_rhogs_batch(rhogid_batch_list, inferhog_concurrent_on, pickles_rhog_folder, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    infer subHOGs for a list of rootHOGs
    """
    logger_hog.debug("Inferring subHOGs for  "+str(len(rhogid_batch_list))+"rootHOGs started.")
    logger_hog.debug("we are not reporting single tone hogs in the output xml. You may check this _config.inferhog_min_hog_size_xml.")
    hogs_rhog_xml_len_batch = []
    for rhogid_num in rhogid_batch_list:
        hogs_rhogs_xml_len = read_infer_xml_rhog(rhogid_num, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder)
        hogs_rhog_xml_len_batch.extend(hogs_rhogs_xml_len)

    return hogs_rhog_xml_len_batch


def read_infer_xml_rhog(rhogid_num, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder):
    """
    infer subHOGs for a  rootHOGs
    """

    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    if not os.path.exists(pickles_subhog_folder):
        os.makedirs(pickles_subhog_folder)

    logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
    rhog_i_prot_address = rhogs_fa_folder + "/HOG_" + str(rhogid_num).zfill(7) + ".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    logger_hog.debug("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")
    (species_tree) = _utils_subhog.read_species_tree_add_internal(_config.species_tree_address)
    # todo double check the issue with the root node name.
    # todo check all the input and needed folders as a validation with nextflow
    (species_tree, species_names_rhog, prot_names_rhog) = _utils_subhog.prepare_species_tree(rhog_i, species_tree, rhogid_num)
    species_names_rhog = list(set(species_names_rhog))
    logger_hog.debug("Number of unique species in rHOG " + str(rhogid_num) + "is " + str(len(species_names_rhog)) + ".")

    if inferhog_concurrent_on:  # for big HOG we use paralelization at the level taxanomic level using concurrent
        hogs_a_rhog_num = infer_hogs_concurrent(species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
    else:
        hogs_a_rhog_num = infer_hogs_for_rhog_levels_recursively(species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
    # Output value hogs_a_rhog_num  is an integer= length. We save the output as pickle file at each taxanomic level.

    #####  Now read the final pickle file for this rootHOG
    root_node_name = species_tree.name
    pickle_subhog_file = pickles_subhog_folder + str(root_node_name) + ".pickle"
    with open(pickle_subhog_file, 'rb') as handle:
        hogs_a_rhog = pickle.load(handle)
    if not _config.keep_subhog_each_pickle:
        shutil.rmtree(pickles_subhog_folder)

    hogs_rhogs_xml = []
    for hog_i in hogs_a_rhog:
        if len(hog_i._members) >= _config.inferhog_min_hog_size_xml:
            # could be improved   # hogs_a_rhog_xml = hog_i.to_orthoxml(**gene_id_name)
            hogs_a_rhog_xml = hog_i.to_orthoxml()
            hogs_rhogs_xml.append(hogs_a_rhog_xml)

    pickles_rhog_file = pickles_rhog_folder + '/file_' + str(rhogid_num) + '.pickle'
    with open(pickles_rhog_file, 'wb') as handle:
        # dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("All subHOGs for the rootHOG as orthoxml format is written in " + pickles_rhog_file)

    del hogs_a_rhog  # to be memory efficient
    gc.collect()
    hogs_rhogs_xml_len = len(hogs_rhogs_xml)
    return hogs_rhogs_xml_len


def infer_hogs_concurrent(species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    infer subHOGs for a rootHOG using multi-threading (in parallel) on different taxanomic levels of species tree
    """

    pending_futures = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=_config.inferhog_max_workers_num) as executor:
        for node in species_tree.traverse(strategy="preorder"):
            node.dependencies_fulfilled = set()  # a set
            # node.infer_submitted = False
            if node.is_leaf():
                future_id = executor.submit(singletone_hog_, node, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
                # singletone_hog_(sub_species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
                # {<Future at 0x7f1b48d9afa0 state=finished raised TypeError>: 'KORCO_'}
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
                    parent_node = species_node.up
                    if not parent_node:  # we reach the root
                        # assert len(pending_futures) == 0, str(species_node_name)+" "+str(rhogid_num)
                        assert species_node.name == species_tree.name
                        break
                    parent_node.dependencies_fulfilled.add(species_node_name)  # a set

                    childrend_parent_nodes = set(node.name for node in parent_node.get_children())
                    if parent_node.dependencies_fulfilled == childrend_parent_nodes:
                        #  if not parent_node.infer_submitted:
                        future_id_parent = executor.submit(infer_hogs_this_level, parent_node, rhogid_num, pickles_subhog_folder_all, prots_to_remove)
                        # parent_node.infer_submitted = True
                        # future_id_parent= parent_node.name+"aaa"
                        pending_futures[future_id_parent] = parent_node.name
                        # for future_id:  del pending_futures[future_id] i need another dictionary the other way arround to removes this futures

    return len(pending_futures) + 1


def infer_hogs_for_rhog_levels_recursively(sub_species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    infer subHOGs for a rootHOG using recursive function to traverse species tree (different taxanomic levels)
    """

    if sub_species_tree.is_leaf():
        singletone_hog_out = singletone_hog_(sub_species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
        #  out 1 =  succesful
        return singletone_hog_out
    children_nodes = sub_species_tree.children

    for node_species_tree_child in children_nodes:
        hogs_chrdn = infer_hogs_for_rhog_levels_recursively(node_species_tree_child, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder)
        # hogs_chrdn should be 1 hogs_chrdn_list.extend(hogs_chrdn)
    infer_hogs_this_level_out = infer_hogs_this_level(sub_species_tree, rhogid_num, pickles_subhog_folder_all)
    # hogs_this_level_list should be one
    return infer_hogs_this_level_out


def singletone_hog_(node_species_tree, rhogid_num, pickles_subhog_folder_all, rhogs_fa_folder):
    """
    create subHOGs for leaves of tree (species level). Each protein is a subHOG.
    """

    node_species_name = node_species_tree.name  # there is only one species (for the one protein)
    this_level_node_name = node_species_name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    if _config.inferhog_resume_subhog:
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
        # open already calculated subhogs , but not completed till root in previous run
        if os.path.exists(pickle_subhog_file):
            if os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
                with open(pickle_subhog_file, 'rb') as handle:
                    # i don't even need to open this even
                    hogs_children_level_list = pickle.load(handle) #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                    if hogs_children_level_list:
                        return len(hogs_children_level_list)
    # logger_hog.debug("reading protien / singletone HOG of  " + str(this_level_node_name))
    rhog_i_prot_address = rhogs_fa_folder +"/HOG_"+str(rhogid_num).zfill(7)+".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    species_names_rhog_nonuniq = [seq.id.split("||")[1] for seq in rhog_i]
    prot_idx_interest_in_rhog = [idx for idx in range(len(species_names_rhog_nonuniq)) if
                                 species_names_rhog_nonuniq[idx] == node_species_name]
    rhog_part = [rhog_i[i] for i in prot_idx_interest_in_rhog]
    hogs_this_level_list = []
    for prot in rhog_part:
        hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
        hogs_this_level_list.append(hog_leaf)
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)+".pickle"
    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs was written as pickle file.")

    return len(hogs_this_level_list)


def read_children_hogs(node_species_tree, rhogid_num, pickles_subhog_folder_all):
    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
    # TODO arrage resume with nextflow
    # if _config.inferhog_resume_subhog: # open already calculated subhogs , but not completed till root in previous run
    #     if os.path.exists(pickle_subhog_file):
    #         if os.path.getsize(pickle_subhog_file) > 3: # bytes
    #             with open(pickle_subhog_file, 'rb') as handle:
    #                 # i don't even need to open this even
    #                 hogs_children_level_list = pickle.load(handle)
    #                 if hogs_children_level_list:
    #                     return len(hogs_children_level_list)
    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list = []
    for child_name in children_name:
        pickle_subhog_file = pickles_subhog_folder + str(child_name) + ".pickle"
        with open(pickle_subhog_file, 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))  # todo, check file exist, how to handle if not
    logger_hog.debug("Finding hogs for rhogid_num: " + str(rhogid_num) + ", for taxonomic level:" + str(
        node_species_tree.name) + " for species sub-tree:\n  " + str(node_species_tree.write(format=1)) + "\n")
    return hogs_children_level_list


def genetree_SD_inference(node_species_tree, rhogid_num, msa_filt_row_col, genetree_msa_file_addr):

    gene_tree_raw = _wrappers.infer_gene_tree(msa_filt_row_col, genetree_msa_file_addr)
    gene_tree = Tree(gene_tree_raw + ";", format=0)
    logger_hog.debug("Gene tree is inferred len"+str(len(gene_tree))+" rhog:"+str(rhogid_num)+", level: "+str(node_species_tree.name))

    if _config.rooting_method == "midpoint":
        r_outgroup = gene_tree.get_midpoint_outgroup()
        gene_tree.set_outgroup(r_outgroup)  # print("Midpoint rooting is done for gene tree.")
    elif _config.rooting_method == "mad":
        gene_tree = _wrappers.mad_rooting(genetree_msa_file_addr)
    # elif _config.rooting_method == "outlier":
    #     gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
    #     outliers = find_outlier_leaves(gene_tree)
    #     r_outgroup = midpoint_rooting_outgroup(gene_tree, leaves_to_exclude=outliers)
    #     gene_tree.set_outgroup(r_outgroup)

    species_suspicious_sd_list = []
    if _config.lable_SD_internal == "species_overlap":
        (gene_tree, species_suspicious_sd_list) = _utils_subhog.lable_sd_internal_nodes(gene_tree)

    elif _config.lable_SD_internal == "reconcilation":
        node_species_tree_nwk_string = node_species_tree.write(format=1)
        node_species_tree_PhyloTree = PhyloTree(node_species_tree_nwk_string, format=1)
        gene_tree_nwk_string = gene_tree.write(format=1)
        gene_tree_PhyloTree = PhyloTree(gene_tree_nwk_string, format=1)
        gene_tree = _utils_subhog.lable_SD_internal_nodes_reconcilation(gene_tree_PhyloTree,node_species_tree_PhyloTree)

    if _config.gene_trees_write:
        tree_nwk_SD_labeled = str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;"
        file_gene_tree = open(genetree_msa_file_addr + "_SD_labeled.nwk", "w")
        file_gene_tree.write(tree_nwk_SD_labeled)
        file_gene_tree.close()

    return (gene_tree, species_suspicious_sd_list)



def infer_hogs_this_level(node_species_tree, rhogid_num, pickles_subhog_folder_all):
    """
    infer subHOGs for a rootHOG at a taxanomic level
    """

    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    assert not node_species_tree.is_leaf(), "issue 1235,singletone hogs are treated elsewhere"+ str(rhogid_num)
    hogs_children_level_list = read_children_hogs(node_species_tree, rhogid_num, pickles_subhog_folder_all)

    if len(hogs_children_level_list) == 1:
        hogs_this_level_list = hogs_children_level_list
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)+".pickle"
        with open(pickle_subhog_file, 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return len(hogs_children_level_list)

    genetree_msa_file_addr = "./genetrees/tree_"+str(rhogid_num)+"_"+str(node_species_tree.name) + ".nwk"
    if len(genetree_msa_file_addr) > 245:
        # there is a limitation on length of file name. I want to  keep it consistent ,msa and gene tree names.
        rand_num = random.randint(1, 10000)
        genetree_msa_file_addr = genetree_msa_file_addr[:245] + str(rand_num)+".nwk"

    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list if len(hog._msa) > 0]
    # sub_msa_list_lowerLevel_ready = [ii for ii in sub_msa_list_lowerLevel_ready_raw if len(ii) > 0]
    logger_hog.debug("Merging "+str(len(sub_msa_list_lowerLevel_ready))+" MSAs for rhog:"+str(rhogid_num)+", level:"+str(node_species_tree.name))
    msa_filt_row_col = []
    if sub_msa_list_lowerLevel_ready:
        if len(sub_msa_list_lowerLevel_ready) > 1:
            merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready, genetree_msa_file_addr)
        else:
            merged_msa = sub_msa_list_lowerLevel_ready #???  when only on  child, the rest msa is empty ?? this could be improved
        logger_hog.debug("All sub-hogs are merged, merged msa is with length of " + str(len(merged_msa)) + " " + str(
        len(merged_msa[0])) + " for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(node_species_tree.name))
        (msa_filt_row_col, msa_filt_col, hogs_children_level_list) = _utils_subhog.filter_msa(merged_msa, genetree_msa_file_addr, hogs_children_level_list)
        # msa_filt_col is used for parent level of HOG. msa_filt_row_col is used for gene tree inference.
    else:
        logger_hog.info("Issue 1455, merged_msa is empty " + str(rhogid_num) + ", for taxonomic level:" + str(node_species_tree.name))

    if len(msa_filt_row_col) > 1 and len(msa_filt_row_col[0]) > 1:

        (gene_tree, species_suspicious_sd_list) = genetree_SD_inference(node_species_tree, rhogid_num, msa_filt_row_col,genetree_msa_file_addr)

        logger_hog.debug("Merging sub-hogs for rhogid_num:"+str(rhogid_num)+", level:"+str(node_species_tree.name))
        hogs_this_level_list = merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, msa_filt_col)
        # for i in hogs_this_level_list: print(i.get_members())
        logger_hog.debug("Hogs of this level is found for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(this_level_node_name))

    else:
        if msa_filt_row_col:
            logger_hog.debug("** hogs_this_level_list is empty. msa_filt_row_col:"+str(len(msa_filt_row_col))+" *"+str(len(msa_filt_row_col[0]))+" !!")
        else:
            logger_hog.debug("** msa_filt_row_col is empty." + str(len(msa_filt_row_col)) +"! .")

        hogs_this_level_list = hogs_children_level_list
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name)+ ".pickle"
    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return len(hogs_children_level_list)


def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa, dubious_list_list):

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

    for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "processed") and n.processed==True):
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

                for dubious_list in dubious_list_list:
                    print("1")

                # for fratgmetns_set in fragmnets_list:
                #     find the subhog for  fratgmetns_set
                #
                #     move the fragments  into subhog 1
                #             update the list of dubios  on the way back
                #     update
                #             fragmnets_list
                #     methods for each subhog
                #
                #             1) add_fragmetns   add to the list dubious all the subhog
                #
                #             2) remove protin     from members hierarcyh  all subhogs excpet subhog
                #
                #             add as a new method
                #             HOG_this_node.mark_prots_dubious(dubious_prots)
                #             result_removal = HOG_this_node.remove_prots_from_hog(prots_to_remove)
                #         if result_removal != 0:
                #         hogs_children_level_list.append(hog_i)
                #
                hogs_this_level_list.append(HOG_this_node)

                # taxnomic_range, rhogid_num, msa = None, fragment_list

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
    # num = random.randint(3, 1000000)
    # plt.savefig("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/hoggraph/" + hogs_children_level_list[0]._hogid[4:] + "file_rndm"+str(num)+".jpg")
    # # plt.show()
    for subHOG in hogs_children_level_list:  # for the single branch  ( D include orthoxml_to_newick.py  subhog and orthoxml_to_newick.py S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)
    # this could be improved
    # we expect to see orthoxml_to_newick.py list of list as ooutput
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    for hog_i in hogs_this_level_list:
        # I recently added this part in 31 march, how did work it before ?
        _utils_subhog.remove_prots_from_hog_hierarchy(hog_i, prots_to_remove)

    for hog_j in hogs_this_level_list:
        hog_j._taxnomic_range = node_species_tree.name

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


    return hogs_this_level_list

#