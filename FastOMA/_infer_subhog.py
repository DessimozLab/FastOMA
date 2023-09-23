

from Bio import SeqIO
import concurrent.futures
import time
import os
import shutil
import pickle
import gc
import random
from ete3 import Tree
import xml.etree.ElementTree as ET

# import networkx as nx
# import matplotlib.pyplot as plt
from . import _wrappers
from . import _utils_subhog
from . import _utils_frag_SO_detection
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
        hogs_rhog_xml_len_batch.append(hogs_rhogs_xml_len)

    return hogs_rhog_xml_len_batch


def read_infer_xml_rhog(rhogid_num, inferhog_concurrent_on, pickles_rhog_folder,  pickles_subhog_folder_all, rhogs_fa_folder):
    """
    infer subHOGs for a  rootHOGs
    """

    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    if not os.path.exists(pickles_subhog_folder):
        os.makedirs(pickles_subhog_folder)

    # if (_config.gene_trees_write or _config.msa_write) and not os.path.exists("./genetrees"):
    #     os.makedirs("./genetrees")

    logger_hog.debug("\n" + "==" * 10 + "\n Start working on root hog: " + str(rhogid_num) + ". \n")
    rhog_i_prot_address = rhogs_fa_folder + "/HOG_" + str(rhogid_num).zfill(7) + ".fa"
    rhog_i = list(SeqIO.parse(rhog_i_prot_address, "fasta"))
    logger_hog.debug("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")
    (species_tree) = _utils_subhog.read_species_tree_add_internal(_config.species_tree_address)

    (species_tree, species_names_rhog, prot_names_rhog) = _utils_subhog.prepare_species_tree(rhog_i, species_tree, rhogid_num)
    species_names_rhog = list(set(species_names_rhog))
    logger_hog.debug("Number of unique species in rHOG " + str(rhogid_num) + " is " + str(len(species_names_rhog)) + ".")

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
            hogs_a_rhog_xml_raw = hog_i.to_orthoxml()    # <generef  >      <paralg object >
            if _config.orthoxml_v03 and 'paralogGroup' in str(hogs_a_rhog_xml_raw) :
                # in version v0.3 of orthoxml, there shouldn't be any paralogGroup at root level. Let's put them inside an orthogroup should be in
                hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(hog_i._hogid)})
                property_element = ET.SubElement(hog_elemnt, "property", attrib={"name": "TaxRange", "value": str(hog_i._tax_now)})
                hog_elemnt.append(hogs_a_rhog_xml_raw)
                hogs_a_rhog_xml = hog_elemnt
            else:
                hogs_a_rhog_xml = hogs_a_rhog_xml_raw
            hogs_rhogs_xml.append(hogs_a_rhog_xml)
        else:
            logger_hog.debug("single tone hog "+str(hog_i._members)+" is not reported")

    pickles_rhog_file = pickles_rhog_folder + '/file_' + str(rhogid_num) + '.pickle'
    with open(pickles_rhog_file, 'wb') as handle:
        # dill_pickle.dump(hogs_rhogs_xml, handle, protocol=dill_pickle.HIGHEST_PROTOCOL)
        pickle.dump(hogs_rhogs_xml, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger_hog.debug("All subHOGs for the rootHOG as OrthoXML format is written in " + pickles_rhog_file)
    # to see orthoxml as string, you might need to do it for different idx
    # idx=0; from xml.dom import minidom; import xml.etree.ElementTree as ET; minidom.parseString(ET.tostring(hogs_rhogs_xml[idx])).toprettyxml(indent="   ")
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
                        future_id_parent = executor.submit(infer_hogs_this_level, parent_node, rhogid_num, pickles_subhog_folder_all)
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
    # logger_hog.debug(" ** inferhog_resume_subhog is " + str(_config.inferhog_resume_subhog))
    if _config.inferhog_resume_subhog:
        # logger_hog.debug("inferhog_resume_subhog is " + str(_config.inferhog_resume_subhog) + " so, we are reading from pickles.")
        pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
        # open already calculated subhogs , but not completed till root in previous run
        if os.path.exists(pickle_subhog_file):
            if os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
                with open(pickle_subhog_file, 'rb') as handle:
                    # i don't even need to open this even
                    # is output of pickle.load(handle) is chlired or this level ?
                    # todo I think I don't need to read the pickle file
                    hogs_this_level_list = pickle.load(handle) #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                    if hogs_this_level_list:
                        logger_hog.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                        return len(hogs_this_level_list)
                    else:
                        logger_hog.debug(" Issue  1238510: the pickle file for single tone is empty "+ str(hogs_this_level_list)+" " + str(rhogid_num))

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
    logger_hog.debug("HOGs for  " + str(this_level_node_name)+" including "+str(len(hogs_this_level_list))+ " hogs is written in pickle file.")

    return len(hogs_this_level_list)


def read_children_hogs(node_species_tree, rhogid_num, pickles_subhog_folder_all):
    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    # pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"
    # if _config.inferhog_resume_subhog:  # open already calculated subhogs , but may not be completed in previous run
    #     if os.path.exists(pickle_subhog_file):
    #         if os.path.getsize(pickle_subhog_file) > 3:  # bytes
    #             with open(pickle_subhog_file, 'rb') as handle:
    #                 # i don't even need to open this even
    #                 hogs_children_level_list = pickle.load(handle)
    #                 if hogs_children_level_list:
    #                     return hogs_children_level_list

    children_name = [child.name for child in node_species_tree.children]
    hogs_children_level_list = []
    for child_name in children_name:
        pickle_subhog_file = pickles_subhog_folder + str(child_name) + ".pickle"
        with open(pickle_subhog_file, 'rb') as handle:
            hogs_children_level_list.extend(pickle.load(handle))  # when there is an error somewhere else, probably with paralelization, for big rootHOG, the root reasong of the error won't shown up. but you stope here at worst which is late.
            #  todo, check file exist, how to handle if not
    logger_hog.debug("Finding hogs for rhogid_num: " + str(rhogid_num) + ", for taxonomic level:" + str(
        node_species_tree.name) + " for species sub-tree:\n  " + str(node_species_tree.write(format=1, format_root_node=True)) + "\n")
    return hogs_children_level_list




def infer_hogs_this_level(node_species_tree, rhogid_num, pickles_subhog_folder_all):
    """
    infer subHOGs for a rootHOG at a taxanomic level
    """

    this_level_node_name = node_species_tree.name
    pickles_subhog_folder = pickles_subhog_folder_all + "/rhog_" + str(rhogid_num) + "/"
    assert not node_species_tree.is_leaf(), "issue 1235,singletone hogs are treated elsewhere"+ str(rhogid_num)
    pickle_subhog_file = pickles_subhog_folder + str(this_level_node_name) + ".pickle"

    # TODO arrage resume with nextflow and also for when read single_tone pickles
    if _config.inferhog_resume_subhog:
        if os.path.exists(pickle_subhog_file) and os.path.getsize(pickle_subhog_file) > 3:  # 3 bytes
            with open(pickle_subhog_file, 'rb') as handle:
                # todo : do I really need to read the pickle file
                hogs_this_level_list = pickle.load(handle)  #[object class HOG HOG:4027_sub1,len=1,taxono=PSETE]
                if hogs_this_level_list:
                    logger_hog.debug("Level " + str(this_level_node_name) + " with " + str(len(hogs_this_level_list)) + " hogs is read from pickle.")
                    return len(hogs_this_level_list)
                else:
                    logger_hog.debug(" Issue  1238510: the pickle file for single tone is empty " + str(hogs_this_level_list) + " " + str(rhogid_num))

    hogs_children_level_list = read_children_hogs(node_species_tree, rhogid_num, pickles_subhog_folder_all)

    if len(hogs_children_level_list) == 1:
        hogs_this_level_list = hogs_children_level_list
        with open(pickle_subhog_file, 'wb') as handle:
            pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return len(hogs_children_level_list)


    genetree_msa_file_addr = "tree_"+str(rhogid_num)+"_"+str(node_species_tree.name) + ".nwk" # genetrees
    if len(genetree_msa_file_addr) > 245:
        # there is a limitation on length of file name. I want to  keep it consistent ,msa and gene tree names.
        rand_num = random.randint(1, 10000)
        genetree_msa_file_addr = genetree_msa_file_addr[:245] + str(rand_num)+".nwk"

    sub_msa_list_lowerLevel_ready = [hog._msa for hog in hogs_children_level_list if len(hog._msa) > 0]
    # sub_msa_list_lowerLevel_ready = [ii for ii in sub_msa_list_lowerLevel_ready_raw if len(ii) > 0]
    logger_hog.debug("Merging "+str(len(sub_msa_list_lowerLevel_ready))+" MSAs for rhog:"+str(rhogid_num)+", level:"+str(node_species_tree.name))
    msa_filt_row_col = []
    prot_dubious_msa_list = []
    if sub_msa_list_lowerLevel_ready:
        if len(sub_msa_list_lowerLevel_ready) > 1:
            merged_msa = _wrappers.merge_msa(sub_msa_list_lowerLevel_ready, genetree_msa_file_addr)
            if _config.fragment_detection:
                prot_dubious_msa_list, seq_dubious_msa_list = _utils_frag_SO_detection.find_prot_dubious_msa(merged_msa)
        else:
            merged_msa = sub_msa_list_lowerLevel_ready   #  when only on  child, the rest msa is empty.
        logger_hog.debug("All sub-hogs are merged, merged_msa "+str(len(merged_msa))+" "+str(len(merged_msa[0]))+" for rhog: "+str(rhogid_num)+", taxonomic level:"+str(node_species_tree.name))
        (msa_filt_row_col, msa_filt_col, hogs_children_level_list) = _utils_subhog.filter_msa(merged_msa, genetree_msa_file_addr, hogs_children_level_list)
        # msa_filt_col is used for parent level of HOG. msa_filt_row_col is used for gene tree inference.
    else:
        logger_hog.info("Issue 1455, merged_msa is empty " + str(rhogid_num) + ", for taxonomic level:" + str(node_species_tree.name))

    if len(msa_filt_row_col) > 1 and len(msa_filt_row_col[0]) > 1:

        gene_tree_raw = _wrappers.infer_gene_tree(msa_filt_row_col, genetree_msa_file_addr)
        gene_tree = Tree(gene_tree_raw + ";", format=0)   # ,quoted_node_names=True
        logger_hog.debug("Gene tree is inferred len "+str(len(gene_tree))+" rhog:"+str(rhogid_num)+", level: "+str(node_species_tree.name))

        if _config.fragment_detection and len(gene_tree) > 2 and prot_dubious_msa_list:
            (gene_tree, hogs_children_level_list, merged_msa_new) = _utils_frag_SO_detection.handle_fragment_msa(prot_dubious_msa_list, seq_dubious_msa_list, gene_tree, node_species_tree, genetree_msa_file_addr, hogs_children_level_list, merged_msa)
        else:
            merged_msa_new = merged_msa

        # when the prot dubious is removed during trimming
        if len(gene_tree) > 1:
            (gene_tree, all_species_dubious_sd_dic) = _utils_subhog.genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr, hogs_children_level_list)
            if _config.low_so_detection and all_species_dubious_sd_dic:
                (gene_tree, hogs_children_level_list) = _utils_frag_SO_detection.handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic, hogs_children_level_list)

            logger_hog.debug("Merging sub-hogs for rhogid_num:"+str(rhogid_num)+", level:"+str(node_species_tree.name))
            # the last element should be merged_msa not the trimmed msa, as we create new hog based on this msa
            hogs_this_level_list = merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa_new)
            # for i in hogs_this_level_list: print(i.get_members())
            logger_hog.debug("After merging subhogs of childrens, "+str(len(hogs_this_level_list))+" subhogs are found for rhogid_num: "+str(rhogid_num)+", for taxonomic level:"+str(this_level_node_name))

        else:
            hogs_this_level_list = hogs_children_level_list
    else:
        if msa_filt_row_col:
            logger_hog.debug("warning id 13805: hogs_this_level_list is empty. msa_filt_row_col:"+str(len(msa_filt_row_col))+"*"+str(len(msa_filt_row_col[0]))+" !!")
        else:
            logger_hog.debug("warning id 13806: msa_filt_row_col is empty." + str(len(msa_filt_row_col)) +"! ")

        hogs_this_level_list = hogs_children_level_list

    with open(pickle_subhog_file, 'wb') as handle:
        pickle.dump(hogs_this_level_list, handle, protocol=pickle.HIGHEST_PROTOCOL)


    # for hog in hogs_this_level_list:
    #     msa_rec_ids = [i.id for i in hog._msa]
    #     if set(hog._members) != set(msa_rec_ids):
    #         logger_hog.debug("issue 123601, the members are not matching with msa"+str(set(hog._members))+"  "+str(msa_rec_ids))
    # this could be becuase of sub-sampling


    return len(hogs_this_level_list)


def merge_subhogs(gene_tree, hogs_children_level_list, node_species_tree, rhogid_num, merged_msa):
    """
    merge subhogs based on the gene tree specieciaton node of gene tree by creating inter-HOG graph (implicitley )
    """

    subhogs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
    hogs_this_level_list = []
    subHOG_to_be_merged_set_other_Snodes = []
    subHOG_to_be_merged_set_other_Snodes_flattned_temp = []
    ##  the following if fore debugging and visualisation of connected component of inter-HOG graph
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

        if not node.is_leaf() and node.name[0] == "S":
            node_leaves_name_raw = [i.name for i in node.get_leaves()]
            # leaves names  with subhog id  'HALSEN_R15425||HALSEN|_|1352015793_sub10149',
            node_leaves_name = [i.split("|_|")[0] for i in node_leaves_name_raw]
            # node_leaves_name =[]
            # for name_i in node_leaves_name_:
            #     node_leaves_name += name_i.split("_|_")

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
            # todo this piece of code could be neater, for extracting connected component of inter-hog graph
            subHOG_to_be_merged = []
            for node_leave_name in node_leaves_name:  # print(node_leave_name)
                for subHOG in hogs_children_level_list:
                    subHOG_members = subHOG._members
                    if node_leave_name in subHOG_members:  # could be improved
                        if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp:
                            subHOG_to_be_merged.append(subHOG)
                            subhogs_id_children_assigned.append(subHOG._hogid)
                        else:  # this hog is already decided to be merged  print(node.name, subHOG._hogid, node_leave_name)
                            if "processed" in node:
                                logger_hog.info("issue 1863"+ str(node.name)+str(subHOG._hogid)+ str(node_leave_name)) # print("processed", node.name) #else: #    print("processed not in ", node.name)  # print(node_leave_name,"is in ",subHOG._hogid)
            if subHOG_to_be_merged:
                if len(subHOG_to_be_merged) == 1:
                    logger_hog.info("issue 125568313"+str(subHOG_to_be_merged)+" "+node.name)

                subHOG_to_be_merged_set = set(subHOG_to_be_merged)
                taxnomic_range = node_species_tree.name
                num_species_tax_speciestree = len(node_species_tree.get_leaves())
                # num_species_tax   is the number of species exist in the species tree at this clade
                HOG_this_node = HOG(subHOG_to_be_merged_set, taxnomic_range, rhogid_num, merged_msa, num_species_tax_speciestree)
                if len(HOG_this_node._msa) == 1:
                    logger_hog.info("issue 1258313"+str(HOG_this_node)+str(HOG_this_node._msa)+" "+node.name  )
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
    ##  the following if fore debugging and visualisation of connected component of inter-HOG graph
    # fig = plt.figure(figsize=(300, 200), dpi=60)
    # pos = nx.spring_layout(hoggraph, k=0.25, iterations=30)  # For better example looking  # smaller k, biger space between
    # nx.draw(hoggraph, pos, with_labels=True, node_color='y', node_size=500, font_size=16) # , alpha=0.4
    # # nx.draw(G, pos,, edge_color="r", font_size=16, with_labels=True)
    # labels = {e: hoggraph.edges[e]['weight'] for e in hoggraph.edges}
    # nx.draw_networkx_edge_labels(hoggraph, pos, edge_labels=labels, font_size=16)
    # num = random.randint(3, 1000000)
    # plt.savefig("./hoggraph/" + hogs_children_level_list[0]._hogid[4:] + "file_rndm"+str(num)+".jpg")
    # plt.show()
    for subHOG in hogs_children_level_list:  # for the single branch  ( D include subhog and S node. )
        if subHOG._hogid not in subhogs_id_children_assigned:  # print("here", subHOG)
            hogs_this_level_list.append(subHOG)
    # if len(hogs_this_level_list)==1:  hogs_this_level_list = [hogs_this_level_list]

    for hog_j in hogs_this_level_list:
        hog_j._tax_now = node_species_tree.name
    ##  the following if fore debugging of connected component of inter-HOG graph
    ## check for conflicts in merging
    #     for i in range(subHOG_to_be_merged_set_other_Snodes):  if
    #         for i in range(subHOG_to_be_merged_set_other_Snodes):  print("*&*& ",node_species_tree.name)
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
