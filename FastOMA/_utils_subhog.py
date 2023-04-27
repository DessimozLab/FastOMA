
from Bio import SeqIO
from ete3 import Phyloxml
from ete3 import Tree
from ete3 import PhyloTree
from Bio.SeqRecord import SeqRecord

from Bio.Align import MultipleSeqAlignment

from Bio.Seq import Seq  # , UnknownSeq
from collections import defaultdict
from typing import List, Tuple
import random
from itertools import combinations
import numpy as np

import sys

from ._config import logger_hog

from os import listdir
# import pickle
# from xml.dom import minidom
# import xml.etree.ElementTree as ET

from . import _config
from . import _wrappers


"""
fragments in the code mean poorly annotated genes that should be one gene.
"""


def list_rhog_fastas(address_rhogs_folder):
    """
     create  list of rootHOG IDs  stored in the folder of rHOG .
     input: folder address
     output: list of rhog Id (integer)
    """
    rhog_files = listdir(address_rhogs_folder)
    rhogid_num_list = []
    for rhog_file in rhog_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)

    return rhogid_num_list


def read_species_tree_add_internal(species_tree_address):
    """
    reading  species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    # logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
    format_tree = species_tree_address.split(".")[-1]
    # print("there shouldnt be any space in the tree name internal node name as well")
    if format_tree == "phyloxml":
        project = Phyloxml()
        project.build_from_file(species_tree_address)
        # Each tree contains the same methods as orthoxml_to_newick.py PhyloTree object
        for species_tree in project.get_phylogeny():
            species_tree = species_tree
        for node_species_tree in species_tree.traverse(strategy="postorder"):
            temp1 = node_species_tree.phyloxml_clade.get_taxonomy()[0]
            if temp1.get_code():
                node_species_tree.name = temp1.get_code()
            else:
                node_species_tree.name = temp1.get_scientific_name()
        # print(len(species_tree)); print(species_tree)
    elif format_tree == "nwk":
        try:
            species_tree = Tree(species_tree_address, format=1)
        except:
            try:
                species_tree = Tree(species_tree_address)
            except:
                logger_hog.error("format of species tree is not known or the file doesn't exist"+species_tree_address )
                sys.exit()
    else:
        logger_hog.error("for now we accept phyloxml or nwk format for input species tree.or the file doesn't exist "+species_tree_address)
        sys.exit()

    # add name for the internal or leaf, if no name is provided
    num_leaves_no_name = 0
    counter_internal = 0
    for node in species_tree.traverse(strategy="postorder"):
        node_name = node.name
        if len(node_name) < 1:
            if node.is_leaf():
                node.name = "leaf_" + str(num_leaves_no_name)
                num_leaves_no_name += 1
            else:
                node.name = "internal_ad_" + str(counter_internal)
                counter_internal += 1

    return species_tree


def genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr, hogs_children_level_list=[]):

    if _config.rooting_method == "midpoint":
        r_outgroup = gene_tree.get_midpoint_outgroup()
        try:
            gene_tree.set_outgroup(r_outgroup)  # print("Midpoint rooting is done for gene tree.")
        except:
            pass

    elif _config.rooting_method == "mad":
        gene_tree = _wrappers.mad_rooting(genetree_msa_file_addr)
    # elif _config.rooting_method == "outlier":
    #     gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
    #     outliers = find_outlier_leaves(gene_tree)
    #     r_outgroup = midpoint_rooting_outgroup(gene_tree, leaves_to_exclude=outliers)
    #     gene_tree.set_outgroup(r_outgroup)

    all_species_dubious_sd_dic = {}
    if _config.label_SD_internal == "species_overlap":
        (gene_tree, all_species_dubious_sd_dic) = label_sd_internal_nodes(gene_tree)

    elif _config.label_SD_internal == "reconcilation":
        node_species_tree_nwk_string = node_species_tree.write(format=1)
        node_species_tree_PhyloTree = PhyloTree(node_species_tree_nwk_string, format=1)
        gene_tree_nwk_string = gene_tree.write(format=1)
        gene_tree_PhyloTree = PhyloTree(gene_tree_nwk_string, format=1)
        gene_tree = label_SD_internal_nodes_reconcilation(gene_tree_PhyloTree, node_species_tree_PhyloTree)

    if hogs_children_level_list:
        for node in gene_tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node_name_old = node.name
                for hog_child in hogs_children_level_list:
                    if node_name_old in hog_child._members:
                        #node_name_new = node_name_old.split("||")[0]+" "+ hog_child._hogid.split("_")[-1]
                        node_name_new = node_name_old + "|_|" + hog_child._hogid.split("_")[-1]
                        # BUCABY_R15453||BUCABY||1286015722_sub10216
                        node.name = node_name_new
                        break

    if _config.gene_trees_write:
        gene_tree.write(format=1, format_root_node=True, outfile=genetree_msa_file_addr+"_SD_labeled.nwk")

    return gene_tree, all_species_dubious_sd_dic


def read_msa(input_msa):
    ids = []
    seqs = []
    species = []
    coords = []
    for rec in input_msa:
        ids.append(rec.id)
        seqs.append(np.frombuffer(str(rec.seq).upper().encode('ascii'), dtype='S1'))
        species.append(rec.id.split('||')[1][:-1])

        # compute coordinates in matrix
        # todo  if all the element is --, this will be empty and causing error
        ii = np.argwhere(seqs[-1] != b'-')
        coords.append((np.min(ii), np.max(ii)))

    ids = np.array(ids)
    seqs = np.vstack(seqs)
    species = np.array(species)
    coords = np.array(coords)
    input_msa_np = (ids, seqs, species, coords)
    return input_msa_np


def compute_identity(s1, s2):
    def compute(s1, s2):
        # ignore gap columns in first sequence
        f = (s1 != b'-')
        n = f.sum()
        return ((s1[f] == s2[f]).sum() / n) if n > 0 else 1

    return compute(s1, s2), compute(s2, s1)


def split_candidates(input_msa_np, margin=0):


    (ids_all, seqs_all, species_all, coords_all) = input_msa_np
    '''
    Count candidates where it is unambiguous which fragments to merge.
    by Alex Warwick Vesztrocy
    '''
    assert 0 <= margin <= 1, "Margin must be in [0,1]"
    candidates = []

    # get species counts
    c = np.unique(species_all, return_counts=True)
    total_seqs = c[1].sum()

    # filter to those with more than 1 seq
    f = (c[1] > 1)
    for sp in c[0][f]:
        ii = np.argwhere(species_all == sp).flatten()
        # print('species', sp, margin, len(ii), total_seqs)
        x = coords_all[ii]

        order = np.argsort(x[:, 0])
        coords = x[order]
        ids = ids_all[ii[order]]
        seqs = seqs_all[ii[order]]

        #  find whether to merge for unambiguous cases
        merge = True
        ident = []
        for i in range(len(coords) - 1):
            (s1, e1) = coords[i]
            (s2, e2) = coords[i + 1]
            i_margin = min(margin * (e1 - s1),
                           margin * (e2 - s2))
            # s1-------e1
            #       s2------e2
            #  (e1-s2) < 0 no overlap; (e1-s2) > 0 overlap.
            overlap = (e1 - s2)
            # print(overlap)
            if overlap < 0:
                ident.append((1, 1))
            elif overlap <= i_margin:
                # minimal overlap, compute identities
                ident.append(
                    compute_identity(
                        seqs[i][s2:e1],
                        seqs[i + 1][s2:e1],
                    )
                )
            else:
                # overlapping too much
                merge = False
                break
        if merge:
            candidates.append((tuple(ids), ident))
    return candidates


def find_prot_dubious_msa(input_msa):

    input_msa_np = read_msa(input_msa)
    candidates = split_candidates(input_msa_np, _config.overlap_fragments)

    prot_dubious_msa_list = [(str(i[0][0]), str(i[0][1])) for i in candidates]  # [i[0] for i in candidates]
    # todo are there only two always? few fragments
    seq_dubious_msa_list = []
    for prot_pair in prot_dubious_msa_list:
        seq_dubious_msa_pair = [i for i in input_msa if i.name in prot_pair]
        seq_dubious_msa_list.append(seq_dubious_msa_pair)

    return prot_dubious_msa_list,  seq_dubious_msa_list  # list of pairs


def insert_dubious_prots_hog_hierarchy_toleaves(hog_host, fragment_host, fragments_list_nothost):
    for subhog in hog_host._subhogs:

        insert_dubious_prots_hog_hierarchy_toleaves(subhog, fragment_host, fragments_list_nothost)

    result_insersion = hog_host.insert_dubious_prots(fragment_host, fragments_list_nothost)

    return 1


def handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic, hogs_children_level_list):
    #  prot_dubious_sd_list, node_species_tree, genetree_msa_file_addr, hogs_children_level_list

    logger_hog.debug("These are  found after removing with msa , all_species_dubious_sd_dic " + str(all_species_dubious_sd_dic)+" which are now being handled.")
    prot_dubious_sd_remove_list = find_prot_dubious_sd_remove(gene_tree, all_species_dubious_sd_dic)

    if prot_dubious_sd_remove_list:
        rest_leaves = set([i.name for i in gene_tree.get_leaves()]) - set(prot_dubious_sd_remove_list)
        gene_tree.prune(rest_leaves, preserve_branch_length=True)
        (gene_tree, all_species_dubious_sd_dic2) = genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr + "_dubious_sd")
        if all_species_dubious_sd_dic2:
            logger_hog.debug("these are found after removing with sd, all_species_dubious_sd_dic2 " + str(all_species_dubious_sd_dic2))
            prot_dubious_sd_remove_list2 = find_prot_dubious_sd_remove(gene_tree, all_species_dubious_sd_dic2)
            if prot_dubious_sd_remove_list2:
                rest_leaves2 = set([i.name for i in gene_tree.get_leaves()]) - set(prot_dubious_sd_remove_list2)
                gene_tree.prune(rest_leaves2, preserve_branch_length=True)
                (gene_tree, all_species_dubious_sd_dic3) = genetree_sd(node_species_tree, gene_tree,genetree_msa_file_addr + "_dubious_sd_2")
                if all_species_dubious_sd_dic3:
                    # todo make it as while to do it for all possible iteration, but there won't many cases for this at least in QFO dataset
                    logger_hog.debug( "issue 13954,these are found after removing with sd two times , all_species_dubious_sd_dic3 " + str(all_species_dubious_sd_dic2))

        hogs_children_level_list_raw = hogs_children_level_list
        logger_hog.debug("** we removed theses sequences"+str(prot_dubious_sd_remove_list))
        for prot_dubious_sd_remove in prot_dubious_sd_remove_list:
            for hog in hogs_children_level_list_raw:
                if prot_dubious_sd_remove in hog._members:
                    result_removing = remove_prots_hog_hierarchy_toleaves(hog, [prot_dubious_sd_remove])
                    if result_removing == 0:  # the hog is empty
                        hogs_children_level_list.remove(hog)

    return (gene_tree, hogs_children_level_list)


def merge_prots_name_hierarchy_toleaves(hog_host, fragment_name_host, fragment_name_remove):

    # todo do I really need this ?
    for subhog in hog_host._subhogs:
        if fragment_name_host in subhog._members:
            merge_prots_name_hierarchy_toleaves(subhog, fragment_name_host, fragment_name_remove)

    result_merging = hog_host.merge_prots_name_hog(fragment_name_host, fragment_name_remove)

    return hog_host


def merge_fragments_hogclass(fragments_set, seq_dubious_msa, hogs_children_level_list, merged_msa):
    # a bit of redundency fragments_set is the names which also are stored in  seq_dubious_msa

    fragments_list = list(fragments_set)
    fragment_name_host = fragments_list[0]
    for hog in hogs_children_level_list:
        if fragment_name_host in hog._members:
            hog_host = hog

    fragments_list_remove = fragments_list[1:]
    for fragment_idx in range(1, len(fragments_list)):  # the 0 element is the host

        fragment_name_remove = fragments_list[fragment_idx]
        # fragment_seq_remove = seq_dubious_msa[fragment_idx]
        # remove
        for hog in hogs_children_level_list:
            if fragment_name_remove in hog._members:
                result_removing = remove_prots_hog_hierarchy_toleaves(hog, [fragment_name_remove])
                if result_removing == 0:  # the hog is empty
                    hogs_children_level_list.remove(hog)
        # merge
        assert len(seq_dubious_msa) == 2


        seq0 = str(seq_dubious_msa[0].seq)
        seq1 = str(seq_dubious_msa[1].seq)
        assert len(seq0) == len(seq1)
        merged_sequence = ""
        for aa_idx in range(len(seq1)):
            if seq0[aa_idx] == '-':
                merged_sequence += seq1[aa_idx]
            elif seq1[aa_idx] == '-':
                merged_sequence += seq0[aa_idx]
            elif seq0[aa_idx] == seq1[aa_idx]:   #  and seq0 != '-' and seq1 != '-'
                merged_sequence += seq0[aa_idx]
            else:
                merged_sequence += "X"
        assert len(merged_sequence) == len(seq0)

        result_merging1 = merge_prots_name_hierarchy_toleaves(hog_host, fragment_name_host, fragment_name_remove)
        hog_host.merge_prots_msa(fragment_name_host, fragment_name_remove, merged_sequence)


        merged_fragment_name = fragment_name_host + "_|_" + fragment_name_remove
        msa_filt_row_col_fragmerged_list =[]
        for seq_rec in merged_msa:
            if seq_rec.id == fragment_name_host:
                seq_rec_edited = SeqRecord(Seq(merged_sequence), id=merged_fragment_name, name=merged_fragment_name)
                msa_filt_row_col_fragmerged_list.append(seq_rec_edited)
            elif seq_rec.id == fragment_name_remove:
                pass
            else:
                msa_filt_row_col_fragmerged_list.append(seq_rec)

        msa_filt_row_col_fragmerged= MultipleSeqAlignment(msa_filt_row_col_fragmerged_list)
        logger_hog.debug("these proteins fragments are merged  into one "+str(merged_fragment_name)+"but reported seperately in orthoxml")

    return fragments_list_remove, hogs_children_level_list, msa_filt_row_col_fragmerged



def handle_fragment_msa(prot_dubious_msa_list, seq_dubious_msa_list, gene_tree, node_species_tree, genetree_msa_file_addr, hogs_children_level_list, merged_msa):
    merged_msa_new = merged_msa
    if not prot_dubious_msa_list: # empty list
        return gene_tree, hogs_children_level_list, merged_msa_new

    logger_hog.debug("** these are found prot_dubious_msa_list " + str(prot_dubious_msa_list))
    fragments_set_list = check_prot_dubious_msa(prot_dubious_msa_list, gene_tree)
    fragments_remove_list = [] # this list include only one of the fragments for each set, the other one is merged version afterwards
    if fragments_set_list and len(fragments_set_list[0]) > 1:
        # logger_hog.debug("** these are found fragments_set_list " + str(fragments_set_list))
        # remove fragments from gene tree
        if _config.merge_fragments_detected:
            for fragments_set_idx, fragments_set in enumerate(fragments_set_list):
                fragments_set = fragments_set_list[fragments_set_idx]
                seq_dubious_msa = seq_dubious_msa_list[fragments_set_idx]
                fragments_list_remove, hogs_children_level_list, merged_msa_new  = merge_fragments_hogclass(fragments_set, seq_dubious_msa, hogs_children_level_list, merged_msa)
                # merged_msa_new  contains the merged_seq
                # note that merged_msa is the merging of different subhog childrend and is not related to merge fragments  of two seqeuncs with bad annotation

                # todo I may need to trim the msa

                #
                gene_tree_raw = _wrappers.infer_gene_tree(merged_msa_new, genetree_msa_file_addr+"_merged")
                gene_tree = Tree(gene_tree_raw + ";", format=0)
                # fragments_remove_list += fragments_list_remove # for now fragments_list_remove include 1 prots
        else:
            fragments_remove_set = set.union(*fragments_set_list)
            rest_leaves = set([i.name for i in gene_tree.get_leaves()]) - fragments_remove_set
            if len(rest_leaves) < 2:
                # todo
                print("** issue 86194")
                hogs_children_level_list = []
                gene_tree =""
                return gene_tree, hogs_children_level_list, merged_msa_new
            else:
                gene_tree.prune(rest_leaves, preserve_branch_length=True)

        (gene_tree, all_species_dubious_sd_dic2) = genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr+"_dubiousMSA")
        if all_species_dubious_sd_dic2:
            # logger_hog.debug("these are  found after removing with msa , all_species_dubious_sd_dic2 "+str(all_species_dubious_sd_dic2))
            (gene_tree, hogs_children_level_list) = handle_fragment_sd(node_species_tree, gene_tree, genetree_msa_file_addr, all_species_dubious_sd_dic2, hogs_children_level_list)

            for fragments_set in fragments_set_list:
                fragments_list = list(fragments_set)
                fragment_host = fragments_list[0]  # host of the new small hog consisting few fragments
                for hog in hogs_children_level_list:
                    if fragment_host in hog._members:
                        hog_host = hog
                        break  # once we found, we don't need to continue searching in hog
                fragments_list_nothost = fragments_list[1:]
                for fragment in fragments_list_nothost:
                    hogs_children_level_list_raw = hogs_children_level_list
                    for hog in hogs_children_level_list_raw:
                        if fragment in hog._members:
                            result_removing = remove_prots_hog_hierarchy_toleaves(hog, [fragment])
                            if result_removing == 0:  # the hog is empty
                                hogs_children_level_list.remove(hog)
                                # print(hogs_children_level_list)
                insert_dubious_prots_hog_hierarchy_toleaves(hog_host, fragment_host, fragments_list_nothost)

    return gene_tree, hogs_children_level_list, merged_msa_new

def check_prot_dubious_msa(prot_dubious_msa_list, gene_tree):

    farthest, max_dist_numNodes = gene_tree.get_farthest_node(topology_only=True)  # furthest from the node
    farthest, max_dist_length = gene_tree.get_farthest_node()  # furthest from the node
    print("max_dist_numNodes, max_dist_length ", max_dist_numNodes, max_dist_length)
    fragments_set_list = []
    gene_tree_leaves_name = set([i.name for i in gene_tree.get_leaves()])
    for prot_dubious_msa_set in prot_dubious_msa_list:
        print(prot_dubious_msa_set)
        fragments = []

        for prot in prot_dubious_msa_set[1:]:
            # todo following could be imporved, during filtering row/col msa, a fragments could be removed and not in gene tree anymore,
            if prot_dubious_msa_set[0] in gene_tree_leaves_name and prot in gene_tree_leaves_name:
                # there might be few fragments, checking the distance of the first one with the rest # todo check all vs all

                # todo substract two terminal branches
                assert len(gene_tree.get_leaves_by_name(prot_dubious_msa_set[0])) == 1,\
                    "the prot name is not in gene tree or there are more than one" + str(prot_dubious_msa_set[0])
                assert len(gene_tree.get_leaves_by_name(prot)) == 1,\
                    "the prot name is not in gene tree or there are more than one" + str(prot)

                dist_numNodes = gene_tree.get_distance(prot_dubious_msa_set[0], prot, topology_only=True)
                dist_length = gene_tree.get_distance(prot_dubious_msa_set[0], prot)
                node_prot = gene_tree.get_leaves_by_name(prot)[0]
                node_prot_dubious = gene_tree.get_leaves_by_name(prot_dubious_msa_set[0])[0]


                dist_length_corrected = dist_length - abs(node_prot.dist- node_prot_dubious.dist)



                print("check_prot_dubious_msa dist_numNodes, dist_length ",dist_numNodes, dist_length)
                if dist_length_corrected < max(0.005, max_dist_length / 5):   # or (dist_length_corrected - 2*max_dist_length)< 0.001
                    # dist_numNodes < max(max_dist_numNodes * 1 / 5, 3) or
                    fragments += [prot_dubious_msa_set[0], prot]
        if fragments:
            fragments_set_list.append(set(fragments))

    return fragments_set_list


def prepare_species_tree(rhog_i, species_tree, rhogid_num):
    """
    orthoxml_to_newick.py function for extracting orthoxml_to_newick.py subtree from the input species tree  orthoxml_to_newick.py.k.orthoxml_to_newick.py pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    assert len(rhog_i) > 0, 'input hog_i is empty, probably previous step find_rhog has issue, rhogs/HOG_B0'+str(rhogid_num)+'is empty?'
    species_names_rhog = []
    prot_names_rhog = []
    for rec in rhog_i:
        # qfo : >tr|A0A0N7KF21|A0A0N7KF21_ORYSJ||ORYSJ_||1000000344 tr|A0A0N7KF21|A0A0N7KF21_ORYSJ Os02g0264501
        # protein OS=Oryza sativa subsp. japonica (Rice) OX=39947 GN=Os02g0264501 PE=4 SV=1
        prot_id = rec.id.split("||")
        prot_name = prot_id[2]   # for debugging  prot_id[0] readable prot name,  for xml prot_id[2]
        species_name = prot_id[1]
        bird_dataset = True
        if species_name.endswith("_") and not bird_dataset:
           species_name = prot_id[1][:-1]
        # if species_name == 'RAT_': species_name = "RATNO_"
        # gene_id = prot_id[2]
        species_names_rhog.append(species_name)
        prot_names_rhog.append(prot_name)
    assert len(species_names_rhog) > 0, "species names list is empty in rhog, probably issue in formating with || in previous step find rhog"

    species_names_uniqe = set(species_names_rhog)

    first_common_ancestor_name = species_tree.get_common_ancestor(species_names_uniqe).name
    species_tree.prune(species_names_uniqe, preserve_branch_length=True)
    species_tree.name = first_common_ancestor_name
    # add internal node name to the tree
    # this has an issue with root name, cannot add the root name
    # print(species_tree.write(format=1, format_root_node=True))
    # counter_internal = 0
    # for node in species_tree.traverse(strategy="postorder"):
    #     node_name = node.name
    #     num_leaves_no_name = 0
    #     if len(node_name) < 1:
    #         if node.is_leaf():
    #             node.name = "leaf_" + str(num_leaves_no_name)
    #         else:
    #             node_children = node.children
    #             # list_children_names = [str(node_child.name) for node_child in node_children]
    #             # node.name = '_'.join(list_children_names)
    #             # ?? to imrpove, if the species tree has internal node name, keep it,
    #             # then checn condition in  _infer_subhog.py, where logger_hog.info("Finding hogs for rhogid_num: "+str(rh
    #             node.name = "internal_" + str(counter_internal)  #  +"_rhg"+str(rhogid_num)  #  for debuging
    #             counter_internal += 1
    # print("Working on the following species tree.")
    # print(species_tree.write(format=1, format_root_node=True))

    return species_tree, species_names_rhog, prot_names_rhog


def find_prot_dubious_sd_remove(gene_tree, all_species_dubious_sd_dic):
    # todo this function need to double check with cases of with and without dubious

    #prot_dubious_sd_allspecies = []
    prot_dubious_sd_remove_list = []
    # todo not sure postorder or preorder
    for node in gene_tree.traverse(strategy="postorder"):
        # print("** now working on node ",node.name) # node_children
        if not node.is_leaf() and 'D' in node.name:
            node_name = node.name #d, intersection, union = node_name.split("_")  # if int(intersection) / int(union) < _config.threshold_dubious_sd:
            if node_name in all_species_dubious_sd_dic: # a duplication node with low score,
                node_children = node.children
                all_species_dubious_sd = all_species_dubious_sd_dic[node_name]
                # prot_dubious_sd_ = []
                for species_dubious_sd in all_species_dubious_sd:
                    child_size = []  # gene tree is binary for fasttree
                    prot_dubious_list = []
                    for node_child in node_children:
                        list_leaves = [i.name for i in node_child.get_leaves()]
                        child_size.append(len(list_leaves))
                        for prot_name in list_leaves:
                            if prot_name.split("||")[1] == species_dubious_sd:
                                prot_dubious_list.append(prot_name)
                    subhogs_list = [i.split("|_|")[1] for i in prot_dubious_list] # subhog id at child level
                    if len(set(subhogs_list)) > 1:
                        # we are removing all sequences of this species on the the side of internal node (gene tree), with least leaves
                        child_size_min_indx = child_size.index(min(child_size))
                        prot_dubious_sd_remove_list.append(prot_dubious_list[child_size_min_indx])

                    else:
                        logger_hog.debug( "This species (protein from the same subhog) is safe to keep "+ str(node_name)+" "+str(species_dubious_sd))
                        #all of them are from the same subhog, so it doesn't matter, a duplication event doesn't affect when all are from the same subhog at children level

    return prot_dubious_sd_remove_list



def label_sd_internal_nodes(tree_out):
    """
    for the input gene tree, run the species overlap method
    and label internal nodes of the gene tree

    output: labeled gene tree
    """
    species_name_dic = {}
    counter_S = 0
    counter_D = 0

    all_species_dubious_sd_dic = {}
    for node in tree_out.traverse(strategy="postorder"):
        # print("** now working on node ",node.name) # node_children
        if node.is_leaf():
            prot_i = node.name
            # species_name_dic[node] = {str(prot_i).split("|")[-1].split("_")[-1]}
            species_name_dic[node] = {str(prot_i).split("||")[1]}
        else:
            node.name = "S/D"
            leaves_list = node.get_leaves()  # print("leaves_list", leaves_list)
            # species_name_set = set([str(prot_i).split("|")[-1].split("_")[-1] for prot_i in leaves_list])
            species_name_set = set([str(prot_i).split("||")[1] for prot_i in leaves_list])
            species_name_dic[node] = species_name_set
            node_children = node.children  # print(node_children)
            node_children_species_list = [species_name_dic[node_child] for node_child in node_children]  # list of sets
            node_children_species_intersection = set.intersection(*node_children_species_list)  # * is for handling list of sets
            node_children_species_union = set.union(*node_children_species_list)

            if node_children_species_intersection:  # print("node_children_species_list",node_children_species_list)
                counter_D += 1
                node.name = "D" + str(counter_D) + "_"+str(len(node_children_species_intersection))+"_"+str(len(node_children_species_union))
                if len(node_children_species_intersection)/ len(node_children_species_union) < _config.threshold_dubious_sd:
                    all_species_dubious_sd_dic[node.name] = list(node_children_species_intersection)
            else:
                counter_S += 1
                node.name = "S" + str(counter_S)
    return tree_out, all_species_dubious_sd_dic






def label_SD_internal_nodes_reconcilation(gene_tree, species_tree):
    """
    for the input gene tree, run the gene/species tree reconciliation method
    and label internal nodes of the gene tree

    output: labeled gene tree
    """
    try:
        gene_tree_reconciled = get_reconciled_tree_zmasek(gene_tree, species_tree, inplace=False)
    except ValueError:
        print("### Algorithm can only work with binary trees. Force resolved polytomies.")
        gene_tree.resolve_polytomy()
        gene_tree_reconciled = get_reconciled_tree_zmasek(gene_tree, species_tree, inplace=False)
    for node in gene_tree_reconciled.traverse(strategy="postorder"):
        if node.is_leaf():
            pass
        else:
            node.name = node.evoltype
    return gene_tree_reconciled


def get_reconciled_tree_zmasek(gtree, sptree, inplace=False):
    """
    from ete3
    https://github.com/etetoolkit/ete/blob/1f587a315f3c61140e3bdbe697e3e86eda6d2eca/ete3/phylo/reconciliation.py

    Reconciles the gene tree with the species tree
    using Zmasek and Eddy's algorithm. Details can be
    found in the paper:
    Christian M. Zmasek, Sean R. Eddy: A simple algorithm
    to infer gene duplication and speciation events on a
    gene tree. Bioinformatics 17(9): 821-828 (2001)
    :argument gtree: gene tree (PhyloTree instance)
    :argument sptree: species tree (PhyloTree instance)
    :argument False inplace: if True, the provided gene tree instance is
       modified. Otherwise a reconciled copy of the gene tree is returned.
    :returns: reconciled gene tree
    """
    # some cleanup operations
    def cleanup(tree):
        for node in tree.traverse():
            node.del_feature("M")

    if not inplace:
        gtree = gtree.copy('deepcopy')

    # check for missing species
    g_node_species_all = []
    for g_node in gtree.get_leaves():
        # # leaves names  with subhog id  'HALSEN_R15425||HALSEN||1352015793||sub10149'
        g_node_species_all.append(g_node.name.split("||")[1])
    species_sptree_all = [i.name for i in sptree.get_leaves()]
    missing_sp = set(g_node_species_all) - set(species_sptree_all)
    if missing_sp:
        raise KeyError("* The following species are not contained in the species tree: "+ ', '.join(missing_sp))

    # initialization
    sp2node = dict()
    for sp_node in sptree.get_leaves():
        sp2node[sp_node.name] = sp_node

    # set/compute the mapping function M(g) for the
    # leaf nodes in the gene tree (see paper for details)
    species = [i.name for i in sptree.get_leaves()]   #sptree.get_species()
    for g_node in gtree.get_leaves():
        g_node_species = g_node.name.split("||")[1]
        g_node.add_feature("M", sp2node[g_node_species])

    # visit each internal node in the gene tree
    # and detect its event (duplication or speciation)
    for node in gtree.traverse(strategy="postorder"):
        if len(node.children) == 0:
            continue # nothing to do for leaf nodes

        if len(node.children) != 2:
            cleanup(gtree)
            raise ValueError("Algorithm can only work with binary trees.")

        lca = node.children[0].M.get_common_ancestor(node.children[1].M) # LCA in the species tree
        node.add_feature("M", lca)

        node.add_feature("evoltype", "S")
        #node.name = "S"
        if id(node.children[0].M) == id(node.M) or id(node.children[1].M) == id(node.M):
                node.evoltype = "D"
                #node.name = "D"

    cleanup(gtree)
    return gtree



def remove_prots_hog_hierarchy_toleaves(hog_ii, prots_to_remove):

    for subhog in hog_ii._subhogs:
        remove_prots_hog_hierarchy_toleaves(subhog, prots_to_remove)

    result_removing = hog_ii.remove_prots_from_hog(prots_to_remove)

    # todo not sure following is needed, as we are removing in next line
    # if result_removing == 0:  # the hog is empty
    #     del hog_ii   # to remove an object of python class

    # remove inside subhog if is empty
    hog_ii._subhogs = [i for i in hog_ii._subhogs if len(i._members) > 0]
    # print(hog._taxnomic_range)
    #if list(prots_to_remove)[0] in hog._members:
    #    print(hog._members, hog._taxnomic_range)
    return result_removing



def msa_filter_col(msa, tresh_ratio_gap_col, gene_tree_file_addr=""):
    # gene_tree_file_addr contains roothog numebr
    # note this is used in hog class as well

    ratio_col_all = []
    length_record = len(msa[0])
    num_records = len(msa)
    keep_cols = []
    for col_i in range(length_record):
        col_values = [record.seq[col_i] for record in msa]
        gap_count=col_values.count("-") + col_values.count("?") + col_values.count(".") +col_values.count("~")
        ratio_col_nongap = 1- gap_count/num_records
        ratio_col_all.append(ratio_col_nongap)
        if ratio_col_nongap > tresh_ratio_gap_col:
            keep_cols.append(col_i)
    #plt.hist(ratio_col_all,bins=100) # , bins=10
    #plt.show()
    #plt.savefig(gene_tree_file_addr+ "filtered_row_"+"_col_"+str(tresh_ratio_gap_col)+".txt.pdf")
    #print("- Columns indecis extracted. Out of ", length_record,"columns,",len(keep_cols),"is remained.")
    msa_filtered_col = []
    for record in msa:
        record_seq = str(record.seq)
        record_seq_edited = ''.join([record_seq[i] for i in keep_cols  ])
        record_edited = SeqRecord(Seq(record_seq_edited), record.id, '', '')
        msa_filtered_col.append(record_edited)

    if _config.msa_write_all and gene_tree_file_addr:
        out_name_msa=gene_tree_file_addr+"filtered_"+"col_"+str(tresh_ratio_gap_col)+".msa.fa"
        handle_msa_fasta = open(out_name_msa, "w")
        SeqIO.write(msa_filtered_col, handle_msa_fasta, "fasta")
        handle_msa_fasta.close()
    # print("- Column-wise filtering of MSA is finished",len(msa_filtered_col),len(msa_filtered_col[0]))
    return msa_filtered_col


def msa_filter_row(msa, tresh_ratio_gap_row, gene_tree_file_addr=""):
    msa_filtered_row = []
    ratio_records=[]
    for record in msa:
        seq = record.seq
        seqLen = len(record)
        gap_count = seq.count("-") + seq.count("?") + seq.count(".") +seq.count("~")
        if seqLen:
            ratio_record_nongap = 1-gap_count/seqLen
            ratio_records.append(round(ratio_record_nongap, 3))
            if ratio_record_nongap > tresh_ratio_gap_row:
                msa_filtered_row.append(record)
        else:
            print("issue 12788 : error , seq len is zero when msa_filter_row")
    if _config.msa_write_all and gene_tree_file_addr:
        out_name_msa = gene_tree_file_addr +"_filtered_row_"+str(tresh_ratio_gap_row)+".msa.fa"
        handle_msa_fasta = open(out_name_msa, "w")
        SeqIO.write(msa_filtered_row, handle_msa_fasta, "fasta")
        handle_msa_fasta.close()
    return msa_filtered_row


# Fragment detection using MSA
#
# def fragment_detector_candidate(merged_msa):
#
#     rec_group_species = {}
#     for seq in merged_msa:
#         species = seq.id.split("||")[1]
#         if species in rec_group_species:
#             rec_group_species[species].append(seq)
#         else:
#             rec_group_species[species]= [seq]
#
#     rec_candidate = {}
#
#     len_aligned = len(seq)
#
#     for species, list_seq in rec_group_species.items():
#         if len(list_seq)>1:
#             num_nongap_list = []
#             for i in range(len(list_seq)):
#                 seq_i= list_seq[i] # seq_i is biopython record
#                 num_nongap_i= len_aligned - seq_i.count("-")
#                 num_nongap_list.append(num_nongap_i)
#                 if num_nongap_i > len_aligned *0.25 and  num_nongap_i < len_aligned *0.75:
#                     for j in range(i):
#                         seq_j= list_seq[j]
#                         num_nongap_j= num_nongap_list[j]
#                         if num_nongap_j > len_aligned *0.25 and  num_nongap_j < len_aligned *0.75:
#                             count_gap_aa = 0
#                             for (chr_i, chr_j) in zip(seq_i, seq_j):
#                                 if (chr_i=='-' and chr_j!='-')  or (chr_i!='-' and chr_j=='-'):
#                                     count_gap_aa +=1
#                             # print(count_gap_aa)
#                             if count_gap_aa  > len_aligned * 0.25: # two seq complment each other
#                                 # TODO the downside is sth like this: seq1=-A-A seq2= A-A-  not fragments
#                                 if species in rec_candidate:
#                                     seq_i_id = seq_i.id
#                                     seq_j_id = seq_j.id
#
#                                     if seq_i_id not in rec_candidate[species]:
#                                         rec_candidate[species] += seq_i_id
#                                     if seq_j_id not in rec_candidate[species]:
#                                         rec_candidate[species] += seq_j_id
#                                 else:
#                                     rec_candidate[species] = [seq_i_id, seq_j_id]  # seq_i is biopython record
#
#     return rec_candidate



def filter_msa(merged_msa, gene_tree_file_addr, hogs_children_level_list):

    msa_filt_row_1 = merged_msa
    # if _config.inferhog_filter_all_msas_row:
    #    msa_filt_row_1 = _utils_subhog.msa_filter_row(merged_msa, _config.inferhog_tresh_ratio_gap_row, gene_tree_file_addr)
    # if   msa_filt_row_1 and len(msa_filt_row_1[0]) >=
    if len(msa_filt_row_1[0]) >= _config.inferhog_min_cols_msa_to_filter:
        # (len(merged_msa) > 10000 and len(merged_msa[0]) > 3000) or (len(merged_msa) > 500 and len(merged_msa[0]) > 5000) or (len(merged_msa) > 200 and len(merged_msa[0]) > 9000):
        # for very big MSA, gene tree is slow. if it is full of gaps, let's trim the msa.
        # logger_hog.debug( "We are doing MSA trimming " + str(rhogid_num) + ", for taxonomic level:" + str(node_species_tree.name))
        # print(len(merged_msa), len(merged_msa[0]))

        if _config.automated_trimAL:
            msa_filt_col = msa_filt_row_1
            msa_filt_row_col = _wrappers.trim_msa(msa_filt_row_1)
        else:
            msa_filt_col = msa_filter_col(msa_filt_row_1, _config.inferhog_tresh_ratio_gap_col, gene_tree_file_addr)
            msa_filt_row_col = msa_filt_col
            if msa_filt_col and msa_filt_col[0] and len(msa_filt_col[0]):
                msa_filt_row_col = msa_filter_row(msa_filt_col, _config.inferhog_tresh_ratio_gap_row, gene_tree_file_addr)

        # compare msa_filt_row_col and msa_filt_col,
        if len(msa_filt_row_col) != len(msa_filt_col):  # some sequences are removed
            set_prot_before = set([i.id for i in msa_filt_col])
            set_prot_after = set([i.id for i in msa_filt_row_col])
            prots_to_remove_level = set_prot_before - set_prot_after
            assert len(prots_to_remove_level), "issue 31235"
            #prots_to_remove |= prots_to_remove_level
            # remove prot from all subhogs
            # hogs_children_level_list_raw = hogs_children_level_list
            # hogs_children_level_list = []
            # for hog_i in hogs_children_level_list_raw:
            #     result_removal = hog_i.remove_prots_from_hog(prots_to_remove)
            #     if result_removal != 0:
            #         hogs_children_level_list.append(hog_i)

        else:
            msa_filt_row_col = msa_filt_col
        merged_msa_filt = msa_filt_row_col
    else:
        msa_filt_row_col = msa_filt_row_1
        msa_filt_col = msa_filt_row_1
        # the msa may be empty
    # if len(msa_filt_row_col) < 2:
    #    msa_filt_row_col = msa_filt_col[:2]

    return (msa_filt_row_col, msa_filt_col, hogs_children_level_list)


class PhyloTree:
    pass

def get_farthest_leaf(tree: PhyloTree, target_leaf: PhyloTree, leaves_list: List[PhyloTree]) -> Tuple[float, PhyloTree]:
    """
    Given a target leaf it returns the farthest leave to it from a given list of leaves
    by Ali Yazdizadeh Kharrazi
    """
    max_dist = 0
    max_leaf: PhyloTree = None
    for leaf in leaves_list:
        if leaf == target_leaf:
            continue

        dist = tree.get_distance(target_leaf, leaf)
        if dist >= max_dist:
            max_dist = dist
            max_leaf = leaf

    print('&', max_leaf, max_dist)
    return max_dist, max_leaf


def midpoint_rooting_longest_path(tree: PhyloTree, leaves_to_exclude=None) -> Tuple[float, PhyloTree, PhyloTree]:
    """
    given a gene tree and optionally a list of leaves to exclude it find the two
    furthest leave in the tree to be used for midpoint rooting
    by Ali Yazdizadeh Kharrazi
    """
    leaves_list = tree.get_leaves()

    if leaves_to_exclude:
        leaves_set = set(leaves_list)
        leaves_set -= set(leaves_to_exclude)
        leaves_list = list(leaves_set)

    random_leaf = random.choice(leaves_list)
    print('&&2', tree, leaves_list, leaves_to_exclude, random_leaf)
    _, first_leaf = get_farthest_leaf(tree=tree, target_leaf=random_leaf, leaves_list=leaves_list)
    _, second_leaf = get_farthest_leaf(tree=tree, target_leaf=first_leaf, leaves_list=leaves_list)
    longest_dist = tree.get_distance(first_leaf, second_leaf)

    print('&&', random_leaf, first_leaf, second_leaf, longest_dist)
    return longest_dist, first_leaf, second_leaf


def midpoint_rooting_outgroup(tree: PhyloTree, leaves_to_exclude=None) -> PhyloTree:
    """
    Using midpoint rooting algorithm find the outgroup to be used to root the tree.
    you can provide a list of leaves to be ignored for example because of long branch
    by Ali Yazdizadeh Kharrazi
    """
    distance, first_leaf, second_leaf = midpoint_rooting_longest_path(tree, leaves_to_exclude)
    distance_first = tree.get_distance(first_leaf)
    distance_second = tree.get_distance(second_leaf)

    farther_node = first_leaf if distance_first > distance_second else second_leaf

    current_distance = 0
    current_node = farther_node
    while current_distance + current_node.dist < distance / 2:
        current_distance += current_node.dist
        current_node = current_node.up

    print('&&&', distance, first_leaf, second_leaf, distance_first, distance_second, farther_node, current_node,
          current_distance)
    return current_node


def find_outlier_leaves(tree: PhyloTree):
    """
    by Ali Yazdizadeh Kharrazi
    """
    distances = defaultdict(list)

    leaves_name = tree.get_leaves()
    for i, j in combinations(leaves_name, 2):
        distances[i].append(tree.get_distance(i, j))
        distances[j].append(tree.get_distance(j, i))

    distances_agg = []
    for leaf in distances.keys():
        distances_agg.append(sum(distances[leaf]))
    q3, q1 = np.percentile(distances_agg, [75, 25])
    iqr = q3 - q1
    threshold = q3 + (1.5 * iqr)

    outliers = []
    for leaf in distances.keys():
        if sum(distances[leaf]) > threshold:
            outliers.append(leaf)
    print('+++', distances_agg, q3, q1, iqr, threshold, outliers)
    return outliers