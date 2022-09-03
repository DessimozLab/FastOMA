import dill as pickle
import logging
from datetime import datetime
from os import listdir
import os
import sys
from typing import Tuple, List
from random import sample

from ete3 import Phyloxml
from ete3 import PhyloTree

import zoo.wrappers.aligners.mafft as mafft
import zoo.wrappers.treebuilders.fasttree as fasttree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import xml.etree.ElementTree as ET
import itertools

import gc
import numpy as np
from ete3.coretype.tree import TreeError

# from utils.midpoint import find_outlier_leaves, midpoint_rooting_outgroup


## the following are needed when we start from orthoxml_to_newick.py rootHOG fasta file.

def read_species_tree(species_tree_address):
    """
    reading orthoxml_to_newick.py species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
    project = Phyloxml()
    project.build_from_file(species_tree_address)
    # Each tree contains the same methods as orthoxml_to_newick.py PhyloTree object
    for species_tree in project.get_phylogeny():
        species_tree = species_tree

    for node_species_tree in species_tree.traverse(strategy="postorder"):
        if node_species_tree.is_leaf():
            temp1 = node_species_tree.phyloxml_clade.get_taxonomy()[0]
            # print(temp1.get_code())
            node_species_tree.name = temp1.get_code()
    # print(len(species_tree)); print(species_tree)
    return (species_tree)


def prepare_species_tree(rhog_i, species_tree):
    """
    orthoxml_to_newick.py function for extracting orthoxml_to_newick.py subtree from the input species tree  orthoxml_to_newick.py.k.orthoxml_to_newick.py pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    species_names_rhog = []
    prot_names_rhog = []
    for rec in rhog_i:
        prot_name = rec.name  # 'tr|E3JPS4|E3JPS4_PUCGT
        # prot_name = prot_name_full.split("|")[1].strip() # # 'tr|E3JPS4|E3JPS4_PUCGT
        species_name = prot_name.split("||")[1][:-1]  # .split("|")[-1].split("_")[-1]
        if species_name == 'RAT': species_name = "RATNO"
        species_names_rhog.append(species_name)
        prot_names_rhog.append(prot_name)

    species_names_uniqe = set(species_names_rhog)
    logger_hog.info("The number of unique species in the rHOG is " + str(len(species_names_uniqe)) + ".")
    species_tree.prune(species_names_uniqe, preserve_branch_length=True)
    # species_tree.write()
    for node in species_tree.traverse(strategy="postorder"):
        node_name = node.name
        if len(node_name) < 1:
            if node.is_leaf():
                num_leaves_no_name = 0
                node.name = "leaf_" + str(num_leaves_no_name)
            else:
                node_children = node.children
                list_children_names = [node_child.name for node_child in node_children]
                node.name = '_'.join(list_children_names)
    print("Working on the following species tree.")
    print(species_tree)

    return (species_tree, species_names_rhog, prot_names_rhog)


def merge_msa(list_msas):
    """
    merge orthoxml_to_newick.py list of MSAs (multiple sequnce aligmnet)
    by run mafft on them.
    Each element of msa should be orthoxml_to_newick.py MultipleSeqAlignment object.

    output: merged (msa)
    """
    logging.debug(list_msas)
    logging.debug(str(list_msas[0][0].id) + "\n")
    wrapper_mafft_merge = mafft.Mafft(list_msas, datatype="PROTEIN")
    wrapper_mafft_merge.options['--merge'].active = True
    merged = wrapper_mafft_merge()
    logger_hog.info(
        str(len(list_msas)) + " msas are merged into one with the length of " + str(len(merged)) + " " + str(
            len(merged[0])))
    return merged


def infer_gene_tree(msa, gene_tree_file_addr):
    """
    infere gene tree using fastTree for the input msa
    and write it as orthoxml_to_newick.py file


    output: gene tree in nwk format
    """
    wrapper_tree = fasttree.Fasttree(msa, datatype="PROTEIN")
    wrapper_tree.options.options['-fastest']
    result_tree1 = wrapper_tree()

    time_taken_tree = wrapper_tree.elapsed_time
    result_tree2 = wrapper_tree.result
    tree_nwk = str(result_tree2["tree"])
    current_time = datetime.now().strftime("%H:%M:%S")
    # for development we write the gene tree, the name of file should be limit in size in linux.
    # danger of overwriting
    # instead -> hash thing
    # ??? hashlib.md5(original_name).hexdig..it()

    # as the debug==True
    if len(gene_tree_file_addr) > 255: gene_tree_file_addr = gene_tree_file_addr[:250] + ".nwk"
    file_gene_tree = open(gene_tree_file_addr, "w")
    file_gene_tree.write(tree_nwk)
    file_gene_tree.write(";\n")
    file_gene_tree.close()

    return tree_nwk


def lable_SD_internal_nodes(tree_out):
    """
    for the input gene tree, run the species overlap method
    and label internal nodes of the gene tree

    output: labeled gene tree
    """
    species_name_dic = {}
    counter_S = 0
    counter_D = 0

    for node in tree_out.traverse(strategy="postorder"):
        # print("** now working on node ",node.name) # node_children
        if node.is_leaf():
            prot_i = node.name
            species_name_dic[node] = {str(prot_i).split("||")[1][:-1]}
        else:
            node.name = "S/D"
            leaves_list = node.get_leaves()  # print("leaves_list", leaves_list)
            species_name_set = set([str(prot_i).split("||")[1][:-1] for prot_i in leaves_list])
            # print("species_name_set", species_name_set)
            species_name_dic[node] = species_name_set

            node_children = node.children  # print(node_children)
            node_children_species_list = [species_name_dic[node_child] for node_child in node_children]  # list of sets
            # print("node_children_species_list", node_children_species_list)
            node_children_species_intersection = set.intersection(*node_children_species_list)

            if node_children_species_intersection:  # print("node_children_species_list",node_children_species_list)
                counter_D += 1
                node.name = "D" + str(counter_D)
            else:
                counter_S += 1
                node.name = "S" + str(counter_S)
    return tree_out


def get_species_name_gene_tree(node_name: str):
    # Species code is the last part of leaf name sep by _
    sp_code = node_name.split("_")[-1]
    if sp_code == 'RAT':
        return 'RATNO'
    return sp_code


def get_species_name_species_tree_qfo(node_name: str):
    # Species code is the last part of leaf name sep by _
    sp_code = node_name
    return sp_code


def get_reconciled_tree_zmasek(gtree, sptree, inplace=False):
    """
    Reconciles the gene tree with the species tree
    using Zmasek and Eddy's algorithm. Details can be
    found in the paper:
    Christian M. Zmasek, Sean R. Eddy: A simple algorithm
    to infer gene duplication and speciation events on orthoxml_to_newick.py
    gene tree. Bioinformatics 17(9): 821-828 (2001)
    :argument gtree: gene tree (PhyloTree instance)
    :argument sptree: species tree (PhyloTree instance)
    :argument False inplace: if True, the provided gene tree instance is
       modified. Otherwise orthoxml_to_newick.py reconciled copy of the gene tree is returned.
    :returns: reconciled gene tree
    """

    # some cleanup operations
    def cleanup(tree):
        for node in tree.traverse(): node.del_feature("M")

    if not inplace:
        gtree = gtree.copy('deepcopy')

    # check for missing species
    missing_sp = gtree.get_species() - sptree.get_species()
    if missing_sp:
        raise KeyError("* The following species are not contained in the species tree: " + ', '.join(missing_sp))

    # initialization
    sp2node = dict()
    for node in sptree.get_leaves(): sp2node[node.species] = node

    # set/compute the mapping function M(g) for the
    # leaf nodes in the gene tree (see paper for details)
    species = sptree.get_species()
    for node in gtree.get_leaves():
        node.add_feature("M", sp2node[node.species])

    # visit each internal node in the gene tree
    # and detect its event (duplication or speciation)
    for node in gtree.traverse(strategy="postorder"):
        if len(node.children) == 0:
            continue  # nothing to do for leaf nodes

        if len(node.children) != 2:
            cleanup(gtree)
            raise ValueError("Algorithm can only work with binary trees.")

        lca = node.children[0].M.get_common_ancestor(node.children[1].M)  # LCA in the species tree
        node.add_feature("M", lca)

        node.add_feature("evoltype", "S")
        if id(node.children[0].M) == id(node.M) or id(node.children[1].M) == id(node.M):
            node.evoltype = "D"

    cleanup(gtree)
    return gtree


def lable_SD_internal_nodes_reconcilation(gene_tree, species_tree):
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


def infer_HOG_rhog3(rhogid_num_list, gene_id_name):  # , address_rhogs_folder, species_tree_address):
    """
    The prot sequences of orthoxml_to_newick.py rootHOG are located in the fasta file address_rhogs_folder+"HOG_rhogid_num.fa,
    we want to infer all subHOGs of this rootHOG for different taxanomic levels.

    output: orthoxml_to_newick.py python dict (HOG_thisLevel):  key=taxanomic level, value= orthoxml_to_newick.py list of subHOGs.
    """

    HOG_thisLevel_list = []
    len_HOG_thisLevel_list = []
    HOG_thisLevel_xml_all = []
    rhogid_num = 0
    for rhogid_num in rhogid_num_list:
        # rhogid_num = rhogid_num_list[rhogid_num_i]
        logger_hog.info(
            "\n" + "=" * 50 + "\n" + "Working on root hog: " + str(rhogid_num) + ". \n")  # +", ",rhogid_num_i,"-th. \n"
        prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        logger_hog.info("number of proteins in the rHOG is " + str(len(rhog_i)) + ".")

        (species_tree) = read_species_tree(species_tree_address)
        (species_tree, species_names_rhog, prot_names_rhog) = prepare_species_tree(rhog_i, species_tree)
        # species_tree.write()  print(species_tree.write())

        dic_sub_hogs = {}
        # finding hogs at each level of species tree (from leaves to root, bottom up)
        for node_species_tree in species_tree.traverse(strategy="postorder"):
            if node_species_tree.is_leaf():
                # each leaf itself is orthoxml_to_newick.py subhog
                continue
            logger_hog.info("\n" + "*" * 15 + "\n" + "Finding hogs for the taxonomic level:" + str(
                node_species_tree.name) + "\n" + str(node_species_tree.write()) + "\n")
            dic_sub_msas = []
            (dic_sub_hogs) = infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs,
                                                 rhogid_num)
            HOG_thisLevel = dic_sub_hogs[node_species_tree.name]
            logger_hog.info("subHOGs in thisLevel are " + ' '.join(["[" + str(i) + "]" for i in HOG_thisLevel]) + " .")

        for hog_i in HOG_thisLevel:
            print(hog_i)
            if len(hog_i._members) > 1:
                # could be improved
                HOG_thisLevel_xml = hog_i.to_orthoxml(**gene_id_name)
                HOG_thisLevel_xml_all.append(HOG_thisLevel_xml)
                # groups_xml.append(HOG_thisLevel_xml)
                # print(hog_i._members)
        # HOG_thisLevel_list.append(HOG_thisLevel)
        del dic_sub_hogs
        del HOG_thisLevel
        gc.collect()

    with open(f'{address_pickles_folder}file_' + str(rhogid_num) + '.pickle', 'wb') as handle:
        pickle.dump(HOG_thisLevel_xml_all, handle, protocol=pickle.HIGHEST_PROTOCOL)

    num_hog = len(HOG_thisLevel_xml_all)

    del HOG_thisLevel_xml_all
    gc.collect()

    return (num_hog)


def infer_HOG_thisLevel(node_species_tree, rhog_i, species_names_rhog, dic_sub_hogs, rhogid_num):
    # during parralleization, there will be orthoxml_to_newick.py problem, few times it wants to creat the folder
    # if not os.path.exists(gene_trees_folder) :
    #    os.mkdir(gene_trees_folder)
    #  File "code7d_4.py", line 470, in infer_HOG_thisLevel
    # os.mkdir(gene_trees_folder)
    # FileExistsError: [Errno 17] File exists: '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2//gene_trees_test_7d_3/'

    if len(rhog_i) == 0:
        logger_hog.warning('There is no protein in the rHOG: ' + str(rhogid_num))
        dic_sub_hogs[node_species_tree.name] = []
        return (dic_sub_hogs)

    elif len(rhog_i) == 1:
        logger_hog.warning('There is only one protein in the rHOG: ' + str(rhogid_num))
        node_species_name = node_species_tree.children[0].name  # there is only one species (for the one protein)
        prot = rhog_i[0]
        sub_hog_leaf = HOG(prot, node_species_name, rhogid_num)
        subHOGs_children = [sub_hog_leaf]
        HOG_this_level = subHOGs_children
        dic_sub_hogs[node_species_tree.name] = HOG_this_level
        return (dic_sub_hogs)

    sub_msa_list_lowerLevel = []  # including subHOGS of lower level
    subHOGs_children = []

    # print("working on node", node_species_tree.name,"with",len(node_species_tree.children),"children.")
    for node_child in node_species_tree.children:
        if node_child.is_leaf():
            node_species_name = node_child.name
            # extracting those proteins of the rHOG that belongs to this species (child node of species tree)
            interest_list = [idx for idx in range(len(species_names_rhog)) if
                             species_names_rhog[idx] == node_species_name]
            rhog_part = [rhog_i[i] for i in interest_list]
            # sub_msa = [MultipleSeqAlignment([i]) for i in rhog_part]             #print("len",len(rhog_part))

            for prot in rhog_part:
                sub_hog_leaf = HOG(prot, node_species_name, rhogid_num)  # node_species_tree.name
                # list_all_hogs_ever.append(sub_hog_leaf)
                subHOGs_children.append(sub_hog_leaf)
        else:  # the child node is an internal node, subHOGs are inferred till now during traversing.
            # print("sub msa for internal node", node_child.name,"is read from dic.")
            if node_child.name in dic_sub_hogs:
                sub_hogs_child = dic_sub_hogs[node_child.name]
                subHOGs_children += sub_hogs_child
            else:
                logger_hog.error("Error 131, no sub msa for the internal node ", node_child.name, node_child, "\n",
                                 dic_sub_hogs)
                assert 2 == 1
    temp11 = []
    for temp in [i._members for i in subHOGs_children]:
        temp11.append([prot.split('|')[2] for prot in temp])
    # print("there are ",len(subHOGs_children),"subHOGs lower of this level:",[i._hogid for i in subHOGs_children],temp11)
    # print("We want to infer subHOGs at this level,i.e. merge few of them.")
    subHOG_to_be_merged_set_other_Snodes = []

    if len(subHOGs_children) == 0:
        logger_hog.error('Error 139, There is no protein in this subhog, for rhog' + str(rhogid_num))

    elif len(subHOGs_children) == 1:
        HOG_this_level = subHOGs_children
        # print("**** error 134 *** ", len(subHOGs_children),subHOGs_children) #return (-1,-1,-1)

    else:

        sub_msa_list_lowerLevel_ready = [hog._msa for hog in subHOGs_children]
        merged_msa = merge_msa(sub_msa_list_lowerLevel_ready)
        logger_hog.info("All subHOGs are merged, merged msa is with length of" + str(len(merged_msa)) + " " + str(
            len(merged_msa[0])) + ".")

        gene_tree_file_addr = gene_trees_folder + "/tree_" + str(rhogid_num) + "_" + str(
            node_species_tree.name) + ".nwk"
        gene_tree_raw = infer_gene_tree(merged_msa, gene_tree_file_addr)
        gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
        logger_hog.info("Gene tree is infered with length of " + str(len(gene_tree)) + ".")
        # gene_tree_i +=1

        R = gene_tree.get_midpoint_outgroup()
        gene_tree.set_outgroup(R)

        gene_tree = lable_SD_internal_nodes(gene_tree)
        # print("Overlap speciation is done for internal nodes of gene tree, as following:")
        print(str(gene_tree.write(format=1))[:-1] + str(gene_tree.name) + ":0;")

        # assigned_leaves_to_hog = []        #sub_msas_list_this_level = []
        subHOGs_id_children_assigned = []  # the same as  subHOG_to_be_merged_all_id
        HOG_this_level = []
        subHOG_to_be_merged_set_other_Snodes = []
        subHOG_to_be_merged_set_other_Snodes_flattned_temp = []
        for node in gene_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n,
                                                                                         "processed") and n.processed == True):  # start from root
            # print("Leaves assigned to hog are ", assigned_leaves_to_hog)   #print("Traversing gene tree. Now at node", node.name)
            if not node.is_leaf():
                node_leaves_name = [i.name for i in node.get_leaves()]
                # print(node_leaves_name)

                if node.name[0] == "S":  # this is orthoxml_to_newick.py sub-hog.
                    subHOG_to_be_merged = []
                    for node_leave_name in node_leaves_name:  # print(node_leave_name)
                        for subHOG in subHOGs_children:
                            subHOG_members = subHOG._members
                            if node_leave_name in subHOG_members:  # could be improved
                                if subHOG._hogid not in subHOG_to_be_merged_set_other_Snodes_flattned_temp:
                                    subHOG_to_be_merged.append(subHOG)
                                    subHOGs_id_children_assigned.append(subHOG._hogid)
                                else:
                                    print("issue 184", node.name, subHOG._hogid, node_leave_name)
                                    if "processed" in node:
                                        print(node.name)
                                    else:
                                        print("processed not in ",
                                              node.name)  # print(node_leave_name,"is in ",subHOG._hogid)
                    if subHOG_to_be_merged:
                        subHOG_to_be_merged_set = set(subHOG_to_be_merged)
                        taxnomic_range = node_species_tree.name
                        HOG_this_node = HOG(subHOG_to_be_merged_set, taxnomic_range, rhogid_num, msa=merged_msa)
                        HOG_this_level.append(HOG_this_node)
                        subHOG_to_be_merged_set_other_Snodes.append([i._hogid for i in subHOG_to_be_merged_set])
                        subHOG_to_be_merged_set_other_Snodes_flattned_temp = [item for items in
                                                                              subHOG_to_be_merged_set_other_Snodes for
                                                                              item in items]
                        #  I don't need to traverse deeper in this clade
                    node.processed = True  # print("?*?*  ", node.name)

            subHOG_to_be_merged_set_other_Snodes_flattned = [item for items in subHOG_to_be_merged_set_other_Snodes for
                                                             item in items]
            if [i._hogid for i in subHOGs_children] == subHOG_to_be_merged_set_other_Snodes_flattned:
                break
        for subHOG in subHOGs_children:  # for the single branch  ( D include orthoxml_to_newick.py  subhog and orthoxml_to_newick.py S node. )
            if subHOG._hogid not in subHOGs_id_children_assigned:  # print("here", subHOG)
                HOG_this_level.append(subHOG)
        prot_list_sbuhog = [i._members for i in HOG_this_level]
        prot_list_sbuhog_short = []
        for prot_sub_list_sbuhog in prot_list_sbuhog:
            prot_list_sbuhog_short.append([prot.split('|')[2] for prot in prot_sub_list_sbuhog])
        logger_hog.info("- " + str(
            len(prot_list_sbuhog_short)) + "HOGs are inferred at the level " + node_species_tree.name + ": " + " ".join(
            [str(i) for i in prot_list_sbuhog_short]))
    # print("By merging ",subHOG_to_be_merged_set_other_Snodes)

    # check for conflicts in merging
    #     for i in range(subHOG_to_be_merged_set_other_Snodes):
    #         if
    #         for i in range(subHOG_to_be_merged_set_other_Snodes):
    # print("*&*& ",node_species_tree.name)
    dic_sub_hogs[node_species_tree.name] = HOG_this_level
    return (dic_sub_hogs)


class HOG:
    _hogid_iter = 1000

    def __init__(self, input_instantiate, taxnomic_range, rhogid_num, msa=None):  # _prot_names
        # the input_instantiate could be either
        #     1) orthoxml_to_newick.py protein as the biopython seq record  SeqRecord(seq=Seq('MAPSSRSPSPRT. ]
        # or  2) orthoxml_to_newick.py set of intances of class HOG   wit orthoxml_to_newick.py big msa
        # those variable starting with _ are local to the class, should not access directly  (although it is possbile)
        self._rhogid_num = rhogid_num
        self.__class__._hogid_iter += 1
        # 0070124
        self._hogid = "HOG:B" + str(self._rhogid_num).zfill(7) + "_sub" + str(self.__class__._hogid_iter)
        self._taxnomic_range = taxnomic_range  # print("**** orthoxml_to_newick.py new HOG is instantiated with id", self._hogid)

        if isinstance(input_instantiate, SeqRecord):  # if len(sub_hogs)==1:
            only_protein = input_instantiate  # only one seq, only on child, leaf
            self._members = set([only_protein.id])
            self._msa = MultipleSeqAlignment([only_protein])
            self._subhogs = []
            # <<class 'Bio.Align.MultipleSeqAlignment'> instance (1 records of length 314) at 7f0e86c713d0>

        elif msa and all(isinstance(x, HOG) for x in input_instantiate):
            # here we want to merge few subHOGs and creat orthoxml_to_newick.py new HOG.   #the n
            sub_hogs = input_instantiate
            hog_members = set()
            for sub_hog in sub_hogs:
                hog_members |= sub_hog.get_members()  # union
            self._members = hog_members  # set.union(*tup)
            self._subhogs = list(input_instantiate)  # full members

            max_num_seq = 30  # subsampling in msa
            records_full = [record for record in msa if record.id in self._members]
            if len(records_full) > max_num_seq:
                records_sub_sampled = sample(records_full, max_num_seq)  # without replacement.
                logger_hog.info(
                    "we are doing subsamping now from " + str(len(records_full)) + " to " + str(max_num_seq) + "seqs.")
            else:
                records_sub_sampled = records_full
            # removing some columns completely gap -  (not x   )
            # now select those proteins
            self._msa = MultipleSeqAlignment(records_sub_sampled)
            # without replacement sampling ,  # self._children = sub_hogs # as legacy  ?
        else:
            logger_hog.error("Error 169,  check the input format to instantiate orthoxml_to_newick.py HOG class")
            assert False

    def __repr__(self):
        return "an object of class HOG of hogID=" + self._hogid + ", length=" + str(
            len(self._members)) + ", taxonomy= " + str(self._taxnomic_range)

    def get_members(self):
        return set(self._members)
        # merge, gene tree, midpoint, lable_SD_internal_nodes, traverse_geneTree_assign_hog

    def to_orthoxml(self, **gene_id_name):  # , indent=0):
        indent = 0
        hog_elemnt = ET.Element('orthologGroup', attrib={"id": str(self._hogid)})
        property_element = ET.SubElement(hog_elemnt, "property",
                                         attrib={"name": "TaxRange", "value": str(self._taxnomic_range)})
        # the following could be improved ???   without this if it will be like, one property is enough
        # <orthologGroup>
        #    <property name="TaxRange" value="GORGO_HUMAN_PANTR"/>
        #    <property name="TaxRange" value="GORGO_HUMAN_PANTR"/>
        # if property_element not in hog_elemnt:
        #    hog_elemnt.append(property_element)
        #    print("*")
        # gene = ET.SubElement(species, "gene", attrib={"id":str(gene_counter), "protId":query_prot_record.id})
        # hog_elemnt = ET.SubElement(species,

        if len(self._subhogs) == 0:
            # print(self._members)

            # print("we are here   ********???--??? ",self._hogid)
            list_member_first = list(self._members)[0]
            geneRef_elemnt = ET.Element('geneRef', attrib={
                'id': str(gene_id_name[list_member_first])})  # # gene_id_name[query_prot_record.id]
            # hog_elemnt.append(geneRef_elemnt)
            # could be improved when the rhog contains only one protein
            return geneRef_elemnt  # hog_elemnt

        def _sorter_key(sh):
            return sh._taxnomic_range

        self._subhogs.sort(key=_sorter_key)  # print(f'{" "*indent}subhog: {self._taxnomic_range}:')
        for sub_clade, sub_hogs in itertools.groupby(self._subhogs, key=_sorter_key):
            list_of_subhogs_of_same_clade = list(
                sub_hogs)  # print(f'{" "*(indent+1)} clade: {sub_clade} with {str(len(list_of_subhogs_of_same_clade))}')
            if len(list_of_subhogs_of_same_clade) > 1:
                paralog_element = ET.Element('paralogGroup')
                for sh in list_of_subhogs_of_same_clade:
                    paralog_element.append(sh.to_orthoxml(**gene_id_name))  # ,**gene_id_name  indent+2
                hog_elemnt.append(paralog_element)
            else:
                hog_elemnt.append(list_of_subhogs_of_same_clade[0].to_orthoxml(**gene_id_name))  # indent+2
        return hog_elemnt

#
# def prepare_xml(rhogid_num_list_temp):
#     species_prot_dic = {}
#     # all_prot_temp_list= []
#     for rhogid_num in rhogid_num_list_temp:
#         prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
#         rhog_i = list(SeqIO.parse(prot_address, "fasta"))
#         for prot_i in rhog_i:
#             species_i = prot_i.id.split("||")[1][:-1]  # .split("|")[-1].split("_")[-1]
#             if species_i in species_prot_dic:
#                 species_prot_dic[species_i].append(prot_i.id)
#             else:
#                 species_prot_dic[species_i] = [prot_i.id]
#             # all_prot_temp_list.append(prot_i.id)
#
#     print("there are species ", len(species_prot_dic))
#     orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
#                                                    "originVersion": "Nov 2021", "version": "0.3"})  #
#     gene_counter = 100000
#     gene_id_name = {}
#     query_species_names_rHOGs = list(species_prot_dic.keys())
#     for species_name in query_species_names_rHOGs:
#         no_gene_species = True  # for code develop ment
#         species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
#         database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
#         genes_xml = ET.SubElement(database_xml, "genes")
#
#         prot_list = species_prot_dic[species_name]
#         for prot_itr in range(len(prot_list)):  # [12:15]
#             prot_i_name = prot_list[prot_itr]
#             gene_id_name[prot_i_name] = gene_counter
#             prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
#             gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
#             gene_counter += 1
#
#     groups_xml = ET.SubElement(orthoxml_file, "groups")
#
#     return (groups_xml, gene_id_name, orthoxml_file)


def distribute_rhogs(rhogs: List[Tuple[str, int]], start_index: int, n_workers: int) -> List[Tuple[str, int]]:
    """Distribute rhogs to workers in orthoxml_to_newick.py way that each worker gets both big and small rhogs.
    inputs:
        rhogs: orthoxml_to_newick.py list of tuples of (rhogs name, rhogs size in bytes)
        start_index: the index of the worker
        n_workers: the number of workers
    outputs:
        orthoxml_to_newick.py list of tuples of (rhogs name, rhogs size in bytes)"""
    if start_index >= n_workers:
        raise ValueError("index out of range")
    else:
        temp = []
        for i in range(start_index, len(rhogs), n_workers):
            temp.append(rhogs[i])
        return temp


if __name__ == "__main__":

    logging.basicConfig()
    logger_hog = logging.getLogger("hog")
    logger_hog.setLevel(logging.INFO)  # WARN
    # make sure addresses end with "/"

    address_working_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/ali_code_31aug/"

    address_rhogs_folder = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/rhog_all_v3_g2_s500/"

    address_pickles_folder = address_working_folder + "pickle/"

    species_tree_address = "/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo2/archive/lineage_tree_qfo.phyloxml"

    gene_trees_folder = address_working_folder + "/gene_trees/"

    address_logs_folder = address_working_folder + "logs/"
    address_group_xml_ortho = address_working_folder + "group_xml_ortho_adjusted_family_40_2sep5pm.pickle"

    this_worker_index = int(sys.argv[1])
    n_workers = int(sys.argv[2])

    with open(address_group_xml_ortho, 'rb') as handle:
        (groups_xml, gene_id_name, orthoxml_file) = pickle.load(handle)

    ## create orthoxml_to_newick.py list of rootHOG IDs  stored in the folder of rHOG .
    rhog_files = listdir(address_rhogs_folder)[:]
    print("#", rhog_files[:4], len(rhog_files))

    file_sizes = [os.path.getsize(f'{address_rhogs_folder}{i}') for i in rhog_files]
    file_name_size_dict = dict(zip(rhog_files, file_sizes))
    file_name_size_list_sorted = sorted(file_name_size_dict.items(), key=lambda x: x[1])
    this_worker_files = distribute_rhogs(file_name_size_list_sorted, this_worker_index, n_workers)
    print("# This worker total rhog size (bytes):", sum([i[1] for i in this_worker_files]))
    this_worker_files = [i[0] for i in this_worker_files]

    # this_worker_files = rhog_files
    rhogid_num_list = []
    for rhog_file in this_worker_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)
    print("##", rhogid_num_list[:4], len(rhogid_num_list))

    print(infer_HOG_rhog3(rhogid_num_list, gene_id_name))

