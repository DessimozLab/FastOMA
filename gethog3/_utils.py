

from os import listdir
from Bio import SeqIO
from ete3 import Phyloxml
from ete3 import Tree
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq  # , UnknownSeq
import logging

import pickle
from os import listdir
from xml.dom import minidom
import xml.etree.ElementTree as ET
import sys
import _config

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger_hog = logging.getLogger("hog")

if _config.logger_level == "INFO":
    logger_hog.setLevel(logging.INFO)  # DEBUG WARN  INFO
if _config.logger_level == "DEBUG" :
    logger_hog.setLevel(logging.DEBUG)  # DEBUG WARN  INFO



#
# TRACE
# DEBUG
# INFO
# WARN
# ERROR
# FATAL


def list_rhog_fastas(address_rhogs_folder):
    """
     create orthoxml_to_newick.py list of rootHOG IDs  stored in the folder of rHOG .
     input: folder address
     output: list of rhog Id (integer)
    """
    rhog_files = listdir(address_rhogs_folder)
    rhogid_num_list= []
    for rhog_file in rhog_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid_num = int(rhog_file.split(".")[0].split("_")[1][1:])
            rhogid_num_list.append(rhogid_num)

    return rhogid_num_list


def read_species_tree_add_internal(species_tree_address):
    """
    reading orthoxml_to_newick.py species tree in Phyloxml format using ete3 package .

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
                logger_hog.error("format of species tree is not known or the file doesn't exist" )
                sys.exit()


    else:
        logger_hog.error("for now we accept phyloxml or nwk format for input species tree.or the file doesn't exist")
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

    return (species_tree)


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
        # qfo : >tr|A0A0N7KF21|A0A0N7KF21_ORYSJ||ORYSJ_||1000000344 tr|A0A0N7KF21|A0A0N7KF21_ORYSJ Os02g0264501 protein OS=Oryza sativa subsp. japonica (Rice) OX=39947 GN=Os02g0264501 PE=4 SV=1
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

    #   print(species_tree.write(format=1))
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
    #
    #             # ?? to imrpove, if the species tree has internal node name, keep it,
    #             # then checn condition in  _inferhog.py, where logger_hog.info("Finding hogs for rhogid_num: "+str(rh
    #
    #             node.name = "internal_" + str(counter_internal)  #  +"_rhg"+str(rhogid_num)  #  for debuging
    #             counter_internal += 1
    # print("Working on the following species tree.")
    # print(species_tree.write(format=1))
    #species_tree.write()

    return (species_tree, species_names_rhog, prot_names_rhog)


def lable_sd_internal_nodes(tree_out):
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
            # species_name_dic[node] = {str(prot_i).split("|")[-1].split("_")[-1]}
            #print(prot_i)
            species_name_dic[node] = {str(prot_i).split("||")[1][:-1]}
        else:
            node.name = "S/D"
            leaves_list = node.get_leaves()  # print("leaves_list", leaves_list)
            # species_name_set = set([str(prot_i).split("|")[-1].split("_")[-1] for prot_i in leaves_list])
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
        node.add_feature("M",lca)

        node.add_feature("evoltype", "S")
        #node.name = "S"
        if id(node.children[0].M) == id(node.M) or id(node.children[1].M) == id(node.M):
                node.evoltype = "D"
                #node.name = "D"

    cleanup(gtree)
    return gtree


def msa_filter_col(msa, tresh_ratio_gap_col, gene_tree_file_addr=""):
    # gene_tree_file_addr contains roothog numebr
    # note this is used in hog class as well

    ratio_col_all = []
    length_record= len(msa[0])
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
    for record in msa :
        record_seq = str(record.seq)
        record_seq_edited  = ''.join([record_seq[i] for i in keep_cols  ])
        record_edited= SeqRecord(Seq(record_seq_edited), record.id, '', '')
        msa_filtered_col.append(record_edited)

    if _config.gene_trees_write and gene_tree_file_addr:
        out_name_msa=gene_tree_file_addr+"filtered_"+"_col_"+str(tresh_ratio_gap_col)+".msa.fa"
        handle_msa_fasta = open(out_name_msa,"w")
        SeqIO.write(msa_filtered_col, handle_msa_fasta,"fasta")
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
    if _config.gene_trees_write and gene_tree_file_addr:
        out_name_msa = gene_tree_file_addr +"filtered_row_"+str(tresh_ratio_gap_row)+".msa.fa"
        handle_msa_fasta = open(out_name_msa, "w")
        SeqIO.write(msa_filtered_row, handle_msa_fasta, "fasta")
        handle_msa_fasta.close()
    return msa_filtered_row


def collect_write_xml():

    gene_id_pickle_file = _config.working_folder + "gene_id_dic_xml.pickle"
    pickles_rhog_folder = _config.working_folder + "pickles_rhog/"
    output_xml_file = _config.working_folder + "hogs.orthoxml"

    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        #gene_id_name = dill_pickle.load(handle)
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)

    #if not _config.write_all_prots_in_header:

    for query_species_name, list_prots in gene_id_name.items():

        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")

        if _config.protein_format_qfo_dataset:
            for (gene_idx_integer, query_prot_name) in list_prots:
                # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ   for qfo benchamrk, the middle should be wirtten in the file
                query_prot_name_pure = query_prot_name.split("|")[1]
                gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})
        else:
            for (gene_idx_integer, query_prot_name) in list_prots:
                query_prot_name_pure = query_prot_name
                gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": query_prot_name_pure})


    pickle_files_adress = listdir(pickles_rhog_folder)

    hogs_a_rhog_xml_all = []
    for pickle_file_adress in pickle_files_adress:
        with open(pickles_rhog_folder + pickle_file_adress, 'rb') as handle:
            hogs_a_rhog_xml_batch = pickle.load(handle)  # hogs_a_rhog_xml_batch is orthoxml_to_newick.py list of hog object.
            hogs_a_rhog_xml_all.extend(hogs_a_rhog_xml_batch)
    print("number of hogs in all batches is ", len(hogs_a_rhog_xml_all))
    groups_xml = ET.SubElement(orthoxml_file, "groups")

    for hogs_a_rhog_xml in hogs_a_rhog_xml_all:
        groups_xml.append(hogs_a_rhog_xml)

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str[:-1000])

    with open(output_xml_file, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    print("orthoxml is written in "+ output_xml_file)
    return 1





from collections import defaultdict
from typing import List, Tuple
import random
from itertools import combinations

import numpy as np
from ete3 import PhyloTree


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