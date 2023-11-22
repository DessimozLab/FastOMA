
from Bio import SeqIO
from ete3 import Phyloxml
from ete3 import Tree
from ete3 import PhyloTree
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from collections import defaultdict
from typing import List, Tuple
import random
from itertools import combinations
import numpy as np
import sys
from os import listdir


from ._config import logger_hog
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
    rhogid_list = []
    for rhog_file in rhog_files:
        if rhog_file.split(".")[-1] == "fa":
            rhogid = rhog_file.split(".")[0].split("_")[1] #  [1:]
            rhogid_list.append(rhogid)

    return rhogid_list


def read_species_tree(species_tree_address):
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
                logger_hog.error("Format of species tree is not known or the file doesn't exist "+species_tree_address )
                sys.exit()
    else:
        logger_hog.error("For now we accept phyloxml or nwk format for input species tree.or the file doesn't exist "+species_tree_address)
        sys.exit()

    return species_tree

## Nevers rooting
def get_ancestors_set(tree):
    full_lineage = {}
    if tree.is_leaf():
        full_lineage[tree.name] = [tree.name]
    else:
        for c in tree.children:
            full_lineage = full_lineage | get_ancestors_set(c)
        for k, v in full_lineage.items():
            full_lineage[k] = v + [tree.name]

    return full_lineage


def get_tax_score_from_species_tree(species_tree):
    pair_score = dict()
    species_lineage = get_ancestors_set(species_tree)
    dict_set = {k: set(v) for k, v in species_lineage.items()}

    sorted_keys = sorted(dict_set.keys())
    for i in range(len(sorted_keys)):
        key_a = sorted_keys[i]
        for j in range(i, len(sorted_keys)):
            key_b = sorted_keys[j]
            pair_score[(key_a, key_b)] = len(dict_set[key_a].intersection(dict_set[key_b]))
    return pair_score


def convert_gene_tree(gene_tree):
    for l in gene_tree.get_leaves():
        l.species = l.name.split('||')[1]


def get_full_score_tree(gene_tree, full_lineage):
    lineage = None
    score = 0
    if gene_tree.is_leaf():
        lineage = full_lineage[gene_tree.species]

    else:
        for c in gene_tree.children:
            cscore, clineage = get_full_score_tree(c, full_lineage)

            score += cscore

            if lineage:
                lineage = lineage.intersection(clineage)
            else:
                lineage = clineage

        score += len(lineage)
    return score, lineage


def get_score_all_root(gtree, stree):
    max_score = 0
    best_tree = [None]

    lineage = get_ancestors_set(stree)
    lineage = {k: set(v) for k, v in lineage.items()}
    convert_gene_tree(gtree)

    edge = 0
    ancestor = []
    default_root = gtree.get_leaves()[0].name
    gtree.set_outgroup(gtree & default_root)
    for node in gtree.traverse():
        if not node.is_root() and node.name != default_root:
            if not node.is_leaf():
                node.name = str(edge)

            ancestor.append(node.name)
            edge += 1
    for i in ancestor:
        gtree.set_outgroup(gtree & default_root)
        gtree.set_outgroup(gtree & i)
        score, _ = get_full_score_tree(gtree, lineage)

        if score >= max_score:
            if score > max_score:
                max_score = score
                best_tree = []
            children = gtree.children

            fleaf1 = max([gtree.get_distance(l) for l in children[0].get_leaves()])
            fleaf2 = max([gtree.get_distance(l) for l in children[1].get_leaves()])

            if max(fleaf1, fleaf2) > 0:
                ratio = min(fleaf1, fleaf2) / max(fleaf1, fleaf2)
                best_tree.append((i, ratio))
            else:
                ratio=1
                best_tree.append((i, ratio))

            #diameter = abs(1 - (max([gtree.get_distance(l) for l in children[0].get_leaves()]) / max(
            #    max([gtree.get_distance(l) for l in children[1].get_leaves()]), 1e-20)))
            #best_tree.append((i, diameter))

    best_tree.sort(key=lambda x: x[1])
    gtree.set_outgroup(gtree & default_root)
    gtree.set_outgroup(gtree & best_tree[0][0]) # is it the min?

    return gtree

def genetree_sd(node_species_tree, gene_tree, genetree_msa_file_addr, hogs_children_level_list=[]):

    if _config.rooting_method == "midpoint":
        r_outgroup = gene_tree.get_midpoint_outgroup()
        try:
            gene_tree.set_outgroup(r_outgroup)  # print("Midpoint rooting is done for gene tree.")
        except:
            pass
    elif  _config.rooting_method == "Nevers_rooting":
        logger_hog.info("Nevers_rooting started for " +str(gene_tree.write(format=1, format_root_node=True)))
        species= Tree("/work/FAC/FBM/DBC/cdessim2/default/smajidi1/qfo/qfo_run_raw/in_folder/species_tree.nwk",format=1)
        gene_tree  = get_score_all_root(gene_tree, species)
        logger_hog.info("Nevers_rooting finished for " + str(gene_tree.write(format=1, format_root_node=True)))

    elif _config.rooting_method == "mad":
        gene_tree = _wrappers.mad_rooting(genetree_msa_file_addr) # todo check with qouted gene tree
    # elif _config.rooting_method == "outlier":
    #     gene_tree = PhyloTree(gene_tree_raw + ";", format=0)
    #     outliers = find_outlier_leaves(gene_tree)
    #     r_outgroup = midpoint_rooting_outgroup(gene_tree, leaves_to_exclude=outliers)
    #     gene_tree.set_outgroup(r_outgroup)
    else:
        logger_hog.warning("rooting method not found !!   * * * * *  *")

    all_species_dubious_sd_dic = {}
    if _config.label_SD_internal == "species_overlap":
        (gene_tree, all_species_dubious_sd_dic) = label_sd_internal_nodes(gene_tree)

    elif _config.label_SD_internal == "reconcilation":
        node_species_tree_nwk_string = node_species_tree.write(format=1)
        node_species_tree_PhyloTree = PhyloTree(node_species_tree_nwk_string, format=1)
        gene_tree_nwk_string = gene_tree.write(format=1, format_root_node=True)
        gene_tree_PhyloTree = PhyloTree(gene_tree_nwk_string, format=1)
        gene_tree = label_SD_internal_nodes_reconcilation(gene_tree_PhyloTree, node_species_tree_PhyloTree)

    # for better viz, the subhog ID is added to leaves of gene tree
    if hogs_children_level_list:
        for node in gene_tree.traverse(strategy="postorder"):
            if node.is_leaf():
                node_name_old_raw = node.name
                if node_name_old_raw.startswith("'"): # quated gene tree
                    node_name_old = node_name_old_raw[1:-1]  # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005'", gene tree is quoted, there are both ' and " !
                else:
                    node_name_old = node_name_old_raw
                # node_name_old = node.name
                for hog_child in hogs_children_level_list:
                    if node_name_old in hog_child._members:
                        #node_name_new = node_name_old.split("||")[0]+" "+ hog_child._hogid.split("_")[-1]
                        node_name_new = "'"+node_name_old + "|_|" + hog_child._hogid.split("_")[-1]+"'"
                        # BUCABY_R15453||BUCABY||1286015722_sub10216
                        node.name = node_name_new
                        break

    if _config.gene_trees_write:
        gene_tree.write(format=1, format_root_node=True, outfile=genetree_msa_file_addr+"_SD_labeled.nwk")

    return gene_tree, all_species_dubious_sd_dic



def prepare_species_tree(rhog_i, species_tree, rhogid):
    """
    orthoxml_to_newick.py function for extracting orthoxml_to_newick.py subtree from the input species tree  orthoxml_to_newick.py.k.orthoxml_to_newick.py pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    assert len(rhog_i) > 0, 'input hog_i is empty, probably previous step find_rhog has issue, rhogs/HOG_B0'+rhogid+'is empty?'
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
    # todo check internal node with one child, we need to report for i or not ?
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
    #             # then checn condition in  _infer_subhog.py, where logger_hog.info("Finding hogs for rhogid: "+str(rh
    #             node.name = "internal_" + str(counter_internal)  #  +"_rhg"+rhogid  #  for debuging
    #             counter_internal += 1
    # print("Working on the following species tree.")
    # print(species_tree.write(format=1, format_root_node=True))

    return species_tree, species_names_rhog, prot_names_rhog



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
            prot_i_raw = node.name
            if prot_i_raw.startswith("'"):
                prot_i = prot_i_raw[1:-1]  # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005'", gene tree is quoted, there are both ' and " !
            else:
                prot_i = prot_i_raw

            #prot_i = node.name # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005'", gene tree is quoted, there are both ' and " !
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
                if len(node_children_species_intersection)/ len(node_children_species_union) <= _config.threshold_dubious_sd:
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
    # todo all over the functions
    #  todo node_leaves_name = [i[1:-1] for i in node_leaves_name_raw1] # "'sp|O67547|SUCD_AQUAE||AQUAE||1002000005'", gene tree is quoted, there are both ' and " !

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
        g_node_species_all.append(g_node.name.split("||")[1]) # todo check quoted
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
        if ratio_col_nongap >= tresh_ratio_gap_col:
            keep_cols.append(col_i)
    #plt.hist(ratio_col_all,bins=100) # , bins=10
    #plt.show()
    #plt.savefig(gene_tree_file_addr+ "filtered_row_"+"_col_"+str(tresh_ratio_gap_col)+".txt.pdf")
    #print("- Columns indecis extracted. Out of ", length_record,"columns,",len(keep_cols),"is remained.")
    msa_filtered_col = []
    for record in msa:
        record_seq = str(record.seq)
        record_seq_edited = ''.join([record_seq[i] for i in keep_cols])
        if len(record_seq_edited) > 1:
            record_edited = SeqRecord(Seq(record_seq_edited), record.id, '', '')
            msa_filtered_col.append(record_edited)

    if _config.msa_write_all and gene_tree_file_addr:
        out_name_msa=gene_tree_file_addr+"_filtered_"+"col_"+str(tresh_ratio_gap_col)+".msa.fa"
        handle_msa_fasta = open(out_name_msa, "w")
        SeqIO.write(msa_filtered_col, handle_msa_fasta, "fasta")
        handle_msa_fasta.close()
    # print("- Column-wise filtering of MSA is finished",len(msa_filtered_col),len(msa_filtered_col[0]))
    return MultipleSeqAlignment(msa_filtered_col)


def msa_filter_row(msa, inferhog_tresh_ratio_gap_row, gene_tree_file_addr=""):
    msa_filtered_row = []
    ratio_records = []
    for record in msa:
        seq = record.seq
        seqLen = len(record)
        gap_count = seq.count("-") + seq.count("?") + seq.count(".") +seq.count("~")
        if seqLen:
            ratio_record_nongap = 1-gap_count/seqLen
            ratio_records.append(round(ratio_record_nongap, 3))
            if ratio_record_nongap >= inferhog_tresh_ratio_gap_row:
                msa_filtered_row.append(record)
        else:
            logger_hog.warning("issue 12788 : error , seq len is zero when msa_filter_row")
    if _config.msa_write_all and gene_tree_file_addr:
        out_name_msa = gene_tree_file_addr +"_filtered_row_"+str(inferhog_tresh_ratio_gap_row)+".msa.fa"
        handle_msa_fasta = open(out_name_msa, "w")
        SeqIO.write(msa_filtered_row, handle_msa_fasta, "fasta")
        handle_msa_fasta.close()

    return MultipleSeqAlignment(msa_filtered_row)



def filter_msa(merged_msa, gene_tree_file_addr, hogs_children_level_list):

    msa_filt_row_1 = merged_msa
    # if _config.inferhog_filter_all_msas_row:
    #  msa_filt_row_1 = _utils_subhog.msa_filter_row(merged_msa, _config.inferhog_tresh_ratio_gap_row, gene_tree_file_addr)
    # if   msa_filt_row_1 and len(msa_filt_row_1[0]) >=
    if len(msa_filt_row_1[0]) >= _config.inferhog_min_cols_msa_to_filter:
        # (len(merged_msa) > 10000 and len(merged_msa[0]) > 3000) or (len(merged_msa) > 500 and len(merged_msa[0]) > 5000) or (len(merged_msa) > 200 and len(merged_msa[0]) > 9000):
        # for very big MSA, gene tree is slow. if it is full of gaps, let's trim the msa.
        # logger_hog.debug( "We are doing MSA trimming " + rhogid + ", for taxonomic level:" + str(node_species_tree.name))
        # print(len(merged_msa), len(merged_msa[0]))

        if _config.automated_trimAL:
            msa_filt_col = msa_filt_row_1
            msa_filt_row_col = _wrappers.trim_msa(msa_filt_row_1)
        else:
            msa_filt_col = msa_filter_col(msa_filt_row_1, _config.inferhog_tresh_ratio_gap_col, gene_tree_file_addr+"_0_")
            msa_filt_row_col = msa_filt_col
            if msa_filt_col and msa_filt_col[0] and len(msa_filt_col[0]):
                msa_filt_row_col_raw = msa_filter_row(msa_filt_col, _config.inferhog_tresh_ratio_gap_row, gene_tree_file_addr+"_1_")
                msa_filt_row_col = msa_filter_col(msa_filt_row_col_raw, _config.inferhog_tresh_ratio_gap_col, gene_tree_file_addr+"_2_")


        # compare msa_filt_row_col and msa_filt_col,
        if len(msa_filt_row_col) != len(msa_filt_col):  # some sequences are removed
            set_prot_before = set([i.id for i in msa_filt_col])
            set_prot_after = set([i.id for i in msa_filt_row_col])
            prots_to_remove_level = set_prot_before - set_prot_after
            for prots_to_remove in prots_to_remove_level:
                logger_hog.debug("** we are removing the sequence "+str(prots_to_remove)+" due to trimming")
                # we may want to tag it in the hog object
            #assert len(prots_to_remove_level), "issue 31235"
            #prots_to_remove |= prots_to_remove_level
            # todo: we may want to remove prot from all subhogs
            # should I remove them from subhog._members ? we are doing so for merging fragments or low species overlap I guess
            # hogs_children_level_list_raw = hogs_children_level_list
            # hogs_children_level_list = []
            # for hog_i in hogs_children_level_list_raw:
            #   for prot_to_remove in prots_to_remove:
            #         result_removal = hog_i.remove_prot_from_hog(prot_to_remove)
            #       if result_removal != 0:
            #             hogs_children_level_list.append(hog_i)
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