

from os import listdir
# from Bio import SeqIO
# import xml.etree.ElementTree as ET
from ete3 import Phyloxml
from ete3 import Tree

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq

import logging
#logging.basicConfig()
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logger_hog = logging.getLogger("hog")
logger_hog.setLevel(logging.INFO)  # WARN  INFO


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


def read_species_tree(species_tree_address):
    """
    reading orthoxml_to_newick.py species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
    format_tree = species_tree_address.split(".")[-1]

    if format_tree == "phyloxml":
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
    elif format_tree == "nwk":
        species_tree = Tree(species_tree_address)
    else:
        print("for now we accept phyloxml or nwk format for input species tree.")

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
        # qfo : >tr|A0A0N7KF21|A0A0N7KF21_ORYSJ||ORYSJ_||1000000344 tr|A0A0N7KF21|A0A0N7KF21_ORYSJ Os02g0264501 protein OS=Oryza sativa subsp. japonica (Rice) OX=39947 GN=Os02g0264501 PE=4 SV=1
        prot_id = rec.id.split("||")
        prot_name = prot_id[2]   # for debugging  prot_id[0] readable prot name,  for xml prot_id[2]
        species_name = prot_id[1][:-1]
        if species_name == 'RAT': species_name = "RATNO"
        # gene_id = prot_id[2]
        species_names_rhog.append(species_name)
        prot_names_rhog.append(prot_name)

    species_names_uniqe = set(species_names_rhog)

    species_tree.prune(species_names_uniqe, preserve_branch_length=True)
    # species_tree.write()
    for node in species_tree.traverse(strategy="postorder"):
        node_name = node.name
        num_leaves_no_name = 0
        if len(node_name) < 1:
            if node.is_leaf():
                node.name = "leaf_" + str(num_leaves_no_name)
            else:
                node_children = node.children
                list_children_names = [node_child.name for node_child in node_children]
                node.name = '_'.join(list_children_names)
    print("Working on the following species tree.")
    # print(species_tree)
    species_tree.write()

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


def msa_filter_row(msa, tresh_ratio_gap_row):  # gene_tree_file_addr

    msa_filtered_row = []
    ratio_records=[]
    for record in msa:
        seq = record.seq
        seqLen = len(record)
        gap_count = seq.count("-") + seq.count("?") + seq.count(".") +seq.count("~")
        ratio_record_nongap= 1-gap_count/seqLen
        ratio_records.append(round(ratio_record_nongap, 3))
        if ratio_record_nongap > tresh_ratio_gap_row:
            msa_filtered_row.append(record)
    # out_name_msa=gene_tree_file_addr +"filtered_row_"+str(tresh_ratio_gap_row)+".txt"
    # handle_msa_fasta = open(out_name_msa,"w")
    # SeqIO.write(msa_filtered_row, handle_msa_fasta,"fasta")
    # handle_msa_fasta.close()
    return msa_filtered_row

def msa_filter_col(msa, tresh_ratio_gap_col):

    ratio_col_all = []
    length_record= len(msa[1])
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

    # out_name_msa=gene_tree_file_addr+ "filtered_row_"+"_col_"+str(tresh_ratio_gap_col)+".txt"
    # handle_msa_fasta = open(out_name_msa,"w")
    # SeqIO.write(msa_filtered_col, handle_msa_fasta,"fasta")
    # handle_msa_fasta.close()
    # print("- Column-wise filtering of MSA is finished",len(msa_filtered_col),len(msa_filtered_col[0]))
    return msa_filtered_col
