

from os import listdir
from Bio import SeqIO
import xml.etree.ElementTree as ET


import logging
logging.basicConfig()
logger_hog = logging.getLogger("hog")
logger_hog.setLevel(logging.INFO)  # WARN




def list_rhog_fastas(address_rhogs_folder):
    """
     create a list of rootHOG IDs  stored in the folder of rHOG .
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




from ete3 import Phyloxml
from ete3 import Tree

def read_species_tree(species_tree_address):
    """
    reading a species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
    format_tree = species_tree_address.split(".")[-1]

    if format_tree == "phyloxml":
        project = Phyloxml()
        project.build_from_file(species_tree_address)
        # Each tree contains the same methods as a PhyloTree object
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


def prepare_species_tree(rhog_i, species_tree, format_prot_name=1):

    """
    a function for extracting a subtree from the input species tree  a.k.a pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    species_names_rhog = []
    prot_names_rhog = []
    for rec in rhog_i:

        if format_prot_name == 1:   # qfo dataset
            prot_name = rec.name  # 'tr|E3JPS4|E3JPS4_PUCGT
            species_name = prot_name.split("|")[-1].split("_")[-1].strip()
            if species_name == 'RAT': species_name = "RATNO"

        elif format_prot_name == 0:  # bird dataset
            # rec.name  CLIRXF_R07389
            prot_name = rec.name
            prot_descrip  = rec.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
            species_name = prot_descrip.split(" ")[1].split("|")[-1]
            # species_name = prot_name.split("_")[0].strip()

        species_names_rhog.append(species_name)
        prot_names_rhog.append(prot_name)

    species_names_uniqe = set(species_names_rhog)
    logger_hog.info("The number of unique species in the rHOG is " + str(len(species_names_uniqe)) + ".")
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
    print(species_tree)

    return (species_tree, species_names_rhog, prot_names_rhog)


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
            species_name_dic[node] = {str(prot_i).split("|")[-1].split("_")[-1]}
        else:
            node.name = "S/D"
            leaves_list = node.get_leaves()  # print("leaves_list", leaves_list)
            species_name_set = set([str(prot_i).split("|")[-1].split("_")[-1] for prot_i in leaves_list])
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


def prepare_xml(rhogid_num_list_input, address_rhogs_folder, format_prot_name, rhogid_batch = 1):
    species_prot_dic = {}
    # all_prot_temp_list= []
    rhogid_len_list = [ ]
    for rhogid_num in rhogid_num_list_input:
        prot_address = address_rhogs_folder + "HOG_B" + str(rhogid_num).zfill(7) + ".fa"
        rhog_i = list(SeqIO.parse(prot_address, "fasta"))
        rhogid_len_list.append(len(rhog_i))

        for prot_i in rhog_i:

            if format_prot_name == 1:  # qfo dataset
                prot_name = prot_i.name  # 'tr|E3JPS4|E3JPS4_PUCGT
                species_i = prot_name.split("|")[-1].split("_")[-1].strip()
                if species_i == 'RAT': species_i = "RATNO"

            elif format_prot_name == 0:  # bird dataset
                # rec.name  CLIRXF_R07389
                # prot_name = prot_i.name
                prot_descrip = prot_i.description  # >CLIRXF_R07389 CLIRXF_R07389|species|CLIRUF
                species_i = prot_descrip.split(" ")[1].split("|")[-1]
                # species_name = prot_name.split("_")[0].strip()

            # species_i = prot_i.id.split("|")[-1].split("_")[-1]
            if species_i in species_prot_dic:
                species_prot_dic[species_i].append(prot_i.id)
            else:
                species_prot_dic[species_i] = [prot_i.id]
            # all_prot_temp_list.append(prot_i.id)

    print("there are species ", len(species_prot_dic))
    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #
    gene_counter = 1000000 + rhogid_batch * 10000
    gene_id_name = {}
    query_species_names_rhogs = list(species_prot_dic.keys())
    for species_name in query_species_names_rhogs:
        no_gene_species = True  # for code develop ment
        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "QFO database ", "version": "2020"})
        genes_xml = ET.SubElement(database_xml, "genes")

        prot_list = species_prot_dic[species_name]
        for prot_itr in range(len(prot_list)):  # [12:15]
            prot_i_name = prot_list[prot_itr]
            gene_id_name[prot_i_name] = gene_counter
            if "|" in prot_i_name:
                prot_i_name_short = prot_i_name.split("|")[1].strip()  # tr|E3JPS4|E3JPS4_PUCGT
            else:
                prot_i_name_short = prot_i_name

            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_counter), "protId": prot_i_name_short})
            gene_counter += 1

    groups_xml = ET.SubElement(orthoxml_file, "groups")


    return (groups_xml, gene_id_name, orthoxml_file, rhogid_len_list)
