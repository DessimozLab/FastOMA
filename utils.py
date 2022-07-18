

from os import listdir


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



import logging
from ete3 import Phyloxml


def read_species_tree(species_tree_address):
    """
    reading a species tree in Phyloxml format using ete3 package .

    output (species_tree)
    """
    logger_hog.info(species_tree_address)
    # print(round(os.path.getsize(species_tree_address)/1000),"kb")
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


    return (species_tree)


def prepare_species_tree(rhog_i, species_tree):
    """
    a function for extracting a subtree from the input species tree  a.k.a pruning,
    based on the names of species in the rootHOG.

    output: species_tree (pruned), species_names_rhog, prot_names_rhog
    """
    species_names_rhog = []
    prot_names_rhog = []
    for rec in rhog_i:
        prot_name = rec.name  # 'tr|E3JPS4|E3JPS4_PUCGT
        # prot_name = prot_name_full.split("|")[1].strip() # # 'tr|E3JPS4|E3JPS4_PUCGT
        species_name = prot_name.split("|")[-1].split("_")[-1]
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
                node.name = "leaf_" + str(num_leaves_no_name)
            else:
                node_children = node.children
                list_children_names = [node_child.name for node_child in node_children]
                node.name = '_'.join(list_children_names)
    print("Working on the following species tree.")
    print(species_tree)

    return (species_tree, species_names_rhog, prot_names_rhog)
