
import sys
from ete3 import Tree
import os
import collections
from os import listdir

from . import _config

def check_species_tree_internalnode(species_tree):

    # All the internal node of the  input species tree should have a name
    for node in species_tree.traverse(strategy="postorder"):
        if not node.is_leaf():
            internal_node_name = []
            if (not node.name) or len(node.name) < 1:
                print("one of the internal node doesn't have a name in species tree", species_tree)
                sys.exit()
            else:
                internal_node_name.append( node.name)

    if len(internal_node_name) != len(set(internal_node_name)):
        print("All the internal node names should be unique. One of the internal node is repeated:")
        print([item for item, count in collections.Counter(internal_node_name).items() if count > 1])
        sys.exit()

    return 1

def check_species_tree_proteome(species_tree):  # list_oma_species
    project_files = listdir("./proteome/")
    print("there are " + str(
        len(project_files)) + "files in the proteome folder. Better to have only fasta/fa format.")

    query_species_names = []
    for file in project_files:
        fasta_format = file.split(".")[-1]
        if fasta_format == "fa" or fasta_format == "fasta":
            file_name_split = file.split(".")[:-1]
            query_species_names.append('.'.join(file_name_split))
    # logger_hog.info("The are " + str(query_species_num) + " species in the proteome folder.")
    query_species_names_set = set(query_species_names)
    if not len(query_species_names) != len(query_species_names_set):
        print("the species name should be unique in the proteome folder.")
        sys.exit()

    try:
        species_tree.prune(set(query_species_names_set))

    except:
        print(
            "there is an mismatch between fasta file name and species tree. We expect to see such tree :  ((AQUAE,CHLTR)inter1,MYCGE)inter2; for proteome folder with these three files:  AQUAE.fa  CHLTR.fa  MYCGE.fa. ")
        sys.exit()

    return 1

def check_omamer_db():
    omamerdb_adress = "omamerdb.h5"
    if _config.inferhog_resume_subhog:
        if not os.path.exists(omamerdb_adress):
            print("file doesnt not exist",omamerdb_adress)
            if os.path.getsize(omamerdb_adress) > 1000:  # 3 bytes
                print("WArning: the omamer db looks very small. are you sure it is correct? ", omamerdb_adress)

    return 1


def check_input_fastoma():
    _config.set_configs()

    # print(_config.in_folder)
    print(_config.logger_level)

    print(_config.species_tree_address)  # nwk format

    species_tree = Tree(_config.species_tree_address, format=1, format_root_node=True)
    check_species_tree_internalnode(species_tree)

    check_species_tree_proteome(species_tree)  # with proteme folder
    # all the species are in the tree

    check_omamer_db()


    # todo check all the input and needed folders as a validation with nextflow

"""

print("there shouldnt be any space in the tree name internal node name as well")
  '/work/FAC/FBM/DBC/cdessim2/default/smajidi1/fastget/qfo3/working_nf//pickles_subhog/rhog_833762/delta/epsilon subdivisions.pickle'

check if all files in proteome are represented in the newick tree.
throw a warning wit few  proteins 


"""