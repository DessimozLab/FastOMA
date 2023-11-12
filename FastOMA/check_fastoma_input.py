

import sys
from ete3 import Tree
import os
import collections

from . import _config
from . import _utils_roothog
from ._utils_subhog import logger_hog



"""
cd in_folder
python check_fastoma_input.py
"""


def check_proteome_files():

    proteome_files = os.listdir("./proteome/")
    logger_hog.info("There are " + str(len(proteome_files)) + " files in the proteome folder. ")

    not_fa = [i for i in proteome_files if not (i.endswith(".fa") or i.endswith(".fasta") )]
    if not_fa:
        logger_hog.warning("We expect that only fa/fasta files are in the proteome folder. Better to remove these " + str(not_fa) )

    fa_fasta = [i.split(".")[-1] for i in proteome_files if  i.endswith(".fa") or i.endswith(".fasta")]
    if len(fa_fasta)<=2:
        logger_hog.error("There are not enough proteomes in the folder ")
        logger_hog.error("Check input failed. FastOMA halted!")
        sys.exit(1)

    if len(set(fa_fasta))>1:
        logger_hog.warning("We expect that all fasta files are with the same format either fa or fasta, we won't include some of them " + str(proteome_files) )


    return 1


def check_proteome(species_names,prot_recs_lists):

    for species_name in species_names:
        num_prots = len(prot_recs_lists[species_name])
        if num_prots <= 2:
            logger_hog.error("The input proteome looks empty or too few proteins or Bio.SeqIO couldn't read it,  in_folder/proteome/" + species_name + "."+fasta_format_keep)

    return 1



def check_hogmap_files():

    hogmap_files = os.listdir("./hogmap_in/")
    logger_hog.info("There are " + str(len(hogmap_files)) + " files in the hogmap folder. ")
    not_hogmap = [i for i in hogmap_files if not i.endswith(".hogmap") ]
    if not_hogmap:
        logger_hog.warning("We expect that only hogmap files are in the hogmap_in folder. Better to remove these " + str(not_hogmap) )

    species_hogmaps = ["".join(i.split(".")[:-2]) for i in hogmap_files if  i.endswith(".hogmap")]
    logger_hog.info("There are "+str(len(species_hogmaps))+" hogmaps.")

    return species_hogmaps



def check_speciestree_internalnode(species_tree):

    # All the internal node of the input species tree should have a name
    for node in species_tree.traverse(strategy="postorder"):
        if not node.is_leaf():
            internal_node_name = []
            if (not node.name) or len(node.name) < 1:
                logger_hog.warning("one of the internal node in species tree doesn't have a name. we'll update the species tree.")
            else:
                internal_node_name.append(node.name)

    if len(internal_node_name) != len(set(internal_node_name)):
        logger_hog.warning("All the internal node names should be unique. One of the internal node is repeated:")
        logger_hog.warning([item for item, count in collections.Counter(internal_node_name).items() if count > 1])
        logger_hog.warning("We'll change the internal node names.")



    return 1

def check_speciestree_leaves(species_tree,species_names):

    leaves_name = [i.name for i in species_tree.get_leaves()]
    logger_hog.info("The species tree has " + str(len(leaves_name))+" leaves.")
    if len(set(leaves_name))!=len(leaves_name):
        logger_hog.error("Leaves name should be unique in the species tree. You may try ete3 prune")
        logger_hog.error("Check input failed. FastOMA halted!")
        sys.exit(1)

    species_names_not_intree = [i for i in species_names if i not in  leaves_name]
    if species_names_not_intree:
        logger_hog.error("These species are not in the species tree "+ str(species_names_not_intree))
        logger_hog.error("This is a tree template:  ((AQUAE,CHLTR)inter1,MYCGE)inter2; for proteome folder with these three files:  AQUAE.fa  CHLTR.fa  MYCGE.fa. ")
        logger_hog.error("Check input failed. FastOMA halted!")
        sys.exit(1)
    leaves_tree_not_proteome = [i for i in leaves_name  if i not in species_names]
    if leaves_tree_not_proteome:
        logger_hog.warning("there are "+str(len(leaves_tree_not_proteome))+" leaves in the species tree that there is no proteome. So we will discard them.")

    return 1

def check_omamer_db(omamerdb_adress="omamerdb.h5"):

    if  os.path.exists(omamerdb_adress):
        if os.path.getsize(omamerdb_adress) > 10000:  # 3 bytes
            omamerdb = True
            # todo we can do some checks on version omamer v2
        else:
            logger_hog.warning("The omamer db looks very small. are you sure it is correct?"+omamerdb_adress)
            omamerdb = False
    else:
        omamerdb =False
        logger_hog.info("OMAmer db does not exist.")

    return omamerdb


def add_internal_node_prune(species_tree,species_names):

    # add name for the internal, if no name is provided, removintg special chars
    counter_internal = 0
    node_names = set()
    for node in species_tree.traverse(strategy="postorder"):
        if not node.is_leaf() :
            node_name = node.name

            if len(node_name) < 3 or node_name in node_names:
                if not node.is_leaf():
                    node.name = "internal_" + str(counter_internal)
                    counter_internal += 1
                    node_names.add(node.name)
                    logger_hog.debug("The internal node name was too small or repeated "+node_name+" which is changed to "+node.name)
            elif any(not c.isalnum() for c in node_name):
                node_name_new = ''.join(e for e in node_name if e.isalnum() or e=="_" or e=="-") # removing special chars
                if node_name_new in node_names:
                    node.name =node_name_new+"_"+str(counter_internal)
                    counter_internal += 1
                else:
                    node.name=node_name_new

                node_names.add(node.name)
                logger_hog.debug("The internal node name has special chars " + node_name + " which is changed to " + node.name)
            else:
                node_names.add(node_name)

    logger_hog.info("The species tree has " + str(len(species_tree)) + " leaves")
    species_tree.prune(species_names) # , preserve_branch_length=True)
    logger_hog.info("After pruning, the species tree has " + str(len(species_tree)) + " leaves")

    species_tree_output= "./species_tree_checked.nwk"
    species_tree.write(format=1, format_root_node=True, outfile=species_tree_output)
    logger_hog.debug("The new species tree is written "+species_tree_output)

    return species_tree


def check_splice(isoform_by_gene_all):
    total_genes = 0
    total_isoforms = 0
    spliced_species_num = 0
    for sp, isoforms in isoform_by_gene_all.items():
        if len(isoforms)>1:
            spliced_species_num+=1

            for isoform in isoforms:
                if len(isoform)>1:
                    total_isoforms+=len(isoform)
                    total_genes+=1

    logger_hog.debug("For "+str(spliced_species_num)+"  species, out of "+str(len(isoform_by_gene_all))+" , we have splice files.")
    if total_genes ==0:
        logger_hog.debug("It seems that in all of the splice files, each line has only one item. Make sure that the splitter in each line is semicolon ; ")
        sys.exit(1)

    logger_hog.debug("In total, for "+ str(total_genes)+" genes, we have " + str(total_isoforms)+"  splices.")

    # make sys back
    return 1






def check_fastoma_input():
    _config.set_configs()

    # print(_config.in_folder)
    print(_config.logger_level)

    print(_config.species_tree_address)  # nwk format

    result = check_proteome_files()


    species_names, prot_recs_lists, fasta_format_keep = _utils_roothog.parse_proteomes()  # optional input folder
    check_proteome(species_names, prot_recs_lists)
    try:
        species_tree = Tree(_config.species_tree_address, format=1)
    except:
        try:
            species_tree = Tree(_config.species_tree_address)
            # todo add check fro Phyloxml
        except:
            logger_hog.error("We have problem with parsing species tree "+str(_config.species_tree_address) + " using ete3 Tree. Maybe there are some special chars.")


    check_speciestree_internalnode(species_tree)
    check_speciestree_leaves(species_tree,species_names)

    add_internal_node_prune(species_tree,species_names)

    hogmap_files = os.path.exists("./hogmap_in/")
    species_hogmaps =[]
    if hogmap_files:
        species_hogmaps = check_hogmap_files()

    #print("species_hogmaps",species_hogmaps)
    #print("species_names",species_names)
    hogmap_complete = False
    if set(species_hogmaps) == set(species_names):
        hogmap_complete = True

    # print("hogmap_complete", hogmap_complete)
    omamerdb = check_omamer_db()

    if not(omamerdb or hogmap_complete):
        logger_hog.error("OMAmer db does not exit and no hogmap provided. ")
        logger_hog.error("Check input failed. FastOMA halted!")
        sys.exit(1)
    else:
        logger_hog.info("OMAmer db or hogmap exist. It looks ok. ")

    #todo  check splice file format . if the name matches with the proteome files.
    splice_files = os.path.exists("./splice/")
    if splice_files:
        logger_hog.debug("splice folder exist. Let's see its inside.")
        isoform_by_gene_all = _utils_roothog.parse_isoform_file(species_names)
        check_splice(isoform_by_gene_all)
    else:
        logger_hog.info("Splice folder doesn't exist and that's ok.")


    logger_hog.info("Input check finished ! ")

