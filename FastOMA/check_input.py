import sys
from ete3 import Tree
import os
import collections

from . import _utils_roothog
from ._wrappers import logger
from . import __version__ as fastoma_version

"""

fastoma-check-input --proteomes proteome --species-tree species_tree.nwk --out-tree species_tree_checked.nwk --hogmap hogmap -vv

"""

def check_proteome_files(proteome_folder):
    proteome_files = os.listdir(proteome_folder)
    logger.info("There are %d files in the proteome folder.", len(proteome_files))

    not_fa = [i for i in proteome_files if not (i.endswith(".fa") or i.endswith(".fasta"))]
    if not_fa:
        logger.warning("We expect that only fa/fasta files are in the proteome folder. Better to remove these %s", not_fa)

    fa_fasta = [i.split(".")[-1] for i in proteome_files if i.endswith(".fa") or i.endswith(".fasta")]
    if len(fa_fasta) <= 2:
        logger.error("There are not enough proteomes in the folder ")
        logger.error("Check input failed. FastOMA halted!")
        return False

    if len(set(fa_fasta)) > 1:
        logger.warning("We expect that all fasta files are with the same format either fa or fasta, we won't include some of them %s.", proteome_files)
        return False

    species_names = [os.path.splitext(i)[0] for i in proteome_files if i.endswith(".fa") or i.endswith(".fasta")]
    return species_names


def check_proteome(species_names, prot_recs_lists, proteome_folder):
    # ids_set=set()
    isOk = True
    num_prots_all = 0
    for species_name in species_names:
        prot_recs_list = prot_recs_lists[species_name]
        prot_recs_ids = [prot_rec.id for prot_rec in prot_recs_list]
        if len(prot_recs_ids)!= len(set(prot_recs_ids)):
            logger.error("It seems that at least of the protein id in fasta file %s is repeated. Note that only part of a fasta record before space is considered.",species_name)
            logger.error("Check input failed. FastOMA halted!")
            return False
        for prot_rec_id in prot_recs_ids:
            if len(prot_rec_id)>60:   # todo 85 is the limit without considering ||s.
                logger.error("The protein ID %s is too long in species %s.fa, which should be changed. Please make sure it will be still unique. ",prot_rec_id,species_name) # (root cause issue due to truncatation by fastTree)
                logger.error("Check input failed. FastOMA halted!")
                return False


        num_prots = len(prot_recs_list) # >EP00159_Fonticula_alba_P004948_XP_009497064.1_small_nuclear_ribonucleoprotein_B_and_B'_Fonticula_alba||691883||1155004948
        num_prots_all += num_prots
        # todo , check duplicated Ids should it be done on all species ?

        if num_prots <= 2:
            logger.error("The input proteome looks empty or too few proteins or Bio.SeqIO couldn't read it, %s/%s.fa", proteome_folder, species_name)
            isOk = False
    # todo write new protoems with cleaned record ids, keep the mapping, to be used in orthoxml writing
    # use the mapping back for orhtoxml
    logger.info("There are %d proteins in total in the input proteome folder.", num_prots_all)
    return isOk


def check_hogmap_files(hogmap_folder):
    hogmap_files = os.listdir(hogmap_folder)
    logger.info("There are %d files in the hogmap folder.", len(hogmap_files))
    not_hogmap = [i for i in hogmap_files if not i.endswith(".hogmap")]
    if not_hogmap:
        logger.warning("We expect that only hogmap files are in the %s folder. Better to remove these %s",
                       hogmap_folder, not_hogmap)

    species_hogmaps = ["".join(i.split(".")[:-2]) for i in hogmap_files if i.endswith(".hogmap")]
    logger.info("There are "+str(len(species_hogmaps))+" hogmaps.")
    return species_hogmaps



def check_speciestree_internalnode(species_tree):

    # All the internal node of the input species tree should have a name
    for node in species_tree.traverse(strategy="postorder"):
        if not node.is_leaf():
            internal_node_name = []
            if (not node.name) or len(node.name) < 1:
                logger.warning("one of the internal node in species tree doesn't have a name. we'll update the species tree.")
            else:
                internal_node_name.append(node.name)

    if len(internal_node_name) != len(set(internal_node_name)):
        logger.warning("All the internal node names should be unique. One of the internal node is repeated:")
        logger.warning([item for item, count in collections.Counter(internal_node_name).items() if count > 1])
        logger.warning("We'll change the internal node names.")
    return 1


def check_speciestree_leaves(species_tree, species_names):
    leaves_name = [i.name for i in species_tree.get_leaves()]
    logger.info("The species tree has %d leaves.", len(leaves_name))
    if len(set(leaves_name)) != len(leaves_name):
        logger.error("Leaves name should be unique in the species tree. You may try ete3 prune")
        logger.error("Check input failed. FastOMA halted!")
        return False

    species_names_not_intree = [i for i in species_names if i not in leaves_name]
    if species_names_not_intree:
        logger.error("These species are not in the species tree:\n %s", "\n ".join(species_names_not_intree))
        logger.error("This is a tree template:  ((AQUAE,CHLTR)inter1,MYCGE)inter2; for proteome folder with these three files:  AQUAE.fa  CHLTR.fa  MYCGE.fa. ")
        logger.error("Check input failed. FastOMA halted!")
        return False

    leaves_tree_not_proteome = [i for i in leaves_name if i not in species_names]
    if leaves_tree_not_proteome:
        logger.warning("there are %d leaves in the species tree that there is no proteome. So we will discard them.", len(leaves_tree_not_proteome))
    return True


def check_omamer_db(omamerdb_adress=None):
    if omamerdb_adress is None:
        logger.warning("omamer_db not passed. will assume you know what you are doing.")
        return True

    if os.path.exists(omamerdb_adress):
        if os.path.getsize(omamerdb_adress) > 10000:  # 3 bytes
            omamerdb = True
            # todo we can do some checks on version omamer v2
        else:
            logger.warning("The omamer db looks very small. are you sure it is correct?"+omamerdb_adress)
            omamerdb = False
    else:
        omamerdb =False
        logger.info("OMAmer db does not exist.")

    return omamerdb


def add_internal_node_prune(species_tree, species_names, out_fname):
    # add name for the internal, if no name is provided, removintg special chars
    counter_internal = 0
    node_names = set()
    for node in species_tree.traverse(strategy="postorder"):
        if not node.is_leaf():
            node_name = node.name

            if len(node_name) < 3 or node_name in node_names:
                if not node.is_leaf():
                    node.name = "internal_" + str(counter_internal)
                    counter_internal += 1
                    node_names.add(node.name)
                    logger.debug("The internal node name was too small or repeated "+node_name+" which is changed to "+node.name)
            elif not all(e.isalnum() or e in ("_", "-") for e in node_name):  #any(not c.isalnum() for c in node_name):
                node_name_new = ''.join(e for e in node_name if e.isalnum() or e in ("_", "-"))  # removing special chars
                if node_name_new in node_names:
                    node_name_new += "_" + str(counter_internal)
                    counter_internal += 1
                logger.debug("The internal node name has special chars or repeated '%s' which is changed to '%s'",
                             node.name, node_name_new)
                node.name = node_name_new
                node_names.add(node.name)
            else:
                node_names.add(node_name)

    logger.info("The species tree has %d leaves", len(species_tree))
    species_tree.prune(species_names)  # , preserve_branch_length=True)
    logger.info("After pruning, the species tree has %d leaves", len(species_tree))

    species_tree_output = "./species_tree_checked.nwk"
    species_tree.write(format=1, format_root_node=True, outfile=out_fname)
    logger.debug("The new species tree is written %s", out_fname)

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

    logger.debug("For "+str(spliced_species_num)+"  species, out of "+str(len(isoform_by_gene_all))+" , we have splice files.")
    if total_genes ==0:
        logger.debug("It seems that in all of the splice files, each line has only one item. Make sure that the splitter in each line is semicolon ; ")
        sys.exit(1)

    logger.debug("In total, for "+ str(total_genes)+" genes, we have " + str(total_isoforms)+"  splices.")

    # make sys back
    return 1


def fastoma_check_input():
    import argparse
    parser = argparse.ArgumentParser(description="checking parameters for FastOMA")
    parser.add_argument("--version", action="version", version="FastOMA v"+fastoma_version)
    parser.add_argument("--proteomes", required=True, help="Path to the folder containing the input proteomes")
    parser.add_argument("--species-tree", required=True, help="Path to the input species tree file in newick format")
    parser.add_argument("--out-tree", required=True, help="Path to output file for sanitised species tree. ")
    parser.add_argument("--splice", help="Path to the folder containing the splice information files")
    parser.add_argument("--hogmap", help="Path to the folder containing the hogmap files")
    parser.add_argument("--omamer_db", help="Path to the omamer database")
    parser.add_argument('-v', action="count", default=0, help="Increase verbosity to info/debug")
    conf = parser.parse_args() # conf_check_input

    logger.setLevel(level=30 - 10 * min(conf.v, 2))
    logger.debug("Arguments: %s", conf)


    species_names = check_proteome_files(conf.proteomes)
    if not species_names:
        logger.error("Halting FastOMA because of invalid proteome input data")
        sys.exit(1)


    try:
        species_tree = Tree(conf.species_tree, format=1)
    except:
        try:
            species_tree = Tree(conf.species_tree)
            # todo add check fro Phyloxml
        except:
            logger.error("We have problem with parsing species tree %s using ete3 Tree. Maybe there are some special chars.", conf.species_tree)
            sys.exit(1)

    check_speciestree_internalnode(species_tree)
    if not check_speciestree_leaves(species_tree, species_names):
        logger.error("Halting FastOMA because of invalid species tree")
        sys.exit(1)
    add_internal_node_prune(species_tree, species_names, conf.out_tree)

    species_names, prot_recs_lists, fasta_format_keep = _utils_roothog.parse_proteomes(conf.proteomes)
    if not check_proteome(species_names, prot_recs_lists, conf.proteomes):
        logger.error("Halting FastOMA because of invalid proteome input data")
        sys.exit(1)

    hogmap_files = conf.hogmap is not None and os.path.exists(conf.hogmap)
    species_hogmaps = []
    if hogmap_files:
        species_hogmaps = check_hogmap_files(conf.hogmap)

    hogmap_complete = False
    if set(species_hogmaps) == set(species_names):
        hogmap_complete = True

    omamerdb = check_omamer_db(conf.omamer_db)

    if not(omamerdb or hogmap_complete):
        logger.error("OMAmer db does not exist and no hogmap provided.")
        logger.error("Check input failed. FastOMA halted!")
        sys.exit(1)
    else:
        logger.info("OMAmer db or hogmap exist. It looks ok.")

    # todo  check splice file format . if the name matches with the proteome files.
    splice_files = conf.splice is not None and os.path.exists(conf.splice)
    if splice_files:
        logger.debug("splice folder exist. Let's see its inside.")
        isoform_by_gene_all = _utils_roothog.parse_isoform_file(species_names, conf.splice)
        check_splice(isoform_by_gene_all)
    else:
        logger.info("Splice folder doesn't exist and that's ok.")
    logger.info("Input check finished ! ")

