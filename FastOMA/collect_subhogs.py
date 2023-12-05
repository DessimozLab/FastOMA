import os
import pickle
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom

from Bio import SeqIO
from ete3 import Tree

from . import _config
from ._config import logger_hog
from .transformer import header_transformer
from .zoo.hog import extract_flat_groups_at_level
from .zoo.hog.convert import orthoxml_to_newick

# This code collect subhogs and writes outputs.


def load_hogs(pickle_folder: Path):
    hogs_all = []
    cnt = 0
    logger_hog.info("reading pickle files from %s", pickle_folder)
    for root, dirs, files in os.walk(pickle_folder, followlinks=True):
        for f in files:
            if f.startswith('file_') and f.endswith('.pickle'):
                with open(os.path.join(root, f), 'rb') as handle:
                    cur = pickle.load(handle)
                hogs_all.extend(cur)
                cnt += 1
                if cnt % 500 == 0:
                    logger_hog.info("read %d batch pickle files so far, resulting in %d roothogs so far",
                                    cnt, len(hogs_all))
    logger_hog.info("number of pickle files is %d.", cnt)
    logger_hog.debug("number of hogs in all batches is %d", len(hogs_all))
    return hogs_all


def max_og_tree(tree, species_dic):
    for node in tree.traverse("preorder"):
        # for node in xml_tree.traverse(strategy="preorder", is_leaf_fn=lambda n: hasattr(n, "attriremoved") and n.attriremoved==True):
        if not node.is_leaf() and hasattr(node, "Ev") and node.Ev == 'duplication':  # node.name[:3] == "dup"
            dup_node = node
            children = dup_node.get_children()
            list_num_species = []
            for child in children:
                child_name_leaves = child.get_leaves()
                species_list = []
                for leaf in child_name_leaves:
                    name = leaf.name
                    if name[:3] == "no_":
                        name = leaf.name.split("_")[-1]
                    if name in species_dic:
                        species_name = species_dic[name]
                        species_list.append(species_name)
                    else:
                        print("species not in the dic ", name)
                species_set = set(species_list)
                list_num_species.append(len(species_set))
            index_max_species = list_num_species.index(max(list_num_species))
            # if there are few children with identical number of species, the case would be not a polytomi but two children with one species
            # num_occurence = [1 for i in list_num_species if i == max(list_num_species)]
            # if len(num_occurence) > 1:
            #    print("please check this case with the developer the tool. The tree has polytomy.")
            child_max_species = children[index_max_species]
            children_to_remove = [i for i in children if i != child_max_species]
            for child_to_remove in children_to_remove:
                for i in child_to_remove.get_leaves():
                    i.in_og = "no"

    og_prot_list = []
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            if hasattr(node, "in_og") and node.in_og == "no":
                pass  # print(node.name)
            else:
                og_prot_list.append(node.name)

    return og_prot_list


def fastoma_collect_subhogs():
    import argparse
    parser = argparse.ArgumentParser(description="collecting all computed HOGs and combine into a single orthoxml")
    parser.add_argument('--pickle-folder', required=True,
                        help="folder containing the pickle files. Will be searched recursively")
    parser.add_argument('--roothogs-folder', required=True, help="folder containing the omamer roothogs")
    parser.add_argument('--gene-id-pickle-file', required=True,
                        help="file containing the gene-id dictionary in pickle format")
    parser.add_argument('--out', required=False, default="output_hog.orthoxml", help="output filename in orthoxml")
    parser.add_argument('-v', action="count", default=0)
    parser.add_argument('--roothog-tsv', default=None, required=False,
                        help="If specified, a tsv file with the given path will be produced containing the roothog asignments as TSV file.")
    parser.add_argument('--marker-groups-fasta', default=None, required=False,
                        help="If specified, a folder named OrthologousFasta and a TSV file with the name provided in "
                             "this argument will be generated that contains single copy groups, i.e. groups which have " 
                             "at most one gene per species. Useful as phylogenetic marker genes to reconstruct species "
                             "trees.")
    parser.add_argument('--id-transform', choices=("UniProt", "noop"), default="noop",
                        help="""ID transformer from fasta files to orthoxml / OrthologGroup
                             protein IDs. By default, no transformation will be done. 
                             Existing values are:
                               noop:      No transformation - entire ID of fasta header
                               UniProt:   '>sp|P68250|1433B_BOVIN' --> P68250""")
    conf = parser.parse_args()
    logger_hog.setLevel(level=30 - 10 * min(conf.v, 2))
    logger_hog.debug(conf)
    id_transformer = header_transformer(conf.id_transform)

    write_hog_orthoxml(conf.pickle_folder, conf.out, conf.gene_id_pickle_file, id_transformer)
    if conf.roothog_tsv is not None:
        write_roothogs(conf.out, conf.roothog_tsv)
    if conf.marker_groups_fasta is not None:
        write_group_files(Path(conf.out), Path(conf.roothogs_folder),
                          conf.marker_groups_fasta,
                          id_transformer=id_transformer)


def write_hog_orthoxml(pickle_folder, output_xml_name, gene_id_pickle_file, id_transformer):
    # todo as input argument/option in nextflow

    # in benchamrk dataset the output prot names should be short
    # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ
    # for qfo benchmark, the middle should be written in the file

    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/", "origin": "OMA",
                                                   "originVersion": "Nov 2021", "version": "0.3"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
    logger_hog.debug("We read the gene_id_name dictionary with %d items", len(gene_id_name))
    logger_hog.debug("Now creating the header of orthoxml")

   #  #### create the header of orthoxml ####
    for query_species_name, list_prots in gene_id_name.items():
        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "NCBITaxId": "1"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "database ", "version": "2023"})
        genes_xml = ET.SubElement(database_xml, "genes")
        for (gene_idx_integer, query_prot_name) in list_prots:
            prot_id = id_transformer.transform(query_prot_name)
            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": prot_id})
    logger_hog.debug("gene_xml is created.")

    #  #### create the groups of orthoxml   ####
    groups_xml = ET.SubElement(orthoxml_file, "groups")
    for hogs_a_rhog_xml in load_hogs(Path(pickle_folder)):
        groups_xml.append(hogs_a_rhog_xml)
    logger_hog.debug("converting the xml object to string.")

    xml_str = minidom.parseString(ET.tostring(orthoxml_file)).toprettyxml(indent="   ")
    # print(xml_str[:-1000])
    logger_hog.debug("writing orthoxml started")
    with open(output_xml_name, "w") as file_xml:
        file_xml.write(xml_str)
    file_xml.close()

    logger_hog.info("orthoxml is written in %s", output_xml_name)


def write_group_files(orthoxml: Path, roothog_folder: Path, output_file_og_tsv=None, output_fasta_groups=None, id_transformer=None):
    if output_file_og_tsv is None:
        output_file_og_tsv = "OrthologousGroups.tsv"
    if output_fasta_groups is None:
        output_fasta_groups = "OrthologousGroupsFasta"
    output_fasta_groups = Path(output_fasta_groups)


    logger_hog.info("Start writing OG fasta files ")

    fasta_format = "fa"  # of the rhogs

    trees, species_dic = orthoxml_to_newick(orthoxml, return_gene_to_species=True)
    logger_hog.info("We extracted %d trees in NHX format from the input HOG orthoxml %s", len(trees),  orthoxml)

    OGs = {}
    for hog_id, tree_string in trees.items():
        try:
            tree = Tree(tree_string, format=1)
        except:
            try:
                tree = Tree(tree_string, format=1, quoted_node_names=True)
            except:
                logger_hog.error("Error loading tree", tree_string)
                raise

        og_prot_list = max_og_tree(tree, species_dic)
        if len(og_prot_list) >= 2:  # a group should have at least 1 member/protein
            OGs[hog_id] = og_prot_list

    logger_hog.info("done extracting trees")

    with open(output_file_og_tsv, 'w') as handle:
        handle.write("Group\tProtein\n")
        for hog_id, og_prot_list in OGs.items():
            for prot in og_prot_list:
                handle.write(f"OG_{hog_id}\t{prot}\n")

    logger_hog.info("We wrote the protein families information in the file %s", output_file_og_tsv)

    output_fasta_groups.mkdir(parents=True, exist_ok=True)
    logger_hog.info("start writing %d OGs as fasta files in folder %s.", len(OGs), output_fasta_groups)
    for hog_id, og_prot_list in OGs.items():  # hog_id="HOG_0667494_sub10524"
        rhog_id = "_".join(hog_id.split("_")[:2])

        rhog_fasta = roothog_folder / (rhog_id + "." + fasta_format)
        omamer_rhogs_all_prots = list(SeqIO.parse(rhog_fasta, "fasta"))

        og_prots = []
        og_prot_list = OGs[hog_id]
        for rhogs_prot in omamer_rhogs_all_prots:
            orig_id, sp, *rest = rhogs_prot.id.split("||")
            prot_id = id_transformer.transform(orig_id)
            if prot_id in og_prot_list:
                rhogs_prot.id = prot_id
                rhogs_prot.description += " [" + sp + "]"
                og_prots.append(rhogs_prot)

        og_id = "OG_" + hog_id  # one OG per rootHOG      # "/HOG_"+ rhogid
        SeqIO.write(og_prots, output_fasta_groups / (og_id + ".fa"), "fasta")
    logger_hog.info("writing done")


def write_roothogs(orthoxml, out):
    if out is None:
        out = "rootHOGs.tsv"

    nr_prot_in_groups = 0
    with open(out, 'wt') as fh:
        fh.write("Group\tProtein\n")
        for nr, grp in enumerate(extract_flat_groups_at_level(orthoxml)):
            # todo instead of integer nr, use HOGid
            for gene in grp:
                fh.write(f"{nr+1}\t{gene.xref}\n")
                nr_prot_in_groups += 1
    logger_hog.info("Nr roothogs: %d; nr proteins in all roothogs: %d", nr+1, nr_prot_in_groups)
    logger_hog.info("We wrote the protein families information in the file %s", out)


if __name__ == "__main__":
    fastoma_collect_subhogs()