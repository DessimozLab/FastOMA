import collections
import gzip
import os
import pickle
import string
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from ete3 import Tree

from ._wrappers import logger
from ._utils_subhog import read_species_tree
from .transformer import header_transformer
from .zoo.hog import extract_flat_groups_at_level, extract_marker_groups_at_level
# from .zoo.hog.convert import orthoxml_to_newick
from . import __version__ as fastoma_version

"""

fastoma-collect-subhogs --pickle-folder pickle_folders  --roothogs-folder omamer_rhogs \
 --gene-id-pickle-file gene_id_dic_xml.pickle --out FastOMA_HOGs.orthoxml  --marker-groups-fasta OrthologousGroups.tsv \
 --roothog-tsv RootHOGs.tsv  --species-tree  species_tree_checked.nwk  -vv
 
"""

# This code collect subhogs and writes outputs.

def iter_hogs(pickle_folder: Path):
    cnt = 0
    nr_hogs = 0
    logger.info("reading pickle files from %s", pickle_folder)
    for root, dirs, files in os.walk(pickle_folder, followlinks=True):
        for f in files:
            if f.startswith('file_') and f.endswith('.pickle'):
                with open(os.path.join(root, f), 'rb') as handle:
                    cur = pickle.load(handle)
                nr_hogs += len(cur)
                yield from cur
                cnt += 1
                if cnt % 500 == 0:
                    logger.info("read %d batch pickle files so far, resulting in %d roothogs so far",
                                    cnt, nr_hogs)
    logger.info("number of pickle files is %d.", cnt)
    logger.debug("number of hogs in all batches is %d", nr_hogs)


def convert_speciestree_to_orthoxml_taxonomy(tree:Tree):
    cur_id = 1
    name2id = {}
    root = ET.Element("taxonomy")
    for node in tree.traverse(strategy="preorder"):
        try:
            parent_xml = node.up.in_xml
        except AttributeError:
            parent_xml = root
        xml_node = ET.SubElement(parent_xml, "taxon", {'id': str(cur_id), 'name': node.name})
        name2id[node.name] = cur_id
        node.add_feature('in_xml', xml_node)
        cur_id += 1
    return root, name2id


def update_hogids(fam, hog, name2taxid):

    dupCnt = collections.defaultdict(int)
    def _getNextSubId(idx):
        """helper method to return the next number at a given depth of
        duplication (idx)"""
        dupCnt[idx - 1] += 1
        return dupCnt[idx - 1]

    def _encodeParalogClusterId(prefix, nr):
        letters = []
        while nr // 26 > 0:
            letters.append(string.ascii_lowercase[nr % 26])
            nr = nr // 26 - 1
        letters.append(string.ascii_lowercase[nr % 26])
        return prefix + ''.join(letters[::-1])

    def _annotateGroupR(node: ET.ElementTree, og: str, idx: int = 0):
        """create the og attributes at the orthologGroup elements
        according to the naming schema of LOFT. ParalogGroup elements
        do not get own attributes (not possible in the xml schema),
        but propagate their sub-names for the subsequent orthologGroup
        elements."""
        if node.tag == "orthologGroup":
            taxrange = node.find('./property[@name="TaxRange"]')
            taxid = name2taxid[taxrange.get('value')]
            node.set('id', f"{og}_{taxid}")
            node.set('taxonId', str(taxid))
            for child in node:
                _annotateGroupR(child, og, idx)
        elif node.tag == "paralogGroup":
            idx += 1
            nextOG = f"{og}.{_getNextSubId(idx)}"
            for i, child in enumerate(list(node)):
                _annotateGroupR(child, _encodeParalogClusterId(nextOG, i), idx)

    omamer_roothog_id = ":".join(hog.get('id').split("_")[0:2])
    fam_elem = ET.Element("property", {"name": "OMAmerRootHOG", "value": omamer_roothog_id})
    hog.insert(1, fam_elem)
    _annotateGroupR(hog, "HOG:{:07d}".format(fam))
    return hog


def fastoma_collect_subhogs():
    import argparse
    parser = argparse.ArgumentParser(description="collecting all computed HOGs and combine into a single orthoxml")
    parser.add_argument("--version", action="version", version="FastOMA v"+fastoma_version)
    parser.add_argument('--pickle-folder', required=True,
                        help="folder containing the pickle files. Will be searched recursively")
    parser.add_argument('--roothogs-folder', required=True, help="folder containing the omamer roothogs")
    parser.add_argument('--gene-id-pickle-file', required=True,
                        help="file containing the gene-id dictionary in pickle format")
    parser.add_argument('--out', required=False, default="output_hog.orthoxml", help="output filename in orthoxml")
    parser.add_argument('-v', action="count", default=0)
    parser.add_argument('--roothog-tsv', default=None, required=False,
                        help="If specified, a tsv file with the given path will be produced containing the roothog "
                             "assignments as TSV file. In addition, a folder named RootHOGsFasta will be generated"
                             "with one fasta file per inferred RootHOG.")
    parser.add_argument('--marker-groups-fasta', default=None, required=False,
                        help="If specified, a folder named OrthologousFasta and a TSV file with the name provided in "
                             "this argument will be generated that contains single copy groups, i.e. groups which have " 
                             "at most one gene per species. Useful as phylogenetic marker genes to reconstruct species "
                             "trees.")
    parser.add_argument('--species-tree', required=True,
                        help="Path to the species tree used to infer the hogs")
    parser.add_argument('--id-transform', choices=("UniProt", "noop"), default="noop",
                        help="""ID transformer from fasta files to orthoxml / OrthologGroup
                             protein IDs. By default, no transformation will be done. 
                             Existing values are:
                               noop:      No transformation - entire ID of fasta header
                               UniProt:   '>sp|P68250|1433B_BOVIN' --> P68250""")
    conf_collect_subhogs = parser.parse_args() # conf_collect_subhogs
    logger.setLevel(level=30 - 10 * min(conf_collect_subhogs.v, 2))
    logger.debug(conf_collect_subhogs)
    id_transformer = header_transformer(conf_collect_subhogs.id_transform)

    write_hog_orthoxml(conf_collect_subhogs.pickle_folder, conf_collect_subhogs.out, conf_collect_subhogs.gene_id_pickle_file,
                       id_transformer=id_transformer, species_tree=conf_collect_subhogs.species_tree)
    if conf_collect_subhogs.roothog_tsv is not None:
        write_roothogs(Path(conf_collect_subhogs.out), Path(conf_collect_subhogs.roothogs_folder),
                       output_file_roothog_tsv=conf_collect_subhogs.roothog_tsv,
                       id_transformer=id_transformer)
    if conf_collect_subhogs.marker_groups_fasta is not None:
        write_group_files(Path(conf_collect_subhogs.out), Path(conf_collect_subhogs.roothogs_folder),
                          output_file_og_tsv=conf_collect_subhogs.marker_groups_fasta,
                          id_transformer=id_transformer)


def write_hog_orthoxml(pickle_folder, output_xml_name, gene_id_pickle_file, id_transformer, species_tree):
    # todo as input argument/option in nextflow

    # in benchamrk dataset the output prot names should be short
    # tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ
    # for qfo benchmark, the middle should be written in the file

    orthoxml_file = ET.Element("orthoXML", attrib={"xmlns": "http://orthoXML.org/2011/",
                                                   "origin": "FastOMA "+fastoma_version,
                                                   "originVersion": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                   "version": "0.5"})  #

    with open(gene_id_pickle_file, 'rb') as handle:
        gene_id_name = pickle.load(handle)
        # gene_id_name[query_species_name] = (gene_idx_integer, query_prot_name)
    logger.debug("We read the gene_id_name dictionary with %d items", len(gene_id_name))

    speciestree = read_species_tree(species_tree)
    taxonomy, name2taxid = convert_speciestree_to_orthoxml_taxonomy(speciestree)
    logger.debug("Now creating the header of orthoxml")

   #  #### create the header of orthoxml ####
    for query_species_name, list_prots in gene_id_name.items():
        species_xml = ET.SubElement(orthoxml_file, "species", attrib={"name": query_species_name, "taxonId": str(name2taxid[query_species_name]), "NCBITaxId": "0"})
        database_xml = ET.SubElement(species_xml, "database", attrib={"name": "database", "version": "2023"})
        genes_xml = ET.SubElement(database_xml, "genes")
        for (gene_idx_integer, query_prot_name) in list_prots:
            prot_id = id_transformer.transform(query_prot_name)
            gene_xml = ET.SubElement(genes_xml, "gene", attrib={"id": str(gene_idx_integer), "protId": prot_id})
    logger.debug("gene_xml is created.")
    orthoxml_file.append(taxonomy)

    scores = ET.SubElement(orthoxml_file, "scores")
    ET.SubElement(scores, "scoreDef", {"id": "CompletenessScore",
                                       "desc": "Fraction of expected species with genes in the (Sub)HOG"})

    #  #### create the groups of orthoxml   ####
    groups_xml = ET.SubElement(orthoxml_file, "groups")
    for fam, hogs_a_rhog_xml in enumerate(iter_hogs(Path(pickle_folder)), start=1):
        groups_xml.append(update_hogids(fam, hogs_a_rhog_xml, name2taxid))
    logger.debug("converting the xml object to string.")
    with open(output_xml_name, 'wb') as fh:
        ET.indent(orthoxml_file, space='  ', level=0)
        orthoxml = ET.ElementTree(orthoxml_file)
        orthoxml.write(fh, encoding="utf-8", xml_declaration=True, )
    logger.info("orthoxml is written in %s", output_xml_name)


def callback_group_and_omamer(node):
    omamer_node = node.find('./{http://orthoXML.org/2011/}property[@name="OMAmerRootHOG"]')
    omamer = omamer_node.get('value') if omamer_node is not None else 'n/a'
    return {"group_id": node.get('id', 'n/a').split('_')[0],
            "omamer_roothog": omamer}


def write_group_files(orthoxml: Path, roothog_folder: Path, output_file_og_tsv=None, output_fasta_groups=None, id_transformer=None):
    if output_file_og_tsv is None:
        output_file_og_tsv = "OrthologousGroups.tsv"
    if output_fasta_groups is None:
        output_fasta_groups = "OrthologousGroupsFasta"
    output_fasta_groups = Path(output_fasta_groups)
    output_fasta_groups.mkdir(parents=True, exist_ok=True)

    logger.info("Start writing OG tsv and fasta files")
    fasta_format = "fa"  # of the rhogs
    nr_prot_in_groups, nr_groups = 0, 0
    with open(output_file_og_tsv, 'w') as tsv:
        tsv.write("Group\tProtein\n")
        for grp, meta in extract_marker_groups_at_level(orthoxml, protein_attribute="protId", callback=callback_group_and_omamer):
            group_members = {g.xref for g in grp}
            group_name = meta['group_id'].replace("HOG:", "OG_")
            nr_prot_in_groups += len(grp)
            nr_groups += 1
            for gene in group_members:
                tsv.write(f"{group_name}\t{gene}\n")

            _write_group_fasta(fasta_format, group_members, group_name, id_transformer, meta, output_fasta_groups,
                               roothog_folder)
    logger.info("writing of %s done. created %d groups containing %d proteins in total",
                    output_file_og_tsv, nr_groups, nr_prot_in_groups)


def write_roothogs(orthoxml: Path, roothog_folder: Path, output_file_roothog_tsv=None, output_fasta_groups=None, id_transformer=None):
    if output_file_roothog_tsv is None:
        output_file_roothog_tsv = "RootHOGs.tsv"
    if output_fasta_groups is None:
        output_fasta_groups = "RootHOGsFasta"
    output_fasta_groups = Path(output_fasta_groups)
    output_fasta_groups.mkdir(parents=True, exist_ok=True)

    logger.info("Start writing RootHOG tsv and fasta files")
    fasta_format = "fa"  # of the rhogs
    nr_prot_in_groups, nr_groups = 0, 0
    with open(output_file_roothog_tsv, 'wt') as tsv:
        tsv.write("RootHOG\tProtein\tOMAmerRootHOG\n")
        for grp, meta in extract_flat_groups_at_level(orthoxml, callback=callback_group_and_omamer):
            group_members = {g.xref for g in grp}
            group_name = meta['group_id']
            # this is the id of the merged roothogs from the placement step
            omamer_roothog = meta['omamer_roothog']
            nr_prot_in_groups += len(grp)
            nr_groups += 1
            for gene in group_members:
                tsv.write(f"{group_name}\t{gene}\t{omamer_roothog}\n")

            _write_group_fasta(fasta_format, group_members, group_name.replace(":", ""), id_transformer, meta, output_fasta_groups,
                               roothog_folder)

    logger.info("writing of %s done. created %d groups containing %d proteins in total",
                    output_file_roothog_tsv, nr_groups, nr_prot_in_groups)


def _write_group_fasta(fasta_format, group_members, group_name, id_transformer, meta, output_fasta_groups,
                       roothog_folder):
    group_seqs = []
    rhog_fasta = roothog_folder / (meta['omamer_roothog'].replace(':', '_') + "." + fasta_format)
    for rec in SeqIO.parse(rhog_fasta, "fasta"):
        orig_id, sp, *rest = rec.id.split("||")
        protid = id_transformer.transform(orig_id)
        if protid in group_members:
            rec.id = protid
            rec.description += " [" + sp + "]"
            group_seqs.append(rec)
    with gzip.open(output_fasta_groups / (group_name + "." + fasta_format + ".gz"), 'wt') as fa_out:
        SeqIO.write(group_seqs, fa_out, "fasta")


if __name__ == "__main__":
    fastoma_collect_subhogs()