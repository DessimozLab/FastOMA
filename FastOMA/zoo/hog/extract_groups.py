from ..utils import auto_open
import collections
from time import time
import xml.etree.ElementTree as etree
import Bio.Phylo
from pathlib import Path
import logging
logger = logging.getLogger(__name__)


class TaxLevel:
    """TaxLevel class

    This class represents a clade of interest to extract the flat groups.
    You can instantiate it with just a clade parameter, in which case it
    must match exactly the <property name="TaxRange" value="XXX"> value
    of your level of interest.

    A better way to use this class is to pass clade and a phyloxml file
    with the clade definitions (a species tree) and select the clade as
    an internal name of this phylogeny. In this case, all the species
    of the orthoxml must match the extend species in this phylogeny.

    """
    def __init__(self, clade: str, fn_phyloxml=None):
        if fn_phyloxml is not None:
            with open(fn_phyloxml, 'rt') as fh:
                tree = Bio.Phylo.read(fh, 'phyloxml')
            target_clade = next(tree.find_clades(target=clade))
            self.species = frozenset([n.taxonomy.code for n in target_clade.get_terminals()])
            self.taxid = target_clade.taxonomy.id.value
        else:
            logger.warning("TaxLevel only initialized with a clade string. Will only match exact taxrange matches!")
            self.species = set([])
            self.taxid = ""
        self.name = clade

    def all_in_clade(self, genes):
        return set(g.species for g in genes).issubset(self.species)

    def does_match(self, name_or_id):
        return name_or_id == self.name or name_or_id == self.taxid


Gene = collections.namedtuple("Gene", "xref species internal_id")


class GroupExtractor(object):
    def __init__(self, target_clade:TaxLevel, gene_attr="protId"):
        self.processed_stats = {'last': time(), 'processed_toplevel': 0}
        self.target_clade = target_clade
        self.gene_attr = gene_attr
        self.genes = {}

    def add_genome_genes(self, genome_node):
        genome_name = genome_node.get('name', None)
        if genome_name is None:
            genome_name = genome_node.get("NCBITaxId")

        generef_2_xref = {}
        for gene in genome_node.findall('.//{http://orthoXML.org/2011/}gene'):
            gene_id = gene.get('id')
            gene_prot_id = gene.get(self.gene_attr)
            generef_2_xref[gene_id] = Gene(gene_prot_id, genome_name, gene_id)
        self.genes.update(generef_2_xref)
        return True

    def _collect_genes(self, node):
        genes = []
        for child in node.iter():
            if child.tag == "{http://orthoXML.org/2011/}geneRef":
                try:
                    genes.append(self.genes[child.get('id')])
                except KeyError:
                    logger.info(f"ignoring gene(id={child.get('id')}), probably in skip set.")
                    pass
            elif child.tag == "{http://orthoXML.org/2011/}orthologGroup":
                genes.extend((n for n in child.text if isinstance(n, Gene)))
        return genes

    def merge_children(self, node):
        genes = self._collect_genes(node)
        node.clear()
        node.text = genes

    def get_group(self, node):
        return node.text

    def analyze_and_yield_groups(self, node, is_toplevel=False):
        genes = self._collect_genes(node)
        if self.target_clade.all_in_clade(genes):
            if is_toplevel:
                logger.debug("dumping toplevel hog with {} genes (hog_id: {})".format(len(genes), node.get('id')))
                yield genes
            else:
                node.clear()
                node.text = genes
        else:
            for group in node.iter("{http://orthoXML.org/2011/}orthologGroup"):
                if group == node:
                    continue
                if self.target_clade.all_in_clade(group.text):
                    logger.debug("found hog with {} genes".format(len(group.text)))
                    yield group.text
            node.clear()
            node.text = genes


def parse_orthoxml(fh, processor:GroupExtractor):
    nsmap = {}
    og_level = 0
    extract_at_depth = -2
    if processor.target_clade == None:
        extract_at_depth = 0

    def fixtag(tag, ns=""):
        return "{" + nsmap[ns] + "}" + tag

    logger.info("start mapping of orthoxml formatted input file")
    for event, elem in etree.iterparse(fh, events=('start-ns', 'start', 'end')):
        if event == 'start-ns':
            ns, url = elem
            nsmap[ns] = url
        elif event == 'start' and elem.tag == fixtag('orthologGroup'):
            og_level += 1
        elif event == 'start' and elem.tag == fixtag('property') and processor.target_clade is not None:
            if og_level > 0 and elem.get('name').lower() in ("taxid", "taxrange", "taxonomic_range", "taxon_id", "ncbitaxid"):
                if processor.target_clade.does_match(elem.get('value')):
                    extract_at_depth = og_level-1
        elif event == 'end':
            if elem.tag == fixtag('orthologGroup'):
                og_level -= 1
                if extract_at_depth == -2:
                    # no level annotations, we need to check ourself if we are at the right level
                    yield from processor.analyze_and_yield_groups(elem, is_toplevel=(og_level == 0))
                elif extract_at_depth >= 0:
                    # we have taxonomic annotations. combine children and yield group in case
                    # we are at the right level
                    processor.merge_children(elem)
                    if og_level == extract_at_depth:
                        logger.debug("dumping annotated group with {} genes".format(len(elem.text)))
                        yield processor.get_group(elem)
                        elem.clear()
                        extract_at_depth = -1 if processor.target_clade is not None else 0
                if og_level == 0:
                    elem.clear()
            elif elem.tag == fixtag('species'):
                processor.add_genome_genes(elem)
                elem.clear()


def extract_flat_groups_at_level(f, protein_attribute="protId", level=None):
    """Iterates over the groups defined in an orthoxml file at a specific taxonomic level.

    This function yields groups of Genes defined in an orthoxml file for a specific
    taxonomic level. You could for example extract groups at the level of Mammalia, that
    had a duplication event prior to the last common ancestor of the mammalia. For that,
    one needs to pass in a `class:TaxLevel` instance for Mammalia. If the level parameter
    is None (default), all the roothogs are returned.

    :Note:
    If you instantiate the TaxLevel class with only a name, that name needs to match
    the <property name="TaxLevel" value="xxx"> exactly. also, if a taxonomic
    level is not annotated, the parser will *NOT* return the group at this level.
    It is much safer to instantiate the TaxLevel instance with a PhyloXML file and the
    associated clade name - in this case any group that contains only proteins of that
    clade will be returned.

    :Example:
    >>> toplevel_groups = []
    >>> oxml = "tests/coverage_test_files/simpleEx.orthoxml"
    >>> for grp in extract_flat_groups_at_level(oxml):
    >>>     toplevel_groups.append(set(g.xref for g in grp))
    [{"XENTR1", "CANFA1", "HUMAN1", "PANTR1", "MOUSE1", "RATNO1"},
     {"HUMAN2", "PANTR2", "CANFA2", "MOUSE2"},
     {"XENTR3", "MOUSE4", "MOUSE3", "CANFA3", "PANTR4", "PANTR3", "HUMAN3"}]


    :param f: filehandle of filename containing the orthoxml file.
    :param protein_attribute: the protein_attribute that should be read from the orthoxml
                              and that should be returned as xref attribute of the Gene.
                              The default value is the `protId` attribute.
    :param level: a `TaxLevel` instance.
    """
    if level is not None:
        if not isinstance(level, TaxLevel):
            raise ValueError("level argument should be a TaxLevel instance.")
    processor = GroupExtractor(target_clade=level, gene_attr=protein_attribute)
    if isinstance(f, (str, bytes, Path)):
        with auto_open(f, 'rt') as fh:
            yield from parse_orthoxml(fh, processor)
    else:
        yield from parse_orthoxml(f, processor)




