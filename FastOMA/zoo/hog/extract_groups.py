import itertools

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


def _callback_group(node):
    return node.get('id')


class GroupExtractor(object):
    def __init__(self, target_clade:TaxLevel, gene_attr="protId", callback=None):
        self.processed_stats = {'last': time(), 'processed_toplevel': 0}
        self.target_clade = target_clade
        self.gene_attr = gene_attr
        self.genes = {}
        self.callback = callback
        if callback == "group_id":
            self.callback = _callback_group

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

    def _count_genes(self, node):
        count = 0
        for gene in node.iter('{http://orthoXML.org/2011/}geneRef'):
            count += 1
        for og in node.iter('{http://orthoXML.org/2011/}orthologGroup'):
            for n in og.text:
                if isinstance(n, Gene):
                    count += 1
        return count

    def _collect_genes(self, node):
        genes = set([])
        if node.tag != "{http://orthoXML.org/2011/}orthologGroup":
            raise RuntimeError("_collect_genes() only works for ortholog groups")
        for child in node.iter():
            if child == node:
                continue
            if child.tag == "{http://orthoXML.org/2011/}geneRef":
                try:
                    genes.add(self.genes[child.get('id')])
                except KeyError:
                    logger.info(f"ignoring gene(id={child.get('id')}), probably in skip set.")
                    pass
            elif child.tag == "{http://orthoXML.org/2011/}orthologGroup":
                genes.update((n for n in child.text if isinstance(n, Gene)))
        return genes

    def merge_children(self, node):
        genes = self._collect_genes(node)
        for child in reversed(node):
            if child.tag in ("{http://orthoXML.org/2011/}orthologGroup", "{http://orthoXML.org/2011/}geneRef", "{http://orthoXML.org/2011/}paralogGroup"):
                node.remove(child)
        node.text = genes

    def get_group(self, node):
        if self.callback is not None:
            return node.text, self.callback(node)
        return node.text

    def analyze_and_yield_groups(self, node, is_toplevel=False):
        genes = self._collect_genes(node)
        if self.target_clade.all_in_clade(genes):
            if is_toplevel:
                logger.debug("dumping toplevel hog with {} genes (hog_id: {})".format(len(genes), node.get('id')))
                if self.callback is not None:
                    yield genes, self.callback(node)
                else:
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
                    yield self.get_group(group)
            node.clear()
            node.text = genes

    def handle_duplication_node(self, elem):
        pass


class MarkerGroupExtractor(GroupExtractor):
    def handle_duplication_node(self, elem):
        nr_children = [self._count_genes(child) for child in elem]
        max_pos = nr_children.index(max(nr_children))
        to_rem = [c for i, c in enumerate(elem) if i != max_pos]
        for child in to_rem:
            elem.remove(child)


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
                        if len(elem.text) > 1:
                            yield processor.get_group(elem)
                        else:
                            logger.debug("won't return group of less than two proteins")
                        elem.clear()
                        extract_at_depth = -1 if processor.target_clade is not None else 0
                if og_level == 0:
                    elem.clear()
            elif elem.tag == fixtag('paralogGroup'):
                processor.handle_duplication_node(elem)
            elif elem.tag == fixtag('species'):
                processor.add_genome_genes(elem)
                elem.clear()


def extract_group_at_level(f, processor):
    if isinstance(f, (str, bytes, Path)):
        with auto_open(f, 'rt') as fh:
            yield from parse_orthoxml(fh, processor)
    else:
        yield from parse_orthoxml(f, processor)


def extract_flat_groups_at_level(f, protein_attribute="protId", level=None, callback=None):
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
    >>> for grp, grp_id in extract_flat_groups_at_level(oxml, callback="group_id"):
    ...     toplevel_groups.append((set(g.xref for g in grp), grp_id))
    [({"XENTR1", "CANFA1", "HUMAN1", "PANTR1", "MOUSE1", "RATNO1"}, "1"),
     ({"HUMAN2", "PANTR2", "CANFA2", "MOUSE2"}, "2"),
     ({"XENTR3", "MOUSE4", "MOUSE3", "CANFA3", "PANTR4", "PANTR3", "HUMAN3"}, "3")]


    :param f: filehandle of filename containing the orthoxml file.

    :param protein_attribute: the protein_attribute that should be read from the orthoxml
                              and that should be returned as xref attribute of the Gene.
                              The default value is the `protId` attribute.

    :param level: a `TaxLevel` instance.

    :param callback: a callback function that accepts the node element of the group that
                     contains all the members of the desired clade.
                     As callback, the user can also pass the string `group_id` in which
                     case the group_id will be returned.

    :returns: If no callback is provided, the method yields lists of Gene elements. Otherwise
              it yields tuples with (List[Gene], return_value_of_callback)

    """
    if level is not None:
        if not isinstance(level, TaxLevel):
            raise ValueError("level argument should be a TaxLevel instance.")
    processor = GroupExtractor(target_clade=level, gene_attr=protein_attribute, callback=callback)
    yield from extract_group_at_level(f, processor)


def extract_marker_groups_at_level(f, protein_attribute="protId", level=None, callback=None):
    """
    Iterates over the groups defined in an orthoxml file at a specific taxonomic level.

    This function yields flat groups of Genes defined in an orthoxml file for a specific
    taxonomic level, which are all orthologous to each other. This is done by always selecting
    the most heavy child group for paralogGroup nodes. You could for example extract groups
    at the level of Mammalia, that had a duplication event prior to the last common ancestor of the mammalia. For that,
    one needs to pass in a `class:TaxLevel` instance for Mammalia. If the level parameter
    is None (default), markers from all the roothogs are returned.

    :Note:
    If you instantiate the TaxLevel class with only a name, that name needs to match
    the <property name="TaxLevel" value="xxx"> exactly. also, if a taxonomic
    level is not annotated, the parser will *NOT* return the group at this level.
    It is much safer to instantiate the TaxLevel instance with a PhyloXML file and the
    associated clade name - in this case any group that contains only proteins of that
    clade will be returned.

    :Example:
    >>> markers = []
    >>> oxml = "tests/coverage_test_files/simpleEx.orthoxml"
    >>> for grp, grp_id in extract_marker_groups_at_level(oxml, callback="group_id"):
    ...     markers.append((set(g.xref for g in grp), grp_id))
    [({"XENTR1", "CANFA1", "HUMAN1", "PANTR1", "MOUSE1", "RATNO1"}, "1"),
     ({"HUMAN2", "PANTR2", "CANFA2", "MOUSE2"}, "2"),
     ({"XENTR3", "MOUSE3", "CANFA3", "PANTR3", "HUMAN3"}, "3")]


    :param f: filehandle of filename containing the orthoxml file.

    :param protein_attribute: the protein_attribute that should be read from the orthoxml
                              and that should be returned as xref attribute of the Gene.
                              The default value is the `protId` attribute.
    :param level: a `TaxLevel` instance.

    :param callback: a callback function that accepts the node element of the group that
                     contains all the members of the desired clade.
                     As callback, the user can also pass the string `group_id` in which
                     case the group_id will be returned.

    :returns: If no callback is provided, the method yields lists of Gene elements. Otherwise
              it yields tuples with (List[Gene], return_value_of_callback)
    """
    if level is not None:
        if not isinstance(level, TaxLevel):
            raise ValueError("level argument should be a TaxLevel instance.")
    processor = MarkerGroupExtractor(target_clade=level, gene_attr=protein_attribute, callback=callback)
    yield from extract_group_at_level(f, processor)