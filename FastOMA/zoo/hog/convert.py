from xml.etree.ElementTree import XMLParser
__all__ = ["orthoxml_to_newick"]


class TaxonNHXMixin:
    def get_tax_nhx(self):
        tags = []
        if self.level:
            tags.append(":S={}".format(self.level))
        if self.taxid:
            tags.append(":T={}".format(self.taxid))
        return tags


class Speciation:
    type = None

    def __init__(self, parent=None):
        self.level = ""
        self.taxid = None
        self.children = []
        self.parent = parent
        if parent is not None:
            parent.add_child(self)

    def add_child(self, e):
        self.children.append(e)

    def set_level(self, level):
        self.level = level

    def set_taxid(self, taxid):
        self.taxid = taxid

    def get_newick_node_name(self):
        if not hasattr(self, 'get_tax_nhx'):
            return self.level.replace(' ', '_')
        return ""

    def as_nhx(self):
        nhx = "[&&NHX"
        t = ",".join([c.as_nhx() for c in self.children])
        if t != "":
            t = "({})".format(t)
        tags = self.get_newick_node_name()

        if self.type:
            nhx += ":Ev={}".format(self.type)
        if hasattr(self, "get_tax_nhx"):
            nhx += "".join(self.get_tax_nhx())
        nhx += "]"
        if len(nhx) > 7:
            tags += nhx
        return "{}{}".format(t, tags)


class Duplication(Speciation):
    type = "duplication"


class Leaf(Speciation):
    def __init__(self, xref, species, parent=None):
        super().__init__(parent=parent)
        self.name = xref
        self.level = species

    def get_newick_node_name(self):
        return self.name


class NHXSpeciation(Speciation, TaxonNHXMixin):
    pass

class NHXDuplication(Duplication, TaxonNHXMixin):
    pass

class NHXLeaf(Leaf, TaxonNHXMixin):
    pass


class OrthoxmlToNewick:

    def __init__(self, xref_tag="protId", encode_levels_as_nhx=True, return_gene_to_species=False):
        self.xref_tag = xref_tag
        self.gene2xref = {}
        self.trees = {}
        self.depth = 0
        self.famid = None
        self.cur_event = None
        self.cur_species = None
        self._use_nhx = encode_levels_as_nhx
        self._return_gene_to_species= return_gene_to_species

    def start(self, tag, attrib):
        if tag == "{http://orthoXML.org/2011/}species":
            self.cur_species = attrib['name']
        if tag == "{http://orthoXML.org/2011/}gene":
            self.gene2xref[attrib['id']] = (attrib[self.xref_tag], self.cur_species)
        elif tag == "{http://orthoXML.org/2011/}geneRef":
            leaf_cls = NHXLeaf if self._use_nhx else Leaf
            self.cur_event.add_child(leaf_cls(*self.gene2xref[attrib['id']]))
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if self.depth == 0:
                self.famid = attrib['id']
            speciation_cls = NHXSpeciation if self._use_nhx else Speciation
            self.cur_event = speciation_cls(self.cur_event)
            self.depth += 1
        elif tag == "{http://orthoXML.org/2011/}paralogGroup":
            dupl_cls = NHXDuplication if self._use_nhx else Duplication
            self.cur_event = dupl_cls(self.cur_event)
        elif tag == "{http://orthoXML.org/2011/}property":
            if attrib['name'] == "TaxRange":
                self.cur_event.set_level(attrib['value'])
            elif attrib['name'].lower() in ("taxid", "taxonid", "taxon_id", "ncbi_taxon_id"):
                self.cur_event.set_taxid(attrib['value'])

    def end(self, tag):
        if tag == "{http://orthoXML.org/2011/}paralogGroup":
            self.cur_event = self.cur_event.parent
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            self.depth -= 1
            if self.depth == 0:
                assert(self.cur_event.parent is None)
                self.trees[self.famid] = self.cur_event.as_nhx() + ";"
            self.cur_event = self.cur_event.parent

    def close(self):
        if self._return_gene_to_species:
            gene2species = {k[0]: k[1] for k in self.gene2xref.values()}
            return self.trees, gene2species
        return self.trees


def orthoxml_to_newick(filename, xref_tag="protId", encode_levels_as_nhx=False, return_gene_to_species=False):
    """function to convert all HOGs from an orthoxml file into newick trees

    This function converts all toplevel orthologGroups into a dictionary of newick trees.
    Duplication nodes are labeled as such using the nhx tag, e.g. a paralogGroup node
    will be translated into an internal node having the nhx label [&&NHX:Ev=duplication]

    :param filename: the filename of the input orthoxml file

    :param xref_tag: the attribute of the <gene> element that should be used to get as label
                     for the leaves labels.

    :param encode_levels_as_nhx: boolean flag indicating whether or not the species information
                                 of the internal and extend nodes should be returned in NHX format
                                 with the :S=<...> and :T=<...> format. otherwise, the TaxRange
                                 value will be used as newick node label for the internal nodes.

    :param return_gene_to_species: boolean flag indicating if a mapping with the gene to species
                                   should be returned.

    :returns either a dict of {roothogid: tree} where tree is in nhx format or a tuple with the
             first element being the tree dictionary and the second being a mapping from
             {gene: species}.
    """

    target = OrthoxmlToNewick(
        xref_tag=xref_tag,
        encode_levels_as_nhx=encode_levels_as_nhx,
        return_gene_to_species=return_gene_to_species)
    parser = XMLParser(target=target)
    with open(filename, 'rb') as xml:
        for chunk in xml:
            parser.feed(chunk)
    return parser.close()
