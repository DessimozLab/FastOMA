from xml.etree.ElementTree import XMLParser
__all__ = ["orthoxml_to_newick"]


class Speciation:
    type = None

    def __init__(self, parent=None):
        self.level = ""
        self.children = []
        self.parent = parent
        if parent is not None:
            parent.add_child(self)

    def add_child(self, e):
        self.children.append(e)

    def set_level(self, level):
        self.level = level

    def as_nhx(self):
        nhx = ""
        t = ",".join([c.as_nhx() for c in self.children])
        tags = self.level
        if self.type:
            tags +="[&&NHX:Ev={}]".format(self.type)
        return "({}){}".format(t, tags)


class Duplication(Speciation):
    type = "duplication"


class Leaf(Speciation):
    def __init__(self, xref, parent=None):
        super().__init__(parent=parent)
        self.name = xref

    def as_nhx(self):
        return self.name


class OrthoxmlToNewick:

    def __init__(self, xref_tag="protId"):
        self.xref_tag = xref_tag
        self.gene2xref = {}
        self.trees = {}
        self.depth = 0
        self.famid = None
        self.cur_event = None

    def start(self, tag, attrib):
        if tag == "{http://orthoXML.org/2011/}gene":
            self.gene2xref[attrib['id']] = attrib[self.xref_tag]
        elif tag == "{http://orthoXML.org/2011/}geneRef":
            self.cur_event.add_child(Leaf(self.gene2xref[attrib['id']]))
        elif tag == "{http://orthoXML.org/2011/}orthologGroup":
            if self.depth == 0:
                self.famid = attrib['id']
            self.cur_event = Speciation(self.cur_event)
            self.depth += 1
        elif tag == "{http://orthoXML.org/2011/}paralogGroup":
            self.cur_event = Duplication(self.cur_event)

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
        return self.trees


def orthoxml_to_newick(filename, xref_tag="protId"):
    """function to convert all HOGs from an orthoxml file into newick trees

    This function converts all toplevel orthologGroups into a dictionary of newick trees.
    Duplication nodes are labeled as such using the nhx tag, e.g. a paralogGroup node
    will be translated into an internal node having the nhx label [&&NHX:Ev=duplication]

    :param filename: the filename of the input orthoxml file
    :param xref_tag: the attribute of the <gene> element that should be used to get as label
                     for the leaves labels.

    """

    target = OrthoxmlToNewick(xref_tag=xref_tag)
    parser = XMLParser(target=target)
    with open(filename, 'rb') as xml:
        for chunk in xml:
            parser.feed(chunk)
    return parser.close()
