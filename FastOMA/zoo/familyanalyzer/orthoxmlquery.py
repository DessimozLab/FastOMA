from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future.builtins import str
from future import standard_library
standard_library.install_hooks()


class ElementError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return str(self.msg)


class OrthoXMLQuery(object):
    """Helper class with predefined queries on an orthoxml tree."""

    ns = {"ns0": "http://orthoXML.org/2011/"}   # xml namespace

    @classmethod
    def getToplevelOrthologGroups(cls, root):
        """returns a list with the toplevel orthologGroup elements
        of the given root element."""
        xquery = ".//{{{ns0}}}groups/{{{ns0}}}orthologGroup".format(**cls.ns)
        return root.findall(xquery)

    @classmethod
    def getTaxRangeNodes(cls, root, recursively=True):
        xPrefix = ".//" if recursively else "./"
        xquery = '{}{{{}}}property[@name="TaxRange"]'.format(xPrefix,
                                                             cls.ns['ns0'])
        return root.findall(xquery)

    @classmethod
    def getTaxidNodes(cls, root, recursively=True):
        xPrefix = ".//" if recursively else "./"
        xquery = '{}{{{}}}property[@name="taxid"]'.format(xPrefix, cls.ns['ns0'])
        return root.findall(xquery)

    @classmethod
    def getGeneRefNodes(cls, root, recursively=True):
        iterfn = root.iter if recursively else root.iterchildren
        iterator = iterfn('{{{}}}geneRef'.format(cls.ns['ns0']))
        return list(iterator)

    @classmethod
    def getGeneFromId(cls, id_, root):
        xquery = ".*//{{{}}}gene[@id='{}']".format(cls.ns['ns0'], id_)
        genes = root.findall(xquery)
        if len(genes) > 1:
            raise ElementError('several gene nodes with id {} '
                               'exist'.format(id_))
        gene = genes[0] if len(genes)>0 else None
        return gene

    @classmethod
    def getGroupsAtLevel(cls, level, root):
        """returns a list with the orthologGroup elements which have a
        TaxRange property equals to the requested level."""
        xquery = (".//{{{0}}}property[@name='TaxRange'][@value='{1}']/..".
                  format(cls.ns['ns0'], level))
        return root.findall(xquery)

    @classmethod
    def getSubNodes(cls, targetNode, root, recursively=True):
        """method which returns a list of all (if recursively
        is set to true) or only the direct children nodes
        having 'targetNode' as their tagname.
        The namespace is automatically added to the tagname."""
        xPrefix = ".//" if recursively else "./"
        xquery = "{}{{{}}}{}".format(xPrefix, cls.ns['ns0'], targetNode)
        return root.findall(xquery)

    @classmethod
    def is_geneRef_node(cls, element):
        """check whether a given element is an instance of a geneRef
        element."""
        return element.tag == '{{{ns0}}}geneRef'.format(**cls.ns)

    @classmethod
    def getLevels(cls, element):
        """returns a list of the TaxRange levels associated to the
        passed orthologGroup element. If the element does not have
        any TaxRange property tags associated, an empty list is
        returned."""
        propTags = cls.getSubNodes("property", element, recursively=False)
        res = [t.get('value') for t in propTags if t.get('name') == 'TaxRange']
        return res

    @classmethod
    def getInputGenes(cls, root, species=None):
        """returns a list of all gene elements in the orthoxml inside
        <species><database> tags, i.e. the list of genes prior to running
        OMA-HOGS. Optionally filtered by species."""
        filter_ = ('[@name="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:species{}/ns:database/'
                      'ns:genes//ns:gene'.format(filter_))
        else:
            xquery = '//ns:gene'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})

    @classmethod
    def getGroupedGenes(cls, root, species=None):
        """ returns a list of all geneRef elements inside <group> tags, i.e.
        the list of genes clustered into families after running OMA-HOGS.
        Optionally filtered by species."""
        filter_ = ('[@name="TaxRange"and@value="{}"]'.format(species)
                   if species is not None else '')
        if filter_ > '':
            xquery = ('/ns:orthoXML/ns:groups/ns:orthologGroup//ns:property{}/'
                      'following-sibling::ns:geneRef'.format(filter_))
        else:
            xquery = '//ns:geneRef'
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})

    @classmethod
    def getScoreNodes(cls, root, score_id=None):
        """returns the associated score nodes for a certain (orthologGroup) node.
        If score_id is not specified, all scores will be returned"""
        xquery = './ns:score'
        if score_id is not None:
            xquery += "[@id='{}']".format(score_id)
        return root.xpath(xquery, namespaces={'ns': cls.ns['ns0']})
