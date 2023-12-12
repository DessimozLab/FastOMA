import itertools
from abc import ABCMeta, abstractmethod
from future.utils import with_metaclass
import dendropy as dpy
import lxml.etree

__all__ = ["iter_pairwise_relations"]


class BaseRelationExtractor(with_metaclass(ABCMeta, object)):
    def __init__(self, obj):
        self.obj = obj
        self._rel_type = 'ortholog'
        self.relations = []

    def default_node_list(self):
        return self.obj

    def extract_pairwise_relations(self, node=None, rel_type='ortholog'):
        self.check_rel_type(rel_type)

        root_nodes = [node]
        if node is None:
            root_nodes = self.default_node_list()

        for group in root_nodes:
            self.relations = []
            self._extract_pw(group)
            for rel in self.relations:
                yield rel

    @abstractmethod
    def is_leaf(self, node):
        pass

    def is_internal_node(self, node):
        return not self.is_leaf(node)

    @abstractmethod
    def leaf_label(self, leaf):
        pass

    def _extract_pw(self, node):
        if self.is_leaf(node):
            return {self.leaf_label(node)}
        elif self.is_internal_node(node):
            nodes_of_children = [self._extract_pw(child) for child in self.get_children(node)]
            if self.shall_extract_from_node(node):
                for child1, child2 in itertools.combinations(nodes_of_children, 2):
                    for gId1, gId2 in itertools.product(child1, child2):
                        self.relations.append((gId1, gId2))
            nodes = set.union(*nodes_of_children)
            return nodes
        else:
            return set([])

    @abstractmethod
    def shall_extract_from_node(self, node):
        pass

    @abstractmethod
    def get_children(self, node):
        """return an iterable over the children nodes of `node`"""
        pass

    def check_rel_type(self, rel_type):
        if callable(rel_type):
            self.shall_extract_from_node = rel_type
        elif rel_type not in ('ortholog', 'paralog'):
            raise ValueError('rel_type argument needs to be one of "ortholog" or "paralog"')
        else:
            self._rel_type = rel_type


class OrthoXmlRelationExtractor(BaseRelationExtractor):
    def __init__(self, fh, id_attribute=None, **kwargs):
        super(OrthoXmlRelationExtractor, self).__init__(fh)
        if id_attribute is not None:
            self.geneRef_to_id = {}
            ns = "{http://orthoXML.org/2011/}"
            for g in self.obj.findall(f'./{ns}species/{ns}database/{ns}genes/{ns}gene'):
                try:
                    self.geneRef_to_id[g.get('id')] = g.get(id_attribute)
                except AttributeError:
                    raise KeyError(f"gene {g} does not have an attribute '{id_attribute}'")
        else:
            self.geneRef_to_id = None

    def default_node_list(self):
        return self.obj.find(".//{http://orthoXML.org/2011/}groups")

    def shall_extract_from_node(self, node):
        return lxml.etree.QName(node).localname == self._rel_type + "Group"

    def is_internal_node(self, node):
        return lxml.etree.QName(node).localname in ('orthologGroup', 'paralogGroup')

    def leaf_label(self, leaf):
        res = leaf.get('id')
        if self.geneRef_to_id:
            res = self.geneRef_to_id[res]
        return res

    def is_leaf(self, node):
        return lxml.etree.QName(node).localname == "geneRef"

    def get_children(self, node):
        return node


class TreeRelationExtractor(BaseRelationExtractor):
    def shall_extract_from_node(self, node):
        annos = set(map(node.annotations.get_value, ('Ev', 'D')))
        annos.remove(None)
        if self._rel_type == 'ortholog':
            if len(annos) == 0 or 'S' in annos or 'F' in annos:
                return True
        elif self._rel_type == 'paralog':
            if 'D' in annos or 'T' in annos:
                return True
        return False

    def get_children(self, node):
        return node.child_node_iter()

    def is_leaf(self, node):
        return node.is_leaf()

    def leaf_label(self, leaf):
        return leaf.taxon.label

    def default_node_list(self):
        return [self.obj.seed_node]

    def __init__(self, obj):
        super(TreeRelationExtractor, self).__init__(obj)


def iter_pairwise_relations(obj, rel_type=None, node=None, **kwargs):
    """iterate over the induced pairwise relations from `obj`.

    This function extracts all the induced pairwise relations of type
    `rel_type` from a the hierarchical structure passed as parameter
    `obj`. As hierarchical structures the method accepts either an
    orthoxml file handle or a labeled :class:`dendropy.datamodel.treemodel.Tree`

    rel_type can be used to specify which type of relations should be
    returned. by default, it returns the orthologous relations. The
    argument can be set to 'ortholog', 'paralog' or - if obj is a
    dendropy tree - it can also be a function. This function takes
    as input an internal tree node and should return True iff the
    relations induced by this node should be taken.

    :param obj: a reference to an orthoxml file opened with lxml or a
        dendropy tree.

    :param rel_type: either 'ortholog', 'paralog'. if obj is a
        dendropy tree, rel_type can also or a function

    :param node: a reference to the root node of the substructure
        from which the pairwise relations should be extracted.
        Defaults to the root of the tree or the <groups> tag
        respectively.

    :param id_attribute: applies only for orthoxml files. If parameter
        is provided, instead of the internal geneRef-id, the value that
        is stored in the gene under the given attribute is returned. So
        meaningful values can be `protId` or `geneId`.

    :param kwargs: additional arguments that are passed. So far
        none will be used. Foreseen examples would be a callback
        function on leaf nodes.

    :return: an iterator yielding tuples of leaf ids of the requested
        evolutionary type.

    :Example:

        >>> tree = dpy.Tree.get_from_string(
        ...       '(1,(2,(3,4))[&&NHX:Ev=D],(5,6));', schema='newick')
        >>> sorted(list(iter_pairwise_relations(tree, rel_type='paralog')))
        [('2', '3'), ('2', '4')]
    """
    if hasattr(obj, 'docinfo') and obj.docinfo.root_name == "orthoXML":
        parser = OrthoXmlRelationExtractor(obj, **kwargs)
    elif isinstance(obj, dpy.Tree):
        parser = TreeRelationExtractor(obj)

    else:
        raise ValueError('cannot extract pairwise relations from obj')

    for rel in parser.extract_pairwise_relations(node=node,
                                                 rel_type=rel_type if rel_type is not None else 'ortholog'):
        yield rel
