from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numbers

from future.builtins import range
from future.builtins import next
from future.builtins import dict
from future.builtins import zip
from future.builtins import open
from future.builtins import map
from future.builtins import str
from future import standard_library
standard_library.install_hooks()
from past.builtins import basestring

import copy
import io
import itertools
import os
import re
from collections import deque
from functools import reduce

from .orthoxmlquery import OrthoXMLQuery
from .newick import NewickLexer, Streamer
from .tools import PROGRESSBAR, Queue, setup_progressbar


class TaxonomyInconsistencyError(Exception):
    pass


class ParseError(Exception):
    pass


class LinearTaxonomyException(Exception):
    pass


class Taxonomy(object):
    def __init__(self):
        raise NotImplementedError("abstract class")

    def __iter__(self):
        return self.hierarchy[self.root].iter_preorder()

    def __getitem__(self, node_name):
        return self.hierarchy[node_name]

    def finialize_init(self):
        self.extract_descendent_species()
        self.extract_younger_nodes()

    def iterParents(self, node, stopBefore=None):
        """iterates over all the taxonomy nodes towards the root
        which are above 'node' and below 'stopBefore'."""

        if node == stopBefore:
            return
        tn = self.hierarchy[node]
        while tn.up is not None and tn.up.name != stopBefore:
            tn = tn.up
            yield tn.name

    def is_ancestor_of(self, anc, desc):
        """ Returns True if `anc' is an ancestor of `desc'"""
        return anc in self.iterParents(desc)

    def levels_between(self, anc, desc):
        """ Counts the number of levels between anc and desc
        If anc == desc, this is 0. If anc is the immediate parent
        of desc this is 1, grandparent 2 etc. Return -1 if anc is
        not an ancestor of desc."""
        if anc == desc:
            return 0
        if not self.is_ancestor_of(anc, desc): # inefficient, but not *too* inefficient (?)
            return -1
        else:
            parents = list(self.iterParents(desc, stopBefore=anc))
            return len(parents) + 1

    def retain(self, leaves):
        """ Returns a deepcopy of self with all leaves not in `leaves'
        pruned away """
        keep = set(leaves)
        for leaf in leaves:
            keep.update(set(self.iterParents(leaf)))
        lose=set(self.hierarchy.keys()).difference(keep)
        return self._prune_set_of_nodes(lose)

    def prune(self, leaves):
        """ Returns a deepcopy of self with `leaves' pruned away """
        all_leaves = {n.name for n in self.hierarchy[self.root].iter_leaves()}
        keep = all_leaves.difference(leaves)
        return self.retain(keep)

    def _prune_set_of_nodes(self, leaves):
        """ Returns a deepcopy of self with `leaves' pruned away """
        newtax = copy.deepcopy(self)
        for loser in leaves:
            if loser in newtax.hierarchy:
                node = newtax.hierarchy[loser]
                node.up.down.remove(node)
                node.up = None
                del newtax.hierarchy[loser]
        return newtax

    def _countParentAmongLevelSet(self, levels):
        """helper method to count for each level how many levels
        are parent levels. e.g. (arrow: is-parent-of)
          A->B->C
              \>D->E
        will return A=0,B=1,C=2,D=2,E=3
        Cached to speed up multiple calls"""
        if not hasattr(self, '_cpals_cache'):
            self._cpals_cache = {}
        levelSet = frozenset(levels)
        try:
            return self._cpals_cache[levelSet]
        except KeyError:
            counts = dict()
            for lev in levelSet:
                t = set(self.iterParents(lev)).intersection(levelSet)
                counts[lev] = len(t)
            self._cpals_cache[levelSet] = counts
            return counts

    def map_potential_internal_speciesname_to_leaf(self, species_name, code):
        if len(self.descendents[species_name]) == 1:
            return species_name
        elif code is not None:
            return code
        code_candidates = []
        for leaf in self.descendents[species_name]:
            if len(leaf) == 5 and re.match(r'[A-Z][A-Z[0-9]{4}', leaf) is not None:
                code_candidates.append(leaf)
        if len(code_candidates) == 1:
            return code_candidates[0]
        raise TaxonomyInconsistencyError("cannot find a 5letter code species for internal species name")

    def mrca(self, species):
        """Returns most recent common ancestor (MRCA) of a set of species
           This is cached to speed up multiple calls """
        if not hasattr(self, '_mrca_cache'):
            self._mrca_cache = {}
        species = frozenset(species)
        try:
            return self._mrca_cache[species]
        except KeyError:
            if len(species) == 1:
                mrca, = species
                self._mrca_cache[species] = mrca
                return mrca
            ancestors = [set(self.iterParents(s)) for s in species]
            common_ancestors = reduce(lambda x, y: x & y, ancestors)
            mrca = self.mostSpecific(common_ancestors)
            self._mrca_cache[species] = mrca
            return mrca

    def mostSpecific(self, levels):
        """returns the most specific (youngest) level among a set of
        levels. it is required that all levels are on one monophyletic
        lineage, otherwise an Exception is raised."""
        # count how often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        counts = self._countParentAmongLevelSet(levels)
        for lev, count in counts.items():
            if count == len(levels)-1:
                return lev
        raise Exception("Non of the element is subelement of all others")

    def mostGeneralLevel(self, levels):
        """returns the most general (oldest) level among a set of levels."""
        # count who often each element is a child of any other one.
        # the one with len(levels)-1 is the most specific level
        counts = self._countParentAmongLevelSet(levels)
        for lev, count in counts.items():
            if count == 0:
                return lev
        raise Exception("Non of the element is the root of all others")

    def younger_than_filter(self, levels, oldest_permitted):
        """
        Filters a set of levels, removing any that are older than the
        oldest_permitted (oldest_permitted=string: node_name)
        """
        # try:
        #     oldest_permitted_node = self.hierarchy[oldest_permitted]
        # except KeyError:
        #     raise Exception('No node with name {} found in '
        #                     'Taxonomy'.format(oldest_permitted))
        # permitted = [node.name
        #              for node in oldest_permitted_node.iter_preorder()]

        return [lev for lev in levels if lev in self.younger_nodes[oldest_permitted]]

    def printSubTreeR(self, fd, lev=None, indent=0):
        if lev is None:
            lev = self.root
        fd.write("{}{}\n".format(" "*2*indent, lev))
        for child in self.hierarchy[lev].down:
            self.printSubTreeR(fd, child.name, indent+1)

    def get_histories(self, parser, verbosity=0):

        histories = {}

        if PROGRESSBAR and verbosity > 0:
            pbar = setup_progressbar('Getting histories', len(self.hierarchy))
            pbar.start()

        for i, level in enumerate(self.hierarchy, start=1):
            history = parser.getFamHistory()
            history.analyzeLevel(level)
            histories[level] = history
            self.hierarchy[level].attach_fam_history(history)
            if PROGRESSBAR and verbosity > 0:
                pbar.update(i)

        if PROGRESSBAR and verbosity > 0:
            pbar.finish()

        self.histories = histories
        return histories

    def get_comparisons(self, parser, verbosity=0):

        if getattr(self, 'histories', None) is None:
            self.get_histories(parser)

        comparisons = {}
        to_compare = [(node, child) for node in self
                                    for child in node.down]

        if PROGRESSBAR and verbosity > 0:
            pbar = setup_progressbar('Comparing', len(to_compare))
            pbar.start()

        for i, (parent, child) in enumerate(to_compare, start=1):
            parent_history = self.histories[parent.name]
            child_history = self.histories[child.name]
            parent_child_comparison = parent_history.compare(child_history)
            comparisons[(parent.name, child.name)] = parent_child_comparison
            child.attach_level_comparison_result(parent_child_comparison)
            if PROGRESSBAR and verbosity > 0:
                pbar.update(i)

        if PROGRESSBAR and verbosity > 0:
            pbar.finish()

        return comparisons

    def newick(self):
        return str(self.hierarchy[self.root]) + ';'

    def __str__(self):
        fd = io.StringIO()
        self.printSubTreeR(fd)
        res = fd.getvalue()
        fd.close()
        return res

    def extract_descendent_species(self):
        """
        Caches some frequently looked-up information - descendent leaves of every node
        """
        self.descendents = {}
        for k, v in self.hierarchy.items():
            self.descendents[k] = set(l.name for l in v.iter_leaves())

    def extract_younger_nodes(self):
        """
        Caches some frequently looked-up information - descendent nodes of every node
        """
        self.younger_nodes = {}
        for k, v in self.hierarchy.items():
            self.younger_nodes[k] = set(n.name for n in v.iter_preorder())


class LinearTaxonomy(Taxonomy):
    """ Linear taxonomy """

    def __init__(self, taxonomy, histories, comparisons=None):
        # member variable error-checking and setup
        self.check_histories(taxonomy, histories)
        if comparisons is None:
            comparisons = self.generate_comparisons(histories)
        self.check_comparisons(comparisons)
        self.check_comparisons_and_histories(histories, comparisons)

        # build taxonomy
        self.hierarchy = dict()
        self.histories = dict()
        nodes = [TaxNode(h.analyzedLevel) for h in histories]
        self.root = nodes[0].name
        for i in range(len(histories) - 1):
            nodes[i].add_child(nodes[i+1])
            nodes[i+1].add_parent(nodes[i])
            self.hierarchy[nodes[i].name] = nodes[i]
            self.histories[nodes[i].name] = histories[i]
            nodes[i].attach_fam_history(histories[i])
            nodes[i+1].attach_level_comparison_result(comparisons[i])
        nodes[-1].attach_fam_history(histories[-1])
        self.hierarchy[nodes[-1].name] = nodes[-1]
        self.histories[nodes[-1].name] = histories[-1]

    def generate_comparisons(self, histories):
        comparisons = list()
        for i in range(len(histories) - 1):
            comparisons.append(histories[i].compare(histories[i+1]))
        return comparisons

    def check_histories(self, taxonomy, histories):
        labels = [h.analyzedLevel for h in histories]
        for i in range(len(labels) - 1):
            if not taxonomy.is_ancestor_of(labels[i], labels[i+1]):
                raise LinearTaxonomyException('Histories are not in a valid '
                                              'taxonomic sequence')

    def check_comparisons(self, comparisons):
        for i in range(len(comparisons) - 1):
            if not comparisons[i].lev2 == comparisons[i+1].lev1:
                raise LinearTaxonomyException('Comparisons do not overlap')

    def check_comparisons_and_histories(self, histories, comparisons):
        if not (len(comparisons) + 1) == len(histories):
            raise LinearTaxonomyException('Incompatible length lists of '
                                          'histories and comparisons')
        inter = self.interleave(histories, comparisons)
        for i in range(1, len(inter), 2):
            h0 = inter[i-1]
            comp = inter[i]
            h1 = inter[i+1]
            if not (h0.analyzedLevel == comp.lev1 and
                    comp.lev2 == h1.analyzedLevel):
                raise LinearTaxonomyException('Histories and comparisons do '
                                              'not overlap')

    def interleave(self, histories, comparisons):
        return [item for item in itertools.chain(
                    *itertools.zip_longest(histories, comparisons))
                    if item is not None]


# class LinearTaxonomyMultiSpecies(LinearTaxonomy):
#     def __init__(self, taxonomy, histories, species_histories):
#         super().__init__(taxonomy, histories, None)
#         top = histories[-1]
#         species_comparisons = [top.compare(h) for h in species_histories]
#         for i, h in enumerate(species_histories):
#             if not taxonomy.is_ancestor_of(top.analyzedLevel, h.analyzedLevel):
#                 raise LinearTaxonomyException('Species {} is not a descendant '
#                                               'of {}'.format(h.analyzedLevel,
#                                                              top.analyzedLevel))
#             node = TaxNode(h.analyzedLevel)
#             node.attach_fam_history(h)
#             node.attach_level_comparison_result(species_comparisons[i])
#             self.hierarchy[top.analyzedLevel].add_child(node)
#             node.add_parent(self.hierarchy[top.analyzedLevel])
#             self.hierarchy[node.name] = node
#             self.histories[node.name] = h


class XMLTaxonomy(Taxonomy):
    def __init__(self, filename):
        raise NotImplementedError("XML Taxonomies have not "
                                  "yet been implemented")


class TaxRangeOrthoXMLTaxonomy(Taxonomy):
    def __init__(self, parser):
        self.parser = parser
        self.extractAdjacencies()
        self.bloat_all()
        self.extractHierarchy()
        self.finialize_init()

    def _parseParentChildRelsR(self, grp):
        levels = None
        if self.parser.is_ortholog_group(grp):
            levels = [l.get('value') for l in grp.findall(
                './{{{ns0}}}property[@name="TaxRange"]'
                .format(**self.parser.ns))]
        directChildNodes = list(grp)
        children = [child for child in directChildNodes
                    if self.parser.is_evolutionary_node(child)]
        geneRefs = [node for node in directChildNodes
                    if OrthoXMLQuery.is_geneRef_node(node)]
        speciesOfGenes = {self.parser.mapGeneToSpecies(x.get('id'))
                          for x in geneRefs}

        # recursively process childreen nodes
        subLevs = speciesOfGenes
        for child in children:
            subLevs.update(self._parseParentChildRelsR(child))

        if levels is not None:
            for parent in levels:
                for child in subLevs:
                    self.adj.add((parent, child))
            subLevs = set(levels)
        return subLevs

    def extractAdjacencies(self):
        self.adj = set()
        for grp in self.parser.getToplevelGroups():
            self._parseParentChildRelsR(grp)

        self.nodes = set(itertools.chain(*self.adj))

    def bloat_all(self):
        """build transitive closure of all parent - child relations"""
        while(self.bloat()):
            pass

    def bloat(self):
        found = False
        for pair in itertools.product(self.nodes, repeat=2):
            first, second = pair
            for node in self.nodes:
                if ((pair not in self.adj) and ((first, node) in self.adj) and
                        ((node, second) in self.adj)):
                    found = True
                    self.adj.add(pair)
        return found

    def extractHierarchy(self):
        self.hierarchy = dict(zip(self.nodes, map(TaxNode, self.nodes)))
        for pair in itertools.product(self.nodes, repeat=2):
            if pair in self.adj:
                if self.good(pair):
                    first, second = pair
                    #print "%s,%s is good" % pair
                    self.hierarchy[first].add_child(self.hierarchy[second])
                    self.hierarchy[second].add_parent(self.hierarchy[first])
        noParentNodes = [z for z in self.nodes if self.hierarchy[z].up is None]
        if len(noParentNodes) != 1:
            raise TaxonomyInconsistencyError(
                "Warning: several/none TaxonomyNodes are roots: {}"
                .format(noParentNodes))
        self.root = noParentNodes[0]

    def good(self, pair):
        first, second = pair
        for node in self.nodes:
            if (first, node) in self.adj and (node, second) in self.adj:
                return False
        return True


class TaxNode(object):

    reg = re.compile(r'\W') # matches anything that's NOT a-z, A-Z, 0-9 or _

    def __init__(self, name, branch_length=None):
        self.name = name
        self.up = None
        self.down = list()
        self.history = None
        self.comparison = None
        self.branch_length = branch_length

    def __str__(self):
        if self.reg.search(self.name):
            label = '"{}"'.format(self.name)
        else:
            label = self.name
        branch_len_str = ":{}".format(self.branch_length) if self.branch_length else ""

        NHX = self._get_node_NHX() + self._get_edge_NHX()
        if self.is_leaf():
            return '{0}{1}{2}'.format(label, branch_len_str,
                            ('[&&NHX{0}]'.format(NHX) if NHX > '' else ''))
        subtree = ', '.join(str(ch) for ch in self.down)
        return '({0}){1}{2}{3}'.format(subtree, label, branch_len_str,
                            ('[&&NHX{0}]'.format(NHX) if NHX > '' else ''))

    def _get_node_NHX(self):
        return (':Genes={}'.format(len(self.history))
                if self.history else '')

    def _get_edge_NHX(self):
        if self.comparison:
            if not hasattr(self.comparison, 'summary'):
                self.comparison.summarise()

            n_ident     = self.comparison.summary['identical']
            n_dupl      = self.comparison.summary['duplicated']
            n_lost      = self.comparison.summary['lost']
            n_novel     = self.comparison.summary['novel']
            n_singleton = self.comparison.summary['singleton']

            NHX = (':Identical={0}'
                   ':Duplicated={1}'
                   ':Lost={2}'.format(n_ident,
                                      n_dupl,
                                      n_lost))
            if self.is_leaf():
                NHX += ':Novel=0:Singleton={}'.format(n_singleton)
            else:
                NHX += ':Novel={}:Singleton=0'.format(n_novel)

            return NHX
        return ''

    def add_child(self, c):
        if not c in self.down:
            self.down.append(c)

    def add_parent(self, p):
        if self.up is not None and self.up != p:
            raise TaxonomyInconsistencyError(
                "Level {} has several parents, at least two: {}, {}"
                .format(self.name, self.up.name, p.name))
        self.up = p

    def is_leaf(self):
        return len(self.down) == 0

    def is_inner(self):
        return not self.is_leaf()

    def is_root(self):
        return self.up is None

    def iter_preorder(self):
        """ Traverse the tree in preorder ordering: ancestors before
        descendents """
        yield self
        for child in self.down:
            for elem in child.iter_preorder():
                yield elem

    def iter_postorder(self):
        """ Traverse the tree in postorder ordering: descendents before
        ancestors """
        for child in self.down:
            for elem in child.iter_postorder():
                yield elem
        yield self

    def iter_levelorder(self):
        """ Traverse the tree in levelorder ordering: first the root, then
        nodes at distance==1 from the root, then distance==2, etc. """
        q = Queue()
        q.enqueue(self)
        while q:
            node = q.dequeue()
            yield node
            for child in node.down:
                q.enqueue(child)

    def iter_leaves(self):
        for elem in self.iter_preorder():
            if elem.is_leaf():
                yield elem

    def iter_inner_nodes(self):
        for elem in self.iter_preorder():
            if elem.is_inner():
                yield elem

    def attach_fam_history(self, history):
        self.history = history

    def attach_level_comparison_result(self, comparison):
        self.comparison = comparison


class NewickTaxonomy(Taxonomy):

    """ Create a taxonomy from a file or filehandle in newick format. The
    file should contain one tree (further trees are ignored). Only the
    tree topology is used - branch lengths and bootstrap support values
    are thown away. The leaf labels should match those in the orthoXML.
    Inner labels should match too, but for OMA XML will be automatically
    generated if auto_annotate == True """

    def __init__(self, fp):
        if isinstance(fp, basestring):
            if not os.path.exists(fp):
              raise Exception('File not found: {0}'.format(fp))
            fp = open(fp)
        self.lexer = NewickLexer(Streamer(fp))
        self.nodes = set()
        self.hierarchy = {}
        self.stack = []
        self.parse()
        self.finialize_init()

    def _get_label(self, tokens):
        """ Get the node data attributes 'label' and 'length'. Assumes these
        will be the next tokens in the stream. Throws ParseError if they are
        not. """
        label = next(self.lexer)
        if label.typ not in (tokens.LABEL, tokens.SUPPORT):
            raise ParseError(
                'Expected a label or a support value, found {0}'.format(
                    label))

        length = next(self.lexer)
        if length.typ != tokens.LENGTH:
            raise ParseError('Expected a length, found {0}'.format(
                length))

        return (label.val if label.typ == tokens.LABEL else None,
                length.val if length.typ == tokens.LENGTH else None)

    def annotate_from_orthoxml(self, xmlparser):
        """ Transfers internal node names from OMA orthoxml """
        def _annotate(self, node, levels_dict):
            if not node.down:
                self.nodes.add(node.name)
                self.hierarchy[node.name] = node
                return
            key = frozenset((n.name for n in node.iter_leaves()))
            node.name = levels_dict[key]
            self.nodes.add(node.name)
            self.hierarchy[node.name] = node
            for child in node.down:
                _annotate(self, child, levels_dict)

        levels = xmlparser.getLevels()
        levels_dict = dict((frozenset(x.split('/')), x) for x in levels)
        root_node = self.hierarchy[self.root]
        root_node.name = 'LUCA'
        self.hierarchy = {'LUCA': root_node}
        self.nodes = set('LUCA')
        self.root = 'LUCA'

        for child in root_node.down:
            _annotate(self, child, levels_dict)

    def populate(self, root):
        if len(root.down) == 0:
            self.hierarchy[root.name] = root
            self.nodes.add(root.name)
            return

        if not root.up:
            self.root = (root.name or 'LUCA')

        self.hierarchy[root.name] = root
        self.nodes.add(root.name)

        for child in root.down:
            self.populate(child)

    def parse(self):
        tmp_name = 1
        tokens = self.lexer.tokens

        for token in self.lexer:
            if token.typ == tokens.EOF:
                return

            elif token.typ == tokens.TREE:
                # invent a name for the node
                n = TaxNode(str(tmp_name))
                tmp_name += 1

                # push it onto the stack
                self.stack.append(n)

                # set as root
                self.root = n

            elif token.typ == tokens.SUBTREE:
                # invent a name and make a node
                n = TaxNode(str(tmp_name))
                tmp_name += 1

                # get parent
                p = self.stack[-1]
                p.add_child(n)
                n.add_parent(p)

                # push subtree onto the stack
                self.stack.append(n)

            elif token.typ == tokens.LEAF:
                label, branch_len = self._get_label(tokens)
                self.nodes.add(label)
                l = TaxNode(label, branch_len)

                # get parent from stack
                p = self.stack[-1]
                p.add_child(l)
                l.add_parent(p)

            elif token.typ == tokens.ENDSUB:
                label, branch_len = self._get_label(tokens)

                # retrieve node from stack
                subtree = self.stack.pop()

                # update node name
                if isinstance(label, str):
                    subtree.name = label
                if isinstance(branch_len, numbers.Number):
                    subtree.branch_length = branch_len

            elif token.typ == tokens.ENDTREE:  # trigger for tree-finalising functions
                self.populate(self.root)
                del self.lexer
                return

            elif token.typ in (tokens.LABEL, tokens.LENGTH, tokens.SUPPORT):
                raise ParseError('Unexpected token in stream: {0}'.format(token))

            else:
                raise ParseError('Not sure what happened')
