import os
import xml.etree.ElementTree as ET
from unittest import TestCase

from ete3 import Tree

from FastOMA._hog_class import (
    HOG,
    attach_scores,
    _species_mrca,
    _species_name_index,
    _species_lineage_index,
    _member_species,
    _count_implied_losses,
    _hog_implied_losses,
    _tax_overlap,
    _combine_tax_overlap,
)


def leaf_hog(prot_id, species_node):
    hog = HOG.__new__(HOG)
    hog._rhogid = "test"
    hog._hogid = "HOG_test_" + prot_id
    hog._tax_now = species_node
    hog._members = {prot_id}
    hog._subhogs = []
    hog._dubious_members = set()
    return hog


def merged_hog(subhogs, tax_node):
    hog = HOG.__new__(HOG)
    hog._rhogid = "test"
    hog._hogid = "HOG_test_" + tax_node.name
    hog._tax_now = tax_node
    hog._subhogs = list(subhogs)
    hog._members = set().union(*(sh._members for sh in subhogs))
    hog._dubious_members = set()
    return hog


def score_value(elem, score_id):
    for score in elem.findall("score"):
        if score.get("id") == score_id:
            return score.get("value")
    return None


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_species_tree_from_lineage_tsv(path):
    """Builds a species tree from a lineage TSV (columns: species name, comma-separated
    "<id> (<rank>)" entries from root to species), e.g. tests/data/tcs_taxonomy.tsv."""
    nodes = {}
    root = None
    with open(path) as fh:
        for line in fh:
            species, _, lineage = line.rstrip("\n").partition("\t")
            if not lineage or species == "query":
                continue
            parent = None
            for entry in lineage.split(","):
                node_id = entry.strip().split(" ", 1)[0]
                node = nodes.get(node_id)
                if node is None:
                    node = Tree(name=node_id)
                    nodes[node_id] = node
                    if parent is not None:
                        parent.add_child(node)
                    elif root is None:
                        root = node
                parent = node
            parent.name = species  # rename the lineage's terminal node to the species name
    for n in root.traverse():
        n.add_feature("size", len(n))
    return root


def load_hog_from_genetree(nwk_path, species_tree):
    """Builds a nested HOG mirroring a gene tree's topology, treating that topology as
    already-reconciled sub-HOG grouping: each internal HOG's taxonomic level is the species-tree
    MRCA of the species under it."""
    gene_tree = Tree(nwk_path)

    def build(node):
        if node.is_leaf():
            species = node.name
            return leaf_hog(f"{species}_p||{species}", species_tree.search_nodes(name=species)[0])
        subhogs = [build(child) for child in node.children]
        members = set().union(*(sh._members for sh in subhogs))
        tax_now = _species_mrca(species_tree, _member_species(members))
        return merged_hog(subhogs, tax_now)

    return build(gene_tree)


class SpeciesNameIndexCachingTests(TestCase):
    """_species_name_index()/_species_mrca() replace ete3's search_nodes()/get_common_ancestor-
    by-name (each an O(species tree size) linear scan -- costly once you're doing one lookup per
    gene on a rootHOG with 100k+ genes) with a name->TreeNode dict, built once and cached on the
    tree's root."""

    def setUp(self):
        self.sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        for n in self.sptree.traverse():
            n.add_feature("size", len(n))

    def test_index_maps_every_name_to_the_right_node(self):
        index = _species_name_index(self.sptree)
        for name in ("M", "A", "B", "sp1", "sp2", "sp3", "sp4"):
            node = index[name]
            self.assertEqual(node.name, name)
            self.assertIs(node, self.sptree.search_nodes(name=name)[0])

    def test_index_is_built_once_and_cached_on_the_root(self):
        first = _species_name_index(self.sptree)
        second = _species_name_index(self.sptree)
        self.assertIs(first, second)
        # also reachable (and the same cache) starting from a non-root node
        deep_node = self.sptree.search_nodes(name="A")[0]
        self.assertIs(_species_name_index(deep_node), first)

    def test_mrca_via_index_matches_mrca_via_search_nodes(self):
        expected = self.sptree.get_common_ancestor(
            self.sptree.search_nodes(name="sp1")[0], self.sptree.search_nodes(name="sp3")[0]
        )
        self.assertIs(_species_mrca(self.sptree, {"sp1", "sp3"}), expected)


class TaxOverlapUnitTests(TestCase):
    """Unit tests for the building blocks of the Moi/Kim taxonomy-overlap TCSScore
    (_species_lineage_index, _tax_overlap, _combine_tax_overlap), isolated from the rest of
    attach_scores (CompletenessScore/ImpliedLosses)."""

    def setUp(self):
        self.sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        for n in self.sptree.traverse():
            n.add_feature("size", len(n))

    def test_species_lineage_is_absolute_not_relative_to_any_mrca(self):
        # unlike the old mrca-relative lineage, this includes the species tree's own root (M)
        # and the species' own name, regardless of what node the index happens to be requested
        # from.
        index = _species_lineage_index(self.sptree)
        self.assertEqual(index["sp1"], frozenset({"sp1", "A", "M"}))
        self.assertEqual(index["sp3"], frozenset({"sp3", "B", "M"}))
        deep_node = self.sptree.search_nodes(name="A")[0]
        self.assertEqual(_species_lineage_index(deep_node), index)

    def test_tax_overlap_leaf_size_comes_from_nrmembers_not_a_recount(self):
        # a leaf hog whose _members set happens to hold 2 entries (e.g. a dubious merged
        # fragment) should report leaf_size=2 straight from NrMembers/len(hog), not "1 per leaf
        # node" the way the reference script's minimal tree class assumes.
        leaf = leaf_hog("p1||sp1", self.sptree.search_nodes(name="sp1")[0])
        leaf._members = {"p1||sp1", "p1b||sp1"}
        self.assertEqual(len(leaf), 2)
        nset, leaf_size, leaf_acc, tax_score = _tax_overlap(leaf, _species_lineage_index(self.sptree))
        self.assertEqual(leaf_size, 2)
        self.assertEqual((nset, leaf_acc, tax_score), (frozenset({"sp1", "A", "M"}), 0, 0))

    def test_combine_ignores_an_empty_nset_instead_of_zeroing_the_whole_intersection(self):
        # a part with an empty nset (an already fully-incongruent branch) must not poison an
        # otherwise-matching sibling by naive-intersecting it down to nothing.
        matching_a = (frozenset({"X", "Y"}), 1, 0, 0)
        matching_b = (frozenset({"X", "Y", "Z"}), 1, 0, 0)
        incongruent = (frozenset(), 1, 0, 0)

        nset, leaf_size, leaf_acc, tax_score = _combine_tax_overlap([matching_a, matching_b, incongruent])
        self.assertEqual(nset, frozenset({"X", "Y"}))  # incongruent's empty nset was excluded
        self.assertEqual(leaf_size, 3)
        self.assertEqual(leaf_acc, 3)  # 0+0+0 (children) + leaf_size
        self.assertEqual(tax_score, len({"X", "Y"}) * 3)  # 0+0+0 (children) + len(nset)*leaf_size


class BalancedSpeciesTreeTests(TestCase):
    """Species tree: M(A(sp1,sp2), B(sp3,sp4))"""

    def setUp(self):
        self.sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        for n in self.sptree.traverse():
            n.add_feature("size", len(n))
        self.M = self.sptree
        self.A = self.sptree.search_nodes(name="A")[0]
        self.B = self.sptree.search_nodes(name="B")[0]

    def _make_full_hog(self):
        leaf1 = leaf_hog("p1||sp1", self.sptree.search_nodes(name="sp1")[0])
        leaf2 = leaf_hog("p2||sp2", self.sptree.search_nodes(name="sp2")[0])
        leaf3 = leaf_hog("p3||sp3", self.sptree.search_nodes(name="sp3")[0])
        leaf4 = leaf_hog("p4||sp4", self.sptree.search_nodes(name="sp4")[0])
        hogA = merged_hog([leaf1, leaf2], self.A)
        hogB = merged_hog([leaf3, leaf4], self.B)
        return merged_hog([hogA, hogB], self.M)

    def test_fully_congruent_hog(self):
        hogM = self._make_full_hog()
        species_of_members = _member_species(hogM.get_members())
        self.assertEqual(species_of_members, {"sp1", "sp2", "sp3", "sp4"})

        elem = ET.Element("orthologGroup")
        attach_scores(elem, hogM, self.M, species_of_members)

        self.assertEqual(score_value(elem, "CompletenessScore"), "1.0")
        self.assertEqual(score_value(elem, "ImpliedLosses"), "0")
        self.assertEqual(score_value(elem, "TCSScore"), "1.0")

    def test_clean_loss_lowers_the_tcs_score(self):
        # TCSScore follows the Moi/Kim taxonomy-overlap formula (see attach_scores' docstring):
        # losing sp2's whole branch removes a genuine nested match, so unlike CompletenessScore
        # (species-level presence) it does drop here relative to test_fully_congruent_hog's
        # 1.0: hogA becomes a bare pass-through of leaf1 (nset={sp1,A,M}, tax_score=0), hogB
        # unchanged (nset={B,M}, leaf_size=2, tax_score=4); folded: leaf_size=3, leaf_acc=5,
        # tax_score=4+len({M})*3=7 -> (7 - 5*1) / 3 = 0.6667.
        leaf1 = leaf_hog("p1||sp1", self.sptree.search_nodes(name="sp1")[0])
        leaf3 = leaf_hog("p3||sp3", self.sptree.search_nodes(name="sp3")[0])
        leaf4 = leaf_hog("p4||sp4", self.sptree.search_nodes(name="sp4")[0])
        hogA = merged_hog([leaf1], self.A)  # sp2 lost, single sub-hog passed through
        hogB = merged_hog([leaf3, leaf4], self.B)
        hogM = merged_hog([hogA, hogB], self.M)

        species_of_members = _member_species(hogM.get_members())
        self.assertEqual(species_of_members, {"sp1", "sp3", "sp4"})

        elem = ET.Element("orthologGroup")
        attach_scores(elem, hogM, self.M, species_of_members)

        self.assertEqual(score_value(elem, "CompletenessScore"), "0.75")
        self.assertEqual(score_value(elem, "ImpliedLosses"), "1")
        self.assertEqual(score_value(elem, "TCSScore"), "0.6667")

    def test_direct_sibling_species_score_zero(self):
        # sp1 and sp2's absolute lineages share {A, M} (nset, len 2), giving tax_score =
        # len(nset)*leaf_size = 2*2 = 4 -- but leaf_acc*len(nset) = 2*2 = 4 too (a plain 2-leaf,
        # single-copy family earns exactly enough "self" bonus to cancel its own match), so
        # (4 - 4) / 2 = 0.0.
        leaf1 = leaf_hog("p1||sp1", self.sptree.search_nodes(name="sp1")[0])
        leaf2 = leaf_hog("p2||sp2", self.sptree.search_nodes(name="sp2")[0])
        hogA = merged_hog([leaf1, leaf2], self.A)

        species_of_members = _member_species(hogA.get_members())
        elem = ET.Element("orthologGroup")
        attach_scores(elem, hogA, self.A, species_of_members)

        self.assertEqual(score_value(elem, "TCSScore"), "0.0")
        self.assertIsNotNone(score_value(elem, "CompletenessScore"))
        self.assertIsNotNone(score_value(elem, "ImpliedLosses"))


class ImpliedLossesUnitTests(TestCase):
    def test_whole_missing_clade_counts_as_one_loss(self):
        sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        losses = _count_implied_losses(sptree, {"sp3", "sp4"})
        self.assertEqual(losses, 1)

    def test_no_losses_when_all_present(self):
        sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        losses = _count_implied_losses(sptree, {"sp1", "sp2", "sp3", "sp4"})
        self.assertEqual(losses, 0)

    def test_single_species_loss_within_a_clade(self):
        sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        losses = _count_implied_losses(sptree, {"sp1", "sp3", "sp4"})
        self.assertEqual(losses, 1)


class HierarchicalImpliedLossesUnitTests(TestCase):
    """_hog_implied_losses() recurses through _subhogs (unlike the flat _count_implied_losses,
    which only sees a single flattened species set) so that a loss hidden behind a sibling
    paralog's presence still gets counted -- while not double-counting a plain (non-duplicated)
    nested loss that _count_implied_losses would already catch on its own."""

    def setUp(self):
        self.sptree = Tree("((sp1,sp2)A,(sp3,sp4)B)M;", format=1)
        for n in self.sptree.traverse():
            n.add_feature("size", len(n))
        self.M, self.A, self.B = self.sptree, *(self.sptree.search_nodes(name=n)[0] for n in ("A", "B"))

    def _sp(self, name):
        return self.sptree.search_nodes(name=name)[0]

    def test_matches_flat_count_without_duplication(self):
        # sp2 lost under A, no duplication anywhere -- single nested loss, counted once (not
        # once per ancestor level).
        hogA = merged_hog([leaf_hog("p1||sp1", self._sp("sp1"))], self.A)
        hogB = merged_hog([leaf_hog("p3||sp3", self._sp("sp3")), leaf_hog("p4||sp4", self._sp("sp4"))], self.B)
        hogM = merged_hog([hogA, hogB], self.M)
        self.assertEqual(_hog_implied_losses(hogM, self.M), 1)

    def test_duplication_loss_hidden_behind_sibling_paralog_is_still_counted(self):
        # two paralogous copies under B (both tax_now=B): copy1 complete (sp3,sp4), copy2 only
        # in sp3. Flattened species of B are still {sp3,sp4} (copy1 covers sp4), so the flat
        # _count_implied_losses would report 0 -- but copy2 truly lost its sp4 copy.
        copy1 = merged_hog([leaf_hog("p3a||sp3", self._sp("sp3")), leaf_hog("p4a||sp4", self._sp("sp4"))], self.B)
        copy2 = merged_hog([leaf_hog("p3b||sp3", self._sp("sp3"))], self.B)
        hogB = merged_hog([copy1, copy2], self.B)

        flat_species = _member_species(hogB.get_members())
        self.assertEqual(flat_species, {"sp3", "sp4"})
        self.assertEqual(_count_implied_losses(self.B, flat_species), 0)

        self.assertEqual(_hog_implied_losses(hogB, self.B), 1)


class SpeciesMrcaOnPrunedTreeTests(TestCase):
    """A rootHOG never containing a given species causes prepare_species_tree() to prune that
    species out of the working tree entirely. Scores must still be computed against the real,
    unpruned species tree topology so such species count as losses."""

    def setUp(self):
        self.full_tree = Tree("((AQUAE,CHLTR)inter1,MYCGE)inter2;", format=1)
        for n in self.full_tree.traverse():
            n.add_feature("size", len(n))
        self.full_tree_copy = self.full_tree.copy()

        # simulate prepare_species_tree(): rootHOG only ever contained CHLTR and MYCGE, so AQUAE
        # (and the now-singleton inter1 node) is pruned away from the working tree.
        self.pruned = self.full_tree.get_common_ancestor({"CHLTR", "MYCGE"})
        self.pruned.prune({"CHLTR", "MYCGE"}, preserve_branch_length=True)

    def test_pruned_tree_hides_the_missing_species(self):
        self.assertEqual(sorted(c.name for c in self.pruned.children), ["CHLTR", "MYCGE"])
        self.assertEqual(_count_implied_losses(self.pruned, {"CHLTR", "MYCGE"}), 0)

    def test_mrca_on_the_full_tree_recovers_the_real_topology(self):
        node = _species_mrca(self.full_tree_copy, {"CHLTR", "MYCGE"})
        self.assertEqual(node.name, self.pruned.name)  # same node identity as the pruned mrca ...
        self.assertEqual(sorted(c.name for c in node.children), ["MYCGE", "inter1"])  # ... but real topology
        self.assertEqual(_count_implied_losses(node, {"CHLTR", "MYCGE"}), 1)

    def test_attach_scores_reports_the_loss_when_scored_against_the_full_tree(self):
        elem = ET.Element("orthologGroup")
        species_of_members = {"CHLTR", "MYCGE"}
        node = _species_mrca(self.full_tree_copy, species_of_members)
        hog = merged_hog(
            [leaf_hog("p1||CHLTR", self.full_tree_copy.search_nodes(name="CHLTR")[0]),
             leaf_hog("p2||MYCGE", self.full_tree_copy.search_nodes(name="MYCGE")[0])],
            node,
        )
        attach_scores(elem, hog, node, species_of_members)

        self.assertEqual(score_value(elem, "CompletenessScore"), "0.6667")
        self.assertEqual(score_value(elem, "ImpliedLosses"), "1")
        self.assertEqual(score_value(elem, "TCSScore"), "0.0")


class GeneTreeMathTCSTests(TestCase):
    """Exercises attach_scores()/_tax_overlap() as pure scoring functions against
    tests/data/tcs_tree*.nwk mapped directly onto species-tree MRCAs via
    load_hog_from_genetree(), using the asymmetric, unbalanced species tree built from
    tests/data/tcs_taxonomy.tsv (unlike the balanced synthetic tree in BalancedSpeciesTreeTests).

    NOTE: a real FastOMA-inferred HOG's sub-hog nesting always follows the species tree exactly
    -- infer_hogs_for_rhog_levels_recursively() recurses over the real species tree, and
    reconciliation only decides whether sibling sub-hogs at a given species-tree node get merged
    (speciation) or kept apart (duplication), never reordering across species-tree clades. So
    tcs_tree2/tcs_tree3 below -- whose topologies genuinely cross species-tree clades -- are not
    structures the real pipeline could ever produce as a HOG's sub-hog tree; they're here purely
    to check the score arithmetic in isolation. Only tcs_tree1 (fully congruent) reflects a
    realistic pipeline output."""

    def setUp(self):
        self.sptree = load_species_tree_from_lineage_tsv(os.path.join(DATA_DIR, "tcs_taxonomy.tsv"))

    def _tcs_score_of(self, tree_filename):
        hog = load_hog_from_genetree(os.path.join(DATA_DIR, tree_filename), self.sptree)
        species_of_members = _member_species(hog.get_members())
        mrca = _species_mrca(self.sptree, species_of_members)

        elem = ET.Element("orthologGroup")
        attach_scores(elem, hog, mrca, species_of_members)
        return score_value(elem, "TCSScore")

    def test_gene_tree_congruent_with_species_tree_scores_highest(self):
        # tcs_tree1: ((SP1,SP2),(SP3,SP4)) matches the real topology -- (SP1,SP2) under genus 4,
        # (SP3,SP4) under order 7 -- exactly. The one case realistically reachable by the
        # pipeline, and (see the other two tests below) the highest-scoring of the three.
        self.assertEqual(self._tcs_score_of("tcs_tree1.nwk"), "2.0")

    def test_gene_tree_partially_incongruent_scores_lower(self):
        # tcs_tree2: (((SP1,SP2),SP3),SP4) nests SP3 with the (SP1,SP2) clade instead of with
        # SP4 as the real species tree has it, so it recovers less nested structure than
        # tcs_tree1's perfect match (1.5 < 2.0) -- but more than tcs_tree3's total mismatch.
        self.assertEqual(self._tcs_score_of("tcs_tree2.nwk"), "1.5")

    def test_gene_tree_fully_incongruent_scores_0(self):
        # tcs_tree3: ((SP1,SP3),(SP2,SP4)) crosses both real clades ((SP1,SP2) and (SP3,SP4)),
        # so no split in the gene tree agrees with the species tree at all -- every internal
        # node's nset collapses to just the family's own top-level match, contributing nothing
        # beyond what a single-copy family would already get for free (same cancellation as
        # test_direct_sibling_species_score_zero, just one level deeper): 0.0.
        self.assertEqual(self._tcs_score_of("tcs_tree3.nwk"), "0.0")


class DuplicationAsymmetricRetentionTCSTests(TestCase):
    r"""A duplication scenario that *is* reachable by the real pipeline (unlike tcs_tree2/tcs_tree3
    above): under the order-7 clade, two paralogous copies exist -- copy1 retained cleanly in
    both SP3 and SP4, copy2 retained only in SP3 (lost in SP4).

    Species tree:
                              1 (class)
                             /         \
                     2 (order)           7 (order)
                        |                /        \
                    3 (family)    8 (family)    11 (family)
                        |              |              |
                    4 (genus)    9 (genus)     12 (subfamily)
                     /    \            |              |
                  SP1    SP2    10 (species)   13 (genus)
                                       |              |
                                      SP3       14 (species)
                                                       |
                                                15 (subspecies)
                                                       |
                                                      SP4

    HOG topology built in the test method below -- copy1 and copy2 both get _tax_now=7 (node
    "7"), matching how the real pipeline's merge_subhogs() always bumps every surviving hog's
    _tax_now to the level it was last processed at, whether or not it found a merge partner
    there. That shared _tax_now is what marks them as a recognized duplication (paralogGroup):

                              1  (top HOG)
                             /            \
                     4 (hogA: SP1,SP2)      7  (hog7: duplication, paralogGroup)
                                            /                    \
                              copy1 (tax_now=7)              copy2 (tax_now=7)
                              /          \                        |
                            SP3          SP4                  SP3 (2nd paralog,
                        (protA_SP3)  (protA_SP4)              protB_SP3 -- no SP4
                                                                 counterpart: lost)

    Every species is still present overall (some copy reaches both SP3 and SP4), so
    CompletenessScore reads 1.0 -- it only sees species-level presence, not per-paralog
    retention. But ImpliedLosses is computed hierarchically (_hog_implied_losses), treating
    copy1 and copy2 as independent lineages so copy1's presence in SP4 can't mask copy2's loss
    there: 1 implied loss (copy2 never reaching SP4).

    TCSScore (see attach_scores' docstring for the formula) works out, per _tax_overlap(), to:
      hogA:  nset={4,3,2,1}      leaf_size=2  leaf_acc=2   tax_score=0+4*2=8
      copy1: nset={7,1}          leaf_size=2  leaf_acc=2   tax_score=0+2*2=4
      copy2: nset={SP3,9,8,7,1}  leaf_size=1  leaf_acc=0   tax_score=0            (bare pass-through)
      hog7:  fold(copy1,copy2)   leaf_size=3  leaf_acc=5   tax_score=4+len({7,1})*3=10
      top:   fold(hogA,hog7)     leaf_size=5  leaf_acc=12  tax_score=18+len({1})*5=23
    TCSScore = (23 - 12*len({1})) / 5 = (23 - 12) / 5 = 2.2."""

    def setUp(self):
        self.sptree = load_species_tree_from_lineage_tsv(os.path.join(DATA_DIR, "tcs_taxonomy.tsv"))

    def _node(self, name):
        return self.sptree.search_nodes(name=name)[0]

    def test_duplication_with_asymmetric_paralog_retention(self):
        hogA = merged_hog(
            [leaf_hog("p1||SP1", self._node("SP1")), leaf_hog("p2||SP2", self._node("SP2"))],
            self._node("4"),
        )
        copy1 = merged_hog(
            [leaf_hog("p3a||SP3", self._node("SP3")), leaf_hog("p4a||SP4", self._node("SP4"))],
            self._node("7"),
        )
        # duplicate of copy1, lost in SP4; wrapped at tax_now=7 like copy1 (see docstring) so
        # both are recognized as one duplication group rather than two unrelated lineages.
        copy2 = merged_hog([leaf_hog("p3b||SP3", self._node("SP3"))], self._node("7"))
        hog7 = merged_hog([copy1, copy2], self._node("7"))
        top = merged_hog([hogA, hog7], self._node("1"))

        species_of_members = _member_species(top.get_members())
        self.assertEqual(species_of_members, {"SP1", "SP2", "SP3", "SP4"})
        mrca = _species_mrca(self.sptree, species_of_members)

        elem = ET.Element("orthologGroup")
        attach_scores(elem, top, mrca, species_of_members)

        self.assertEqual(score_value(elem, "CompletenessScore"), "1.0")
        self.assertEqual(score_value(elem, "ImpliedLosses"), "1")
        self.assertEqual(score_value(elem, "TCSScore"), "2.2")
