import xml.etree.ElementTree as ET
from unittest import TestCase

from ete3 import Tree

from FastOMA._hog_class import (
    HOG,
    attach_scores,
    _member_species,
    _count_implied_losses,
    _induced_congruent_score,
    _tcs_partial,
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

    def test_clean_loss_does_not_penalize_tcs(self):
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
        self.assertEqual(score_value(elem, "TCSScore"), "1.0")

    def test_direct_sibling_species_have_no_tcs_score(self):
        leaf1 = leaf_hog("p1||sp1", self.sptree.search_nodes(name="sp1")[0])
        leaf2 = leaf_hog("p2||sp2", self.sptree.search_nodes(name="sp2")[0])
        hogA = merged_hog([leaf1, leaf2], self.A)

        species_of_members = _member_species(hogA.get_members())
        elem = ET.Element("orthologGroup")
        attach_scores(elem, hogA, self.A, species_of_members)

        self.assertIsNone(score_value(elem, "TCSScore"))
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
