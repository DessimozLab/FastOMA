import io
import pickle
import xml.etree.ElementTree as ET

import pytest

from FastOMA._hog_class import HOG, Representative, split_hog, root_orthoxml, ScoreFlags, _species_mrca
from subhog_helpers import pid, rec, msa, species_tree, leaf_hog, merged_hog, leaf_members

NO_SCORES = ScoreFlags(False, False, False)

A1, A2, A3 = pid("a1", "A", 1), pid("a2", "A", 2), pid("a3", "A", 3)
B1, B2 = pid("b1", "B", 11), pid("b2", "B", 12)
C1 = pid("c1", "C", 21)


@pytest.fixture
def sp():
    return species_tree()


def test_leaf_hog(sp):
    h = leaf_hog(A1, sp)
    assert h.get_members() == {A1}
    assert len(h) == 1
    assert h.taxname == "A"
    assert h.hogid == "HOG_r_sub10001"
    assert [r.get_id() for r in h.get_representatives()] == [A1]
    assert [r.id for r in h.get_msa()] == [A1]
    assert h._subhogs == [] and h.get_dubious_members() == set()
    assert h.active


def test_hog_stores_only_the_name_of_its_level(sp):
    # HOGs are pickled after every taxonomic level; they must not drag a copy of the species tree along
    top = merged_hog([leaf_hog(A1, sp), leaf_hog(B1, sp)], sp & "AB")
    assert top.taxname == "AB" and [h.taxname for h in top._subhogs] == ["A", "B"]

    class NoTreePickler(pickle.Pickler):
        # fails on any ete object anywhere in the pickled object graph, whatever the ete version
        def reducer_override(self, obj):
            if type(obj).__module__.split(".")[0] in ("ete3", "ete4"):
                raise AssertionError(f"tree object leaked into pickle: {type(obj)}")
            return NotImplemented

    NoTreePickler(io.BytesIO()).dump(top)


def test_hog_requires_level_name_not_tree_node(sp):
    with pytest.raises(TypeError):
        HOG(rec(A1), sp & "A", "r")


def test_hogids_are_sequential_per_process(sp):
    ids = [leaf_hog(x, sp).hogid for x in (A1, A2)]
    assert ids == ["HOG_r_sub10001", "HOG_r_sub10002"]


def test_merged_hog_members_are_union_of_subhogs(sp):
    subs = [leaf_hog(A1, sp), leaf_hog(B1, sp), leaf_hog(C1, sp)]
    h = merged_hog(subs, sp, rep_ids=[A1, C1])
    assert h.get_members() == {A1, B1, C1} == leaf_members(h)
    assert h._subhogs == subs
    # the HOG keeps only the representatives' rows of the MSA
    assert sorted(r.id for r in h.get_msa()) == [A1, C1]
    assert A1 in h and h.get_representatives()[0] in h and subs[1] in h
    assert "x||A||9" not in h


def test_get_members_returns_a_copy(sp):
    h = leaf_hog(A1, sp)
    h.get_members().add("x")
    assert h.get_members() == {A1}


def test_merged_hog_rejects_duplicate_representatives(sp):
    subs = [leaf_hog(A1, sp), leaf_hog(B1, sp)]
    with pytest.raises(AssertionError, match="Non-unique"):
        HOG(subs, sp.name, "r", msa=msa(rec(A1), rec(B1)), representatives=[Representative(rec(A1)), Representative(rec(A1))])


def test_merged_hog_rejects_representative_not_member(sp):
    subs = [leaf_hog(A1, sp), leaf_hog(B1, sp)]
    with pytest.raises(AssertionError, match="subset"):
        HOG(subs, sp.name, "r", msa=msa(rec(A1), rec(B1), rec(C1)), representatives=[Representative(rec(C1))])


def test_merged_hog_with_empty_msa_is_rejected(sp):
    # characterises current behaviour: an empty MSA is falsy, so the constructor hits its generic error branch
    with pytest.raises(AssertionError):
        HOG([leaf_hog(A1, sp)], sp.name, "r", msa=msa(), representatives=[Representative(rec(A1))])


def test_merged_hog_column_filter_keeps_representatives(sp, subhog_conf):
    s1 = "A" * 300 + "-" * 300
    s2 = "-" * 300 + "A" * 300
    subs = [leaf_hog(A1, sp, seq=s1), leaf_hog(B1, sp, seq=s2)]
    h = merged_hog(subs, sp, conf=subhog_conf())
    assert len(h) == 2
    assert len(h.get_msa()) == 2


def test_merged_hog_column_filter_on_wide_msa(sp, subhog_conf):
    s1 = "A" * 550 + "-" * 50
    s2 = "A" * 600
    subs = [leaf_hog(A1, sp, seq=s1), leaf_hog(B1, sp, seq=s2)]
    h = merged_hog(subs, sp, conf=subhog_conf())
    # columns with gap ratio 0.5 are removed (threshold is strict)
    assert [len(r) for r in h.get_msa()] == [550, 550]


def test_get_subhog_path(sp):
    a1, b1, c1 = leaf_hog(A1, sp), leaf_hog(B1, sp), leaf_hog(C1, sp)
    ab = merged_hog([a1, b1], sp & "AB")
    top = merged_hog([ab, c1], sp)
    assert top.get_subhog_path(B1) == [ab, b1]
    assert top.get_subhog_path(B1, max_depth=1) == [ab]
    assert top.get_subhog_path(C1) == [c1]
    with pytest.raises(KeyError):
        top.get_subhog_path("x||A||9")


# Representative

def test_representative_from_seqrecord():
    r = Representative(rec(A1))
    assert r.get_id() == A1
    assert r.get_subelements() == {A1}
    assert r.get_species() == {"A"}
    assert r.enabled


def test_representative_with_elements_collects_subelements_and_species():
    a, b = Representative(rec(A1)), Representative(rec(B1))
    r = Representative(a, [a, b])
    assert r.get_id() == A1
    assert r.get_subelements() == {A1, B1}
    assert r.get_species() == {"A", "B"}


def test_representative_with_elements_does_not_mutate_source():
    a, b = Representative(rec(A1)), Representative(rec(B1))
    Representative(a, [a, b])
    assert a.get_subelements() == {A1}
    assert a.get_species() == {"A"}


# find_representative

def test_find_representative_by_id(sp):
    h = merged_hog([leaf_hog(A1, sp), leaf_hog(B1, sp)], sp)
    assert h.find_representative(B1).get_id() == B1


def test_merged_hog_needs_representatives(sp):
    with pytest.raises(ValueError, match="representatives"):
        HOG([leaf_hog(A1, sp), leaf_hog(B1, sp)], sp.name, "r", msa=msa(rec(A1), rec(B1)))


def test_find_representative_of_non_representative_member(sp):
    a1, b1 = leaf_hog(A1, sp), leaf_hog(B1, sp)
    rep = Representative(rec(A1), [Representative(rec(B1))])
    h = HOG([a1, b1], sp.name, "r", msa=msa(rec(A1), rec(B1)), representatives=[rep])
    assert h.find_representative(B1) is rep


# split_hog

def _two_level_hog(sp):
    """root HOG = [AB HOG (a1, b1; only a1 is a representative), c1]."""
    a1, b1, c1 = leaf_hog(A1, sp), leaf_hog(B1, sp), leaf_hog(C1, sp)
    ab = merged_hog([a1, b1], sp & "AB", rep_ids=[A1], aln=msa(rec(A1), rec(B1), rec(C1)))
    top = merged_hog([ab, c1], sp, rep_ids=[A1, C1], aln=msa(rec(A1), rec(B1), rec(C1)))
    return top, ab, c1


def test_split_hog_on_disjoint_first_level_subhogs(sp):
    top, ab, c1 = _two_level_hog(sp)
    parts = split_hog(top, "root", [A1], [C1])
    assert len(parts) == 2
    assert [h.get_members() for h in parts] == [{A1, B1}, {C1}]
    assert [h._subhogs for h in parts] == [[ab], [c1]]
    assert all(h.taxname == top.taxname for h in parts)
    assert [[r.get_id() for r in h.get_representatives()] for h in parts] == [[A1], [C1]]


def test_split_hog_requires_two_partitions(sp):
    top, _, _ = _two_level_hog(sp)
    with pytest.raises(ValueError):
        split_hog(top, "root", [A1])


def test_split_hog_rejects_overlapping_partitions(sp):
    top, _, _ = _two_level_hog(sp)
    with pytest.raises(ValueError):
        split_hog(top, "root", [A1], [A1, C1])


def test_split_hog_conflicting_partitions_returns_none(sp):
    a1, b1, b2 = leaf_hog(A1, sp), leaf_hog(B1, sp), leaf_hog(B2, sp)
    ab = merged_hog([a1, b1], sp & "AB")
    top = merged_hog([ab, b2], sp)
    # characterises current behaviour: a1 and b1 were merged below this level, so the split is
    # "not implemented" and the caller (_resolve_conflicts) has to deal with it
    assert split_hog(top, "root", [A1], [B1, B2]) is None


def test_split_hog_keeps_all_members(sp):
    a1, a2, b1 = leaf_hog(A1, sp), leaf_hog(A2, sp), leaf_hog(B1, sp)
    ab = merged_hog([a1, a2, b1], sp & "AB", rep_ids=[A1, A2])
    parts = split_hog(ab, "AB", [A1], [A2])
    assert set().union(*(h.get_members() for h in parts)) == {A1, A2, B1}


# to_orthoxml

def _tags(elem):
    return [c.tag for c in elem]


def _taxrange(elem):
    return [p.get("value") for p in elem.findall("property") if p.get("name") == "TaxRange"]


def test_to_orthoxml_nested(sp):
    top, ab, c1 = _two_level_hog(sp)
    e = top.to_orthoxml(sp, NO_SCORES)
    assert e.tag == "orthologGroup" and e.get("id") == top.hogid
    assert _taxrange(e) == ["root"]
    assert _tags(e) == ["property", "orthologGroup", "geneRef"]
    inner = e.find("orthologGroup")
    assert inner.get("id") == ab.hogid and _taxrange(inner) == ["AB"]
    assert sorted(g.get("id") for g in inner.findall("geneRef")) == ["1", "11"]
    assert e.find("geneRef").get("id") == "21"


def test_to_orthoxml_with_scores(sp):
    top, _, _ = _two_level_hog(sp)
    e = top.to_orthoxml(sp)
    assert [s.get("id") for s in e.findall("score")] == ["CompletenessScore", "ImpliedLosses", "TCSScore"]
    assert _tags(e)[:4] == ["score", "score", "score", "property"]


@pytest.mark.parametrize("pickled", [False, True], ids=["in-memory", "pickled-per-level"])
@pytest.mark.parametrize("present,expected_losses", [(["A", "B", "C", "D"], 0), (["A", "C"], 2)],
                         ids=["complete", "two-clades-lost"])
def test_implied_losses_against_full_species_tree(present, expected_losses, pickled):
    # as in read_infer_xml_rhog: HOGs live on the pruned working tree, scores use an unpruned copy, and
    # infer_hogs_this_level writes each level to its own pickle file
    import pickle
    full = species_tree("((A,B)AB,(C,D)CD)root;")
    work = full.copy()
    work.prune(present)

    def level(h):
        return pickle.loads(pickle.dumps([h]))[0] if pickled else h

    leaves = [level(leaf_hog(pid(s.lower(), s, i), work)) for i, s in enumerate(present, 1)]
    if len(present) == 4:
        top = merged_hog([level(merged_hog(leaves[:2], work & "AB")), level(merged_hog(leaves[2:], work & "CD"))], work)
    else:
        top = merged_hog(leaves, work)
    e = level(top).to_orthoxml(full)
    assert e.find("score[@id='ImpliedLosses']").get("value") == str(expected_losses)


def test_implied_losses_with_duplication_independent_of_tree_copies():
    # ((A,B)AB,(C,D)CD)root with a duplication at AB: copy 1 has A and B, copy 2 only A; C present, D lost.
    # The value must not depend on whether the levels share one tree object.
    import pickle

    def build(tree, level):
        a1, a2 = (level(leaf_hog(pid(n, "A", i), tree)) for n, i in (("a1", 1), ("a2", 2)))
        b1, c1 = level(leaf_hog(pid("b1", "B", 3), tree)), level(leaf_hog(pid("c1", "C", 4), tree))
        copy1 = level(merged_hog([a1, b1], tree & "AB"))
        copy2 = level(merged_hog([a2], tree & "AB"))
        return level(merged_hog([copy1, copy2, c1], tree))

    def losses(top, scoring_tree):
        return top.to_orthoxml(scoring_tree).find("score[@id='ImpliedLosses']").get("value")

    single = species_tree("((A,B)AB,(C,D)CD)root;")
    expected = losses(build(single, lambda h: h), single)
    full = species_tree("((A,B)AB,(C,D)CD)root;")
    work = full.copy()
    pickled = build(work, lambda h: pickle.loads(pickle.dumps([h]))[0])
    assert losses(pickled, full) == expected == "2"   # D lost, and B lost in copy 2


def test_to_orthoxml_inparalogs_give_paralog_group(sp):
    h = merged_hog([leaf_hog(A1, sp), leaf_hog(A2, sp)], sp & "A")
    e = h.to_orthoxml(sp, NO_SCORES)
    assert e.tag == "paralogGroup"
    assert sorted(g.get("id") for g in e) == ["1", "2"]


def test_to_orthoxml_duplication_below_speciation(sp):
    a1, a2, b1 = leaf_hog(A1, sp), leaf_hog(A2, sp), leaf_hog(B1, sp)
    h = merged_hog([a1, a2, b1], sp & "AB")
    e = h.to_orthoxml(sp, NO_SCORES)
    assert e.tag == "orthologGroup" and _taxrange(e) == ["AB"]
    assert _tags(e) == ["property", "paralogGroup", "geneRef"]


def test_to_orthoxml_taxrange_is_mrca_of_members(sp):
    # HOG lives at the root level but only has members from A and B
    h = merged_hog([leaf_hog(A1, sp), leaf_hog(B1, sp)], sp)
    e = h.to_orthoxml(sp, NO_SCORES)
    assert _taxrange(e) == ["AB"]


def test_to_orthoxml_single_subhog_is_not_wrapped(sp):
    a1, b1 = leaf_hog(A1, sp), leaf_hog(B1, sp)
    ab = merged_hog([a1, b1], sp & "AB")
    top = merged_hog([ab], sp)
    e = top.to_orthoxml(sp, NO_SCORES)
    assert e.get("id") == ab.hogid


def test_to_orthoxml_singleton_is_generef(sp):
    e = leaf_hog(A1, sp).to_orthoxml(sp, NO_SCORES)
    assert e.tag == "geneRef" and e.get("id") == "1"


def test_to_orthoxml_merged_fragments_leaf(sp):
    h = leaf_hog(A1 + "_|_" + A2, sp)
    e = h.to_orthoxml(sp, NO_SCORES)
    assert e.tag == "paralogGroup"
    assert [(p.get("name"), p.get("value")) for p in e.findall("property")] == [("Type", "DubiousMergedfragment")]
    assert [g.get("id") for g in e.findall("geneRef")] == ["1", "2"]


def test_to_orthoxml_every_member_once(sp):
    a1, a2, b1, b2, c1 = (leaf_hog(x, sp) for x in (A1, A2, B1, B2, C1))
    ab1 = merged_hog([a1, b1], sp & "AB")
    ab2 = merged_hog([a2, b2], sp & "AB")
    top = merged_hog([ab1, ab2, c1], sp)
    xml = ET.tostring(top.to_orthoxml(sp))
    ids = [g.get("id") for g in ET.fromstring(xml).iter("geneRef")]
    assert sorted(ids) == sorted(m.split("||")[2] for m in top.get_members())


# ---------------------------------------------------------------- top-level groups (root_orthoxml)
# A HOG that nothing joins at a level is still wrapped in a new HOG at that level, so a HOG starting
# with a duplication usually reaches the top of the rootHOG inside several such wrappers (issue #128).

PRIMATES = "(((HUMAN,CHIMP)Hominini,MACAQUE)Primates,MOUSE)Mammals;"
H1, H2 = pid("h1", "HUMAN", 1), pid("h2", "HUMAN", 2)
CH1, CH2 = pid("c1", "CHIMP", 3), pid("c2", "CHIMP", 4)


def _top_scores(element):
    return {s.get("id"): s.get("value") for s in element.findall("score")}


def _top_taxrange(element):
    return next(p.get("value") for p in element.findall("property") if p.get("name") == "TaxRange")


def test_species_mrca_of_a_single_species_is_its_leaf():
    t = species_tree(PRIMATES)
    assert _species_mrca(t, {"HUMAN"}) is t & "HUMAN"
    assert _species_mrca(t, {"HUMAN", "CHIMP"}) is t & "Hominini"


def test_root_group_of_a_duplication_sits_directly_above_the_duplicated_clade():
    # duplication in the Hominini lineage; macaque and mouse have no copy
    t = species_tree(PRIMATES)
    copy1 = merged_hog([leaf_hog(H1, t), leaf_hog(CH1, t)], t & "Hominini")
    copy2 = merged_hog([leaf_hog(H2, t), leaf_hog(CH2, t)], t & "Hominini")
    at_primates = merged_hog([copy1, copy2], t & "Primates")
    at_mammals = merged_hog([at_primates], t)   # nothing joins it at Mammals

    for hog in (at_primates, at_mammals):
        e, = root_orthoxml(hog, t)
        assert e.tag == "orthologGroup" and e.get("id") == hog.hogid
        assert _top_taxrange(e) == "Primates"
        # one loss (macaque); mouse is not a loss, the gene need not be older than Primates
        assert _top_scores(e)["CompletenessScore"] == "0.6667"
        assert _top_scores(e)["ImpliedLosses"] == "1"
        assert [c.tag for c in e if c.tag != "score" and c.tag != "property"] == ["paralogGroup"]
    assert _top_scores(root_orthoxml(at_primates, t)[0]) == _top_scores(root_orthoxml(at_mammals, t)[0])


def test_root_group_of_inparalogs_of_one_species_sits_at_its_parent():
    t = species_tree(PRIMATES)
    dup = merged_hog([leaf_hog(H1, t), leaf_hog(H2, t)], t & "Hominini")
    e, = root_orthoxml(merged_hog([merged_hog([dup], t & "Primates")], t), t)
    assert _top_taxrange(e) == "Hominini"
    assert _top_scores(e)["CompletenessScore"] == "0.5"
    assert _top_scores(e)["ImpliedLosses"] == "1"   # chimp


M1, M2 = pid("m1", "MOUSE", 5), pid("m2", "MOUSE", 6)


def _ortholog_trio(h, c, m, t):
    return merged_hog([merged_hog([leaf_hog(h, t), leaf_hog(c, t)], t & "Hominini"), leaf_hog(m, t)], t)


def test_duplication_older_than_the_tree_root_splits_into_separate_groups():
    # every species has both copies: the duplication is older than the root, there is no level to put
    # a single ancestral gene at, so the copies are separate HOGs
    t = species_tree(PRIMATES)
    copies = [_ortholog_trio(H1, CH1, M1, t), _ortholog_trio(H2, CH2, M2, t)]
    groups = root_orthoxml(merged_hog(copies, t), t)
    assert [g.tag for g in groups] == ["orthologGroup", "orthologGroup"]
    assert sorted(g.get("id") for g in groups) == sorted(c.hogid for c in copies)
    assert [_top_taxrange(g) for g in groups] == ["Mammals", "Mammals"]
    assert sorted(sorted(r.get("id") for r in g.iter("geneRef")) for g in groups) == [["1", "3", "5"], ["2", "4", "6"]]


def test_single_gene_copies_of_a_duplication_older_than_the_root_are_not_written():
    t = species_tree(PRIMATES)
    lone_copy = merged_hog([leaf_hog(M2, t)], t)   # boxed up to Mammals, nothing joined it
    trio = _ortholog_trio(H1, CH1, M1, t)
    groups = root_orthoxml(merged_hog([trio, lone_copy], t), t)
    assert [g.get("id") for g in groups] == [trio.hogid]


def test_copy_of_a_split_duplication_that_starts_with_a_duplication_is_wrapped_again():
    t = species_tree(PRIMATES)
    H3, CH3 = pid("h3", "HUMAN", 7), pid("c3", "CHIMP", 8)
    hominini_dup = merged_hog([merged_hog([leaf_hog(H2, t), leaf_hog(CH2, t)], t & "Hominini"),
                               merged_hog([leaf_hog(H3, t), leaf_hog(CH3, t)], t & "Hominini")], t & "Primates")
    copy2 = merged_hog([hominini_dup], t)
    groups = root_orthoxml(merged_hog([_ortholog_trio(H1, CH1, M1, t), copy2], t), t)
    assert sorted(_top_taxrange(g) for g in groups) == ["Mammals", "Primates"]
    wrapper, = [g for g in groups if _top_taxrange(g) == "Primates"]
    assert [c.tag for c in wrapper if c.tag not in ("score", "property")] == ["paralogGroup"]
    assert _top_scores(wrapper)["CompletenessScore"] == "0.6667"


def test_root_group_without_duplication_is_not_wrapped():
    t = species_tree(PRIMATES)
    hog = merged_hog([merged_hog([leaf_hog(H1, t), leaf_hog(CH1, t)], t & "Hominini")], t & "Primates")
    e, = root_orthoxml(hog, t)
    assert ET.tostring(e) == ET.tostring(hog.to_orthoxml(t))
    assert e.tag == "orthologGroup" and _top_taxrange(e) == "Hominini"
