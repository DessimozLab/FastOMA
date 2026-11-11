import pickle

import pytest
from ete3 import Tree

from FastOMA import _infer_subhog
from FastOMA._infer_subhog import LevelHOGProcessor
from subhog_helpers import pid, rec, msa, species_tree, leaf_hog, merged_hog, members_of

A1, A2, A3 = pid("a1", "A", 1), pid("a2", "A", 2), pid("a3", "A", 3)
B1, B2, B3 = pid("b1", "B", 11), pid("b2", "B", 12), pid("b3", "B", 13)
C1 = pid("c1", "C", 21)


def q(*names):
    """newick leaf names need quoting because of the '|' characters."""
    return tuple(f"'{n}'" for n in names)


@pytest.fixture
def sp():
    return species_tree()


@pytest.fixture
def make_lp(sp, subhog_conf):
    def make(ids, level="AB", frozen=(), **conf):
        hogs = [leaf_hog(i, sp) for i in ids]
        for i in frozen:
            h = leaf_hog(i, sp)
            h.active = False
            hogs.append(h)
        return LevelHOGProcessor(sp & level, hogs, "r", subhog_conf(**conf))
    return make


def genetree(nwk):
    return Tree(nwk, format=0, quoted_node_names=True)


def test_constructor_separates_frozen_hogs(make_lp):
    lp = make_lp([A1, B1], frozen=[A3])
    assert {h.get_members().pop() for h in lp.subhogs.values()} == {A1, B1}
    assert {h.get_members().pop() for h in lp._frozen_hogs_before} == {A3}
    assert set(lp._rep_lookup) == {A1, B1}


def test_invalid_msa_filter_method(make_lp):
    with pytest.raises(ValueError):
        make_lp([A1, B1], msa_filter_method="nope")


# infer_reconciliation

def test_reconciliation_labels(make_lp):
    lp = make_lp([A1, A2, B1, B2])
    t = genetree("((%s,%s),(%s,%s));" % q(A1, B1, A2, B2))
    lp.infer_reconciliation(t, sos_threshold=0.1)
    assert (t.evoltype, t.sos, t.so_tuple, t.name) == ("D", 1.0, (2, 2), "D1_2_2")
    assert sorted((c.name, c.evoltype, c.sos) for c in t.children) == [("S1", "S", 0.0), ("S2", "S", 0.0)]
    leaf = t & A1
    assert leaf.hogid == lp._rep_lookup[A1].hog.hogid
    # species sets are only kept on the leaves
    assert not hasattr(t, "species") and not hasattr(t.children[0], "species")


@pytest.mark.parametrize("threshold,expected", [(1 / 3, "S"), (0.3, "D")])
def test_reconciliation_threshold_is_exclusive(make_lp, threshold, expected):
    lp = make_lp([A1, A2, B1, C1], level="root")
    # root: {A,B} vs {A,C} -> sos = 1/3
    t = genetree("((%s,%s),(%s,%s));" % q(A1, B1, A2, C1))
    lp.infer_reconciliation(t, sos_threshold=threshold)
    assert t.sos == pytest.approx(1 / 3)
    assert t.evoltype == expected


@pytest.mark.xfail(strict=True, reason="BUG: species overlap of a multifurcation intersects all children, "
                                       "so a polytomy with two copies of a species is labelled S (_infer_subhog.py:546)")
def test_reconciliation_polytomy_with_duplicated_species(make_lp):
    lp = make_lp([A1, A2, B1])
    t = genetree("(%s,%s,%s);" % q(A1, A2, B1))
    lp.infer_reconciliation(t, sos_threshold=0.1)
    assert t.evoltype == "D"


# merge_subhogs / process

def _reconciled(lp, nwk):
    t = genetree(nwk)
    lp.infer_reconciliation(t, sos_threshold=lp.conf.threshold_dubious_sd)
    return t


def test_merge_subhogs_two_orthogroups(make_lp, sp):
    lp = make_lp([A1, A2, B1, B2])
    t = _reconciled(lp, "((%s:1,%s:1):1,(%s:1,%s:1):1);" % q(A1, B1, A2, B2))
    new = lp.merge_subhogs(t, msa(rec(A1), rec(A2), rec(B1), rec(B2)))
    assert sorted(sorted(h.get_members()) for h in new) == [[A1, B1], [A2, B2]]
    assert all(h.taxname == "AB" for h in new)
    for h in new:
        assert {r.get_id() for r in h.get_representatives()} <= h.get_members()


def test_merge_subhogs_conserves_members(make_lp):
    ids = [A1, A2, A3, B1, B2]
    lp = make_lp(ids)
    t = _reconciled(lp, "(((%s:1,%s:1):1,%s:2):1,(%s:1,%s:1):1);" % q(A1, B1, A3, A2, B2))
    new = lp.merge_subhogs(t, msa(*(rec(i) for i in ids)))
    assert members_of(new) == set(ids)
    assert sum(len(h) for h in new) == len(ids)


def _patch_pipeline(lp, monkeypatch, aln, tree_nwk=None):
    monkeypatch.setattr(lp, "align_subhogs", lambda: aln)
    if tree_nwk is not None:
        monkeypatch.setattr(lp, "infer_genetree_from_msa", lambda m: genetree(tree_nwk))
    else:
        def fail(m):
            raise AssertionError("gene tree should not be inferred")
        monkeypatch.setattr(lp, "infer_genetree_from_msa", fail)


def test_process_keeps_frozen_hogs(make_lp, monkeypatch):
    lp = make_lp([A1, A2, B1, B2], frozen=[A3])
    aln = msa(rec(A1), rec(A2), rec(B1), rec(B2))
    _patch_pipeline(lp, monkeypatch, aln, "((%s:1,%s:1):1,(%s:1,%s:1):1);" % q(A1, B1, A2, B2))
    new = lp.process()
    assert members_of(new) == {A1, A2, A3, B1, B2}
    assert sorted(sorted(h.get_members()) for h in new) == [[A1, B1], [A2, B2], [A3]]


def _all_gappy_msa():
    # every column has a gap ratio of 0.5 -> column filter removes all columns
    return msa((A1, "AAAA----"), (B1, "----AAAA"), (A2, "AAAA----"), (B2, "----AAAA"))


def test_process_when_filter_empties_msa_returns_unmerged_subhogs(make_lp, monkeypatch):
    lp = make_lp([A1, A2, B1, B2], min_col_trim=0)
    active = list(lp.subhogs.values())
    _patch_pipeline(lp, monkeypatch, _all_gappy_msa())
    # characterises current behaviour: the whole level becomes a no-op (seen on real data, e.g. D0040139)
    assert lp.process() == active


@pytest.mark.xfail(strict=True, reason="BUG: process() returns early without the frozen HOGs when filtering empties "
                                       "the MSA (_infer_subhog.py:690); their proteins are lost")
def test_process_when_filter_empties_msa_keeps_frozen_hogs(make_lp, monkeypatch):
    lp = make_lp([A1, A2, B1, B2], frozen=[A3], min_col_trim=0)
    _patch_pipeline(lp, monkeypatch, _all_gappy_msa())
    assert members_of(lp.process()) == {A1, A2, A3, B1, B2}


def test_filter_msa_freezes_hogs_without_representatives(make_lp):
    lp = make_lp([A1, A2, B1, B2], min_col_trim=0)
    b2_hog = lp._rep_lookup[B2].hog
    aln = msa((A1, "ACDEFGHIKL"), (A2, "ACDEFGHIKL"), (B1, "ACDEFGHIKL"), (B2, "A---------"))
    filtered = lp.filter_msa(aln)
    assert [r.id for r in filtered] == [A1, A2, B1]
    assert b2_hog.hogid not in lp.subhogs and not b2_hog.active
    assert b2_hog in lp._frozen_hogs_before
    assert B2 not in lp._rep_lookup


# find_most_divergent_representatives_from_genetree

DIVERGENT_TREE = "(((%s:1,%s:1.2):0.5,%s:3):0.2,((%s:0.4,%s:0.8):2,%s:1):0.3);"


@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_representatives_invariants(make_lp, n):
    ids = [A1, B1, C1, A2, B2, B3]
    lp = make_lp(ids, level="root", number_of_samples_per_hog=n)
    t = genetree(DIVERGENT_TREE % q(*ids))
    reps = lp.find_most_divergent_representatives_from_genetree(t)
    rep_ids = [r.get_id() for r in reps]
    assert 1 <= len(reps) <= n
    assert len(set(rep_ids)) == len(rep_ids)
    assert set(rep_ids) <= set(ids)
    assert set().union(*(r.get_subelements() for r in reps)) == set(ids)


def test_representatives_small_tree_keeps_all_leaves(make_lp):
    ids = [A1, B1, C1]
    lp = make_lp(ids, level="root", number_of_samples_per_hog=5)
    t = genetree("((%s:1,%s:1):1,%s:1);" % q(*ids))
    reps = lp.find_most_divergent_representatives_from_genetree(t)
    assert sorted(r.get_id() for r in reps) == sorted(ids)


def test_representatives_do_not_modify_input_tree(make_lp):
    ids = [A1, B1, C1, A2, B2, B3]
    lp = make_lp(ids, level="root", number_of_samples_per_hog=2)
    t = genetree(DIVERGENT_TREE % q(*ids))
    before = t.write()
    lp.find_most_divergent_representatives_from_genetree(t)
    assert t.write() == before


# read_infer_xml_rhog: wrapping of root-level paralogGroups

def _run_rhog(in_tmp_dir, subhog_conf, monkeypatch, make_hogs):
    """read_infer_xml_rhog on rootHOG "r" of A1..C1, with the subHOG inference replaced by
    ``make_hogs(sptree, a1, a2, a3, b1, c1)`` giving the HOGs at the root level."""
    ids = [A1, A2, A3, B1, C1]
    with open("HOG_r.fa", "w") as fh:
        for i in ids:
            fh.write(f">{i}\nACDEFGHIKL\n")
    with open("species_tree.nwk", "w") as fh:
        fh.write("((A,B)AB,C)root;\n")

    def fake_inference(sptree, rhogid, pickles_subhog_folder_all, rhogs_fa_folder, conf):
        hogs = make_hogs(sptree, *(leaf_hog(i, sptree) for i in ids))
        with open(f"{pickles_subhog_folder_all}/rhog_{rhogid}/{sptree.name}.pickle", "wb") as fh:
            pickle.dump(hogs, fh)

    monkeypatch.setattr(_infer_subhog, "infer_hogs_for_rhog_levels_recursively", fake_inference)
    (in_tmp_dir / "out").mkdir()
    n = _infer_subhog.read_infer_xml_rhog("r", False, "out", str(in_tmp_dir), str(in_tmp_dir),
                                          subhog_conf(species_tree="species_tree.nwk", v=0))
    with open("out/file_r.pickle", "rb") as fh:
        return n, pickle.load(fh)


@pytest.fixture
def rhog_run(in_tmp_dir, subhog_conf, monkeypatch):
    # the in-paralogs a1, a2 end up at the root level (as when nothing merges with them)
    return _run_rhog(in_tmp_dir, subhog_conf, monkeypatch,
                     lambda t, a1, a2, a3, b1, c1: [merged_hog([a1, a2], t), merged_hog([b1, c1], t), a3])


def test_read_infer_xml_rhog_splits_a_duplication_older_than_the_root(in_tmp_dir, subhog_conf, monkeypatch):
    # the copies span the whole tree (A, B and C): there is no level above the duplication, so each
    # copy is a top-level group of its own; the single-gene copies a2 and a3 are not written
    def make_hogs(t, a1, a2, a3, b1, c1):
        copies = [merged_hog([merged_hog([a1, b1], t & "AB"), c1], t),
                  merged_hog([a2], t), merged_hog([a3], t)]
        return [merged_hog(copies, t)]

    n, groups = _run_rhog(in_tmp_dir, subhog_conf, monkeypatch, make_hogs)
    assert n == len(groups) == 1
    assert groups[0].tag == "orthologGroup"
    assert sorted(r.get("id") for r in groups[0].iter("geneRef")) == ["1", "11", "21"]


def test_read_infer_xml_rhog_skips_singletons(rhog_run, in_tmp_dir):
    n, groups = rhog_run
    assert n == len(groups) == 2
    ids = sorted(g.get("id") for e in groups for g in e.iter("geneRef"))
    assert ids == ["1", "11", "2", "21"]
    # the per-level pickles are removed again
    assert not (in_tmp_dir / "rhog_r").exists()


def test_read_infer_xml_rhog_wraps_root_paralog_group(rhog_run):
    _, groups = rhog_run
    wrapped = groups[0]
    assert wrapped.tag == "orthologGroup"
    assert [c.tag for c in wrapped if c.tag not in ("score", "property")] == ["paralogGroup"]


def test_read_infer_xml_rhog_wrapper_taxrange_is_parent_of_mrca(rhog_run):
    # the in-paralogs of A were duplicated in A's lineage: the gene is older than A, but nothing
    # says it is older than AB (issue #128)
    _, groups = rhog_run
    assert [p.get("value") for p in groups[0].findall("property") if p.get("name") == "TaxRange"] == ["AB"]


@pytest.mark.xfail(strict=True, reason="BUG: a root-level paralogGroup is wrapped in an orthologGroup with a single "
                                       "child (_infer_subhog.py:112-122); orthoXML groups need >= 2 children")
def test_read_infer_xml_rhog_groups_have_two_children(rhog_run):
    _, groups = rhog_run
    for e in groups:
        assert len([c for c in e if c.tag not in ("score", "property")]) >= 2
