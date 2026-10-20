import random

import pytest
from ete3 import Tree
from ete3.parser.newick import NewickError

from FastOMA import _wrappers
from FastOMA._infer_subhog import LevelHOGProcessor
from FastOMA._utils_subhog import find_outlier_leaves, midpoint_rooting_outgroup
from subhog_helpers import species_tree, leaf_hog

P1 = "sp|P1|G3P_MYCGE||MYCGE||1000000005"
Q2 = "tr|Q2|X_Y_AQUAE||AQUAE||1001"
P3 = "sp|P3|Z_CHLTR||CHLTR||1002"
P4 = "sp|P4|W_CHLTR||CHLTR||1003"
IDS = [P1, Q2, P3, P4]
SPECIES_NWK = "(MYCGE,(AQUAE,CHLTR)AC)root;"
# unrooted (trifurcating) FastTree-like tree; longest path P4 - Q2 = 3.55
UNROOTED = "('%s':0.1,'%s':0.5,('%s':0.2,'%s':3.0):0.05);" % tuple(IDS)


@pytest.fixture
def make_lp(subhog_conf):
    def make(method):
        sp = species_tree(SPECIES_NWK)
        hogs = [leaf_hog(i, sp) for i in IDS]
        return LevelHOGProcessor(sp, hogs, "r", subhog_conf(gene_rooting_method=method))
    return make


def gtree(nwk=UNROOTED):
    return Tree(nwk, format=0, quoted_node_names=True)


def root_split(t):
    return sorted(sorted(c.get_leaf_names()) for c in t.children)


@pytest.mark.parametrize("method", ["midpoint", "midpoint-dendropy"])
def test_midpoint_rooting(make_lp, method):
    t = make_lp(method).infer_rooted_genetree(gtree())
    assert sorted(t.get_leaf_names()) == sorted(IDS)
    assert len(t.children) == 2
    assert root_split(t) == [[P1, P3, Q2], [P4]]


def test_midpoint_ete3_splits_root_edge_in_half(make_lp):
    t = make_lp("midpoint").infer_rooted_genetree(gtree())
    # characterises current behaviour: ete3 places the root in the middle of the edge, not at the
    # true midpoint (1.775 from P4); this affects the representative selection downstream
    assert sorted(c.dist for c in t.children) == pytest.approx([1.5, 1.5])


def test_midpoint_dendropy_uses_true_midpoint(make_lp):
    t = make_lp("midpoint-dendropy").infer_rooted_genetree(gtree())
    assert sorted(c.dist for c in t.children) == pytest.approx([1.225, 1.775])


@pytest.mark.parametrize("method", ["midpoint", "midpoint-dendropy"])
def test_midpoint_two_leaves(make_lp, method):
    t = make_lp(method).infer_rooted_genetree(gtree("('%s':1,'%s':2);" % (P1, Q2)))
    assert sorted(t.get_leaf_names()) == sorted([P1, Q2])
    assert sorted(c.dist for c in t.children) == pytest.approx([1.5, 1.5])


def test_midpoint_zero_branch_lengths(make_lp):
    t = make_lp("midpoint").infer_rooted_genetree(gtree("('%s':0,'%s':0,'%s':0);" % (P1, Q2, P3)))
    assert sorted(t.get_leaf_names()) == sorted([P1, Q2, P3])
    assert len(t.children) == 2


def test_outlier_rooting_ignores_long_branch(make_lp):
    random.seed(1)
    t = make_lp("outlier").infer_rooted_genetree(gtree())
    assert sorted(t.get_leaf_names()) == sorted(IDS)
    # P4 is excluded from the midpoint computation, so the root moves away from its branch
    assert root_split(t) == [[P1, P3, P4], [Q2]]


def test_nevers_rooting_reads_species_tree_from_cwd(make_lp, in_tmp_dir):
    # characterises current behaviour: the species tree is read from ./species_tree.nwk (_infer_subhog.py:507)
    with pytest.raises(NewickError):
        make_lp("Nevers_rooting").infer_rooted_genetree(gtree())
    (in_tmp_dir / "species_tree.nwk").write_text(SPECIES_NWK + "\n")
    t = make_lp("Nevers_rooting").infer_rooted_genetree(gtree())
    assert sorted(t.get_leaf_names()) == sorted(IDS)
    assert len(t.children) == 2


@pytest.mark.xfail(strict=True, reason="BUG: mad rooting passes the TreeNode, but _wrappers.mad_rooting expects a "
                                       "file path (_infer_subhog.py:512 vs _wrappers.py:138)")
def test_mad_rooting_gets_a_file(make_lp, monkeypatch):
    seen = []
    monkeypatch.setattr(_wrappers, "mad_rooting", lambda arg: seen.append(arg) or gtree())
    make_lp("mad").infer_rooted_genetree(gtree())
    assert isinstance(seen[0], str)


def test_invalid_rooting_method(make_lp):
    with pytest.raises(ValueError):
        make_lp("bogus").infer_rooted_genetree(gtree())


# outlier helpers

def test_find_outlier_leaves_long_branch():
    t = Tree("((a:1,b:1):1,(c:1,(d:1,e:30):1):1);")
    assert [n.name for n in find_outlier_leaves(t)] == ["e"]


def test_find_outlier_leaves_balanced_tree():
    t = Tree("((a:1,b:1):1,(c:1,d:1):1);")
    assert find_outlier_leaves(t) == []


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_midpoint_rooting_outgroup_with_excluded_leaf(seed):
    random.seed(seed)
    t = Tree("((a:1,b:1):1,(c:1,(d:1,e:30):1):1);")
    og = midpoint_rooting_outgroup(t, leaves_to_exclude=[t & "e"])
    t.set_outgroup(og)
    assert root_split(t) == [["a", "b"], ["c", "d", "e"]]


def test_midpoint_rooting_outgroup_two_leaves():
    random.seed(0)
    t = Tree("(a:1,b:1);")
    og = midpoint_rooting_outgroup(t)
    assert og.name in ("a", "b")
