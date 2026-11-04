import pytest
from ete3 import Tree

from FastOMA._utils_subhog import prepare_species_tree, read_species_tree, list_rhog_fastas
from subhog_helpers import pid, rec

TREE = "((A,B)AB,(C,D)CD)root;"


def recs(*species):
    return [rec(pid(f"p{i}", s, 100 + i)) for i, s in enumerate(species)]


def nwk(t):
    return t.write(format=9)


def test_prepare_species_tree_returns_mrca_subtree():
    sp = Tree(TREE, format=1)
    pruned, species, prots, full = prepare_species_tree(recs("A", "B", "A"), sp, "r")
    assert pruned.name == "AB"
    assert sorted(pruned.get_leaf_names()) == ["A", "B"]
    assert species == ["A", "B", "A"]
    assert prots == ["100", "101", "102"]
    # characterises current behaviour: the subtree stays attached to the input tree
    assert pruned.up is not None and pruned.get_tree_root() is sp


def test_prepare_species_tree_prunes_absent_species():
    sp = Tree(TREE, format=1)
    pruned, _, _, full = prepare_species_tree(recs("A", "B", "C"), sp, "r")
    assert pruned.write(format=8, format_root_node=True) == "((A,B)AB,C)root;"
    # the full tree is an unpruned copy with the original clade sizes
    assert sorted(full.get_leaf_names()) == ["A", "B", "C", "D"]
    assert {n.name: n.size for n in full.traverse()} == {"root": 4, "AB": 2, "CD": 2, "A": 1, "B": 1, "C": 1, "D": 1}
    assert pruned.size == 4 and (pruned & "C").size == 1


def test_prepare_species_tree_unknown_species():
    with pytest.raises(ValueError):
        prepare_species_tree(recs("A", "Q"), Tree(TREE, format=1), "r")


def test_prepare_species_tree_empty_roothog():
    with pytest.raises(AssertionError):
        prepare_species_tree([], Tree(TREE, format=1), "r")


@pytest.mark.xfail(strict=True, reason="BUG: for a single-species rootHOG get_common_ancestor({'A'}) returns the root, "
                                       "giving a unary root instead of the species leaf (_utils_subhog.py:283)")
def test_prepare_species_tree_single_species():
    pruned, _, _, _ = prepare_species_tree(recs("A", "A"), Tree(TREE, format=1), "r")
    assert pruned.is_leaf() and pruned.name == "A"


def test_read_species_tree_newick_with_internal_names(tmp_path):
    p = tmp_path / "sp.nwk"
    p.write_text(TREE + "\n")
    t = read_species_tree(str(p))
    assert [n.name for n in t.traverse("preorder")] == ["root", "AB", "A", "B", "CD", "C", "D"]


def test_read_species_tree_missing_file(tmp_path):
    with pytest.raises(SystemExit):
        read_species_tree(str(tmp_path / "missing.nwk"))


def test_read_species_tree_unknown_extension(tmp_path):
    p = tmp_path / "sp.tree"
    p.write_text(TREE + "\n")
    with pytest.raises(SystemExit):
        read_species_tree(str(p))


def test_list_rhog_fastas(tmp_path):
    for name in ["HOG_D0001.fa", "HOG_E0002.fa", "HOG_D0003.fa.gz", "HOG_D0004.fa.bak", "notes.txt"]:
        (tmp_path / name).write_text("")
    assert sorted(list_rhog_fastas(str(tmp_path))) == ["D0001", "E0002"]
