"""Tests for the UnionFind data structure (zoo/unionfind.py, also used by _utils_roothog.py)."""
import pytest

from FastOMA import _utils_roothog
from FastOMA.zoo import unionfind


@pytest.fixture
def UF():
    return unionfind.UnionFind


def test_roothog_utils_use_the_zoo_unionfind():
    assert _utils_roothog.UnionFind is unionfind.UnionFind


def components(uf):
    return {frozenset(c) for c in uf.get_components()}


def test_union_and_find(UF):
    uf = UF()
    uf.union("a", "b")
    uf.union("c", "d", "e")
    uf.union("b", "e")
    uf.union("x")
    assert components(uf) == {frozenset("abcde"), frozenset("x")}
    assert uf["a"] == uf.find("d")


def test_initial_elements_are_singletons(UF):
    uf = UF([1, 2, 3])
    uf.union(1, 2)
    assert components(uf) == {frozenset({1, 2}), frozenset({3})}


def test_find_unknown_creates_singleton(UF):
    uf = UF()
    assert uf.find("q") == "q"
    assert list(uf) == ["q"]


def test_union_without_arguments(UF):
    # characterises current behaviour
    with pytest.raises(ValueError):
        UF().union()


def test_remove_non_root(UF):
    uf = UF()
    uf.union("a", "b", "c")
    member = next(x for x in "abc" if x != uf["a"])
    uf.remove(member)
    assert components(uf) == {frozenset(set("abc") - {member})}


def test_remove_root(UF):
    uf = UF()
    uf.union("a", "b", "c")
    root = uf["a"]
    uf.remove(root)
    assert components(uf) == {frozenset(set("abc") - {root})}


def test_remove_element_others_point_through(UF):
    # a -> b -> c chain (no path compression yet): removing b must not break find(a)
    uf = UF()
    uf.parents.update({"a": "b", "b": "c", "c": "c"})
    uf.weights["c"] = 3
    uf.remove("b")
    assert components(uf) == {frozenset("ac")}
    assert uf.weights[uf["a"]] == 2


def test_remove_singleton(UF):
    uf = UF()
    uf.union("a", "b")
    uf["x"]
    uf.remove("x")
    assert components(uf) == {frozenset("ab")}


def test_union_with_shared_roots_counts_weight_once(UF):
    uf = UF()
    uf.union("a", "x")
    uf.union("b", "c")
    uf.union("a", "b", "c")
    assert components(uf) == {frozenset("abcx")}
    assert uf.weights[uf["a"]] == 4


def test_union_same_element_twice(UF):
    uf = UF("ab")
    uf.union("a", "a", "b")
    assert uf.weights[uf["a"]] == 2


def test_remove_root_of_chained_set(UF):
    # parents {'a': 'a', 'b': 'a', 'c': 'a', 'd': 'c'}: d still points to c, not to the root
    uf = UF()
    uf.union("a", "b")
    uf.union("c", "d")
    uf.union("a", "c")
    root = uf["a"]
    uf.remove(root)
    others = set("abcd") - {root}
    assert components(uf) == {frozenset(others)}
    assert uf.weights[uf[next(iter(others))]] == 3
    assert root not in uf.weights
