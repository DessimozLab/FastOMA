"""Tests for the UnionFind data structure (zoo/unionfind.py and its copy in _utils_roothog.py)."""
import pytest

from FastOMA import _utils_roothog
from FastOMA.zoo import unionfind


@pytest.fixture(params=[unionfind.UnionFind, _utils_roothog.UnionFind], ids=["zoo", "utils_roothog"])
def UF(request):
    return request.param


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


@pytest.mark.xfail(strict=True, reason="BUG: removing the representative of a set breaks find() for the "
                                       "other members (KeyError, zoo/unionfind.py:78)")
def test_remove_root(UF):
    uf = UF()
    uf.union("a", "b", "c")
    root = uf["a"]
    uf.remove(root)
    assert components(uf) == {frozenset(set("abc") - {root})}
