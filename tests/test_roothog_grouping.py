"""Tests for the in-memory rootHOG grouping logic of fastoma-infer-roothogs (_utils_roothog.py)."""
import os
import subprocess
import sys
from collections import Counter

import networkx as nx
import pytest

from FastOMA import _utils_roothog as U
from roothog_helpers import hit, partition


@pytest.mark.parametrize("hogid, expected", [
    ("HOG:E0990677.124f", "E0990677"),
    ("HOG:E0990677", "E0990677"),
    ("HOG:C0886513.1b.2a", "C0886513"),
])
def test_extract_root_hog_id(hogid, expected):
    assert U.extract_root_hog_id(hogid) == expected


def test_extract_root_hog_id_requires_prefix():
    # characterises current behaviour: ids without the 'HOG:' prefix are not supported
    with pytest.raises(IndexError):
        U.extract_root_hog_id("E0990677.1a")


def test_group_prots_by_roothogs_uses_first_hit_only():
    hogmaps = {"A": {"a": [hit("HOG:X1.1a", 1), hit("HOG:Y1", 900)]},
               "B": {"b": [hit("HOG:X1", 1)]}}
    assert U.group_prots_by_roothogs(hogmaps) == {"X1": [("A", "a"), ("B", "b")]}


class TestResolveSingletons:
    def test_singleton_moves_to_multi_member_alt_rhog(self, roothog_conf):
        hogmaps = {"S1": {"p": [hit("HOG:A1", 500)], "q": [hit("HOG:A1", 500)],
                          "s": [hit("HOG:B1", 500), hit("HOG:A1.3", 300)]}}
        res = U.resolve_singletons(U.group_prots_by_roothogs(hogmaps), hogmaps, roothog_conf())
        assert res == {"A1": [("S1", "p"), ("S1", "q"), ("S1", "s")]}

    @pytest.mark.parametrize("score, moved", [(70, True), (69.9, False)])
    def test_alt_hit_threshold_is_inclusive(self, roothog_conf, score, moved):
        hogmaps = {"S1": {"p": [hit("HOG:A1", 500)], "q": [hit("HOG:A1", 500)],
                          "s": [hit("HOG:B1", 500), hit("HOG:A1", score)]}}
        res = U.resolve_singletons(U.group_prots_by_roothogs(hogmaps), hogmaps,
                                   roothog_conf(mergHOG_fscore_thresh=70))
        assert (("S1", "s") in res["A1"]) == moved
        assert ("B1" in res) != moved

    def test_alt_hit_to_singleton_rhog_does_not_move(self, roothog_conf):
        # characterises current behaviour: q's alternative hit is p's best rootHOG, but since A1 is a
        # singleton too, p and q are not merged (only singletons sharing an *alternative* rootHOG are)
        hogmaps = {"S1": {"p": [hit("HOG:A1.1a", 500)],
                          "q": [hit("HOG:B1", 500), hit("HOG:A1.2b", 300)]}}
        res = U.resolve_singletons(U.group_prots_by_roothogs(hogmaps), hogmaps, roothog_conf())
        assert partition(res) == {frozenset({("S1", "p")}), frozenset({("S1", "q")})}

    def test_singletons_sharing_alt_rhog_are_merged(self, roothog_conf):
        hogmaps = {"S1": {"p": [hit("HOG:A1", 500), hit("HOG:C1", 300)]},
                   "S2": {"q": [hit("HOG:B1", 500), hit("HOG:C1", 300)]},
                   "S3": {"r": [hit("HOG:D1", 500), hit("HOG:C1", 10)]}}
        res = U.resolve_singletons(U.group_prots_by_roothogs(hogmaps), hogmaps, roothog_conf())
        assert partition(res) == {frozenset({("S1", "p"), ("S2", "q")}), frozenset({("S3", "r")})}
        # the merged group is named after one of the members' rootHOGs (which one depends on the hash seed)
        assert (set(res) - {"D1"}) <= {"A1", "B1"}


class TestFindRhogCandidatePairs:
    @staticmethod
    def hogmaps(n_shared=12, alt_score=300):
        # n_shared proteins with best hit A1 and an alternative hit to B1, plus two B1-only proteins
        hm = {"S1": {f"p{i}": [hit("HOG:A1.1a", 500), hit("HOG:B1", alt_score)] for i in range(n_shared)}}
        hm["S1"].update({f"b{i}": [hit("HOG:B1", 500)] for i in range(2)})
        return hm

    def candidates(self, conf, **kw):
        hm = self.hogmaps(**kw)
        return U.find_rhog_candidate_pairs(hm, U.group_prots_by_roothogs(hm), conf)

    def test_pair_found_bigger_first(self, roothog_conf):
        # sizes count all hits above threshold: A1=12, B1=14
        assert self.candidates(roothog_conf()) == [("B1", "A1")]

    def test_fscore_threshold_is_strict(self, roothog_conf):
        # characterises current behaviour: strict '>' here, while resolve_singletons uses '>='
        assert self.candidates(roothog_conf(mergHOG_fscore_thresh=300)) == []

    def test_shared_count_must_exceed_threshold(self, roothog_conf):
        assert self.candidates(roothog_conf(), n_shared=10) == []

    def test_big_rhogs_are_not_merged(self, roothog_conf):
        # rootHOGs must be smaller than big_rhog_size / 6 (here 10)
        assert self.candidates(roothog_conf(big_rhog_size=60)) == []

    def test_pairs_with_non_best_rhog_are_skipped(self, roothog_conf):
        hm = self.hogmaps()
        rhogs_prots = U.group_prots_by_roothogs(hm)
        del rhogs_prots["B1"]
        assert U.find_rhog_candidate_pairs(hm, rhogs_prots, roothog_conf()) == []


def test_cluster_rhogs_connected_components():
    clusters = U.cluster_rhogs([("a", "b"), ("c", "d"), ("b", "e")])
    assert {frozenset(c) for c in clusters} == {frozenset("abe"), frozenset("cd")}
    assert all(isinstance(c, list) for c in clusters)


class TestHCS:
    def test_two_triangles_joined_by_bridge_are_split(self):
        pairs = [("a", "b"), ("b", "c"), ("a", "c"), ("d", "e"), ("e", "f"), ("d", "f"), ("c", "d")]
        clusters = U.cluster_rhogs_nx([list("abcdef")], pairs)
        assert {frozenset(c) for c in clusters} == {frozenset("abc"), frozenset("def")}

    def test_complete_graph_stays_together(self):
        nodes = list("abcd")
        pairs = [(x, y) for i, x in enumerate(nodes) for y in nodes[i + 1:]]
        assert [set(c) for c in U.cluster_rhogs_nx([nodes], pairs)] == [set(nodes)]

    def test_path_of_three_keeps_one_pair(self):
        # characterises current behaviour: HCS cuts a path, the leftover node is not merged at all
        clusters = U.cluster_rhogs_nx([list("abc")], [("a", "b"), ("b", "c")])
        assert len(clusters) == 1 and len(clusters[0]) == 2 and "b" in clusters[0]

    def test_highly_connected(self):
        g = nx.complete_graph(4)
        assert U.highly_connected(g, nx.minimum_edge_cut(g))
        assert not U.highly_connected(nx.path_graph(3), {(0, 1)})

    def test_clusters_of_1000_rhogs_are_not_merged(self):
        nodes = [f"R{i}" for i in range(1000)]
        clusters = U.cluster_rhogs_nx([nodes], [("R0", "R1")])
        # characterises current behaviour: huge clusters are returned as singletons instead of dropped
        assert len(clusters) == 1000 and all(len(c) == 1 for c in clusters)

    @pytest.mark.slow
    @pytest.mark.xfail(strict=True, reason="BUG: merge partition depends on PYTHONHASHSEED via set order "
                                           "fed into HCS (_utils_roothog.py:576/1028)")
    def test_merge_partition_independent_of_hash_seed(self):
        code = ("from FastOMA import _utils_roothog as U\n"
                "p=[('R1','R2'),('R2','R3'),('R3','R4'),('R4','R5')]\n"
                "print(sorted(sorted(c) for c in U.cluster_rhogs_nx(U.cluster_rhogs(p), p)))")
        results = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                                  env=dict(os.environ, PYTHONHASHSEED=str(seed))).stdout
                   for seed in range(8)}
        assert len(results) == 1


class TestMergeRhogs2:
    @staticmethod
    def patch_clusters(monkeypatch, clusters):
        monkeypatch.setattr(U, "find_rhog_candidate_pairs", lambda *a: [])
        monkeypatch.setattr(U, "cluster_rhogs_nx", lambda cl, cp: [list(c) for c in clusters])

    def test_cluster_merged_into_smallest_id(self, in_tmp_dir, monkeypatch, roothog_conf):
        self.patch_clusters(monkeypatch, [["B", "A"]])
        rp = {"A": [("s", "a1")], "B": [("s", "b1"), ("s", "b2")], "C": [("s", "c1")]}
        res = U.merge_rhogs2({}, rp, roothog_conf())
        assert set(res) == {"A", "C"}
        assert sorted(res["A"]) == [("s", "a1"), ("s", "b1"), ("s", "b2")]
        assert (in_tmp_dir / "merging_rhogs.txt").read_text().splitlines()[1] == "A\tB\t"

    def test_merge_with_real_candidates(self, in_tmp_dir, roothog_conf):
        hm = TestFindRhogCandidatePairs.hogmaps()
        res = U.merge_rhogs2(hm, U.group_prots_by_roothogs(hm), roothog_conf())
        assert list(res) == ["A1"] and len(res["A1"]) == 14

    @pytest.mark.xfail(strict=True, reason="BUG: `newhost` is not reset per cluster and the overflow host "
                                           "loses its own proteins (_utils_roothog.py:667-713)")
    @pytest.mark.parametrize("clusters", [[["A", "B", "C"]], [["A", "B", "C", "D"], ["W", "X", "Y", "Z"]]])
    def test_overflowing_clusters_keep_all_proteins(self, in_tmp_dir, monkeypatch, roothog_conf, clusters):
        self.patch_clusters(monkeypatch, clusters)
        rp = {k: [("s", k + str(i)) for i in range(2)] for c in clusters for k in c}
        before = Counter(p for v in rp.values() for p in v)
        res = U.merge_rhogs2({}, rp, roothog_conf(big_rhog_size=3))
        assert Counter(p for v in res.values() for p in v) == before
        # and no rootHOG mixes proteins of different clusters
        cluster_of = {k: i for i, c in enumerate(clusters) for k in c}
        for members in res.values():
            assert len({cluster_of[p[1][0]] for p in members}) == 1


class TestFilterBigRoothogs:
    def test_low_score_proteins_dropped_from_big_rhog(self, roothog_conf):
        hm = {"A": {"a": [hit("HOG:X1.1a", 96)], "b": [hit("HOG:X1", 95)],
                    "c": [hit("HOG:Y1", 50), hit("HOG:X1", 99)], "z": [hit("HOG:Z1", 1)]}}
        rp = {"X1": [("A", "a"), ("A", "b"), ("A", "c")], "Z1": [("A", "z")]}
        res = U.filter_big_roothogs(hm, rp, roothog_conf(big_rhog_size=2, big_fscore_thresh=95))
        # strict '>' on the score; the score of the *first* hit is used even if it belongs to another
        # rootHOG (c was moved here from Y1) -- characterises current behaviour
        assert res == {"X1": [("A", "a")], "Z1": [("A", "z")]}

    def test_big_rhog_removed_if_nothing_passes(self, roothog_conf):
        hm = {"A": {p: [hit("HOG:X1", 10)] for p in "abc"}}
        rp = {"X1": [("A", p) for p in "abc"]}
        assert U.filter_big_roothogs(hm, rp, roothog_conf(big_rhog_size=2)) == {}

    def test_small_rhogs_untouched(self, roothog_conf):
        hm = {"A": {p: [hit("HOG:X1", 10)] for p in "ab"}}
        rp = {"X1": [("A", "a"), ("A", "b")]}
        assert U.filter_big_roothogs(hm, rp, roothog_conf()) == rp

    def test_filter_nonchild_rootHOG(self, roothog_conf, monkeypatch):
        monkeypatch.setattr(U, "filter_nonchild_rootHOG", True)
        hm = {"A": {"a": [hit("HOG:X1.1a", 99)], "b": [hit("HOG:X1", 99)], "c": [hit("HOG:X1.2b", 99)]}}
        rp = {"X1": [("A", p) for p in "abc"]}
        res = U.filter_big_roothogs(hm, rp, roothog_conf(big_rhog_size=2))
        assert res == {"X1": [("A", "a"), ("A", "c")]}


def test_create_root_hogs(in_tmp_dir, roothog_conf):
    hm = {"S1": {"a": [hit("HOG:A1.1a", 500)],
                 "p": [hit("HOG:P1", 500), hit("HOG:Q1", 300)]},
          "S2": {"a": [hit("HOG:A1", 500)],
                 "q": [hit("HOG:R1", 500), hit("HOG:Q1.4c", 300)]},
          "S3": {"s": [hit("HOG:Z1", 500), hit("HOG:A1.2b", 400)],
                 "u": [hit("HOG:U1", 500)]}}
    res = U.create_root_hogs(hm, roothog_conf())
    assert partition(res) == {
        frozenset({("S1", "a"), ("S2", "a"), ("S3", "s")}),   # s resolved via its alternative hit
        frozenset({("S1", "p"), ("S2", "q")}),                # singletons sharing alternative Q1
        frozenset({("S3", "u")}),                             # remains a singleton
    }
    assert set(res) >= {"A1", "U1"}
