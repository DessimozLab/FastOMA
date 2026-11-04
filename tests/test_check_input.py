"""Tests for the input validation of fastoma-check-input (check_input.py)."""
import logging

import pytest
from ete3 import Tree

from FastOMA import check_input as C
from conftest import TESTDATA_DIR
from roothog_helpers import rec

PROT = "MKVLAAGIVLLLAAGCSSHE"


def make_files(folder, names):
    folder.mkdir(exist_ok=True)
    for n in names:
        (folder / n).write_text(">p\nMKV\n")
    return str(folder)


class TestCheckProteomeFiles:
    def test_species_names(self, tmp_path):
        res = C.check_proteome_files(make_files(tmp_path / "p", ["A.fa", "B.fa", "Homo.sapiens.fa"]))
        assert sorted(res) == ["A", "B", "Homo.sapiens"]

    def test_needs_three_proteomes(self, tmp_path):
        assert C.check_proteome_files(make_files(tmp_path / "p", ["A.fa", "B.fa"])) is False

    def test_mixed_extensions_rejected(self, tmp_path):
        assert C.check_proteome_files(make_files(tmp_path / "p", ["A.fa", "B.fa", "C.fasta"])) is False

    @pytest.mark.xfail(strict=True, reason="BUG: hidden files are taken as proteomes, e.g. the ._A.fa AppleDouble "
                                           "files of archives made with macOS tar (check_input.py:39); check-input "
                                           "then reports '._A' as missing from the species tree")
    def test_hidden_files_ignored(self, tmp_path):
        folder = make_files(tmp_path / "p", ["A.fa", "B.fa", "C.fa", "._A.fa", ".DS_Store"])
        assert sorted(C.check_proteome_files(folder)) == ["A", "B", "C"]


class TestCheckProteome:
    @staticmethod
    def check(recs):
        return C.check_proteome(["A"], {"A": recs}, "folder")

    def test_ok(self):
        assert self.check([rec(f"p{i}", PROT) for i in range(3)]) is True

    @pytest.mark.parametrize("recs", [
        pytest.param([rec("p", PROT)] * 3, id="duplicate-ids"),
        pytest.param([rec("x" * 61, PROT), rec("b", PROT), rec("c", PROT)], id="id-longer-than-60"),
        pytest.param([rec(f"p{i}", "ACGTACGTAC") for i in range(3)], id="dna"),
        pytest.param([rec("a", PROT), rec("b", PROT)], id="two-proteins"),
    ])
    def test_rejected(self, recs):
        assert self.check(recs) is False

    @pytest.mark.xfail(strict=True, reason="BUG: empty proteome raises ZeroDivisionError (check_input.py:65) "
                                           "instead of being reported as invalid")
    def test_empty_proteome_rejected(self):
        assert self.check([]) is False


class TestCheckHogmapFiles:
    def test_species_names(self, tmp_path):
        folder = make_files(tmp_path / "h", ["A.fa.hogmap", "B.fasta.hogmap", "notes.txt"])
        assert sorted(C.check_hogmap_files(folder)) == ["A", "B"]

    @pytest.mark.xfail(strict=True, reason="BUG: hidden ._A.fa.hogmap files (macOS tar) count as hogmaps "
                                           "(check_input.py:91), so the hogmap set never equals the proteomes")
    def test_hidden_files_ignored(self, tmp_path):
        assert C.check_hogmap_files(make_files(tmp_path / "h", ["A.fa.hogmap", "._A.fa.hogmap"])) == ["A"]

    @pytest.mark.xfail(strict=True, reason="BUG: dots are dropped from species names (check_input.py:91), "
                                           "so they never match the proteome names")
    def test_dotted_species_name(self, tmp_path):
        assert C.check_hogmap_files(make_files(tmp_path / "h", ["Homo.sapiens.fa.hogmap"])) == ["Homo.sapiens"]


class TestSpeciesTreeChecks:
    def test_leaves_ok(self):
        assert C.check_speciestree_leaves(Tree("((A,B),C);"), ["A", "B"]) is True

    def test_leaves_missing_species(self):
        assert C.check_speciestree_leaves(Tree("((A,B),C);"), ["A", "D"]) is False

    def test_leaves_duplicated(self):
        assert C.check_speciestree_leaves(Tree("((A,B),A);"), ["A", "B"]) is False

    @pytest.mark.xfail(strict=True, reason="BUG: name list is reset inside the loop (check_input.py:100), "
                                           "duplicate internal names are never reported")
    def test_internalnode_duplicates_reported(self, caplog):
        with caplog.at_level(logging.WARNING, logger="FastOMA"):
            C.check_speciestree_internalnode(Tree("((A,B)x,(C,D)x)r;", format=1))
        assert any("unique" in r.message for r in caplog.records)

    @pytest.mark.xfail(strict=True, reason="BUG: UnboundLocalError for trees without internal nodes "
                                           "(check_input.py:106)")
    def test_internalnode_single_leaf(self):
        C.check_speciestree_internalnode(Tree("A;"))


class TestAddInternalNodePrune:
    @staticmethod
    def prune(tmp_path, nwk, species):
        out = tmp_path / "tree.nwk"
        C.add_internal_node_prune(Tree(nwk, format=1), species, str(out))
        return out.read_text().strip()

    @pytest.mark.parametrize("nwk, species, expected", [
        ("((A,B),(C,D));", "ABCD", "((A:1,B:1)internal_0:1,(C:1,D:1)internal_1:1)internal_2:0;"),
        ("((A,B)inner,(C,D)inner)root;", "ABCD", "((A:1,B:1)inner:1,(C:1,D:1)internal_0:1)root:0;"),
        ("((A,B)abc,(C,D)'a.bc')root;", "ABCD", "((A:1,B:1)abc:1,(C:1,D:1)abc_0:1)root:0;"),
        # names shorter than 3 characters are replaced
        ("((A,B)ab,(C,D)cd)root;", "ABC", "((A:1,B:1)internal_0:1,C:1)root:0;"),
    ])
    def test_rename_and_prune(self, tmp_path, nwk, species, expected):
        assert self.prune(tmp_path, nwk, list(species)) == expected

    @pytest.mark.xfail(strict=True, reason="BUG: a name made only of special chars becomes empty "
                                           "(check_input.py:168)")
    def test_special_char_only_name(self, tmp_path):
        tree = Tree(self.prune(tmp_path, "((A,B)'!!!',(C,D)cde)root;", list("ABCD")), format=1)
        assert all(n.name for n in tree.traverse())

    @pytest.mark.xfail(strict=True, reason="BUG: internal node named like a leaf makes ete3 prune raise "
                                           "'Ambiguous node name' (check_input.py:180)")
    def test_internal_node_named_like_leaf(self, tmp_path):
        res = self.prune(tmp_path, "((AQUAE,CHLTR)MYCGE,MYCGE)root;", ["AQUAE", "CHLTR", "MYCGE"])
        assert {n.name for n in Tree(res, format=1).get_leaves()} == {"AQUAE", "CHLTR", "MYCGE"}


def test_check_omamer_db(tmp_path):
    small = tmp_path / "small.h5"
    small.write_bytes(b"x" * 10)
    assert C.check_omamer_db(None) is True
    assert C.check_omamer_db(str(tmp_path / "missing.h5")) is False
    assert C.check_omamer_db(str(small)) is False
    assert C.check_omamer_db(str(TESTDATA_DIR / "test.h5")) is True


class TestCheckSplice:
    def test_ok(self):
        assert C.check_splice({"A": [["a1", "a2"], ["b1"]], "B": []}) is None

    def test_no_isoforms_exits(self):
        with pytest.raises(SystemExit):
            C.check_splice({"A": [["a1"], ["b1"]]})

    @pytest.mark.xfail(strict=True, reason="BUG: species with a single splice line are ignored "
                                           "(len(isoforms)>1 at check_input.py:195)")
    def test_single_gene_with_isoforms(self):
        assert C.check_splice({"A": [["a1", "a2"]]}) is None
