"""Tests for the file parsing/writing helpers of fastoma-infer-roothogs (_utils_roothog.py)."""
import pickle

import pytest
from Bio import SeqIO

from FastOMA import _utils_roothog as U
from roothog_helpers import fasta_ids, hit, rec, write_fasta, write_hogmap

SEQ = "MKVLAAGIVLLLAAGCSSHE"


def test_add_species_name_prot_id_numbering():
    recs = {"AQUAE": [rec("a1"), rec("a2")], "CHLTR": [rec("c1")]}
    res = U.add_species_name_prot_id(recs)
    assert {sp: {k: r.id for k, r in d.items()} for sp, d in res.items()} == {
        "AQUAE": {"a1": "a1||AQUAE||1000000001", "a2": "a2||AQUAE||1000000002"},
        "CHLTR": {"c1": "c1||CHLTR||1001000001"},
    }
    # the input records are modified in place
    assert recs["CHLTR"][0].id == "c1||CHLTR||1001000001"


def test_add_species_name_prot_id_truncates_long_ids():
    res = U.add_species_name_prot_id({"A": [rec("L" * 240)]})
    assert list(res["A"]) == ["L" * 230]


class TestParseProteomes:
    def test_sorted_species_and_min_length(self, tmp_path):
        write_fasta(tmp_path / "B.fa", [("b1", SEQ), ("b_short", "MKV")])
        write_fasta(tmp_path / "A.fa", [("a1", SEQ)])
        species, recs, ext = U.parse_proteomes(str(tmp_path), min_sequence_length=10)
        assert species == ["A", "B"]
        assert ext == "fa"
        assert {sp: [r.id for r in rs] for sp, rs in recs.items()} == {"A": ["a1"], "B": ["b1"]}

    def test_fasta_extension(self, tmp_path):
        write_fasta(tmp_path / "A.fasta", [("a1", SEQ)])
        assert U.parse_proteomes(str(tmp_path))[2] == "fasta"

    @pytest.mark.xfail(strict=True, reason="BUG: files without an extension crash parse_proteomes "
                                           "(file.rsplit at _utils_roothog.py:166); check_input only warns")
    def test_non_fasta_files_ignored(self, tmp_path):
        write_fasta(tmp_path / "A.fa", [("a1", SEQ)])
        (tmp_path / "README").write_text("x")
        assert U.parse_proteomes(str(tmp_path))[0] == ["A"]

    @pytest.mark.xfail(strict=True, reason="BUG: hidden files such as the ._A.fa AppleDouble files of macOS tar "
                                           "archives are parsed as proteomes (_utils_roothog.py:165-170)")
    def test_hidden_files_ignored(self, tmp_path):
        write_fasta(tmp_path / "A.fa", [("a1", SEQ)])
        (tmp_path / "._A.fa").write_bytes(b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X        ")
        assert U.parse_proteomes(str(tmp_path))[0] == ["A"]

    def test_mixed_extensions_fail(self, tmp_path):
        # characterises current behaviour: only the last seen extension is used; check_input rejects this
        write_fasta(tmp_path / "A.fa", [("a1", SEQ)])
        write_fasta(tmp_path / "B.fasta", [("b1", SEQ)])
        with pytest.raises(FileNotFoundError):
            U.parse_proteomes(str(tmp_path))


class TestParseHogmapOmamer:
    @staticmethod
    def parse(tmp_path, proteome_ids, rows):
        write_hogmap(tmp_path / "A.fa.hogmap", rows)
        recs = {"A": [rec(i) for i in proteome_ids]}
        U.add_species_name_prot_id(recs)  # parse_hogmap_omamer expects the renamed records
        return U.parse_hogmap_omamer(recs, "fa", folder=str(tmp_path))

    def test_hits_unmapped_and_unknown(self, tmp_path):
        hogmaps, unmapped = self.parse(tmp_path, ["p1", "p2", "p3"], [
            ("p1", "HOG:A1.1a", 500.5), ("p1", "HOG:B1", 30),
            ("p2", "N/A", None),
            ("other", "HOG:C1", 100),
        ])
        assert dict(hogmaps["A"]) == {"p1": [hit("HOG:A1.1a", 500.5), hit("HOG:B1", 30)]}
        assert dict(unmapped) == {"A": ["p2"]}

    def test_hogmap_values_are_strings(self, tmp_path):
        hogmaps, _ = self.parse(tmp_path, ["p1"], [("p1", "HOG:A1", 12.5)])
        assert hogmaps["A"]["p1"][0].score == "12.5"

    def test_protein_missing_from_hogmap(self, tmp_path):
        # characterises current behaviour: a protein without any hogmap row is neither mapped nor
        # unmapped, so it is silently dropped from the pipeline (see test_roothog_stage.py)
        hogmaps, unmapped = self.parse(tmp_path, ["p1", "p2"], [("p1", "HOG:A1", 100)])
        assert "p2" not in hogmaps["A"] and "p2" not in unmapped["A"]

    @pytest.mark.xfail(strict=True, reason="BUG: ids >230 chars are truncated (_utils_roothog.py:256) and then "
                                           "no longer match the hogmap qseqid (:289), the protein vanishes")
    def test_long_ids_are_kept(self, tmp_path):
        long_id = "L" * 240
        hogmaps, unmapped = self.parse(tmp_path, [long_id], [(long_id, "HOG:A1", 100)])
        assert len(hogmaps["A"]) + len(unmapped["A"]) == 1


def test_write_rhog(tmp_path):
    prot_recs_all = U.add_species_name_prot_id({"A": [rec("a1"), rec("a2")], "B": [rec("b1")]})
    rhogs = {"X1": [("A", "a1"), ("B", "b1")], "Y1": [("A", "a2")],
             "Z1": [("A", "a2"), ("B", "too_short")]}  # too_short is not in prot_recs_all
    out = tmp_path / "rhogs"
    written = U.write_rhog(rhogs, prot_recs_all, str(out), min_rhog_size=2)
    assert written == ["X1"]
    assert sorted(p.name for p in out.iterdir()) == ["HOG_X1.fa"]
    assert [r.id for r in SeqIO.parse(out / "HOG_X1.fa", "fasta")] == ["a1||A||1000000001", "b1||B||1001000001"]


def test_collect_unmapped_singleton(tmp_path):
    prot_recs_all = U.add_species_name_prot_id({"A": [rec("a1"), rec("a2"), rec("a3")], "B": [rec("b1")]})
    rhogs = {"X1": [("A", "a1"), ("B", "b1")], "Y1": [("A", "a2")]}
    unmapped = {"A": ["a3", "short_and_removed"]}
    out = tmp_path / "singletons.fa"
    assert U.collect_unmapped_singleton(rhogs, unmapped, prot_recs_all, str(out)) == 2
    assert fasta_ids(out) == ["a3", "a2"]


class TestWriteClusters:
    MMSEQS_OUT = (">r1\n>a||S1||1 desc\nMKV\n>b||S2||2\nMKA\n"   # 2 species -> written
                  ">r2\n>c||S1||3\nMKK\n"                        # singleton -> skipped
                  ">r3\n>d||S1||4\nMAA\n>e||S1||5\nMCC\n"        # single species -> discarded
                  ">r4\n>f||S1||6\nMAA\n>g||S2||7\nMCC\n")       # last cluster -> written

    def test_write_clusters(self, tmp_path):
        clustered = tmp_path / "all_seqs.fasta"
        clustered.write_text(self.MMSEQS_OUT)
        out = tmp_path / "rhogs"
        out.mkdir()
        assert U.write_clusters(str(clustered), str(out), 2) == 2
        assert {p.name: fasta_ids(p) for p in out.iterdir()} == {
            "HOG_clust10000.fa": ["a", "b"], "HOG_clust10001.fa": ["f", "g"]}

    def test_empty_file(self, tmp_path):
        clustered = tmp_path / "all_seqs.fasta"
        clustered.write_text("")
        assert U.write_clusters(str(clustered), str(tmp_path), 2) == 0


@pytest.mark.tools("mmseqs")
def test_run_linclust_and_write_clusters(in_tmp_dir, proteome_recs):
    prot_recs_all = U.add_species_name_prot_id(proteome_recs)
    unmapped = {sp: list(d) for sp, d in prot_recs_all.items()}
    assert U.collect_unmapped_singleton({}, unmapped, prot_recs_all) == 28
    clustered = U.run_linclust()
    (in_tmp_dir / "rhogs").mkdir()
    n = U.write_clusters(clustered, "rhogs", 2)
    files = list((in_tmp_dir / "rhogs").glob("HOG_clust*.fa"))
    assert n == len(files) > 0
    for f in files:  # every written cluster spans at least two species
        with open(f) as fh:
            species = {line.split()[0].split("||")[1] for line in fh if line.startswith(">")}
        assert len(species) >= 2


@pytest.fixture
def proteome_recs():
    from conftest import TESTDATA_DIR
    return U.parse_proteomes(str(TESTDATA_DIR / "in_folder" / "proteome"), 50)[1]


class TestIsoforms:
    def test_parse_isoform_files(self, tmp_path):
        (tmp_path / "A.splice").write_text("i1;i2\ni3\n\n")
        assert U.parse_isoform_files(["A", "B"], str(tmp_path)) == {
            "A": [["i1", "i2"], ["i3"], [""]],  # blank lines become [''] -- characterises current behaviour
            "B": []}

    def test_select_best_isoform_uses_max_score_over_all_hits(self):
        hm = {"A": {"i1": [hit("HOG:X1", 100), hit("HOG:Y1", 900)], "i2": [hit("HOG:Z1", 500)]}}
        assert U.select_best_isoforms(["A"], {"A": [["i1", "i2"]]}, hm) == ({"A": ["i1"]}, {"A": ["i2"]})

    def test_select_best_isoform_tie_takes_last(self):
        hm = {"A": {"i1": [hit("HOG:X1", 5)], "i2": [hit("HOG:X1", 5)]}}
        assert U.select_best_isoforms(["A"], {"A": [["i1", "i2"]]}, hm) == ({"A": ["i2"]}, {"A": ["i1"]})

    def test_select_best_isoform_none_mapped_takes_first(self):
        assert U.select_best_isoforms(["A"], {"A": [["i1", "i2"]]}, {"A": {}}) == ({"A": ["i1"]}, {"A": ["i2"]})

    def test_filter_hogmaps_by_isoforms(self):
        hm = {"A": {"i1": 1, "i2": 2, "x": 3}, "B": {"b": 4}}
        assert U.filter_hogmaps_by_isoforms(hm, {"A": ["i2"], "B": []}) == {"A": {"i1": 1, "x": 3}, "B": {"b": 4}}

    def test_handle_splice_variants(self, tmp_path):
        (tmp_path / "A.splice").write_text("i1;i2\n")
        hm = {"A": {"i1": [hit("HOG:X1", 10)], "i2": [hit("HOG:X1", 20)], "x": [hit("HOG:X1", 1)]}}
        res = U.handle_splice_variants(["A"], hm, str(tmp_path))
        assert res["selected_isoforms"] == {"A": ["i2"]}
        assert res["non_selected"] == {"A": ["i1"]}
        assert set(res["filtered_hogmaps"]["A"]) == {"i2", "x"}

    def test_write_selected_isoforms(self, in_tmp_dir):
        recs = {"A": [rec("i1"), rec("i2"), rec("x")]}
        U.add_species_name_prot_id(recs)
        data = {"selected_isoforms": {"A": ["i2"]}, "isoform_by_gene": {"A": [["i1", "i2"]]}}
        U.write_selected_isoforms(data, recs)
        assert (in_tmp_dir / "selected_isoforms" / "A_selected_isoforms.fa").read_text() == \
            f">i2 minor_isoforms=i1\n{SEQ}\n>x \n{SEQ}\n"


class TestSaveGeneIdMapping:
    def test_with_isoforms(self, tmp_path):
        prot_recs_all = U.add_species_name_prot_id({"A": [rec("g1a"), rec("g1b"), rec("g2")]})
        iso = {"selected_isoforms": {"A": ["g1b"]}, "isoform_by_gene": {"A": [["g1a", "g1b"]]}}
        out = tmp_path / "map.pickle"
        res = U.save_gene_id_mapping(prot_recs_all, iso, out_path=str(out))
        assert res == {"A": [U.Gene(1000000001, "g1a", 1000000002), U.Gene(1000000002, "g1b", 1000000002),
                             U.Gene(1000000003, "g2", None)]}
        with open(out, "rb") as fh:
            assert pickle.load(fh) == res

    def test_without_isoforms(self, tmp_path):
        prot_recs_all = U.add_species_name_prot_id({"A": [rec("g1")]})
        res = U.save_gene_id_mapping(prot_recs_all, None, out_path=str(tmp_path / "m.pickle"))
        assert res == {"A": [U.Gene(1000000001, "g1", None)]}
