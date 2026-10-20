"""Tests for the output stage, FastOMA/collect_subhogs.py (fastoma-collect-subhogs)."""
import csv
import gzip
import os
import pickle
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from ete3 import Tree
from lxml import etree

from FastOMA import collect_subhogs as cs
from FastOMA.transformer import header_transformer
from fastoma_testutils import NS, geneRef_counts, orthoxml_validation_errors
import output_fixtures as fx
from output_fixtures import accession, prot_id


_REAL_WALK = os.walk


def _walk_in_order(order):
    """os.walk replacement that visits the pickle sub folders in the given order."""
    def walk(top, **kwargs):
        entries = list(_REAL_WALK(top, **kwargs))
        root, rest = entries[0], entries[1:]
        rest.sort(key=lambda e: order(Path(e[0]).name))
        yield root
        yield from rest
    return walk


@pytest.fixture
def sorted_walk(monkeypatch):
    # HOG numbering follows the os.walk order (see TestHogNumbering); fix it to rhog_1, rhog_2
    monkeypatch.setattr(cs.os, "walk", _walk_in_order(lambda name: name))


def _collect(inputs, transform="noop", out="FastOMA_HOGs.orthoxml"):
    tr = header_transformer(transform)
    cs.write_hog_orthoxml(str(inputs["pickle_folder"]), out, inputs["gene_id_pickle"], tr,
                          str(inputs["species_tree"]))
    cs.write_roothogs(Path(out), inputs["rhogs"], output_file_roothog_tsv="RootHOGs.tsv", id_transformer=tr)
    cs.write_group_files(Path(out), inputs["rhogs"], output_file_og_tsv="OrthologousGroups.tsv", id_transformer=tr)
    return Path(out)


@pytest.fixture
def inputs(in_tmp_dir):
    return fx.build_collect_inputs(in_tmp_dir / "in")


@pytest.fixture
def collected(inputs, sorted_walk):
    return _collect(inputs)


def _taxonomy():
    return cs.convert_speciestree_to_orthoxml_taxonomy(Tree(fx.SPECIES_TREE, format=1))


def _ids(elem):
    return [e.get("id") for e in elem.iter("orthologGroup")]


def _tsv_rows(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        return next(reader), {tuple(r) for r in reader}


class TestTaxonomy:
    def test_preorder_ids(self):
        root, name2taxid = _taxonomy()
        assert name2taxid == {"root": 1, "ABC": 2, "AB": 3, "A": 4, "B": 5, "C": 6, "D": 7}
        assert [(t.get("id"), t.get("name")) for t in root.iter("taxon")][:3] == [("1", "root"), ("2", "ABC"), ("3", "AB")]
        ab = root.find("./taxon/taxon/taxon")
        assert [c.get("name") for c in ab] == ["A", "B"]

    def test_unnamed_internal_nodes_collide(self):
        # characterises current behaviour: relies on check_input naming all internal nodes
        root, name2taxid = cs.convert_speciestree_to_orthoxml_taxonomy(Tree("((A,B),C);", format=1))
        assert len(list(root.iter("taxon"))) == 5
        assert name2taxid == {"": 2, "A": 3, "B": 4, "C": 5}


class TestUpdateHogids:
    def test_nested_duplications(self):
        _, name2taxid = _taxonomy()
        hog = cs.update_hogids(1, fx.rhog1_hog(), name2taxid)
        assert _ids(hog) == [
            "HOG:0000001_1",
            "HOG:0000001.1a_2", "HOG:0000001.1a.1a_3", "HOG:0000001.1a.1b_3",
            # the duplication counter is per depth over the whole rootHOG, not per parent
            "HOG:0000001.1b_2", "HOG:0000001.1b.2a_3", "HOG:0000001.1b.2b_3",
        ]
        assert [e.get("taxonId") for e in hog.iter("orthologGroup")] == ["1", "2", "3", "3", "2", "3", "3"]

    def test_root_level_duplication(self):
        _, name2taxid = _taxonomy()
        hog = cs.update_hogids(42, fx.rhog2_hog(), name2taxid)
        assert _ids(hog) == ["HOG:0000042_3", "HOG:0000042.1a_3", "HOG:0000042.1b_3"]
        assert all(p.get("id") is None for p in hog.iter("paralogGroup"))

    def test_letter_overflow_after_z(self):
        _, name2taxid = _taxonomy()
        pg = fx._pg(*[fx._og("AB", fx._g("a1"), fx._g("b1"), score=False) for _ in range(30)])
        hog = cs.update_hogids(7, fx._og("root", pg, fx._g("c1"), id_="HOG_E1_sub1", score=False), name2taxid)
        ids = [o.get("id") for o in pg]
        assert ids[0] == "HOG:0000007.1a_3"
        assert ids[25] == "HOG:0000007.1z_3"
        assert ids[26:] == ["HOG:0000007.1aa_3", "HOG:0000007.1ab_3", "HOG:0000007.1ac_3", "HOG:0000007.1ad_3"]
        assert len(set(ids)) == 30

    def test_property_in_paralog_group_shifts_letters(self):
        # characterises current behaviour: the <property> child of a fragment paralogGroup
        # (_hog_class.py:521) takes index 0, so its orthologGroup children start at 'b'.
        # _hog_class only puts geneRefs into such groups today, so no id is affected yet.
        _, name2taxid = _taxonomy()
        pg = fx._pg(fx._og("AB", fx._g("a1"), fx._g("b1")), fx._og("AB", fx._g("a2"), fx._g("b2")),
                    prop="Dubiousfragment")
        cs.update_hogids(1, fx._og("root", pg, fx._g("c1"), id_="HOG_E1_sub1"), name2taxid)
        assert [o.get("id") for o in pg.iter("orthologGroup")] == ["HOG:0000001.1b_3", "HOG:0000001.1c_3"]

    def test_omamer_roothog_property_after_scores(self):
        _, name2taxid = _taxonomy()
        hog = cs.update_hogids(1, fx.rhog1_hog(), name2taxid)
        assert [(c.tag, c.get("id") or c.get("name")) for c in hog][:3] == [
            ("score", "CompletenessScore"), ("property", "OMAmerRootHOG"), ("property", "TaxRange")]
        assert hog[1].get("value") == "HOG:E0000001"

    def test_omamer_roothog_property_without_scores(self):
        _, name2taxid = _taxonomy()
        hog = cs.update_hogids(1, fx.rhog1_hog(score=False), name2taxid)
        assert (hog[0].get("name"), hog[0].get("value")) == ("OMAmerRootHOG", "HOG:E0000001")

    def test_omamer_roothog_from_unexpected_id(self):
        # characterises current behaviour: the first two '_' separated fields joined by ':'
        _, name2taxid = _taxonomy()
        hog = cs.update_hogids(1, fx._og("AB", fx._g("a1"), fx._g("b1"), id_="weird"), name2taxid)
        assert hog.find("./property[@name='OMAmerRootHOG']").get("value") == "weird"

    def test_missing_taxrange_raises(self):
        hog = ET.Element("orthologGroup", {"id": "HOG_E1_sub1"})
        with pytest.raises(AttributeError):
            cs.update_hogids(1, hog, {"root": 1})


class TestWriteHogOrthoxml:
    def test_header(self, collected):
        root = etree.parse(str(collected)).getroot()
        assert root.tag == "{http://orthoXML.org/2011/}orthoXML"
        assert root.get("version") == "0.5"
        assert root.get("origin").startswith("FastOMA ")
        assert root.get("originVersion")  # a timestamp
        species = {s.get("name"): (s.get("taxonId"), s.get("NCBITaxId")) for s in root.iterfind("o:species", NS)}
        assert species == {"A": ("4", "0"), "B": ("5", "0"), "C": ("6", "0"), "D": ("7", "0")}
        genes = {g.get("id"): g.get("protId") for g in root.iterfind(".//o:gene", NS)}
        assert len(genes) == 28 and genes["2003"] == "sp|P02003|B3_B"
        assert all(g.get("main_isoform") is None for g in root.iterfind(".//o:gene", NS))

    def test_groups_reference_every_placed_protein_once(self, collected):
        counts = geneRef_counts(collected)
        assert set(counts) == {prot_id(n) for n in fx.RHOG1_MEMBERS + fx.RHOG2_MEMBERS}
        assert max(counts.values()) == 1

    def test_hog_ids_in_file(self, collected):
        tree = etree.parse(str(collected))
        tops = tree.findall("./o:groups/o:orthologGroup", NS)
        assert [t.get("id") for t in tops] == ["HOG:0000001_1", "HOG:0000002_3"]
        assert [t.find("o:property[@name='OMAmerRootHOG']", NS).get("value") for t in tops] == [
            "HOG:E0000001", "HOG:E0000002"]
        assert len(tree.findall(".//o:orthologGroup", NS)) == 10

    def test_score_defs_only_for_present_scores(self, collected):
        root = etree.parse(str(collected)).getroot()
        assert [d.get("id") for d in root.iterfind("o:scores/o:scoreDef", NS)] == ["CompletenessScore"]

    def test_no_scores_element_without_scores(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in", scores=False)
        out = _collect(inputs)
        assert etree.parse(str(out)).find("o:scores", NS) is None

    def test_main_isoform_attribute(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in", isoforms=True)
        out = _collect(inputs)
        genes = {g.get("id"): g.get("main_isoform") for g in etree.parse(str(out)).iterfind(".//o:gene", NS)}
        assert genes["1007"] == "1001" and genes["1001"] == "1001" and genes["1002"] is None

    def test_nested_duplications_are_schema_valid(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in")
        (inputs["pickle_folder"] / "rhog_2" / "file_HOG_E0000002.pickle").unlink()
        assert orthoxml_validation_errors(_collect(inputs)) == []

    @pytest.mark.xfail(strict=True, reason="BUG: root-level duplication is wrapped in an orthologGroup with a "
                                           "single paralogGroup child, schema needs >= 2 (_infer_subhog.py:112-121)")
    def test_root_level_duplication_is_schema_valid(self, collected):
        assert orthoxml_validation_errors(collected) == []

    @pytest.mark.xfail(strict=True, reason="BUG: gene@main_isoform is not allowed by the OrthoXML 0.5 schema "
                                           "(collect_subhogs.py:208-209)")
    def test_main_isoform_output_is_schema_valid(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in", isoforms=True)
        (inputs["pickle_folder"] / "rhog_2" / "file_HOG_E0000002.pickle").unlink()
        assert orthoxml_validation_errors(_collect(inputs)) == []

    def test_only_file_pickles_are_read(self, inputs, sorted_walk):
        with open(inputs["pickle_folder"] / "rhog_1" / "other.pickle", "wb") as fh:
            pickle.dump([fx.rhog2_hog()], fh)
        assert [h.get("id") for h in cs.iter_hogs(inputs["pickle_folder"])] == [
            "HOG_E0000001_sub10001", "HOG_E0000002_sub10002"]

    def test_empty_pickle_folder(self, inputs, sorted_walk, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        inputs["pickle_folder"] = empty
        out = _collect(inputs, out="empty.orthoxml")
        tree = etree.parse(str(out))
        assert tree.find("o:groups", NS) is not None and len(tree.find("o:groups", NS)) == 0
        assert _tsv_rows("RootHOGs.tsv")[1] == set()


class TestHogNumbering:
    def _top_ids(self, out):
        return {t.find("o:property[@name='OMAmerRootHOG']", NS).get("value"): t.get("id")
                for t in etree.parse(str(out)).findall("./o:groups/o:orthologGroup", NS)}

    def test_numbering_follows_pickle_listing_order(self, inputs, monkeypatch):
        # characterises current behaviour: rootHOGs are numbered in os.walk order (collect_subhogs.py:223)
        monkeypatch.setattr(cs.os, "walk", _walk_in_order(lambda name: name), raising=True)
        assert self._top_ids(_collect(inputs, out="a.orthoxml")) == {"HOG:E0000001": "HOG:0000001_1",
                                                                      "HOG:E0000002": "HOG:0000002_3"}
        monkeypatch.setattr(cs.os, "walk", _walk_in_order(lambda name: [-ord(c) for c in name]))
        assert self._top_ids(_collect(inputs, out="b.orthoxml")) == {"HOG:E0000001": "HOG:0000002_1",
                                                                      "HOG:E0000002": "HOG:0000001_3"}

    @pytest.mark.xfail(strict=True, reason="BUG?/contract undecided: HOG ids depend on the directory listing "
                                           "order of the pickle folders, i.e. on Nextflow task completion "
                                           "(collect_subhogs.py:51,223, FastOMA.nf collect_subhogs stageAs)")
    def test_numbering_independent_of_listing_order(self, inputs, monkeypatch):
        monkeypatch.setattr(cs.os, "walk", _walk_in_order(lambda name: name))
        first = self._top_ids(_collect(inputs, out="a.orthoxml"))
        monkeypatch.setattr(cs.os, "walk", _walk_in_order(lambda name: [-ord(c) for c in name]))
        assert self._top_ids(_collect(inputs, out="b.orthoxml")) == first


class TestRootHogFiles:
    def test_tsv(self, collected):
        header, rows = _tsv_rows("RootHOGs.tsv")
        assert header == ["RootHOG", "Protein", "OMAmerRootHOG"]
        assert rows == ({("HOG:0000001", prot_id(n), "HOG:E0000001") for n in fx.RHOG1_MEMBERS}
                        | {("HOG:0000002", prot_id(n), "HOG:E0000002") for n in fx.RHOG2_MEMBERS})

    def test_fasta_files(self, collected):
        assert sorted(os.listdir("RootHOGsFasta")) == ["HOG0000001.fa.gz", "HOG0000002.fa.gz"]
        with gzip.open("RootHOGsFasta/HOG0000002.fa.gz", "rt") as fh:
            recs = {r.id: str(r.seq) for r in SeqIO.parse(fh, "fasta")}
        assert recs == {prot_id(n): fx.sequence(n) for n in fx.RHOG2_MEMBERS}

    def test_fasta_header_format(self, collected):
        # characterises current behaviour: the internal 'id||species||numericId' stays in the
        # description because only rec.id is replaced (collect_subhogs.py:305-306)
        with gzip.open("RootHOGsFasta/HOG0000001.fa.gz", "rt") as fh:
            headers = {line.rstrip("\n") for line in fh if line.startswith(">")}
        assert f">{prot_id('a1')} {prot_id('a1')}||A||1001 desc of a1 [A]" in headers
        assert len(headers) == len(fx.RHOG1_MEMBERS)


class TestOrthologousGroupFiles:
    def test_tsv(self, collected):
        header, rows = _tsv_rows("OrthologousGroups.tsv")
        assert header == ["Group", "Protein"]
        # heaviest paralogGroup child wins, ties go to the first child (zoo/hog/extract_groups.py:143-148)
        assert rows == ({("OG_0000001", prot_id(n)) for n in ("a1", "b1", "c1", "d1")}
                        | {("OG_0000002", prot_id(n)) for n in ("a5", "b5")})

    def test_heavier_paralog_child_wins(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in")
        hog = fx._og("root", fx._pg(fx._og("AB", fx._g("a1"), fx._g("b1")),
                                    fx._og("ABC", fx._g("a2"), fx._g("b2"), fx._g("c2"))),
                     fx._g("d1"), id_="HOG_E0000001_sub1")
        with open(inputs["pickle_folder"] / "rhog_1" / "file_HOG_E0000001.pickle", "wb") as fh:
            pickle.dump([hog], fh)
        _collect(inputs)
        rows = _tsv_rows("OrthologousGroups.tsv")[1]
        assert {p for g, p in rows if g == "OG_0000001"} == {prot_id(n) for n in ("a2", "b2", "c2", "d1")}

    def test_at_most_one_protein_per_species(self, collected):
        groups = {}
        for g, p in _tsv_rows("OrthologousGroups.tsv")[1]:
            groups.setdefault(g, []).append(p.rsplit("_", 1)[1])
        assert all(len(sp) == len(set(sp)) >= 2 for sp in groups.values())

    def test_groups_smaller_than_two_are_dropped(self, in_tmp_dir, sorted_walk):
        inputs = fx.build_collect_inputs(in_tmp_dir / "in")
        # a duplication of A only, as wrapped by _infer_subhog: the marker group would be {a5}
        hog = fx._og("A", fx._pg(fx._g("a5"), fx._g("a6")), id_="HOG_E0000002_sub1")
        with open(inputs["pickle_folder"] / "rhog_2" / "file_HOG_E0000002.pickle", "wb") as fh:
            pickle.dump([hog], fh)
        _collect(inputs)
        assert {g for g, _ in _tsv_rows("OrthologousGroups.tsv")[1]} == {"OG_0000001"}
        assert {g for g, _, _ in _tsv_rows("RootHOGs.tsv")[1]} == {"HOG:0000001", "HOG:0000002"}
        assert sorted(os.listdir("OrthologousGroupsFasta")) == ["OG_0000001.fa.gz"]

    def test_fasta_files(self, collected):
        assert sorted(os.listdir("OrthologousGroupsFasta")) == ["OG_0000001.fa.gz", "OG_0000002.fa.gz"]
        with gzip.open("OrthologousGroupsFasta/OG_0000001.fa.gz", "rt") as fh:
            recs = list(SeqIO.parse(fh, "fasta"))
        assert {r.id for r in recs} == {prot_id(n) for n in ("a1", "b1", "c1", "d1")}
        assert {r.description.rsplit(" ", 1)[1] for r in recs} == {"[A]", "[B]", "[C]", "[D]"}

    def test_protein_missing_from_roothog_fasta_is_silently_dropped(self, in_tmp_dir, sorted_walk):
        # characterises current behaviour: no error if a group member is not in the rootHOG fasta
        inputs = fx.build_collect_inputs(in_tmp_dir / "in")
        fa = inputs["rhogs"] / "HOG_E0000002.fa"
        fa.write_text("".join(line for line in fa.read_text().splitlines(keepends=True)
                              if "A5_A" not in line and line != "MKVA5\n"))
        _collect(inputs)
        with gzip.open("OrthologousGroupsFasta/OG_0000002.fa.gz", "rt") as fh:
            assert [r.id for r in SeqIO.parse(fh, "fasta")] == [prot_id("b5")]
        assert ("OG_0000002", prot_id("a5")) in _tsv_rows("OrthologousGroups.tsv")[1]


class TestUniProtIdTransform:
    @pytest.fixture
    def collected_up(self, inputs, sorted_walk):
        return _collect(inputs, transform="UniProt")

    def test_orthoxml_protids(self, collected_up):
        prot_ids = {g.get("protId") for g in etree.parse(str(collected_up)).iterfind(".//o:gene", NS)}
        assert prot_ids == {accession(n) for n in fx.PROTEIN_NAMES}

    def test_tsvs_and_fastas_use_the_same_ids(self, collected_up):
        roothog_prots = {p for _, p, _ in _tsv_rows("RootHOGs.tsv")[1]}
        assert roothog_prots == {accession(n) for n in fx.RHOG1_MEMBERS + fx.RHOG2_MEMBERS}
        og_prots = {p for g, p in _tsv_rows("OrthologousGroups.tsv")[1] if g == "OG_0000002"}
        assert og_prots == {accession("a5"), accession("b5")}
        with gzip.open("OrthologousGroupsFasta/OG_0000002.fa.gz", "rt") as fh:
            assert {r.id for r in SeqIO.parse(fh, "fasta")} == og_prots
        with gzip.open("RootHOGsFasta/HOG0000001.fa.gz", "rt") as fh:
            assert {r.id for r in SeqIO.parse(fh, "fasta")} == {accession(n) for n in fx.RHOG1_MEMBERS}


class TestCli:
    def test_end_to_end(self, inputs, sorted_walk, monkeypatch):
        import sys
        monkeypatch.setattr(sys, "argv", [
            "fastoma-collect-subhogs", "--pickle-folder", str(inputs["pickle_folder"]),
            "--roothogs-folder", str(inputs["rhogs"]), "--gene-id-pickle-file", str(inputs["gene_id_pickle"]),
            "--out", "FastOMA_HOGs.orthoxml", "--marker-groups-fasta", "OrthologousGroups.tsv",
            "--roothog-tsv", "RootHOGs.tsv", "--species-tree", str(inputs["species_tree"])])
        cs.fastoma_collect_subhogs()
        for name in ("FastOMA_HOGs.orthoxml", "RootHOGs.tsv", "OrthologousGroups.tsv",
                     "RootHOGsFasta", "OrthologousGroupsFasta"):
            assert Path(name).exists(), name
