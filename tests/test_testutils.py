"""Tests for the shared test infrastructure (conftest fixtures and fastoma_testutils)."""
import gzip

import pytest
from lxml import etree

from conftest import TESTDATA_DIR
from fastoma_testutils import (
    canonical_orthoxml,
    geneRef_counts,
    orthoxml_element_counts,
    orthoxml_validation_errors,
    read_group_partition,
    read_pairs,
)

EXPECTED = TESTDATA_DIR / "expected_output"

ORTHOXML_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<orthoXML xmlns="http://orthoXML.org/2011/" version="0.5" origin="test" originVersion="1">
  <species name="A" taxonId="2" NCBITaxId="0"><database name="A" version="1"><genes>
    <gene id="1" protId="a1"/><gene id="2" protId="a2"/></genes></database></species>
  <species name="B" taxonId="3" NCBITaxId="0"><database name="B" version="1"><genes>
    <gene id="3" protId="b1"/></genes></database></species>
  <taxonomy><taxon id="1" name="AB"><taxon id="2" name="A"/><taxon id="3" name="B"/></taxon></taxonomy>
  <groups>
    <orthologGroup id="{gid}"><property name="TaxRange" value="AB"/>
      {children}
    </orthologGroup>
  </groups>
</orthoXML>"""


def _xml(gid, children):
    return etree.fromstring(ORTHOXML_TEMPLATE.format(gid=gid, children=children).encode()).getroottree()


class TestConfFixtures:
    def test_subhog_conf_has_cli_defaults(self, subhog_conf):
        conf = subhog_conf()
        assert conf.number_of_samples_per_hog == 5
        assert conf.msa_filter_method == "col-row-threshold"
        assert conf.gap_ratio_row == 0.3 and conf.gap_ratio_col == 0.5
        assert conf.store_tcs_score is True

    def test_subhog_conf_overrides_and_rejects_unknown_options(self, subhog_conf):
        assert subhog_conf(number_of_samples_per_hog=3).number_of_samples_per_hog == 3
        assert subhog_conf().number_of_samples_per_hog == 5, "overrides must not leak into the base conf"
        with pytest.raises(AttributeError):
            subhog_conf(no_such_option=1)

    def test_roothog_conf_has_cli_defaults(self, roothog_conf):
        conf = roothog_conf()
        assert conf.mergHOG_fscore_thresh == 70
        assert conf.big_rhog_size == 50_000
        assert conf.min_sequence_length == 50


class TestPartitionsAndPairs:
    def test_partition_ignores_group_names_and_row_order(self, tmp_path):
        a = tmp_path / "a.tsv"
        b = tmp_path / "b.tsv"
        a.write_text("RootHOG\tProtein\nHOG:1\tx\nHOG:1\ty\nHOG:2\tz\n")
        b.write_text("RootHOG\tProtein\nHOG:7\tz\nHOG:9\ty\nHOG:9\tx\n")
        assert read_group_partition(a) == read_group_partition(b)
        assert read_group_partition(a) == {frozenset({"x", "y"}), frozenset({"z"})}

    def test_pairs_are_unordered_and_read_gzip(self, tmp_path):
        p = tmp_path / "pairs.tsv.gz"
        with gzip.open(p, "wt") as fh:
            fh.write("a\tb\nc\td\n")
        assert read_pairs(p) == {frozenset("ab"), frozenset("dc")}


class TestCanonicalOrthoxml:
    def test_independent_of_ids_and_child_order(self):
        x1 = _xml("HOG:1", '<geneRef id="1"/><geneRef id="3"/>')
        x2 = _xml("HOG:99", '<geneRef id="3"/><geneRef id="1"/>')
        assert canonical_orthoxml(x1) == canonical_orthoxml(x2)

    def test_detects_changed_nesting(self):
        flat = _xml("HOG:1", '<geneRef id="1"/><geneRef id="2"/><geneRef id="3"/>')
        dup = _xml("HOG:1", '<paralogGroup><geneRef id="1"/><geneRef id="2"/></paralogGroup><geneRef id="3"/>')
        assert canonical_orthoxml(flat) != canonical_orthoxml(dup)

    def test_geneRef_counts_use_protids(self):
        x = _xml("HOG:1", '<geneRef id="1"/><geneRef id="3"/>')
        assert geneRef_counts(x) == {"a1": 1, "b1": 1}

    def test_schema_rejects_invalid_document(self):
        # a paralogGroup with a single child violates the OrthoXML 0.5 schema
        bad = _xml("HOG:1", '<paralogGroup><geneRef id="1"/></paralogGroup><geneRef id="3"/>')
        assert orthoxml_validation_errors(bad)


class TestExpectedOutputConsistency:
    """testdata/expected_output is what the README promises; its files must agree with each other."""

    def test_orthoxml_is_schema_valid(self):
        assert orthoxml_validation_errors(EXPECTED / "FastOMA_HOGs.orthoxml") == []

    def test_every_protein_is_referenced_exactly_once(self):
        counts = geneRef_counts(EXPECTED / "FastOMA_HOGs.orthoxml")
        assert counts and max(counts.values()) == 1
        assert orthoxml_element_counts(EXPECTED / "FastOMA_HOGs.orthoxml")["gene"] == len(counts)

    def test_roothogs_match_orthoxml_top_level_groups(self):
        tree = etree.parse(str(EXPECTED / "FastOMA_HOGs.orthoxml"))
        canonical = canonical_orthoxml(tree)

        def genes(group):
            if group[0] == "gene":
                return {group[1]}
            return set().union(*(genes(c) for c in group[2]))

        from_xml = frozenset(frozenset(genes(g)) for g in canonical)
        assert from_xml == read_group_partition(EXPECTED / "RootHOGs.tsv")

    def test_orthologous_groups_are_within_roothogs(self):
        roothogs = read_group_partition(EXPECTED / "RootHOGs.tsv")
        for og in read_group_partition(EXPECTED / "OrthologousGroups.tsv"):
            assert any(og <= rh for rh in roothogs), og
