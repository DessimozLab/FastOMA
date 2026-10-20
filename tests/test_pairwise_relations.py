"""Tests for pairwise relation extraction (fastoma-helper pw-rel, FastOMA/zoo/hog/transform.py)."""
import gzip
import io
import subprocess
import sys
from itertools import combinations, product

import dendropy
import pytest
from lxml import etree

from FastOMA import collect_subhogs as cs
from FastOMA.zoo.hog import transform
from FastOMA.transformer import header_transformer
from conftest import TESTDATA_DIR
from fastoma_testutils import read_pairs
import output_fixtures as fx
from output_fixtures import prot_id

EXPECTED = TESTDATA_DIR / "expected_output"


def _cross(xs, ys):
    return {frozenset((prot_id(x), prot_id(y))) for x, y in product(xs, ys)}


def _within(xs):
    return {frozenset((prot_id(x), prot_id(y))) for x, y in combinations(xs, 2)}


EXPECTED_ORTHOLOGS = (
    set().union(*(_cross([f"a{i}"], [f"b{i}"]) for i in "1234"))
    | _cross(["c1"], ["a1", "b1", "a2", "b2"])
    | _cross(["c2"], ["a3", "b3", "a4", "b4"])
    | _cross(["d1"], ["a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "c1", "c2"])
    | _cross(["a5"], ["b5"])
    | _cross(["a6", "a7"], ["b6"])
)
EXPECTED_PARALOGS = (
    _cross(["a1", "b1"], ["a2", "b2"])
    | _cross(["a3", "b3"], ["a4", "b4"])
    | _cross(["a1", "b1", "a2", "b2", "c1"], ["a3", "b3", "a4", "b4", "c2"])
    | _within(["a6", "a7"])                      # fragment paralogGroup
    | _cross(["a5", "b5"], ["a6", "a7", "b6"])
)


@pytest.fixture
def fixture_orthoxml(tmp_path):
    inputs = fx.build_collect_inputs(tmp_path / "in")
    out = tmp_path / "hogs.orthoxml"
    cs.write_hog_orthoxml(str(inputs["pickle_folder"]), str(out), inputs["gene_id_pickle"],
                          header_transformer("noop"), str(inputs["species_tree"]))
    return out


def _pairs(orthoxml, rel_type, **kw):
    if not isinstance(orthoxml, etree._ElementTree):
        orthoxml = etree.parse(str(orthoxml))
    return list(transform.iter_pairwise_relations(orthoxml, rel_type=rel_type, **kw))


def _xml(groups, genes='<gene id="1" protId="a1"/><gene id="2" protId="a2"/><gene id="3" protId="b1"/>'):
    return etree.parse(io.BytesIO(
        f'<orthoXML xmlns="http://orthoXML.org/2011/" version="0.5" origin="t" originVersion="1">'
        f'<species name="A"><database name="d" version="1"><genes>{genes}</genes></database></species>'
        f'<groups>{groups}</groups></orthoXML>'.encode()))


class TestOrthoxmlRelations:
    def test_orthologs(self, fixture_orthoxml):
        pairs = _pairs(fixture_orthoxml, "ortholog", id_attribute="protId")
        assert len(pairs) == 25
        assert {frozenset(p) for p in pairs} == EXPECTED_ORTHOLOGS

    def test_paralogs(self, fixture_orthoxml):
        pairs = _pairs(fixture_orthoxml, "paralog", id_attribute="protId")
        assert len(pairs) == 40
        assert {frozenset(p) for p in pairs} == EXPECTED_PARALOGS

    def test_orthologs_and_paralogs_are_disjoint(self):
        assert not EXPECTED_ORTHOLOGS & EXPECTED_PARALOGS

    def test_default_is_ortholog_and_internal_ids(self, fixture_orthoxml):
        pairs = {frozenset(p) for p in _pairs(fixture_orthoxml, None)}
        assert frozenset(("1005", "2005")) in pairs and len(pairs) == 25

    def test_invalid_rel_type(self, fixture_orthoxml):
        with pytest.raises(ValueError):
            _pairs(fixture_orthoxml, "xenolog")

    def test_rejects_other_objects(self):
        with pytest.raises(ValueError):
            list(transform.iter_pairwise_relations("not an orthoxml"))

    def test_properties_and_scores_are_ignored(self):
        doc = _xml('<orthologGroup><score id="s" value="1"/><property name="TaxRange" value="X"/>'
                   '<geneRef id="1"/><geneRef id="3"/></orthologGroup>')
        assert _pairs(doc, "ortholog", id_attribute="protId") == [("a1", "b1")]

    def test_gene_without_attribute_gives_none(self):
        # characterises current behaviour: the KeyError branch (transform.py:82-85) is dead code
        doc = _xml('<orthologGroup><geneRef id="1"/><geneRef id="4"/></orthologGroup>',
                   genes='<gene id="1" protId="a1"/><gene id="4"/>')
        assert _pairs(doc, "ortholog", id_attribute="protId") == [("a1", None)]

    @pytest.mark.xfail(strict=True, reason="BUG: an orthologGroup without children raises TypeError in "
                                           "set.union() (zoo/hog/transform.py:52)")
    def test_empty_group_has_no_pairs(self):
        doc = _xml('<orthologGroup><geneRef id="1"/><geneRef id="3"/></orthologGroup><orthologGroup/>')
        assert _pairs(doc, "ortholog", id_attribute="protId") == [("a1", "b1")]


class TestTreeRelations:
    def test_docstring_example(self):
        tree = dendropy.Tree.get_from_string('(1,(2,(3,4))[&&NHX:Ev=D],(5,6));', schema='newick')
        assert sorted(transform.iter_pairwise_relations(tree, rel_type='paralog')) == [('2', '3'), ('2', '4')]


class TestHelperCli:
    def test_expected_output_orthologs(self, tmp_path, fixed_hashseed):
        out = tmp_path / "orthologs.tsv.gz"
        subprocess.run([sys.executable, "-m", "FastOMA.helper_scripts", "pw-rel", "--orthoxml",
                        str(EXPECTED / "FastOMA_HOGs.orthoxml"), "--out", str(out)],
                       check=True, env=fixed_hashseed)
        with gzip.open(out, "rt") as fh:
            lines = fh.read().splitlines()
        assert all(len(line.split("\t")) == 2 for line in lines)  # no header, two protIds
        assert read_pairs(out) == read_pairs(EXPECTED / "orthologs.tsv")
        assert len(lines) == len(read_pairs(out))

    def test_paralog_type_on_fixture(self, tmp_path, fixture_orthoxml, monkeypatch):
        from FastOMA import helper_scripts
        out = tmp_path / "paralogs.tsv"
        monkeypatch.setattr(sys, "argv", ["fastoma-helper", "pw-rel", "--type", "paralog",
                                          "--orthoxml", str(fixture_orthoxml), "--out", str(out)])
        helper_scripts.main()
        assert read_pairs(out) == EXPECTED_PARALOGS
