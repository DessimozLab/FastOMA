"""Property tests of fastoma-infer-subhogs on real rootHOGs (tests/data/subhog).

small/    D0455823 (45 proteins, 19 species), D0111149 (2 proteins)
conflict/ D0625237 (35 proteins; split_hog loses proteins), D0662268 (14 proteins; conflict resolution)
The species trees are the original trees pruned to the species of the rootHOGs (branch lengths set to 1).
"""
import os
import pickle
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from ete3 import Tree
from lxml import etree

from FastOMA._utils_subhog import prepare_species_tree
from fastoma_testutils import NS, canonical_orthoxml, geneRef_counts, orthoxml_validation_errors

DATA = Path(__file__).parent / "data" / "subhog"

# run the CLI, but keep the per-level pickles so the final HOG objects can be inspected
RUNNER = ("import FastOMA._infer_subhog as m; m.keep_subhog_each_pickle = True; "
          "from FastOMA.infer_subhogs import fastoma_infer_subhogs; fastoma_infer_subhogs()")

_ROOTHOGS = [("small", "D0455823"), ("small", "D0111149"), ("conflict", "D0625237"), ("conflict", "D0662268")]


def roothogs(xfail=None):
    """parametrisation over the rootHOGs; ``xfail`` maps rootHOG ids to the reason of a known bug."""
    params = []
    for dataset, rhogid in _ROOTHOGS:
        marks = [pytest.mark.slow] if dataset == "conflict" else []
        if xfail and rhogid in xfail:
            marks.append(pytest.mark.xfail(strict=True, reason=xfail[rhogid]))
        params.append(pytest.param(dataset, rhogid, id=rhogid, marks=marks))
    return params


ROOTHOGS = roothogs()

pytestmark = pytest.mark.tools("mafft", "FastTree")


class Run:
    def __init__(self, dataset, workdir, parallel):
        self.dataset = dataset
        self.workdir = workdir
        self.species_tree_path = DATA / f"species_tree_{dataset}.nwk"
        cmd = [sys.executable, "-c", RUNNER,
               "--input-rhog-folder", str(DATA / dataset),
               "--species-tree", str(self.species_tree_path),
               "--output-pickles", "pickles", "-vv"]
        if parallel:
            cmd.append("--parallel")
        # same environment as the fixed_hashseed fixture (which is function scoped)
        env = dict(os.environ, PYTHONHASHSEED="0")
        res = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
        assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]

    def orthoxml(self, rhogid):
        return etree.parse(str(self.workdir / f"ortho_{rhogid}.orthoxml"))

    def input_ids(self, rhogid):
        return [r.id for r in SeqIO.parse(str(DATA / self.dataset / f"HOG_{rhogid}.fa"), "fasta")]

    def final_hogs(self, rhogid):
        recs = list(SeqIO.parse(str(DATA / self.dataset / f"HOG_{rhogid}.fa"), "fasta"))
        root = prepare_species_tree(recs, Tree(str(self.species_tree_path), format=1), rhogid)[0]
        with open(self.workdir / f"rhog_{rhogid}" / f"{root.name}.pickle", "rb") as fh:
            return pickle.load(fh)


_runs = {}


@pytest.fixture(scope="module")
def infer_run(tmp_path_factory):
    def get(dataset, parallel=False):
        key = (dataset, parallel)
        if key not in _runs:
            workdir = tmp_path_factory.mktemp(f"subhogs_{dataset}_{'par' if parallel else 'seq'}")
            _runs[key] = Run(dataset, workdir, parallel)
        return _runs[key]
    yield get
    _runs.clear()


def gene_species(tree):
    return {g.get("id"): sp.get("name") for sp in tree.iterfind(".//o:species", NS) for g in sp.iterfind(".//o:gene", NS)}


def taxrange(elem):
    p = elem.find("o:property[@name='TaxRange']", NS)
    return None if p is None else p.get("value")


def ortholog_groups(tree):
    return tree.findall(".//o:groups//o:orthologGroup", NS)


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_no_protein_twice(infer_run, dataset, rhogid):
    counts = geneRef_counts(infer_run(dataset).orthoxml(rhogid))
    assert counts and max(counts.values()) == 1


@pytest.mark.parametrize("dataset,rhogid", roothogs(xfail={
    "D0625237": "BUG: split_hog drops subhogs not reached by any partition (_hog_class.py:636-643)"}))
def test_all_input_proteins_in_final_hogs(infer_run, dataset, rhogid):
    run = infer_run(dataset)
    members = [p for h in run.final_hogs(rhogid) for m in h.get_members() for p in m.split("_|_")]
    assert sorted(members) == sorted(run.input_ids(rhogid))


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_written_genes_are_final_hogs_with_two_or_more_members(infer_run, dataset, rhogid):
    run = infer_run(dataset)
    hogs = run.final_hogs(rhogid)
    expected = {m.split("||")[2] for h in hogs if len(h) >= 2 for m in h.get_members()}
    written = set(geneRef_counts(run.orthoxml(rhogid), use_protid=False))
    assert written == expected


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_top_level_groups_have_two_or_more_genes(infer_run, dataset, rhogid):
    tree = infer_run(dataset).orthoxml(rhogid)
    for group in tree.find(".//o:groups", NS):
        assert len(group.findall(".//o:geneRef", NS)) >= 2


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_taxrange_nesting_follows_species_tree(infer_run, dataset, rhogid):
    run = infer_run(dataset)
    sptree = Tree(str(run.species_tree_path), format=1)
    nodes = {n.name: n for n in sptree.traverse()}
    for group in ortholog_groups(run.orthoxml(rhogid)):
        parent = group.getparent()
        while parent is not None and etree.QName(parent).localname == "paralogGroup":
            parent = parent.getparent()
        if etree.QName(parent).localname != "orthologGroup":
            continue
        child_node, parent_node = nodes[taxrange(group)], nodes[taxrange(parent)]
        if group.getparent() is parent:
            # directly nested orthologGroups are at strictly lower levels
            assert parent_node in child_node.get_ancestors()
        else:
            # below a duplication the copies can span the same clade
            assert child_node is parent_node or parent_node in child_node.get_ancestors()


@pytest.mark.parametrize("dataset,rhogid", roothogs(xfail={
    "D0625237": "BUG: root-level wrapper uses the HOG's level as TaxRange instead of the MRCA (_infer_subhog.py:120)"}))
def test_taxrange_is_mrca_of_member_species(infer_run, dataset, rhogid):
    run = infer_run(dataset)
    sptree = Tree(str(run.species_tree_path), format=1)
    tree = run.orthoxml(rhogid)
    species = gene_species(tree)
    for group in ortholog_groups(tree):
        member_species = {species[g.get("id")] for g in group.iterfind(".//o:geneRef", NS)}
        mrca = sptree.get_common_ancestor(*member_species) if len(member_species) > 1 else sptree & member_species.pop()
        assert taxrange(group) == mrca.name, group.get("id")


@pytest.mark.parametrize("dataset,rhogid", roothogs(xfail={
    "D0625237": "BUG: root-level paralogGroup wrapped in an orthologGroup with a single child (_infer_subhog.py:112-122)"}))
def test_orthoxml_schema(infer_run, dataset, rhogid):
    tree = infer_run(dataset).orthoxml(rhogid)
    # characterises current behaviour: the per-rootHOG debug file has no version/originVersion attributes
    tree.getroot().set("version", "0.5")
    tree.getroot().set("originVersion", "test")
    assert orthoxml_validation_errors(tree) == []


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_parallel_and_sequential_give_same_hogs(infer_run, dataset, rhogid):
    seq = infer_run(dataset).orthoxml(rhogid)
    par = infer_run(dataset, parallel=True).orthoxml(rhogid)
    assert canonical_orthoxml(seq) == canonical_orthoxml(par)


@pytest.mark.parametrize("dataset,rhogid", ROOTHOGS)
def test_output_pickle_matches_orthoxml(infer_run, dataset, rhogid):
    run = infer_run(dataset)
    with open(run.workdir / "pickles" / f"file_{rhogid}.pickle", "rb") as fh:
        groups = pickle.load(fh)
    assert len(groups) == len(run.orthoxml(rhogid).find(".//o:groups", NS))


def test_cli_defaults(subhog_conf):
    conf = subhog_conf()
    assert conf.threshold_dubious_sd == pytest.approx(0.1)
    assert conf.number_of_samples_per_hog == 5
    assert conf.overlap_fragments == pytest.approx(0.15)
    assert conf.gene_rooting_method == "midpoint"
    assert conf.msa_filter_method == "col-row-threshold"
    assert (conf.gap_ratio_row, conf.gap_ratio_col, conf.min_col_trim) == (0.3, 0.5, 50)
    assert (conf.store_completeness_score, conf.store_implied_losses_score, conf.store_tcs_score) == (True, True, True)
    assert not conf.parallel and not conf.msa_write and not conf.gene_trees_write
    assert conf.output_pickles == "pickle_hogs"
