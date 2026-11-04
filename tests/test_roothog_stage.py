"""Stage-level characterisation tests for fastoma-infer-roothogs, run as a subprocess.

The rootHOG partition depends on set iteration order (see test_roothog_grouping.py), so the
subprocesses run with a fixed PYTHONHASHSEED and we compare partitions, never file names.
"""
import os
import random
import subprocess
import sys

import pytest

from conftest import TESTDATA_DIR
from roothog_helpers import fasta_ids, read_rhog_partition, write_fasta, write_hogmap


def run_infer_roothogs(workdir, proteomes, hogmaps, *extra):
    cmd = [sys.executable, "-m", "FastOMA.infer_roothogs", "--proteomes", str(proteomes),
           "--hogmap", str(hogmaps), "--out-rhog-folder", "omamer_rhogs", *map(str, extra)]
    subprocess.run(cmd, cwd=workdir, check=True, capture_output=True,
                   env=dict(os.environ, PYTHONHASHSEED="0"))
    return workdir / "omamer_rhogs"


def test_testdata_partition(tmp_path):
    out = run_infer_roothogs(tmp_path, TESTDATA_DIR / "in_folder" / "proteome",
                             TESTDATA_DIR / "expected_output" / "hogmap")
    expected = [
        {"sp|O67618|LEPA_AQUAE", "sp|O66429|EFTU_AQUAE", "sp|O84067|LEPA_CHLTR", "sp|P0CD71|EFTU_CHLTR",
         "sp|P13927|EFTU_MYCGE"},
        {"sp|O67118|DNAK_AQUAE", "sp|P17821|DNAK_CHLTR", "sp|P47547|DNAK_MYCGE"},
        {"sp|O67032|RF1_AQUAE", "sp|O84026|RF1_CHLTR", "sp|P47500|RF1_MYCGE"},
        {"sp|O67828|ATPB_AQUAE", "sp|O66907|ATPA_AQUAE", "sp|P47639|ATPB_MYCGE", "sp|P47641|ATPA_MYCGE"},
        {"sp|O67943|CH60_AQUAE", "sp|P0C0Z7|CH60_CHLTR"},
        {"sp|O66778|ENO_AQUAE", "sp|O84591|ENO_CHLTR"},
        {"sp|O67736|FOLD_AQUAE", "sp|O84081|FOLD_CHLTR"},
        {"sp|O67161|G3P_AQUAE", "sp|P0CE13|G3P_CHLTR", "sp|P47543|G3P_MYCGE"},
        {"sp|O67547|SUCD_AQUAE", "tr|O84829|O84829_CHLTR"},
        {"sp|O66686|TPIS_AQUAE", "sp|O84332|TPIS_CHLTR"},
    ]
    assert read_rhog_partition(out) == frozenset(frozenset(g) for g in expected)
    assert len(list(out.glob("HOG_*.fa"))) == len(expected)


# --- synthetic dataset exercising all branches of the stage ------------------------------------------
#
# Building and running takes a few seconds (mmseqs), hence `slow`.
# CLI options: --big-rhog-size 30 (so merge candidates must be < 5 proteins), --mergHOG-shared-thresh 1
SPECIES = ["S1", "S2", "S3", "S4"]
BIG = [f"{sp}_b{i}" for sp in SPECIES for i in range(8)]
EXPECTED_PARTITION = frozenset(frozenset(g) for g in [
    {"S1_a", "S2_a", "S3_a", "S4_x"},      # S4_x is a singleton resolved via its alternative hit to A1
    {"S1_p", "S2_q"},                      # singletons sharing the alternative rootHOG Q1
    {"S1_m", "S2_m", "S3_m", "S4_m"},      # rootHOGs M1 and M2 merged
    set(BIG[0::2]),                        # big rootHOG B1, low-score half filtered out
])
SENT_TO_LINCLUST = {"S3_u", "S4_n", "S2_n1", "S2_n2"}
EXPECTED_TO_VANISH = {"S1_short", "S1_a2"} | set(BIG[1::2])   # too short / minor isoform / big-rhog filter
KNOWN_BUG_VANISH = {"S2_missing"}                               # absent from the hogmap


def build_synthetic(root):
    rng = random.Random(42)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def seq():
        return "M" + "".join(rng.choice(aa) for _ in range(79))

    hits = {sp: [] for sp in SPECIES}
    prots = {sp: [] for sp in SPECIES}

    def add(sp, prot, *rows, sequence=None):
        prots[sp].append((prot, sequence or seq()))
        hits[sp].extend((prot, hogid, score) for hogid, score in rows)

    add("S1", "S1_a", ("HOG:A1.1a", 500))
    add("S1", "S1_a2", ("HOG:A1", 200))             # minor isoform of S1_a
    add("S2", "S2_a", ("HOG:A1", 500))
    add("S3", "S3_a", ("HOG:A1.2b", 500))
    add("S4", "S4_x", ("HOG:Z9", 500), ("HOG:A1.3", 300))
    add("S1", "S1_p", ("HOG:P1", 500), ("HOG:Q1", 300))
    add("S2", "S2_q", ("HOG:R1", 500), ("HOG:Q1.1a", 300))
    add("S1", "S1_m", ("HOG:M1", 500), ("HOG:M2", 300))
    add("S2", "S2_m", ("HOG:M1", 500), ("HOG:M2", 300))
    add("S3", "S3_m", ("HOG:M2", 500), ("HOG:M1", 300))
    add("S4", "S4_m", ("HOG:M2", 500), ("HOG:M1", 300))
    for i, prot in enumerate(BIG):
        add(prot[:2], prot, ("HOG:B1.1a", 99 if i % 2 == 0 else 90))
    shared = seq()                                   # identical sequences: linclust must group them
    add("S3", "S3_u", ("HOG:U1", 500), sequence=shared)
    add("S4", "S4_n", ("N/A", None), sequence=shared)
    add("S2", "S2_n1", ("N/A", None))
    add("S2", "S2_n2", ("N/A", None))                # minor isoform of S2_n1, both unmapped
    add("S1", "S1_short", ("HOG:A1", 500), sequence="MKVLAAGIVL")
    prots["S2"].append(("S2_missing", seq()))

    for sp in SPECIES:
        write_fasta(root / "proteome" / f"{sp}.fa", prots[sp])
        write_hogmap(root / "hogmap" / f"{sp}.fa.hogmap", hits[sp])
    (root / "splice").mkdir()
    (root / "splice" / "S1.splice").write_text("S1_a;S1_a2\n")
    (root / "splice" / "S2.splice").write_text("S2_n1;S2_n2\n")
    return {p for sp in SPECIES for p, _ in prots[sp]}


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("synthetic")
    all_proteins = build_synthetic(root)
    out = run_infer_roothogs(root, root / "proteome", root / "hogmap", "--splice", root / "splice",
                             "--big-rhog-size", 30, "--mergHOG-shared-thresh", 1)
    return root, out, all_proteins


@pytest.mark.slow
def test_synthetic_partition(synthetic_run):
    _, out, _ = synthetic_run
    assert read_rhog_partition(out, "HOG_[!c]*.fa") == EXPECTED_PARTITION


@pytest.mark.slow
def test_synthetic_selected_isoforms_written(synthetic_run):
    root, _, _ = synthetic_run
    assert fasta_ids(root / "selected_isoforms" / "S1_selected_isoforms.fa").count("S1_a2") == 0


@pytest.mark.slow
@pytest.mark.tools("mmseqs")
def test_synthetic_linclust(synthetic_run):
    root, out, _ = synthetic_run
    assert set(fasta_ids(root / "singleton_unmapped.fa")) == SENT_TO_LINCLUST
    assert read_rhog_partition(out, "HOG_clust*.fa") == {frozenset({"S3_u", "S4_n"})}


@pytest.mark.slow
@pytest.mark.tools("mmseqs")
def test_synthetic_every_protein_accounted_for(synthetic_run):
    root, out, all_proteins = synthetic_run
    in_rhogs = {p for g in read_rhog_partition(out) for p in g}
    in_singletons = set(fasta_ids(root / "singleton_unmapped.fa"))
    unaccounted = all_proteins - in_rhogs - in_singletons - EXPECTED_TO_VANISH
    assert unaccounted == KNOWN_BUG_VANISH


@pytest.mark.slow
@pytest.mark.tools("mmseqs")
@pytest.mark.xfail(strict=True, reason="BUG: proteins without any hogmap row are dropped silently "
                                       "(parse_hogmap_omamer, _utils_roothog.py:264)")
def test_synthetic_protein_missing_from_hogmap_not_lost(synthetic_run):
    root, out, _ = synthetic_run
    kept = {p for g in read_rhog_partition(out) for p in g} | set(fasta_ids(root / "singleton_unmapped.fa"))
    assert KNOWN_BUG_VANISH <= kept


@pytest.mark.slow
@pytest.mark.tools("mmseqs")
@pytest.mark.xfail(strict=True, reason="BUG: unmapped proteins are not filtered by isoform selection, minor "
                                       "isoforms reach linclust (infer_roothogs.py:50-60, :93)")
def test_synthetic_minor_isoforms_not_clustered(synthetic_run):
    root, _, _ = synthetic_run
    assert len({"S2_n1", "S2_n2"} & set(fasta_ids(root / "singleton_unmapped.fa"))) == 1
