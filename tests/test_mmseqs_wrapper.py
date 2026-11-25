"""Tests for the mmseqs easy-cluster wrapper (FastOMA/zoo/wrappers/clusterers/mmseqs.py)."""
import os
import stat

import pytest

from FastOMA.zoo.wrappers import WrapperError
from FastOMA.zoo.wrappers.clusterers import cluster_sequences, find_mmseqs
from roothog_helpers import rec

SEQ_A = "MKVLAAGIVLLLAAGCSSHEWQRTYPLKDEFGHIKLMNPQRSTVWY"
SEQ_B = "MTTQAPTFLQQAQAAGLSPRVSVLVRRGDTRQALAAAREAGADYII"


def fake_mmseqs(tmp_path, script):
    """an executable 'mmseqs' shell script; the arguments of each call are appended to args.txt."""
    exe = tmp_path / "mmseqs"
    exe.write_text(f'#!/bin/sh\necho "$@" >> {tmp_path}/args.txt\n{script}\n')
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    return str(exe)


@pytest.mark.tools("mmseqs")
def test_clusters_records_with_uniprot_ids(tmp_path):
    # mmseqs reduces "sp|P1|NAME" to the accession; the clusters must still map back to the records
    recs = [rec("sp|P1|A_HUMAN||HUMAN||1", SEQ_A), rec("sp|P1|A_MOUSE||MOUSE||2", SEQ_A),
            rec("tr|Q9|B_HUMAN||HUMAN||3", SEQ_B)]
    clusters = cluster_sequences(recs, threads=1, tmp_dir=str(tmp_path))
    assert sorted(sorted(r.id for r in c) for c in clusters) == [
        ["sp|P1|A_HUMAN||HUMAN||1", "sp|P1|A_MOUSE||MOUSE||2"], ["tr|Q9|B_HUMAN||HUMAN||3"]]
    assert all(c[0] is r for c in clusters for r in recs if r.id == c[0].id)  # the records themselves
    assert os.listdir(tmp_path) == []  # temporary files removed


def test_no_records_does_not_run_mmseqs(tmp_path):
    assert cluster_sequences([], executable=fake_mmseqs(tmp_path, "exit 1")) == []
    assert not (tmp_path / "args.txt").exists()


def test_options_are_passed_only_when_set(tmp_path):
    exe = fake_mmseqs(tmp_path, 'printf "0\\t0\\n0\\t1\\n" > "$3_cluster.tsv"')
    clusters = cluster_sequences([rec("a", SEQ_A), rec("b", SEQ_A)], executable=exe, threads=3, coverage=0.5,
                                 split_memory_limit="1200M", tmp_dir=str(tmp_path))
    assert [[r.id for r in c] for c in clusters] == [["a", "b"]]
    args = (tmp_path / "args.txt").read_text().split()
    assert args[0] == "easy-cluster" and args[4:] == ["--threads", "3", "-c", "0.5", "--split-memory-limit", "1200M"]


@pytest.mark.parametrize("script, exit_code, message", [
    ("echo 'Input does not exist' >&2; exit 1", 1, "exited with code 1: Input does not exist"),
    ("kill -9 $$", 137, "terminated by signal 9"),   # e.g. killed for running out of memory
])
def test_failure_raises_with_exit_code(tmp_path, script, exit_code, message):
    with pytest.raises(WrapperError, match=message) as e:
        cluster_sequences([rec("a", SEQ_A)], executable=fake_mmseqs(tmp_path, script), tmp_dir=str(tmp_path))
    assert e.value.exit_code == exit_code


def test_find_mmseqs(tmp_path):
    assert find_mmseqs(str(tmp_path / "missing")) is None
    assert find_mmseqs(fake_mmseqs(tmp_path, "exit 1")) is None   # `mmseqs version` fails
    assert find_mmseqs(fake_mmseqs(tmp_path, "exit 0")) == str(tmp_path / "mmseqs")
