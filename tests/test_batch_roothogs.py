"""Tests for splitting rootHOG FASTA files into big/rest batches (batch_roothogs.py)."""
import pytest

from FastOMA import batch_roothogs as B


def make_files(folder, sizes):
    for name, size in sizes.items():
        path = folder / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"A" * size)


def batches(folder):
    return {d.name: sorted(f.name for f in d.iterdir()) for d in sorted(folder.iterdir())}


def run(tmp_path, tag="out"):
    big, rest = tmp_path / tag / "big", tmp_path / tag / "rest"
    B.folder_1h_rhog(tmp_path / "in", big, rest)
    return batches(big), batches(rest)


def test_big_files_get_own_batch_in_descending_size(tmp_path):
    make_files(tmp_path / "in", {"big1.fa": 500_000, "big2.fa": 450_000, "huge.fa": 1_300_000,
                                 "edge.fa": 400_000, "s1.fa": 300_000, "sub/s2.fa": 200, "x.txt": 10})
    big, rest = run(tmp_path)
    assert big == {"0": ["huge.fa"], "1": ["big1.fa"], "2": ["big2.fa"]}
    # threshold is strict; subfolders are searched; non-.fa files ignored
    assert rest == {"0": ["edge.fa", "s1.fa", "s2.fa"]}


def test_rest_batches_close_after_exceeding_max_size(tmp_path):
    make_files(tmp_path / "in", {f"s{i}.fa": 390_000 for i in range(7)})
    _, rest = run(tmp_path)
    # the file that pushes the batch above 2 MB stays in that batch
    assert [len(v) for v in rest.values()] == [6, 1]


def test_empty_input(tmp_path):
    (tmp_path / "in").mkdir()
    assert run(tmp_path) == ({}, {})


def test_rerun_into_existing_output_fails(tmp_path):
    # characterises current behaviour: batch folders are created with mkdir() without exist_ok
    make_files(tmp_path / "in", {"s.fa": 10})
    run(tmp_path)
    with pytest.raises(FileExistsError):
        run(tmp_path)


def test_same_basename_in_subfolders_rejected(tmp_path):
    # the file name is the rootHOG id downstream and batches are flat, so two files with the same
    # name would silently overwrite each other
    make_files(tmp_path / "in", {"a/HOG_1.fa": 10, "b/HOG_1.fa": 20})
    with pytest.raises(ValueError, match="HOG_1.fa"):
        run(tmp_path)
