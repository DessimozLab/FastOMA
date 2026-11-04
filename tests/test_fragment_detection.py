import numpy as np
import pytest

from FastOMA._utils_frag_SO_detection import read_msa, split_candidates, compute_identity, find_prot_dubious_msa
from subhog_helpers import pid, msa

H1, H2, H3 = pid("h1", "HUMAN", 1), pid("h2", "HUMAN", 2), pid("h3", "HUMAN", 3)
M1 = pid("m1", "MOUSE", 4)


def test_read_msa_ids_and_coordinates():
    ids, seqs, species, coords = read_msa(msa((H1, "aa--"), (M1, "-AA-")))
    assert list(ids) == [H1, M1]
    assert seqs.shape == (2, 4)
    # sequences are upper-cased
    assert seqs[0].tobytes() == b"AA--"
    assert coords.tolist() == [[0, 1], [1, 2]]


def test_read_msa_all_gap_row():
    _, _, _, coords = read_msa(msa((H1, "----"), (M1, "AAAA")))
    assert coords.tolist() == [[0, 0], [0, 3]]


@pytest.mark.xfail(strict=True, reason="BUG: read_msa strips the last character of species names "
                                       "(_utils_frag_SO_detection.py:36)")
def test_read_msa_species_names():
    _, _, species, _ = read_msa(msa((pid("a", "AQUAE", 1), "AA"), (pid("b", "AQUAF", 2), "AA")))
    assert list(species) == ["AQUAE", "AQUAF"]


def test_compute_identity():
    s1 = np.frombuffer(b"AC-D", dtype="S1")
    s2 = np.frombuffer(b"AGWD", dtype="S1")
    assert compute_identity(s1, s2) == (pytest.approx(2 / 3), pytest.approx(2 / 4))
    gap = np.frombuffer(b"--", dtype="S1")
    assert compute_identity(gap, gap) == (1, 1)


def test_split_candidates_non_overlapping_fragments():
    cands = split_candidates(read_msa(msa((H2, "-----AAAA"), (H1, "AAAA-----"), (M1, "AAAAAAAAA"))), 0.15)
    # ordered by start coordinate
    assert cands == [((H1, H2), [(1, 1)])]


def test_split_candidates_small_overlap_within_margin():
    cands = split_candidates(read_msa(msa((H1, "AAAAAAAAAAW----------"), (H2, "---------AWAAAAAAAAAA"))), 0.15)
    assert len(cands) == 1
    assert cands[0][0] == (H1, H2)


def test_split_candidates_large_overlap():
    assert split_candidates(read_msa(msa((H1, "AAAAAA----"), (H2, "--AAAAAA--"))), 0.15) == []


def test_split_candidates_third_full_length_copy_blocks_merge():
    aln = msa((H1, "AAAA------"), (H2, "------AAAA"), (H3, "AAAAAAAAAA"))
    assert split_candidates(read_msa(aln), 0.15) == []


def test_split_candidates_identity_is_not_used():
    # characterises current behaviour: identity of the overlap is computed but does not prevent a merge
    cands = split_candidates(read_msa(msa((H1, "AAAAAA----"), (H2, "----WWWWWW"))), 0.5)
    assert cands == [((H1, H2), [(0.0, 0.0)])]


def test_split_candidates_margin_must_be_fraction():
    with pytest.raises(AssertionError):
        split_candidates(read_msa(msa((H1, "AA"), (H2, "AA"))), 2)


def test_find_prot_dubious_msa(subhog_conf):
    aln = msa((H1, "AAAA------"), (H2, "------AAAA"), (M1, "AAAAAAAAAA"))
    for r in aln:
        r.name = r.id  # as for records parsed by SeqIO (mafft output)
    prots, seqs = find_prot_dubious_msa(aln, subhog_conf())
    assert prots == [[H1, H2]]
    assert [[r.id for r in s] for s in seqs] == [[H1, H2]]


def test_find_prot_dubious_msa_matches_sequences_by_name(subhog_conf):
    # characterises current behaviour: sequences are looked up by record.name, not record.id
    # (_utils_frag_SO_detection.py:142)
    aln = msa((H1, "AAAA------"), (H2, "------AAAA"))
    assert find_prot_dubious_msa(aln, subhog_conf()) == ([[H1, H2]], [[]])


def test_find_prot_dubious_msa_none(subhog_conf):
    assert find_prot_dubious_msa(msa((H1, "AAAA"), (M1, "AAAA")), subhog_conf()) == ([], [])


@pytest.mark.xfail(strict=True, reason="BUG: species names differing only in the last character are treated as one "
                                       "species (_utils_frag_SO_detection.py:36)")
def test_find_prot_dubious_msa_distinct_species(subhog_conf):
    aln = msa((pid("a", "HUMAN1", 1), "AAAA------"), (pid("b", "HUMAN2", 2), "------AAAA"))
    assert find_prot_dubious_msa(aln, subhog_conf())[0] == []
