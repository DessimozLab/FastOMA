import pytest

from FastOMA._utils_subhog import MSAFilter, MSAFilterElbow, MSAFilterTrimAL
from subhog_helpers import pid, msa

A1, B1, C1, D1 = pid("a1", "A", 1), pid("b1", "B", 2), pid("c1", "C", 3), pid("d1", "D", 4)


@pytest.fixture
def flt(subhog_conf):
    return MSAFilter(None, subhog_conf(min_col_trim=0))


def seqs(aln):
    return [str(r.seq) for r in aln]


def test_defaults_from_conf(subhog_conf):
    f = MSAFilter(None, subhog_conf())
    assert (f.gap_ratio_row, f.gap_ratio_col, f.min_col_trim) == (0.3, 0.5, 50)


def test_gap_ratios_count_all_gap_symbols(flt):
    assert list(flt._get_gap_ratios(msa((A1, "-.~?XA"), (B1, "AAAAAA")))) == [0.5, 0.5, 0.5, 0.5, 0, 0]


def test_column_threshold_is_exclusive(flt):
    # gap ratio exactly 0.5 is removed
    assert seqs(flt.msa_filter_col(msa((A1, "A-A"), (B1, "AAA")))) == ["AA", "AA"]


def test_column_filter_keeps_msa_without_gappy_columns(flt):
    aln = msa((A1, "ACD"), (B1, "AC-"), (C1, "ACD"))
    assert flt.msa_filter_col(aln) is aln


def test_row_threshold_is_min_nongap_fraction(flt):
    # gap_ratio_row is the minimal fraction of non-gap characters; exactly 0.3 is kept
    aln = msa((A1, "AAA-------"), (B1, "AA--------"))
    assert [r.id for r in flt.msa_filter_row(aln)] == [A1]


def test_filter_msa_reports_removed_rows(flt):
    aln = msa((A1, "ACDEFGHIKL"), (B1, "ACDEFGHIKL"), (C1, "A---------"), (D1, "ACDEFGHIKL"))
    filtered, removed = flt.filter_msa(aln)
    assert [r.id for r in filtered] == [A1, B1, D1]
    assert removed == {C1}


def test_filter_msa_removes_empty_columns_after_row_filter(subhog_conf):
    f = MSAFilter(None, subhog_conf(min_col_trim=0, gap_ratio_col=0.7, gap_ratio_row=0.9))
    E1 = pid("e1", "E", 5)
    aln = msa((A1, "ACDEFGHIKL-"), (B1, "ACDEFGHIKL-"), (C1, "A---------W"), (D1, "A---------W"), (E1, "A---------W"))
    # the last column survives the column filter, but only the removed rows have residues there
    filtered, removed = f.filter_msa(aln)
    assert seqs(filtered) == ["ACDEFGHIKL", "ACDEFGHIKL"]
    assert removed == {C1, D1, E1}


def test_filter_msa_skips_short_alignments(subhog_conf):
    f = MSAFilter(None, subhog_conf(min_col_trim=5))
    aln = msa((A1, "A---"), (B1, "AAAA"))
    filtered, removed = f.filter_msa(aln)
    assert filtered is aln and removed == set()


def test_filter_msa_single_sequence(flt):
    filtered, removed = flt.filter_msa(msa((A1, "AA--AA")))
    assert seqs(filtered) == ["AAAA"] and removed == set()


def test_filter_msa_without_gaps_is_unchanged(flt):
    aln = msa((A1, "ACDE"), (B1, "ACDE"))
    filtered, removed = flt.filter_msa(aln)
    assert seqs(filtered) == ["ACDE", "ACDE"] and removed == set()


def test_column_filter_with_single_column_left_drops_rows(flt):
    # characterises current behaviour: records shorter than 2 columns are dropped (_utils_subhog.py:485)
    assert len(flt.msa_filter_col(msa((A1, "AA-"), (B1, "-AA")))) == 0


@pytest.mark.xfail(strict=True, reason="BUG: rows dropped by the column filter are not reported as removed "
                                       "(_utils_subhog.py:485 and :456)")
def test_filter_msa_all_columns_gappy_reports_rows(flt):
    aln = msa((A1, "AA--"), (B1, "--AA"))
    filtered, removed = flt.filter_msa(aln)
    assert {r.id for r in filtered} | removed == {A1, B1}


def test_filter_msa_empty_msa(flt):
    # characterises current behaviour: an empty MSA is not handled
    with pytest.raises(IndexError):
        flt.filter_msa(msa())


# MSAFilterElbow

def test_elbow_without_gaps_keeps_everything(subhog_conf):
    f = MSAFilterElbow(None, subhog_conf(min_col_trim=0))
    assert seqs(f.msa_filter_col(msa((A1, "ACDE"), (B1, "ACDE")))) == ["ACDE", "ACDE"]


def test_elbow_never_uses_threshold_below_gap_ratio_col(subhog_conf):
    f = MSAFilterElbow(None, subhog_conf(min_col_trim=0, gap_ratio_col=0.5))
    aln = msa((A1, "ACDEFGHI"), (B1, "ACDEFG--"), (C1, "ACDEF---"), (D1, "ACDE----"))
    # gap ratios: 0,0,0,0,0.25,0.5,0.75,0.75 -> at least everything below 0.5 is kept
    assert len(f.msa_filter_col(aln)[0]) >= 5


def test_elbow_removes_all_gap_columns(subhog_conf):
    f = MSAFilterElbow(None, subhog_conf(min_col_trim=0))
    aln = msa((A1, "AAA-"), (B1, "AAA-"), (C1, "AAA-"))
    assert seqs(f.msa_filter_col(aln)) == ["AAA"] * 3


# MSAFilterTrimAL

@pytest.mark.tools("trimal")
def test_trimal_filter_reports_removed_sequences(subhog_conf):
    f = MSAFilterTrimAL(None, subhog_conf())
    aln = msa(*((pid(f"p{i}", "S", i), "MKVLAAGIVGLLLAAPLAQAKTEEAPKAEEGHHHH") for i in range(4)))
    filtered, removed = f.filter_msa(aln)
    assert {r.id for r in filtered} | removed == {r.id for r in aln}
