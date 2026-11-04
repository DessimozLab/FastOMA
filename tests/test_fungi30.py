"""fastoma-infer-subhogs on rootHOGs selected from the fungi-30 run (see tests/fungi30.py).

golden families: the output must equal the stored snapshot (only in the snapshot's tool environment).
drops families:  known to drop proteins; the invariant test is a strict xfail until that is fixed.
All families:    every input protein must end up in a final HOG.

Set FASTOMA_REQUIRE_GOLDEN=1 (as CI does) to fail instead of skip when the local mafft/FastTree differ
from the snapshot environment.
"""
import difflib
import os

import pytest

import fungi30 as F

pytestmark = pytest.mark.tools("mafft", "FastTree")

FAMILIES = F.families()


def _params(kind=None, xfail_drops=False):
    params = []
    for fam in FAMILIES:
        if kind and fam["kind"] != kind:
            continue
        marks = [pytest.mark.slow] if fam["tier"] == "slow" else []
        if xfail_drops and fam["kind"] == "drops":
            marks.append(pytest.mark.xfail(strict=True, reason="BUG: fastoma-infer-subhogs drops proteins of this "
                                                               "rootHOG: " + fam["reason"]))
        params.append(pytest.param(fam["rhog"], id=fam["rhog"], marks=marks))
    return params


_runs = {}


@pytest.fixture(scope="module")
def run_of(tmp_path_factory):
    """work dir of the (cached) fastoma-infer-subhogs run that contains the given rootHOG."""
    by_rhog = {r: key for key, rhogids in F.groups().items() for r in rhogids}

    def get(rhogid):
        key = by_rhog[rhogid]
        if key not in _runs:
            workdir = tmp_path_factory.mktemp("fungi30_{}_{}".format(*key))
            F.prepare_input(F.groups()[key], workdir / "in")
            F.run_infer_subhogs(workdir / "in", workdir, parallel=key[1] == "par")
            _runs[key] = workdir
        return _runs[key]
    yield get
    _runs.clear()


@pytest.fixture(scope="module")
def golden_env():
    expected, actual = F.golden_environment(), F.tool_versions()
    expected_tools = {k: v for k, v in expected.items() if k in actual}
    if expected_tools != actual:
        msg = (f"tool versions {actual} differ from the golden snapshot environment {expected}; "
               "regenerate with `python tests/fungi30.py regenerate` if this is intended")
        if os.environ.get("FASTOMA_REQUIRE_GOLDEN"):
            pytest.fail(msg)
        pytest.skip(msg)


@pytest.mark.parametrize("rhogid", _params(xfail_drops=True))
def test_no_protein_dropped(run_of, rhogid):
    members = F.final_members(run_of(rhogid), rhogid)
    assert sorted(members) == sorted(F.input_records(rhogid))


@pytest.mark.parametrize("rhogid", _params(kind="golden"))
def test_matches_golden_snapshot(golden_env, run_of, rhogid):
    actual = F.golden_text(F.output_groups(run_of(rhogid), rhogid), F.input_records(rhogid))
    expected = F.golden_path(rhogid).read_text()
    if actual != expected:
        diff = "".join(difflib.unified_diff(expected.splitlines(True), actual.splitlines(True),
                                            "golden", "actual", n=2))
        pytest.fail(f"output of {rhogid} differs from tests/data/fungi30/golden/{rhogid}.txt:\n{diff[:6000]}")


def test_every_family_has_input_and_golden_snapshot():
    for fam in FAMILIES:
        assert (F.DATA / "rhogs" / f"HOG_{fam['rhog']}.fa.gz").exists(), fam["rhog"]
        assert F.golden_path(fam["rhog"]).exists() == (fam["kind"] == "golden"), fam["rhog"]


def test_parse_tool_versions():
    v = F.parse_tool_versions("v7.505 (2022/Apr/10)",
                              "Detailed usage for FastTree 2.1.11 Double precision (No SSE3), OpenMP (10 threads):")
    assert v == {"mafft": "v7.505 (2022/Apr/10)", "FastTree": "2.1.11 Double precision (No SSE3), OpenMP"}
