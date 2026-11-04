"""What the report notebook (FastOMA/fastoma_notebook_stat.ipynb) expects from the pipeline outputs.

The notebook is not executed here (papermill is an optional extra); these tests pin the input
formats it parses via FastOMA/zoo/hog/extract_hog_info.py and the names it reads.
"""
import io
import json
import re

import pytest

from FastOMA import collect_subhogs as cs
from FastOMA.transformer import header_transformer
from FastOMA.zoo.hog import extract_hog_info
from conftest import REPO_DIR, TESTDATA_DIR
import output_fixtures as fx

NOTEBOOK = REPO_DIR / "FastOMA" / "fastoma_notebook_stat.ipynb"


def _code_cells():
    nb = json.loads(NOTEBOOK.read_text())
    return [c for c in nb["cells"] if c["cell_type"] == "code"]


@pytest.fixture
def fixture_orthoxml(tmp_path):
    inputs = fx.build_collect_inputs(tmp_path / "in", isoforms=True)
    out = tmp_path / "FastOMA_HOGs.orthoxml"
    cs.write_hog_orthoxml(str(inputs["pickle_folder"]), str(out), inputs["gene_id_pickle"],
                          header_transformer("noop"), str(inputs["species_tree"]))
    return out


def _parse(path_or_text):
    analyser = extract_hog_info.SpeciesAnalyser()
    fh = io.StringIO(path_or_text) if isinstance(path_or_text, str) else open(path_or_text)
    with fh:
        records = list(extract_hog_info.parse_orthoxml(fh, analyser))
    return records, analyser


class TestHogInfoParsing:
    def test_records(self, fixture_orthoxml):
        records, _ = _parse(fixture_orthoxml)
        assert len(records) == 10
        assert all(set(r) == {"id", "level", "CompletenessScore", "nr_members", "is_roothog"} for r in records)
        roots = {r["id"]: r for r in records if r["is_roothog"]}
        assert {(r["level"], r["nr_members"]) for r in roots.values()} == {("root", 11), ("AB", 5)}
        assert all(isinstance(r["CompletenessScore"], float) for r in records)

    def test_species_summary_counts_main_isoforms(self, fixture_orthoxml):
        _, analyser = _parse(fixture_orthoxml)
        summary = {s["species"]: s for s in analyser.summary()}
        # a7 is a minor isoform of a1; A has 7 proteins, 6 genes, a7 is in a group
        assert summary["A"] == {"species": "A", "genes": 6, "not_in_group": 0, "proteins": 7, "minor_splice": 1}
        assert summary["B"] == {"species": "B", "genes": 7, "not_in_group": 1, "proteins": 7, "minor_splice": 0}
        assert summary["D"]["not_in_group"] == 6

    def test_expected_output_parses(self):
        records, analyser = _parse(TESTDATA_DIR / "expected_output" / "FastOMA_HOGs.orthoxml")
        assert sum(r["is_roothog"] for r in records) == 12
        assert {s["species"] for s in analyser.summary()} == {"AQUAE", "CHLTR", "MYCGE"}

    def test_requires_version_0_5(self, fixture_orthoxml):
        text = fixture_orthoxml.read_text().replace('version="0.5"', 'version="0.4"', 1)
        with pytest.raises(RuntimeError, match="0.5"):
            _parse(text)

    def test_without_scores_there_is_no_completeness_column(self, fixture_orthoxml):
        text = re.sub(r"\s*<score [^>]*/>", "", fixture_orthoxml.read_text())
        records, _ = _parse(text)
        assert all("CompletenessScore" not in r for r in records)


class TestNotebook:
    def test_parameters_cell(self):
        # FastOMA.nf (fastoma_report) passes -p output_folder, min_sequence_length, proteome_folder
        params = [c for c in _code_cells() if "parameters" in c.get("metadata", {}).get("tags", [])]
        assert len(params) == 1
        source = "".join(params[0]["source"])
        for name in ("output_folder", "min_sequence_length", "proteome_folder"):
            assert re.search(rf"^{name}\s*=", source, re.M), name

    def test_reads_these_outputs(self):
        source = "\n".join("".join(c["source"]) for c in _code_cells())
        assert '"FastOMA_HOGs.orthoxml"' in source
        assert '"species_tree_checked.nwk"' in source
        assert '"phylostratigraphy.html"' in source

    def test_completeness_plots_are_guarded(self):
        for cell in _code_cells():
            source = "".join(cell["source"])
            if "x=\"CompletenessScore\"" in source or "y='CompletenessScore'" in source:
                assert "CompletenessScore\" in" in source or "'CompletenessScore' in" in source, source
