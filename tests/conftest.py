"""Shared pytest fixtures for the FastOMA unit tests.

Markers (registered in pyproject.toml):
  tools -- needs external binaries (mafft, FastTree, mmseqs, omamer); skipped if missing
  slow  -- takes more than a few seconds; deselect with ``-m "not slow"``
"""
import os
import shutil
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).parent
REPO_DIR = TESTS_DIR.parent
DATA_DIR = TESTS_DIR / "data"
TESTDATA_DIR = REPO_DIR / "testdata"

# make tests/fastoma_testutils.py importable from every test module
sys.path.insert(0, str(TESTS_DIR))

# LevelHOGProcessor.process() raises instead of repairing a lost or duplicated protein, so the tests
# see bugs the safety net would hide (inherited by the infer-subhogs subprocesses of the stage tests)
os.environ.setdefault("FASTOMA_STRICT_INVARIANTS", "1")

REQUIRED_TOOLS = {"mafft": "mafft", "FastTree": "FastTree", "mmseqs": "mmseqs", "omamer": "omamer"}


def pytest_collection_modifyitems(config, items):
    """Skip tests marked ``tools`` when the binaries they name are missing.

    Usage: ``@pytest.mark.tools("mafft", "FastTree")``; without arguments all tools are required.
    """
    for item in items:
        marker = item.get_closest_marker("tools")
        if marker is None:
            continue
        needed = marker.args or tuple(REQUIRED_TOOLS)
        missing = [t for t in needed if shutil.which(REQUIRED_TOOLS.get(t, t)) is None]
        if missing:
            item.add_marker(pytest.mark.skip(reason="missing external tool(s): " + ", ".join(missing)))


class _ParsedArgs(Exception):
    def __init__(self, conf):
        self.conf = conf


def _capture_cli_conf(monkeypatch, module, entry_point, argv):
    """Run a fastoma CLI entry point only up to argument parsing and return its Namespace.

    All entry points call ``setup_logging(conf.v, ...)`` right after ``parse_args()``; replacing it
    with a function that raises lets us get the real argparse defaults without running anything.
    """
    def stop(*args, **kwargs):
        raise _ParsedArgs(stop.conf)

    import argparse
    orig_parse_args = argparse.ArgumentParser.parse_args

    def parse_args(self, *args, **kwargs):
        stop.conf = orig_parse_args(self, *args, **kwargs)
        return stop.conf

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", parse_args)
    monkeypatch.setattr(module, "setup_logging", stop)
    monkeypatch.setattr(sys, "argv", [entry_point] + list(argv))
    try:
        getattr(module, entry_point.replace("-", "_"))()
    except _ParsedArgs as e:
        return e.conf
    raise AssertionError(f"{entry_point} did not call setup_logging after parsing arguments")


@pytest.fixture
def subhog_conf(monkeypatch):
    """Factory for the configuration Namespace of fastoma-infer-subhogs with the real CLI defaults.

    Keyword arguments override single fields, e.g. ``subhog_conf(number_of_samples_per_hog=3)``.
    """
    from FastOMA import infer_subhogs

    with monkeypatch.context() as m:
        base = _capture_cli_conf(m, infer_subhogs, "fastoma-infer-subhogs",
                                 ["--input-rhog-folder", "rhogs", "--species-tree", "species_tree.nwk"])

    def make(**overrides):
        conf = type(base)(**vars(base))
        for key, value in overrides.items():
            if not hasattr(conf, key):
                raise AttributeError(f"infer-subhogs has no option {key!r}")
            setattr(conf, key, value)
        return conf
    return make


@pytest.fixture
def roothog_conf(monkeypatch):
    """Factory for the configuration Namespace of fastoma-infer-roothogs with the real CLI defaults."""
    from FastOMA import infer_roothogs

    with monkeypatch.context() as m:
        base = _capture_cli_conf(m, infer_roothogs, "fastoma-infer-roothogs",
                                 ["--proteomes", "proteome", "--out-rhog-folder", "omamer_rhogs"])

    def make(**overrides):
        conf = type(base)(**vars(base))
        for key, value in overrides.items():
            if not hasattr(conf, key):
                raise AttributeError(f"infer-roothogs has no option {key!r}")
            setattr(conf, key, value)
        return conf
    return make


@pytest.fixture(autouse=True)
def reset_hogid_counter():
    """HOG ids come from a class-level counter; reset it so every test sees the same ids."""
    from FastOMA._hog_class import HOG
    saved = HOG._hogid_iter
    HOG._hogid_iter = 10000
    yield
    HOG._hogid_iter = saved


@pytest.fixture
def in_tmp_dir(tmp_path, monkeypatch):
    """Run the test inside a fresh temporary directory (many functions write into the CWD)."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def fixed_hashseed(monkeypatch):
    """Environment for subprocesses with a fixed PYTHONHASHSEED.

    Several FastOMA steps iterate over sets, so their output depends on the hash seed
    (see test plan). Subprocess-based tests use this to stay reproducible.
    """
    env = dict(os.environ, PYTHONHASHSEED="0")
    return env
