"""Tests for the fasta header -> protein id transformers (FastOMA/transformer.py)."""
import pytest

from FastOMA.transformer import (
    ExtractUniProtAccessionFastaHeaderTransformer,
    NoOpFastaHeaderTransformer,
    header_transformer,
)


@pytest.mark.parametrize("name, cls", [("noop", NoOpFastaHeaderTransformer), ("NOOP", NoOpFastaHeaderTransformer),
                                       ("UniProt", ExtractUniProtAccessionFastaHeaderTransformer),
                                       ("uniprot", ExtractUniProtAccessionFastaHeaderTransformer)])
def test_factory(name, cls):
    assert isinstance(header_transformer(name), cls)


def test_factory_unknown_name_returns_none():
    # characterises current behaviour: argparse choices guard the CLI, the factory itself does not
    assert header_transformer("foo") is None


def test_noop():
    assert header_transformer("noop").transform("sp|P68250|1433B_BOVIN") == "sp|P68250|1433B_BOVIN"


@pytest.mark.parametrize("header, expected", [
    ("sp|P68250|1433B_BOVIN", "P68250"),
    ("tr|A0A0N7KCI6|A0A0N7KCI6_ORYSJ", "A0A0N7KCI6"),
    ("sp|Q9XYZ1|", "Q9XYZ1"),
    ("ts|P68250|x", "P68250"),  # prefix is any two of [sptr]
])
def test_uniprot_extracts_accession(header, expected):
    assert header_transformer("UniProt").transform(header) == expected


@pytest.mark.parametrize("header", [
    "sp|P68250-2|1433B_BOVIN",  # isoform accessions are not recognised
    "P68250",
    "sp|p68250|x",
    "ENSP00000354587",
])
def test_uniprot_keeps_unrecognised_headers(header, caplog):
    assert header_transformer("UniProt").transform(header) == header
    assert "cannot extract uniprot accession" in caplog.text
