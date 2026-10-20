"""Helpers to compare FastOMA outputs by meaning rather than by bytes.

HOG/OG numbering, TSV row order and the order of elements in the orthoXML can change
between runs (task completion order, set iteration order) without any change in the
inferred orthology. These helpers reduce outputs to order- and id-independent forms.
"""
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

from lxml import etree

ORTHOXML_NS = "http://orthoXML.org/2011/"
NS = {"o": ORTHOXML_NS}
ORTHOXML_XSD = Path(__file__).parent / "data" / "orthoxml-0.5.xsd"


def _open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def read_group_partition(tsv_path, group_col=0, member_col=1, header=True):
    """Partition of members into groups from a 'group<TAB>member' TSV, ignoring group names and row order.

    Works for RootHOGs.tsv and OrthologousGroups.tsv.
    Returns a frozenset of frozensets.
    """
    groups = defaultdict(set)
    with _open_text(tsv_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        if header:
            next(reader, None)
        for row in reader:
            if row:
                groups[row[group_col]].add(row[member_col])
    return frozenset(frozenset(m) for m in groups.values())


def read_pairs(tsv_path):
    """Set of unordered pairs from a two-column TSV (e.g. orthologs.tsv.gz)."""
    with _open_text(tsv_path) as fh:
        return {frozenset(line.rstrip("\n").split("\t")[:2]) for line in fh if line.strip()}


def _parse(orthoxml):
    if isinstance(orthoxml, etree._ElementTree):
        return orthoxml
    if isinstance(orthoxml, etree._Element):
        return orthoxml.getroottree()
    return etree.parse(str(orthoxml))


def gene_id_to_protid(orthoxml):
    tree = _parse(orthoxml)
    return {g.get("id"): g.get("protId") for g in tree.iterfind(".//o:gene", NS)}


def geneRef_counts(orthoxml, use_protid=True):
    """Counter of how often each gene is referenced from groups (a protein must appear at most once)."""
    tree = _parse(orthoxml)
    id2prot = gene_id_to_protid(tree) if use_protid else {}
    return Counter(id2prot.get(r.get("id"), r.get("id")) for r in tree.iterfind(".//o:groups//o:geneRef", NS))


def canonical_group(elem, id2prot):
    """Order- and id-independent nested representation of an orthologGroup/paralogGroup.

    Keeps: group type, TaxRange property, and the nesting of genes (as protIds).
    Drops: group ids, scores, other properties, taxonIds.
    """
    tag = etree.QName(elem).localname
    if tag == "geneRef":
        return ("gene", id2prot.get(elem.get("id"), elem.get("id")))
    taxrange = None
    children = []
    for child in elem:
        ctag = etree.QName(child).localname
        if ctag == "property" and child.get("name") == "TaxRange":
            taxrange = child.get("value")
        elif ctag in ("orthologGroup", "paralogGroup", "geneRef"):
            children.append(canonical_group(child, id2prot))
    return (tag, taxrange, tuple(sorted(children, key=repr)))


def canonical_orthoxml(orthoxml):
    """Canonical form of a whole orthoXML file: sorted tuple of canonical top-level groups."""
    tree = _parse(orthoxml)
    id2prot = gene_id_to_protid(tree)
    groups = tree.find(".//o:groups", NS)
    if groups is None:
        return ()
    return tuple(sorted((canonical_group(g, id2prot) for g in groups), key=repr))


def orthoxml_element_counts(orthoxml):
    tree = _parse(orthoxml)
    return {
        "orthologGroup": len(tree.findall(".//o:orthologGroup", NS)),
        "paralogGroup": len(tree.findall(".//o:paralogGroup", NS)),
        "gene": len(tree.findall(".//o:gene", NS)),
        "geneRef": len(tree.findall(".//o:geneRef", NS)),
    }


_schema = None


def orthoxml_schema():
    global _schema
    if _schema is None:
        _schema = etree.XMLSchema(etree.parse(str(ORTHOXML_XSD)))
    return _schema


def orthoxml_validation_errors(orthoxml):
    """List of schema violations (empty if valid against OrthoXML 0.5)."""
    schema = orthoxml_schema()
    tree = _parse(orthoxml)
    if schema.validate(tree):
        return []
    return [f"line {e.line}: {e.message}" for e in schema.error_log]
