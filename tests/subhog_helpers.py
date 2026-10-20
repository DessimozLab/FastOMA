"""Builders for small HOG / species tree / MSA fixtures used by the subHOG-stage tests.

Protein ids follow the FastOMA convention ``<protein>||<species>||<numeric id>``.
"""
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree

from FastOMA._hog_class import HOG, Representative


def pid(name, species, num):
    return f"{name}||{species}||{num}"


def rec(seq_id, seq="ACDEFGHIKL"):
    return SeqRecord(Seq(seq), id=seq_id)


def msa(*records):
    """MSA from SeqRecords or (id, seq) tuples."""
    return MultipleSeqAlignment([r if isinstance(r, SeqRecord) else rec(*r) for r in records])


def species_tree(nwk="((A,B)AB,C)root;", with_size=True):
    """ete3 species tree with named internal nodes; optionally annotated with ``size``
    (number of leaves), as prepare_species_tree does for the scores."""
    t = Tree(nwk, format=1)
    if with_size:
        for n in t.traverse():
            n.add_feature("size", len(n))
    return t


def species_of(seq_id):
    return seq_id.split("||")[1]


def leaf_hog(seq_id, sptree, rhogid="r", seq="ACDEFGHIKL"):
    """Singleton HOG of one protein at its species' leaf node of ``sptree``."""
    return HOG(rec(seq_id, seq), sptree & species_of(seq_id), rhogid)


def merged_hog(subhogs, taxnode, rep_ids=None, rhogid="r", conf=None, aln=None):
    """HOG at ``taxnode`` merging ``subhogs``.

    The MSA defaults to the stacked rows of the subhogs' MSAs; representatives default to
    all representatives of the subhogs (or the given ``rep_ids``)."""
    if aln is None:
        aln = MultipleSeqAlignment([r for h in subhogs for r in h.get_msa()])
    by_id = {r.id: r for r in aln}
    if rep_ids is None:
        rep_ids = [r.get_id() for h in subhogs for r in h.get_representatives()]
    reps = [Representative(by_id[i]) for i in rep_ids]
    return HOG(list(subhogs), taxnode, rhogid, msa=aln, representatives=reps, conf_infer_subhhogs=conf)


def leaf_members(hog):
    """Members collected from the leaves of the subhog hierarchy (should equal hog.get_members())."""
    if not hog._subhogs:
        return set(hog.get_members())
    return set().union(*(leaf_members(sh) for sh in hog._subhogs))


def members_of(hogs):
    return set().union(set(), *(h.get_members() for h in hogs))
