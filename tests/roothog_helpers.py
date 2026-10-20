"""Small builders for rootHOG-stage tests: hogmap hits, hogmap/FASTA files and rootHOG partitions."""
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from FastOMA._utils_roothog import HOGMapData

HOGMAP_HEADER = ("qseqid", "hogid", "hoglevel", "family_p", "family_count", "family_normcount",
                 "subfamily_score", "subfamily_count", "qseqlen", "subfamily_medianseqlen", "qseq_overlap")


def hit(hogid, score, seqlen=100, medlen=100):
    """One omamer hit as stored by parse_hogmap_omamer (all fields are strings)."""
    return HOGMapData(hogid, str(score), str(seqlen), str(medlen))


def rec(prot_id, seq="MKVLAAGIVLLLAAGCSSHE"):
    return SeqRecord(Seq(seq), id=prot_id, description="")


def write_fasta(path, records):
    """records: iterable of (id, sequence)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{pid}\n{seq}\n" for pid, seq in records))
    return path


def write_hogmap(path, rows, comment=True):
    """rows: iterable of (qseqid, hogid, family_p); hogid 'N/A' marks an unmapped protein (omamer v2)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["!omamer-version: 2.0.0"] if comment else []
    lines.append("\t".join(HOGMAP_HEADER))
    for qseqid, hogid, score in rows:
        if hogid == "N/A":
            vals = [qseqid, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "100", "N/A", "N/A"]
        else:
            vals = [qseqid, hogid, "LUCA", str(score), "10", "1.0", "1.0", "10", "100", "100", "1.0"]
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n")
    return path


def fasta_ids(path):
    """Protein ids (header up to the first space and '||') of a FASTA file."""
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].split("||")[0])
    return ids


def read_rhog_partition(folder, pattern="HOG_*.fa"):
    """Partition of protein ids over the rootHOG FASTA files in folder, ignoring file names and order."""
    return frozenset(frozenset(fasta_ids(f)) for f in Path(folder).glob(pattern))


def partition(rhogs_prots):
    """Order-independent form of a rhogs_prots dict ({rhogid: [(species, prot), ...]})."""
    return frozenset(frozenset(v) for v in rhogs_prots.values())
