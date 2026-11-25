"""Clustering of protein sequences with MMseqs2 (``mmseqs easy-cluster``).

Unlike the aligner and tree builder wrappers, this is a plain function: FastOMA calls exactly one
mmseqs command, so there is nothing to swap behind a common interface.

    clusters = cluster_sequences(records, threads=4, coverage=0.8)

mmseqs reads only the first word of a FASTA header as id, and for UniProt-style ids
(``sp|P12345|NAME``) it even reduces that to the accession. The records are therefore passed to
mmseqs under their index, and the clusters are mapped back to the records, so any record id works.
"""
import logging
import os
import shutil
import subprocess
import tempfile

from .. import WrapperError

logger = logging.getLogger(__name__)

DEFAULT_EXECUTABLE = "mmseqs"


def find_mmseqs(executable=DEFAULT_EXECUTABLE):
    """the path of a working mmseqs binary, or None."""
    path = shutil.which(executable)
    if path is None:
        return None
    try:
        subprocess.run([path, "version"], check=True, capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError) as e:
        logger.error("mmseqs (%s) is not working: %s", path, e)
        return None
    return path


def cluster_sequences(records, executable=DEFAULT_EXECUTABLE, threads=None, coverage=None, sensitivity=None,
                      min_seq_id=None, cov_mode=None, tmp_dir=None):
    """cluster the sequences of `records` (Bio.SeqRecord) with ``mmseqs easy-cluster``.

    Returns the clusters as lists of records, in the order mmseqs reports them, each starting with
    its representative; singletons included. Options left at None are not passed, so mmseqs uses its
    defaults (threads: all cores). The temporary files go to a new folder in `tmp_dir` (default: the
    current directory, as /tmp is often small on clusters) and are removed afterwards.

    Raises WrapperError if mmseqs fails; its exit_code is that of mmseqs, or 128 + signal number if
    mmseqs was killed (e.g. out of memory), so the caller can exit with it and Nextflow can retry
    the task with more resources."""
    records = list(records)
    if not records:
        return []
    path = shutil.which(executable)
    if path is None:
        raise WrapperError(f"{executable} not found on the PATH")
    options = {"--threads": threads, "-c": coverage, "-s": sensitivity, "--min-seq-id": min_seq_id,
               "--cov-mode": cov_mode}
    with tempfile.TemporaryDirectory(prefix="mmseqs_", dir=tmp_dir or os.getcwd()) as work:
        fasta = os.path.join(work, "input.fa")
        with open(fasta, "wt") as fh:
            for idx, rec in enumerate(records):
                fh.write(f">{idx}\n{rec.seq}\n")
        prefix = os.path.join(work, "result")
        cmd = [path, "easy-cluster", fasta, prefix, os.path.join(work, "tmp")]
        cmd += [str(arg) for opt, value in options.items() if value is not None for arg in (opt, value)]
        logger.info("running %s", " ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            if res.returncode < 0:
                termination, exit_code = f"was terminated by signal {-res.returncode}", 128 - res.returncode
            else:
                termination, exit_code = f"exited with code {res.returncode}", res.returncode
            last_line = (res.stderr.strip() or res.stdout.strip()).splitlines()[-1:] or [""]
            logger.error("mmseqs stdout:\n%s\nmmseqs stderr:\n%s", res.stdout[-3000:], res.stderr[-3000:])
            raise WrapperError(f"mmseqs easy-cluster {termination}: {last_line[0]}", exit_code=exit_code)
        return _read_clusters(prefix + "_cluster.tsv", records)


def _read_clusters(cluster_tsv, records):
    """clusters of the "representative <TAB> member" file of mmseqs, with the index ids mapped back."""
    clusters = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            rep, member = line.split()
            clusters.setdefault(rep, []).append(records[int(member)])
    return list(clusters.values())
