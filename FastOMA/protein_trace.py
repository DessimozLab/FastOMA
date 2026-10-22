"""Merge the per-step protein event files into one fate per input protein.

fastoma-protein-trace --proteomes proteome --events events/ --roothog-tsv RootHOGs.tsv

Outputs (in the current directory):
  protein_fates.tsv          one row per input protein, with exactly one final fate, sorted by species and protein
  protein_events.tsv.gz      all recorded events, sorted by species and protein
  protein_trace_summary.txt  counts per fate, rootHOGs with most losses, consistency checks

Fates, in the order they are decided:
  too_short                         shorter than --min-sequence-length, never used
  in_orthoxml                       member of a HOG of size >= 2 in the orthoXML
  lost_in_split_hog                 dropped when split_hog did not cover all child subHOGs
  lost_at_level                     dropped at a taxonomic level for another reason
  final_singleton_not_written       ended alone in its HOG; orthoXML only keeps HOGs of size >= 2
  missing_at_end                    in a rootHOG file, but in no final HOG and no loss recorded
  roothog_not_processed             in a rootHOG file, but no subHOG result (e.g. failed task)
  isoform_not_selected              minor splicing variant, another isoform of the gene was selected
  singleton_roothog_dropped         alone in its rootHOG and not placed in a cluster
  no_omamer_hit_dropped             no OMAmer hit and not placed in a cluster
  removed_by_filter_big_roothogs    removed from a very big rootHOG (low family score)
  unaccounted                       no event at all (should not happen)

The event files are sorted on disk (GNU sort), and the proteins are then processed
one at a time, so memory stays at a few hundred bytes per input protein no matter
how many events were recorded.
"""

import argparse
import collections
import csv
import gzip
import itertools
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

from . import __version__ as fastoma_version
from ._protein_trace import COLUMNS, EVENT_FILE_SUFFIX
from .transformer import header_transformer

LOSS_EVENTS = ("lost_in_split_hog", "lost_at_level")


def read_proteomes(folder):
    """(species, protein) -> length."""
    proteins = {}
    for f in sorted(os.listdir(folder)):
        species, ext = f.rsplit(".", 1)
        if ext not in ("fa", "fasta"):
            continue
        species = sys.intern(species)
        for rec in SeqIO.parse(os.path.join(folder, f), "fasta"):
            proteins[(species, rec.id)] = len(rec)
    return proteins


def _parse_event_lines(lines):
    _i = sys.intern   # species, step, event and roothog repeat a lot, share one copy of each
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != len(COLUMNS):
            continue
        species, protein, step, event, roothog, detail = fields
        yield _i(species), protein, _i(step), _i(event), _i(roothog), detail


def sorted_events(paths, tmp_dir, sort_memory):
    """Yield ((species, protein), [(step, event, roothog, detail), ...]), ordered by (species, protein).

    The event files are sorted on disk with GNU sort, so only the events of one protein
    are held in memory at a time. The sort is stable and the files are given in name order
    (infer_roothogs before infer_subhogs), so the events of a protein keep the order in
    which they were written. LC_ALL=C makes sort compare bytes, which for UTF-8 is the
    same order as comparing Python str tuples.
    """
    if not paths:
        return
    with tempfile.NamedTemporaryFile("wb", dir=tmp_dir, suffix=".files0") as file_list:
        file_list.write(b"".join(os.fsencode(p) + b"\0" for p in paths))
        file_list.flush()
        cmd = ["sort", "-s", "-t", "\t", "-k1,1", "-k2,2", "-S", sort_memory, "-T", str(tmp_dir),
               "--files0-from=" + file_list.name]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, env=dict(os.environ, LC_ALL="C"),
                              encoding="utf-8") as proc:
            for key, group in itertools.groupby(_parse_event_lines(proc.stdout), key=lambda f: (f[0], f[1])):
                yield key, [f[2:] for f in group]
        if proc.returncode:
            raise RuntimeError(f"sorting the event files failed: {' '.join(cmd)} exited with {proc.returncode}")


def join_proteins_events(protein_keys, event_groups):
    """Merge the sorted input proteins with the sorted event groups.

    Yields (key, events, is_input_protein); input proteins without any event get [].
    """
    keys = iter(protein_keys)
    key = next(keys, None)
    for ev_key, evs in event_groups:
        while key is not None and key < ev_key:
            yield key, [], True
            key = next(keys, None)
        if key == ev_key:
            yield key, evs, True
            key = next(keys, None)
        else:
            yield ev_key, evs, False
    while key is not None:
        yield key, [], True
        key = next(keys, None)


def event_files(inputs):
    for p in inputs:
        p = Path(p)
        if p.is_dir():
            yield from sorted(p.rglob("*" + EVENT_FILE_SUFFIX))
        elif p.suffix == EVENT_FILE_SUFFIX:
            yield p


def read_roothog_tsv(path):
    """protein -> tuple of final RootHOG ids (a protein id could in principle occur in two species)."""
    final = {}
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                final[row["Protein"]] = final.get(row["Protein"], ()) + (sys.intern(row["RootHOG"]),)
    return final


def decide(evs):
    """Return (fate, step, roothog, detail) from the events of one protein."""
    by_event = collections.defaultdict(list)
    for step, event, roothog, detail in evs:
        by_event[event].append((step, roothog, detail))

    def first(name):
        step, roothog, detail = by_event[name][0]
        return name, step, roothog, detail

    if "too_short" in by_event:
        return first("too_short")
    if "in_final_hog" in by_event:
        return first("in_final_hog")
    for loss in LOSS_EVENTS:
        if loss in by_event:
            return first(loss)
    if "final_singleton_not_written" in by_event:
        return first("final_singleton_not_written")
    if "missing_at_end" in by_event:
        return first("missing_at_end")
    if "roothog_written" in by_event:
        _, step, roothog, detail = first("roothog_written")
        return "roothog_not_processed", step, roothog, detail
    if "isoform_not_selected" in by_event:
        return first("isoform_not_selected")
    if "filter_big_roothogs_removed" in by_event:
        _, step, roothog, detail = first("filter_big_roothogs_removed")
        return "removed_by_filter_big_roothogs", step, roothog, detail
    reason = next((e for e in ("cluster_single_species_discarded", "cluster_too_small", "sent_to_clustering") if e in by_event), "")
    if "no_omamer_hit" in by_event:
        return "no_omamer_hit_dropped", "infer_roothogs", "", reason
    if "roothog_too_small" in by_event or "omamer_best_hit" in by_event:
        roothog = (by_event.get("roothog_too_small") or by_event["omamer_best_hit"])[0][1]
        return "singleton_roothog_dropped", "infer_roothogs", roothog, reason
    return "unaccounted", "", "", ""


def size_bucket(n):
    if n is None:
        return "-"
    return "2-9" if n < 10 else "10-99" if n < 100 else "100-999" if n < 1000 else ">=1000"


def fastoma_protein_trace():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version="FastOMA v" + fastoma_version)
    parser.add_argument("--proteomes", required=True, help="input proteome folder")
    parser.add_argument("--events", required=True, nargs="+", help="event files or folders containing *.events")
    parser.add_argument("--roothog-tsv", help="RootHOGs.tsv from fastoma-collect-subhogs (final RootHOG ids)")
    parser.add_argument("--id-transform", choices=("UniProt", "noop"), default="noop",
                        help="ID transformer that fastoma-collect-subhogs used for the protein IDs in RootHOGs.tsv")
    parser.add_argument("--out-dir", default=".", help="where to write the three output files")
    parser.add_argument("--strict", action="store_true", help="exit with an error if any protein is unaccounted")
    parser.add_argument("--tmp-dir", help="folder for the temporary files of the on-disk sort (default: --out-dir)")
    parser.add_argument("--sort-memory", default="1G", help="main memory buffer of the on-disk sort (sort -S)")
    conf = parser.parse_args()

    proteins = read_proteomes(conf.proteomes)
    files = list(event_files(conf.events))
    final_roothog = read_roothog_tsv(conf.roothog_tsv)
    id_transformer = header_transformer(conf.id_transform)
    out = Path(conf.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(conf.tmp_dir) if conf.tmp_dir else out
    tmp_dir.mkdir(parents=True, exist_ok=True)
    n_events = 0
    n_unknown = 0

    fates = collections.Counter()
    fates_by_size = collections.Counter()
    losses_by_rhog = collections.Counter()
    rhog_size = {}
    issues = collections.Counter()

    with open(out / "protein_fates.tsv", "w", newline="") as fh, \
            gzip.open(out / "protein_events.tsv.gz", "wt") as events_fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["species", "protein", "length", "omamer_roothog", "omamer_n_hits", "roothog_changes",
                    "roothog", "roothog_size", "fate", "fate_step", "fate_detail", "final_roothog"])
        events_fh.write("\t".join(COLUMNS) + "\n")
        event_groups = sorted_events(files, tmp_dir, conf.sort_memory)
        for key, evs, is_input in join_proteins_events(sorted(proteins), event_groups):
            for ev in evs:
                events_fh.write("\t".join(key + ev) + "\n")
            n_events += len(evs)
            if not is_input:
                n_unknown += 1
                continue
            length = proteins[key]
            fate, step, rhog, detail = decide(evs)
            if fate == "in_final_hog":
                fate = "in_orthoxml"
            omamer = next(((r, d) for s, e, r, d in evs if e == "omamer_best_hit"), ("", ""))
            changes = ";".join(f"{e.rsplit('_', 1)[0]}:{d.replace('from=', '')}->{r}" for s, e, r, d in evs if e.endswith("_moved"))
            written = [(r, d) for s, e, r, d in evs if e == "roothog_written"]
            if len(written) > 1:
                issues["in_more_than_one_roothog_file"] += 1
            roothog, size = (written[0][0], int(written[0][1].split("size=")[1].split()[0])) if written else (rhog, None)
            # RootHOGs.tsv has the protein IDs as written to the orthoXML, i.e. after the id transform
            finals = final_roothog.get(id_transformer.transform(key[1]), ()) if fate == "in_orthoxml" else ()
            if fate == "in_orthoxml" and not finals:
                issues["in_final_hog_but_not_in_RootHOGs.tsv"] += 1
            if len(finals) > 1:
                issues["protein_id_in_several_RootHOGs.tsv_rows"] += 1
            w.writerow([key[0], key[1], length, omamer[0], omamer[1].replace("n_hits=", ""), changes,
                        roothog, "" if size is None else size, fate, step, detail, ",".join(finals)])
            fates[fate] += 1
            fates_by_size[(fate, size_bucket(size))] += 1
            if written:
                rhog_size[roothog] = size
                if fate in LOSS_EVENTS or fate in ("missing_at_end", "roothog_not_processed"):
                    losses_by_rhog[roothog] += 1

    if n_unknown:
        issues["events_for_proteins_not_in_proteomes"] = n_unknown

    total = len(proteins)
    buckets = ["2-9", "10-99", "100-999", ">=1000", "-"]
    lines = [f"FastOMA protein trace ({len(files)} event files, {n_events:,} events)", "",
             f"input proteins: {total:,}", "", "fate                               proteins        %"]
    for fate, n in fates.most_common():
        lines.append(f"{fate:32} {n:10,} {100 * n / total:8.2f}")
    lines += ["", "fate by rootHOG size (proteins)", f"{'fate':32}" + "".join(f"{b:>10}" for b in buckets)]
    for fate, _ in fates.most_common():
        lines.append(f"{fate:32}" + "".join(f"{fates_by_size[(fate, b)]:10,}" for b in buckets))
    lines += ["", "rootHOGs losing most proteins (lost, rootHOG size, rootHOG)"]
    for rhog, n in losses_by_rhog.most_common(20):
        lines.append(f"  {n:7,} {rhog_size.get(rhog, 0):7,}  {rhog}")
    lines += ["", "consistency checks"]
    lines.append(f"  every input protein has exactly one fate: yes ({sum(fates.values()):,} rows)")
    lines.append(f"  unaccounted proteins: {fates.get('unaccounted', 0):,}")
    for issue, n in sorted(issues.items()):
        lines.append(f"  {issue}: {n:,}")
    summary = "\n".join(lines) + "\n"
    (out / "protein_trace_summary.txt").write_text(summary)
    print(summary)

    if conf.strict and fates.get("unaccounted", 0):
        sys.exit(f"{fates['unaccounted']} proteins are unaccounted for")


if __name__ == "__main__":
    fastoma_protein_trace()
