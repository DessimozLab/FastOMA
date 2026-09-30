"""Per-protein event log, so that the fate of every input protein can be traced.

Each FastOMA step appends events to its own file in the task directory
(protein_events_<step>_<uuid>.events, one per process so that parallel
Nextflow tasks never collide). The columns are, tab separated and without header:

    species  protein  step  event  roothog  detail

`fastoma-protein-trace` (FastOMA/protein_trace.py) merges all event files at the
end of the run into one row per input protein.

Only changes are recorded (e.g. a protein being moved to another rootHOG, or
getting lost), plus the final state of every protein, so the files stay small.
"""

import atexit
import threading
import uuid

EVENT_FILE_SUFFIX = ".events"
COLUMNS = ["species", "protein", "step", "event", "roothog", "detail"]

_lock = threading.Lock()
_handles = {}


def _handle(step):
    fh = _handles.get(step)
    if fh is None:
        fh = open(f"protein_events_{step}_{uuid.uuid4().hex[:12]}{EVENT_FILE_SUFFIX}", "wt")
        _handles[step] = fh
    return fh


def event(species, protein, step, event_name, roothog="", detail=""):
    line = "\t".join(str(x).replace("\t", " ") for x in (species, protein, step, event_name, roothog, detail))
    with _lock:
        _handle(step).write(line + "\n")


def split_member(member):
    """Internal FastOMA member ids look like 'protId||species||numericId', or several of
    those joined by '_|_' for merged fragments. Returns a list of (species, protId)."""
    out = []
    for part in str(member).split("_|_"):
        fields = part.split("||")
        if len(fields) >= 2:
            out.append((fields[1], fields[0]))
    return out


def member_event(member, step, event_name, roothog="", detail=""):
    for species, protein in split_member(member):
        event(species, protein, step, event_name, roothog, detail)


def rhog_assignment(rhogs_prots):
    """(species, protId) -> rootHOG id, for diffing the rootHOG grouping steps."""
    return {sp_prot: rhogid for rhogid, sp_prots in rhogs_prots.items() for sp_prot in sp_prots}


def record_rhog_changes(before, after, step, event_prefix):
    """Record proteins that changed rootHOG, or disappeared, between two assignments."""
    for sp_prot, old in before.items():
        new = after.get(sp_prot)
        if new is None:
            event(*sp_prot, step, f"{event_prefix}_removed", old)
        elif new != old:
            event(*sp_prot, step, f"{event_prefix}_moved", new, f"from={old}")
    for sp_prot, new in after.items():
        if sp_prot not in before:
            event(*sp_prot, step, f"{event_prefix}_added", new)


def close():
    with _lock:
        for fh in _handles.values():
            fh.close()
        _handles.clear()


atexit.register(close)
