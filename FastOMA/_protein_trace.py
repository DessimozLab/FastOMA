"""Per-protein event log, so that the fate of every input protein can be traced.

Each FastOMA step appends events to its own file in the task directory
(protein_events_<step>_<uuid>.events, one per process so that parallel
Nextflow tasks never collide). The columns are, tab separated and without header:

    species  protein  step  event  roothog  detail

`fastoma-protein-trace` (FastOMA/protein_trace.py) merges all event files at the
end of the run into one row per input protein.

Nothing is recorded unless enable() was called (the --protein-trace option of
fastoma-infer-roothogs and fastoma-infer-subhogs). Expect a few events per protein,
i.e. roughly 250-400 bytes per input protein in total.
"""

import atexit
import threading
import uuid

EVENT_FILE_SUFFIX = ".events"
COLUMNS = ["species", "protein", "step", "event", "roothog", "detail"]

_lock = threading.Lock()
_handles = {}
_enabled = False


def enable():
    global _enabled
    _enabled = True


def is_enabled():
    return _enabled


def _handle(step):
    fh = _handles.get(step)
    if fh is None:
        fh = open(f"protein_events_{step}_{uuid.uuid4().hex[:12]}{EVENT_FILE_SUFFIX}", "wt")
        _handles[step] = fh
    return fh


def event(species, protein, step, event_name, roothog="", detail=""):
    if not _enabled:
        return
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
    if not _enabled:
        return
    for species, protein in split_member(member):
        event(species, protein, step, event_name, roothog, detail)


def rhog_assignment(rhogs_prots):
    """(species, protId) -> rootHOG id, for diffing the rootHOG grouping steps (empty if tracing is off)."""
    if not _enabled:
        return {}
    return {sp_prot: rhogid for rhogid, sp_prots in rhogs_prots.items() for sp_prot in sp_prots}


def record_rhog_changes(before, after, step, event_prefix):
    """Record proteins that changed rootHOG, or disappeared, between two assignments."""
    if not _enabled:
        return
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
