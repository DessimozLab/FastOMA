"""Golden-snapshot fixture built from rootHOGs of the fungi-30 run (tests/data/fungi30).

The rootHOGs listed in families.tsv were selected from a full fungi-30 run because they reach specific
code paths of fastoma-infer-subhogs (conflict resolution, split_hog, MSA trimming, fragments, deep
nesting, ...). They come in two kinds:

golden  the per-rootHOG output (structure and scores) is compared with a stored snapshot in golden/.
        The snapshot captures current behaviour, it is not a statement of correctness.
drops   fastoma-infer-subhogs drops proteins of these rootHOGs (a known bug). How many it drops depends
        on small differences in the MSA, so they are only checked for the invariant "every input protein
        is in a final HOG" (strict xfail until the bug is fixed) and are not snapshotted.

The output depends on the mafft/FastTree versions, so the snapshots are only compared in the
environment they were generated in (golden/ENVIRONMENT.json, i.e. the dessimozlab/fastoma image and the
apt packages of Ubuntu 24.04 used in CI).

Regenerate the snapshots after an intended change (needs docker; runs the working-tree code with the
tools of the container image):

    python tests/fungi30.py regenerate [--tier fast|slow|all] [--image dessimozlab/fastoma:dev]
"""
import argparse
import csv
import gzip
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).parent
REPO_DIR = TESTS_DIR.parent
DATA = TESTS_DIR / "data" / "fungi30"
SPECIES_TREE = DATA / "species_tree.nwk"
GOLDEN = DATA / "golden"
ENVIRONMENT_FILE = GOLDEN / "ENVIRONMENT.json"

# keep the per-level pickles so the final HOG objects can be inspected for the invariant
RUNNER = ("import FastOMA._infer_subhog as m; m.keep_subhog_each_pickle = True; "
          "from FastOMA.infer_subhogs import fastoma_infer_subhogs; fastoma_infer_subhogs()")


def families():
    with open(DATA / "families.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def groups():
    """families grouped by (tier, mode); each group is one fastoma-infer-subhogs run."""
    out = {}
    for fam in families():
        out.setdefault((fam["tier"], fam["mode"]), []).append(fam["rhog"])
    return out


# ---------------------------------------------------------------- environment

def _first_line(cmd):
    res = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=30)
    return (res.stdout + res.stderr).strip().splitlines()[0]


def parse_tool_versions(mafft_line, fasttree_line):
    """normalise `mafft --version` and the first line of `FastTree -expert`.

    The number of OpenMP threads is dropped, but whether FastTree was built with OpenMP is kept."""
    fasttree = fasttree_line.replace("Detailed usage for FastTree", "").strip().rstrip(":")
    fasttree = re.sub(r"OpenMP \(\d+ threads\)", "OpenMP", fasttree)
    return {"mafft": mafft_line.strip(), "FastTree": fasttree}


def tool_versions():
    return parse_tool_versions(_first_line(["mafft", "--version"]), _first_line(["FastTree", "-expert"]))


def golden_environment():
    with open(ENVIRONMENT_FILE) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- running

def prepare_input(rhogids, dest):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    for rhogid in rhogids:
        with gzip.open(DATA / "rhogs" / f"HOG_{rhogid}.fa.gz", "rb") as src, open(dest / f"HOG_{rhogid}.fa", "wb") as out:
            shutil.copyfileobj(src, out)
    return dest


def infer_subhogs_args(input_dir, species_tree, parallel):
    args = ["--input-rhog-folder", str(input_dir), "--species-tree", str(species_tree),
            "--output-pickles", "pickles"]
    return args + (["--parallel"] if parallel else [])


def run_infer_subhogs(input_dir, workdir, parallel):
    """runs the working-tree fastoma-infer-subhogs on input_dir inside workdir."""
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=str(REPO_DIR))
    cmd = [sys.executable, "-c", RUNNER] + infer_subhogs_args(input_dir, SPECIES_TREE, parallel)
    res = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"fastoma-infer-subhogs failed:\n{res.stdout[-3000:]}\n{res.stderr[-3000:]}")


# ---------------------------------------------------------------- reading results

def input_records(rhogid):
    with gzip.open(DATA / "rhogs" / f"HOG_{rhogid}.fa.gz", "rt") as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def final_members(workdir, rhogid):
    """all proteins in the final HOG objects of a rootHOG (singletons included)."""
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from ete3 import Tree
    from FastOMA._utils_subhog import prepare_species_tree

    recs = [SeqRecord(Seq(""), id=i) for i in input_records(rhogid)]
    root = prepare_species_tree(recs, Tree(str(SPECIES_TREE), format=1), rhogid)[0]
    with open(Path(workdir) / f"rhog_{rhogid}" / f"{root.name}.pickle", "rb") as fh:
        hogs = pickle.load(fh)
    return [p for h in hogs for m in h.get_members() for p in m.split("_|_")]


def output_groups(workdir, rhogid):
    """the orthoXML group elements that fastoma-infer-subhogs passes on to collect-subhogs."""
    with open(Path(workdir) / "pickles" / f"file_{rhogid}.pickle", "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------- rendering

_SCORE_ORDER = ("CompletenessScore", "ImpliedLosses", "TCSScore")


def _render(elem, id2prot):
    """an indented, order-independent text form of an orthoXML group (lines list)."""
    if elem.tag == "geneRef":
        return [id2prot.get(elem.get("id"), elem.get("id"))]
    props = {p.get("name"): p.get("value") for p in elem.findall("property")}
    scores = {s.get("id"): s.get("value") for s in elem.findall("score")}
    if elem.tag == "orthologGroup":
        head = "orthologGroup " + props.get("TaxRange", "?")
    else:
        head = "paralogGroup" + (f" [{props['Type']}]" if "Type" in props else "")
    score_txt = " ".join(f"{k}={scores[k]}" for k in _SCORE_ORDER if k in scores)
    if score_txt:
        head += "  " + score_txt
    children = [_render(c, id2prot) for c in elem if c.tag in ("orthologGroup", "paralogGroup", "geneRef")]
    children.sort(key=lambda lines: "\n".join(lines))
    return [head] + ["  " + line for child in children for line in child]


def golden_text(groups_, input_ids):
    """snapshot text of one rootHOG: header line, then one block per top-level group."""
    id2prot = {i.split("||")[2]: i.split("||")[0] for i in input_ids}
    blocks = sorted("\n".join(_render(g, id2prot)) for g in groups_)
    n_ref = sum(1 for g in groups_ for _ in g.iter("geneRef"))
    header = f"# {len(input_ids)} input proteins, {len(groups_)} groups, {n_ref} proteins in groups"
    return header + "\n\n" + "\n\n".join(blocks) + "\n"


def golden_path(rhogid):
    return GOLDEN / f"{rhogid}.txt"


# ---------------------------------------------------------------- regeneration

def _docker_run(image, mount, workdir, args):
    # run as the calling user, otherwise a Linux docker daemon leaves root-owned files in the mounted
    # repository that the caller cannot remove; HOME must then point to a writable folder
    user = ["--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp"] if hasattr(os, "getuid") else []
    cmd = ["docker", "run", "--rm", "--platform", "linux/amd64", *user, "-v", f"{mount}:{mount}", "-w", str(workdir),
           "-e", "PYTHONHASHSEED=0", "-e", f"PYTHONPATH={REPO_DIR}", "-e", "PYTHONDONTWRITEBYTECODE=1",
           "--entrypoint", "sh", image, "-c", args]
    return subprocess.run(cmd, capture_output=True, text=True)


def _tools(env):
    return {k: v for k, v in env.items() if k != "image"}


def regenerate(tier, image):
    GOLDEN.mkdir(parents=True, exist_ok=True)
    fams = {f["rhog"]: f for f in families()}
    res = _docker_run(image, REPO_DIR, REPO_DIR, "mafft --version </dev/null 2>&1; FastTree -expert </dev/null 2>&1 | head -1")
    lines = res.stdout.strip().splitlines()
    env = parse_tool_versions(lines[0], lines[-1])
    env["image"] = image
    print("environment:", env)
    # ENVIRONMENT.json is shared by all snapshots: a partial regeneration with other tool versions would
    # leave the other tier's snapshots from the old tools while claiming the new ones
    if tier != "all" and ENVIRONMENT_FILE.exists() and _tools(golden_environment()) != _tools(env):
        sys.exit(f"the tool versions of {image} ({_tools(env)}) differ from those of the existing snapshots "
                 f"({_tools(golden_environment())}); regenerate all tiers with --tier all")
    # docker can only see the repository (mounted at the same path), so work in a folder inside it
    with tempfile.TemporaryDirectory(dir=TESTS_DIR, prefix=".fungi30-") as tmp:
        for (grp_tier, mode), rhogids in groups().items():
            if tier != "all" and grp_tier != tier:
                continue
            golden_ids = [r for r in rhogids if fams[r]["kind"] == "golden"]
            if not golden_ids:
                continue
            workdir = Path(tmp) / f"{grp_tier}_{mode}"
            workdir.mkdir()
            prepare_input(golden_ids, workdir / "in")
            args = " ".join(["/app/bin/python", "-c", repr(RUNNER)] +
                            infer_subhogs_args(workdir / "in", SPECIES_TREE, mode == "par"))
            print(f"running {grp_tier}/{mode}: {', '.join(golden_ids)}", flush=True)
            res = _docker_run(image, REPO_DIR, workdir, args)
            if res.returncode != 0:
                sys.exit(f"fastoma-infer-subhogs failed in {image}:\n{res.stdout[-3000:]}\n{res.stderr[-3000:]}")
            for rhogid in golden_ids:
                golden_path(rhogid).write_text(golden_text(output_groups(workdir, rhogid), input_records(rhogid)))
                print("  wrote", golden_path(rhogid).relative_to(REPO_DIR))
    with open(ENVIRONMENT_FILE, "w") as fh:
        json.dump(env, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    sys.path.insert(0, str(REPO_DIR))
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("regenerate", help="regenerate golden snapshots with docker")
    p.add_argument("--tier", choices=("fast", "slow", "all"), default="all")
    p.add_argument("--image", default="dessimozlab/fastoma:dev")
    conf = parser.parse_args()
    regenerate(conf.tier, conf.image)
