# fungi-30 rootHOG fixture

23 rootHOGs (about 1,700 proteins) taken from a full FastOMA run on the fungi-30 dataset (`-profile fungi`,
FastOMA 0.5.1, `dessimozlab/fastoma:dev` image, November 2026) and used by `tests/test_fungi30.py`.

- `rhogs/HOG_<id>.fa.gz`: the rootHOG FASTA files exactly as `fastoma-infer-roothogs` wrote them in that run
  (protein ids `<UniProt id>||<proteome>||<numeric id>`). The sequences are from UniProt reference proteomes
  (CC BY 4.0).
- `species_tree.nwk`: the `species_tree_checked.nwk` of that run (all 30 proteomes, with polytomies).
- `families.tsv`: why each rootHOG was selected, its tier (`fast` runs on every PR, `slow` nightly), whether
  it was run with `--parallel` (`par`, i.e. a `hog_big` rootHOG) and its kind:
  - `golden`: the output is compared with `golden/<id>.txt`.
  - `drops`: `fastoma-infer-subhogs` drops proteins of this rootHOG (known bug). It is only checked for the
    "no protein dropped" invariant, which is a strict xfail until the bug is fixed. How many proteins are
    dropped changes with small MSA differences (e.g. 367 vs 11 of 496 for F0412257 with mafft 7.505 vs
    7.526), so these rootHOGs are not snapshotted.
- `golden/`: one text snapshot per golden rootHOG (nesting, TaxRange and scores of the groups passed on to
  `collect-subhogs`) and `ENVIRONMENT.json` with the tool versions they were generated with.

## Selection

Every rootHOG of the run was classified from the `hog_big`/`hog_rest` task logs (conflict resolution,
split_hog, MSA trimming, fragment detection, dropped proteins, runtime) and the final orthoXML (nesting
depth, paralogGroups, scores). In that run, 414 of 22,878 rootHOGs dropped 8,113 of 265,464 proteins.

Running these rootHOGs in the `dessimozlab/fastoma:dev` image reproduced the results of the full run
exactly. With mafft 7.526 instead of 7.505, clust12762 and all large rootHOGs gave different results.

## Snapshots are not correctness

The snapshots record current behaviour, including known bugs (e.g. the ImpliedLosses values). When a change
of the output is intended, regenerate them and review the diff:

```bash
python tests/fungi30.py regenerate            # needs docker; --tier fast|slow|all
git diff tests/data/fungi30/golden
```

Regeneration runs the working-tree code inside the container image, so the snapshots match CI (Ubuntu 24.04
ships the same mafft 7.505 / FastTree 2.1.11 builds). With other tool versions, the snapshot tests are skipped
(or fail with `FASTOMA_REQUIRE_GOLDEN=1`, as in CI).
