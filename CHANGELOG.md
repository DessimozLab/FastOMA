# Changelog

All notable changes to FastOMA. Changes that alter the inferred orthology are listed under
"Changed output", so results of different versions can be compared knowingly.

## [Unreleased]

### Requirements
- **Nextflow >= 25.10 is now required** (was >= 22.10.4). With Nextflow 26, boolean flags on the
  command line such as `--protein_trace` failed the parameter validation ("Value is [string] but
  should be [boolean]"); the newer nf-schema plugin (2.7.3) that fixes this needs Nextflow 25.10.

### Added
- `--denovo`: infers the rootHOGs by clustering all proteins with mmseqs instead of placing them
  with OMAmer, e.g. for simulated proteomes. OMAmer is not run, and its database is not needed.
- `--mmseqs_coverage`, `--mmseqs_cov_mode`, `--mmseqs_sensitivity` and `--mmseqs_min_seq_id` set
  the mmseqs clustering of the proteins without OMAmer placement (with `--denovo`: of all
  proteins). Unset, the mmseqs defaults are used as before. The memory of the mmseqs clustering
  is capped at 60% of the memory of the infer_roothogs task, but at least 2 GB (mmseqs
  `--split-memory-limit`); with `--denovo`, the task memory grows with the size of the proteomes.
- `--protein_trace`: records what happens to every input protein and writes `protein_fates.tsv`,
  `protein_events.tsv.gz` and `protein_trace_summary.txt` (#106).
- `TCSScore` and `ImpliedLosses` scores for every (sub)HOG in `FastOMA_HOGs.orthoxml`, next to
  `CompletenessScore`. Each score can be turned off with `--disable_completeness_score`,
  `--disable_implied_losses_score` and `--disable_tcs_score`.
- `utils/generate_splice_file_Ensembl.py`, which writes splice files from Ensembl proteomes.

### Changed output
- **Splice genes without any OMAmer hit keep their longest isoform**, no longer the first one listed
  in the splice file; a fragment listed first did not cluster with its orthologs. With `--denovo`
  this applies to every splice gene.
- **No proteins are dropped any more.** fastoma-infer-subhogs lost proteins in several places
  (about 3% of all proteins on a 30-species fungi dataset), fastoma-infer-roothogs in some
  others (proteins missing from the hogmap, ids longer than 230 characters, oversized merges of
  rootHOGs, splice genes whose first isoform is too short). Both steps now check at every level
  that every protein is accounted for (#113, #114, #130, #131). Proteins that were lost before
  now appear in the output; proteins recovered inside infer-subhogs are kept as separate HOGs,
  not merged back into their families.
- **Top-level groups that start with a duplication** (a `paralogGroup` in an extra
  `orthologGroup`) could get a TaxRange up to 6 levels too old. They now get the level directly
  above the duplication, and its CompletenessScore and HOG id (which contains the taxid of the
  TaxRange) change with it. When the duplication is at the root of the species
  tree, the copies become separate top-level groups (#128, #129).
- The proteins of merged rootHOGs are listed in a fixed order. Before, the order depended on
  Python's hash seed, so MSAs and gene trees of these rootHOGs differ slightly from earlier
  runs (#131).

### Fixed
- If the proteome folder (`proteome/` in `--input`, or `--proteome_folder`) was a symbolic link,
  the run stopped after `check_input` and still reported success: no proteome was found in it.
  Symlinked folders work now, and a proteome folder without `.fa`/`.fasta` files stops the run
  with an error (#83).
- With Nextflow 26, `--write_msas false`, `--write_genetrees false`, `--report false`,
  `--force_pairwise_ortholog_generation false` and `--debug_enabled false` on the command line
  turned the option on, and `--memory_multiplier` / `--time_multiplier` made the tasks fail
  (Nextflow 26 passes command-line values as strings).
- mmseqs is limited to the CPUs Nextflow assigns to the task; before, it used all cores of the
  machine. When mmseqs is killed (e.g. out of memory), infer-roothogs exits with the same code, so
  Nextflow retries the task with more resources.
- `--help` printed the help message and then failed with "Specified param ... does not exist in
  JSON schema"; `--help <parameter>` and `--helpFull` work too.
- Input checks: hidden files in the input folders (e.g. the `._<name>` files of archives made on
  macOS, `.DS_Store`) are ignored; empty proteomes, files without an extension, species trees
  without internal nodes and duplicate internal node names are reported instead of crashing
  (#118, #125).
- The rootHOG batch folders are staged as inputs, so they are visible inside the container
  wherever the work directory is (#108, #109).
- When mafft or FastTree is killed (e.g. out of memory on SLURM), FastOMA exits with the same
  return code, so Nextflow retries the task with more resources.
- Several internal bugs in the HOG class, UnionFind, the batching of rootHOGs and the pairwise
  ortholog extraction (#122, #125).
- The docker image contains papermill again, which the report needs.

## [0.5.1]
- Bug with input handling fixed.
- upgraded github actions for CI/CD.

## [0.5.0]
- renamed input_folder parameter to input. input accepts now also (remote) archive tarball files.
- better configuration setup (close to nf-core)
- improved resource allocation for nextflow
- improved handling of alternative splicing variants in reporting
- adding test profile and nf-test based CI checks

## [0.4.1]
- bump omamer to 2.1.2 (fixes a bug in the version 2.1.1, see [#74](https://github.com/DessimozLab/omamer/issues/74)

## [0.4.0]
- Improvements for nextflow: alternative version selection, README updates
- Split HOG and sampling improvements
- ensure orthologGroup at the MRCA of all genes in the group
- Docker improvements, additional labels and tools
- New gene families with mmseqs easy-cluster if mmseqs is installed
- Improved input checking
- Merge rootHOGs, improved handling of singleton using omamer multi-hits
- Various documentation and usability updates

## [0.3.5]
- Fixes an issue with reaching the maximum recursion limit. (#31)
- Fixes a problem with parallel execution for big families. (#44)

## [0.3.4]
- Fixing a bug in marker gene groups extraction. Before, more than one gene per species were possible

## [0.3.3]
- improvements for nextflow (selection of alternative versions) and updates on readme

## [0.3.1]
- spliting HOG and sampling

## [0.1.6]
- adding dynamic resources, additional and improved output

## [0.1.5]
- docker, add help, clean nextflow

## [0.1.4]
- new gene families with linclust if mmseqs is installed, using quoted protein name to handle species chars, check input first

## [0.1.3]
- merge rootHOGs and handle singleton using omamer multi-hits

## [0.1.2]
- improve rootHOG inference, splice, OMAmerv2 with multi-hits

## [0.1.0]
- improve nextflow pipeline and outputs.

## Pre-releases
- v0.0.6: use `--fragment-detection` for `infer-subhogs` and `--low-so-detection --fragment-detection`
- v0.0.6: using input hogmpa
- v0.0.5: adding pip setup.py
- v0.0.4: simple nextflow
