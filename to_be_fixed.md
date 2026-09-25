# Correction Tracker

## Completed

- [x] **PBSPro profile parses on current Nextflow.** Removed invalid top-level
	variable declarations from `conf/unsw_katana.config`. The profile was parsed
	successfully by Nextflow 26.04.6.

- [x] **Missing required inputs fail early.** `main.nf` now errors clearly when
	`input_csv` or both library sections are absent. The Katana profile falls back
	to the repository Conda YAML files when shared environment paths are not
	supplied.

- [x] **`global_libraries` is implemented.** Global libraries apply to every
	sample in `input_csv`; `raw_libraries` adds libraries only to its matching
	sample ID. `config.example.yml` documents the schema.

- [x] **Long-only and short-only runs are supported.** Long-read mapping and
	assembly, short-read mapping and MEGAHIT, and hybrid SPAdes branches are now
	scheduled only when their required libraries are present.

- [x] **Long-read BAM staging collisions are prevented.** `MINIMAP2_MAP` now
	includes `${task.index}` in each BAM filename.

- [x] **BAM merges run independently per sample.** Long- and short-read
	mappings carry a `groupKey` with each sample's expected mapping count, so a
	merge and downstream assembly start once that sample's mappings complete
	without waiting for unrelated samples.

- [x] **Assembly results are published.** Flye, Myloasm, MEGAHIT, and SPAdes
	result directories are copied under `<outdir>/<sample-id>/`.

- [x] **Interrupted assemblies can continue.** Per-assembler state is retained
	under `<outdir>/.resume/<sample-id>/`; Flye, Myloasm, MEGAHIT, and SPAdes use
	their respective continuation command when prior state exists.

- [x] **Empty reference folders fail.** Reference concatenation now verifies
	that at least one FASTA exists and sorts input paths before concatenation for
	reproducible ordering.

## Site Configuration Required

- [ ] **Set PBS queue and project for the target cluster.** Add `pbs_queue` and
	`pbs_project` to the run parameter file when required by Katana. The profile
	forwards them to PBSPro but cannot infer site/account values.

- [ ] **Use a shared high-performance work directory.** Supply `-work-dir`
	pointing to filesystem storage visible to both the Nextflow launcher and PBS
	compute nodes, for example `/srv/scratch/<user>/reassmflow-work`.

- [ ] **Use accessible paths in the run config.** Replace all placeholder paths
	in `config.example.yml`, including the sample CSV paths and shared Conda
	environment locations, before submitting a real run.

## Validation Status

- [x] Nextflow 26.04.6 compiles the workflow and constructs the complete
	process graph with `conf/unsw_katana.config`.
- [ ] Run a `-stub-run` using an accessible small reference and read set.
- [ ] Submit a small real PBSPro test after confirming the queue, project,
	Conda paths, and work directory for the target site.

