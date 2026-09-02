# reassmflow

Initial DSL2 workflow for preparing a concatenated reference and recovering
long and/or short reads by reference mapping, followed by Flye and Myloasm
long-read assembly.

Example:

```bash
nextflow run main.nf -params-file config.yml
```

Pass `config.yml` explicitly because Nextflow automatically loads
`nextflow.config`, but not YAML parameter files.

Define long and short read libraries under `raw_libraries`. A sample can reuse
a library from `global_libraries` by listing its name in `libraries`; see
`config.example.yml`. At least one long- or short-read library is required.
Final assembly directories are copied to `<outdir>/<sample-id>/`.

Required tools are `seqkit`, `minimap2`, `bowtie2`, `samtools`, `pigz`, and
either Flye 2.9.6 or Myloasm >= 0.6.0 with `mylotools`.

Conda environments are assigned by process label: `mapping` uses
`envs/mapping.yaml`, `lr_assm` uses `envs/lr-env.yaml`, and future
short-read assembly processes use `envs/sr-env.yaml`. Enable Conda with
`nextflow run main.nf -params-file config.yml -with-conda`.

For UNSW Katana, set `conda_envs.mapping`, `conda_envs.lr_assm`, and
`conda_envs.sr_assm` in `config.yml` to absolute paths for existing shared
Conda environments. Set `pbs_queue` and `pbs_project` when required by the
site. Run with `-c conf/unsw_katana.config`; choose a shared high-performance
filesystem for Nextflow's `-work-dir`.
