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

Short reads can be supplied with `--short_r1`, `--short_r2`, and optionally
`--short_singletons`. Required tools are `seqkit`, `minimap2`, `bowtie2`,
`samtools`, `pigz`, and either Flye 2.9.6 or Myloasm >= 0.6.0 with
`mylotools`.

Conda environments are assigned by process label: `mapping` uses
`envs/mapping.yaml`, `lr_assm` uses `envs/lr-env.yaml`, and future
short-read assembly processes use `envs/sr-env.yaml`. Enable Conda with
`nextflow run main.nf -params-file config.yml -with-conda`.

For UNSW Katana, use `-c conf/unsw_katana.config` after exporting
`REASSMFLOW_MAPPING_ENV`, `REASSMFLOW_LR_ASSM_ENV`, and
`REASSMFLOW_SR_ASSM_ENV` as absolute paths to the existing shared Conda
environment directories.
