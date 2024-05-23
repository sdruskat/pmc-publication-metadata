<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

# Snakemake workflow: Extract LUTs from PMC OA non-commercial, commercial and other bulk metadata

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for extracting accessible lookup tables for PMC publication dates.
The workflow attempts to extract the earliest publication date for each publication from
the metadata and maps them to publication identifiers. The lookup tables are in JSON,
where identifiers have the format `PMC<numerical identifier` and dates have the format
`YYYY[-MM[-DD]]`.

```json
{
  "PMC1234567": "2020-02-20",
  "PMC1234568": "2020-02",
  "PMC1234569": "2020"
}
```

## Documentation

The technical documentation and description of outputs is in [workflow/documentation.md](workflow/documentation.md).

## Running the workflow

You need to have `conda` installed to create and activate a new environment.

```bash
conda env create -n pmc-metadata --file conda-environment.yaml
conda activate pmc-metadata
```

Run with `-–keep-storage-local-copies` to avoid downloading resources over and over again.
Also run with `--software-deployment-method conda` to use global conda packages.

```shell
snakemake --keep-storage-local-copies --software-deployment-method conda -c <number-of-cores-to-use>
```

You can use [`run.sh`](run.sh) to run the workflow this way, and with 12 cores.

## Running on a cluster with Slurm

Create a Snakemake profile, then run as follows.

```yaml
# ~/.config/snakemake/<profilename>/config.v8+.yaml
executor: slurm
jobs: 100
local-storage-prefix: /scratch/$USER/snakemake-scratch
latency-wait: 300
default-resources:
  slurm_partition: cpu
  nodes: 1
  tasks: 1
  cpus_per_task: 1
```

```shell
nohup snakemake --keep-storage-local-copies --software-deployment-method conda --profile <profilename> --directory <working directory> --configfile config/config.yml &
```

# Citation

If you use this workflow in your work, please cite it using the metadata provided in [`CITATION.cff`](CITATION.cff).

# License

This work is licensed as specific in the [REUSE 3.0 Specification](https://reuse.software/spec/). 
Please consult the single file licenses.
