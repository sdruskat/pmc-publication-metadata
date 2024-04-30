# Snakemake workflow: Extract LUTs from PMC OA non-comm bulk metadata

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for extracting accessible lookup tables for PMC publication metadata.

## Documentation

The technical documentation and description of outputs is in [workflow/documentation.md](workflow/documentation.md).

## Running the workflow

You need to have `conda` installed to create and activate a new environment.

```bash
conda env create -n pmc-metadata --file conda-environment.yaml
conda activate pmc-metadata
```

You also need to get an access token from [Zenodo](https://zenodo.org) and set it to the following
two environment variables: 

```shell
export SNAKEMAKE_STORAGE_ZENODO_ACCESS_TOKEN=<Your Zenodo access token>
export SNAKEMAKE_STORAGE_ZENODO_RESTRICTED_ACCESS_TOKEN=<Your Zenodo access token>
```

Run with `-–keep-storage-local-copies` to avoid downloading resources over and over again.
Also run with `--software-deployment-method conda` to use global conda packages.

```shell
snakemake --keep-storage-local-copies --software-deployment-method conda -c <number-of-cores-to-use>
```

You can use [`run.sh`](run.sh) to run the workflow this way, and with 12 cores.

# Citation

If you use this workflow in your work, please cite it using the metadata provided in [`CITATION.cff`](CITATION.cff).

# License

This work is licensed as specific in the [REUSE 3.0 Specification](https://reuse.software/spec/). 
Please consult the single file licenses.
