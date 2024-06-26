<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

# Workflow configuration

Configure this workflow by making settings in `config.yml`:

- `pmc_baseline_date`: Set the date string for the baseline files, e.g., "2023-12-18".
To see the latest available baseline files, see https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/xml/ or
https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/xml/.
- `subset`: Set the list of strings for PMC ID prefix subsets to be considered when running the workflow, e.g.
    - "001"
    - "002"
    - "003"
    - "004"
    - "005"
    - "006"
    - "007"
    - "008"
    - "009"
    - "010"
- `outdir` (_optional_): Set the output directory, e.g., "/scratch/$USER/somedir" 