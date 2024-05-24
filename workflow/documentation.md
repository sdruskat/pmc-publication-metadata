<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

# Documentation

This workflow takes as input FTP downloads of PMC publication metadata in JATS XML as included in PMC OA bulk downloads.
The output is a `.tar.gz` archive containing one JSON lookup table for each first digit of all encountered identifiers,
named `PMC<first digit>.json`, e.g., `PMC1.json` (which contains all identifiers starting with the digit `1`).

## Steps ("Rules")

# FIXME Update

1. The baseline tar.gz archives for the [configured baseline date](../config/README.md) are downloaded from the
PMC FTP server, and extracted to give access to the JATS XML files containing the publication metadata.
2. From each XML file for a PMC publication, the following metadata is extracted and written into a lookup table that 
maps the PMC publication identifier (e.g., `PMC5123456`) to the value for the version:
    - Publication date in `YYYY[-MM[-DD]]` format.
Each lookup table is saved to a JSON file.
3. The lookup tables are written into JSON files names `PMC<zero-padded PMC identifier prefix>.json`.
4. All lookup tables are merged into a single file `pmc-publication-dates.json`.
5. An attempt is made to patch missing metadata into the merged lookup table.
    - For all PMC identifiers that are present in the *Extract-URLs* dataset [^3], it is checked if they are 
      in the lookup table, attempts are made to retrieve missing metadata from the single OAI-PMH record of
      the publication. For identifiers that are still missing, results from a manual check are patched into
      the table.
6. The patched merged lookup table is split into files. For each digit for which PMC identifiers exist,
one file is written containing the PMC identifier records starting with this digit.
E.g., `PMC1.json` contains the data for all PMC publications whose identifier starts with `1`.
7. The files containing the LUT for a specific PMC identifier prefix are archived in a tarball that is saved in
`results/`

A graphical overview  of the rules is given below:

![DAG of the rules described above, generated via `snakemake --dag | dot -Tsvg > dag.svg` run in the repository root.](../dag.svg)
