<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

# Documentation

This workflow takes as input PMC publication metadata in JATS XML as included in PMC OA bulk downloads.

## Steps ("Rules")

1. The baseline tar.gz archives for the [configured baseline date](../config/README.md) are downloaded from the
PMC FTP server, and extracted to give access to the JATS XML files containing the publication metadata.
2. From each XML file for a PMC publication, the following metadata is extracted and written into a lookup table that 
maps the PMC publication identifier (e.g., `PMC5123456`) to the value for the version:
    - Publication date in `YYYY[-MM[-DD]]` format.
Each lookup table is saved to a JSON file.
3. The lookup tables are written into JSON files names `PMC<zero-padded PMC identifier prefix>.json`.
4. The files containing the LUT for a specific PMC identifier prefix are archived in a tarball that is saved in
`results/`

A graphical overview  of the rules is given below:

![Rulegraph of the rules described above, generated via `snakemake --rulegraph | dot -Tsvg > rulegraph.svg` run in the repository root.](../rulegraph.svg)
