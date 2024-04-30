#!/bin/bash

# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: CC0-1.0

echo "Running Snakemake workflow with 12 cores, keeping local storage."
snakemake --keep-storage-local-copies --software-deployment-method conda -c 12
echo "Done."

