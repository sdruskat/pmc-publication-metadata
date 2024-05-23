# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: MIT

import json
import logging

merged_lut = dict()

log = logging.getLogger(__name__)
file_handler = logging.FileHandler(str(snakemake.log), mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)

def deep_merge_lut(exist: dict, new: dict) -> dict:
    return exist | new


for lut_file in snakemake.input:
    log.debug(f"Merging {lut_file}.")
    with open(lut_file, "r") as fi:
        data = json.load(fi)

    merged_lut = merged_lut | data

with open(snakemake.output[0], "w") as of:
    json.dump(merged_lut, of)
