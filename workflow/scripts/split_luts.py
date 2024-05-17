# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: MIT

import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)

file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)


patched_lut = snakemake.input.patched_lut
out_dir = snakemake.output.outdir

luts = dict()

lut_map = {i: {} for i in range(1, 10)}

with open(patched_lut, "r") as lutf:
    data = json.load(lutf)

    data_len = len(data)

    for i, (k, v) in enumerate(data.items()):
        if i % 100000 == 0:
            log.debug(f"Splitting item {i} of {data_len} ({i / data_len * 100}%)")
        digit = k.lstrip("CMP")[0]
        lut_map[int(digit)][k] = v

od = Path(out_dir)
od.mkdir(parents=True, exist_ok=True)

for dgt, pmc_map in lut_map.items():
    log.debug(f"Writing LUT PMC{dgt}.")
    with open(f"{out_dir}/PMC{dgt}.json", "w") as out:
        json.dump(pmc_map, out)
