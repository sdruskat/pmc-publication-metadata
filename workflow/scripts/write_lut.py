import os
from lxml import etree

in_dir = snakemake.input

for _file in os.listdir(in_dir):
    if _file.endswith(".xml"):
        with open(os.path.join(in_dir, _file), "r") as xml_file:
            d = etree.parse(xml_file)
            e = d.getroot()
            e.xpath(".//abstract")
            e.xpath(".//abstract/p")[0].text

print(f"FULL INPUT: {snakemake.input}")
# for x in snakemake.input:
#     print(f"INPUT:; {x}")
print(f"OUTPUT: {snakemake.output}")
from pathlib import Path
Path(str(snakemake.output)).touch()