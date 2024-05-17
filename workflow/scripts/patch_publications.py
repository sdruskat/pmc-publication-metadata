# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: MIT

import logging
import json
from datetime import datetime, MAXYEAR
import xml.etree.ElementTree as ET
import calendar

import xmlschema
from sickle import Sickle


log = logging.getLogger(__name__)

file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)


class Date:
    """
    Represents a date that can be expressed as a datetime object and a string.
    """

    def __init__(self, year: str, month: str = None, day: str = None):
        self.year = year
        self.month = month
        self.day = day

    def string(self) -> str |None:
        """
        Returns the string representation of the Date, with the pattern
        "<year>[-<month>[-<day>]]".

        :return: the string representation of the Date
        """
        if self.year is not None:
            m_str = f"-{self.month.zfill(2)}" if self.month is not None else ""
            d_str = f"-{self.day.zfill(2)}" if self.month is not None and self.day is not None else ""
            return self.year + m_str + d_str
        else:
            return None

    def datetime(self) -> datetime | None:
        """
        Returns the precise, or the latest possible datetime for the Date,
        where "precise" is the correct datetime for the object when all of
        year, month, day are given, and "latest possible" is the last
        day of the month in a year where month and year are given, and the
        last day of the year when only a year is given.

        :return: the latest datetime object of the Date
        """
        if self.year is not None:
            y = int(self.year)
            m = int(self.month) if self.month is not None else 12
            d = int(self.day) if self.day is not None else calendar.monthrange(y, m)[1]
            return datetime(y, m, d)
        else:
            return None


def get_earliest_date_str(dates: list[Date]) -> str | None:
    """
    Determines the earliest date from a set of string tuples representing the year, month,
    and day parts of a date.

    :param dates: A list of string tuples representing dates
    :return: the string representation iof the earliest date in the given dates
    """
    earliest_date = Date(str(MAXYEAR))

    for date in dates:
        dt = date.datetime()
        if dt is not None:
            if dt < earliest_date.datetime():
                earliest_date = date

    if earliest_date.datetime().year < MAXYEAR:
        return earliest_date.string()
    else:
        return None


def extract_lut_from_xml(xml_source: str, xsd_file: str):
    xsd = xmlschema.XMLSchema(xsd_file)

    prefix = "{http://arxiv.org/OAI/arXivRaw/}"
    oai_prefix = "{http://www.openarchives.org/OAI/2.0/}"

    _lut = dict()

    tree = ET.fromstring(xml_source)

    metadata = tree.find(f"{oai_prefix}GetRecord/{oai_prefix}record/{oai_prefix}metadata")
    if metadata:
        arxiv_data = metadata.find(f"{prefix}arXivRaw")
        data, errors = xsd.to_dict(arxiv_data, validation="lax")
        if len(data) > 0:
            arxiv_id = data[f"{prefix}id"]
            versions = data[f"{prefix}version"]
            if not arxiv_id or not versions:
                log.warning(f"Metadata did not contain both, versions and ID, in {xml_source}.")
            else:
                for version in versions:
                    v = version["@version"]
                    date = version[f"{prefix}date"]
                    if date:
                        date_obj = datetime.strptime(date, "%a, %d %b %Y %H:%M:%S %Z")
                        date_f = datetime.strftime(date_obj, "%Y-%m-%d")
                        _lut[f"{arxiv_id}{v}"] = date_f
                        log.debug(f"Recorded {date_f} for {arxiv_id}{v}")

    return _lut


def patch_manually(lut_data):
    """
    Patch manually missing dates that could not be retrieved from the
    PMC OAI-PMH interface.

    :param lut_data: The LUT to patch
    :return: The patched LUT
    """
    patch_data = {}

    return lut_data | patch_data


def assert_versions() -> list[str]:
    """
    Test if all IDs in all Extract-URLs ArXiv JSONs are in the compiled LUT.
    Return missing IDs.

    Must only run once every 3 seconds.
    """

    with open(snakemake.input.lut, "r") as lut:
        lut_data = json.load(lut)

    missing = []

    for jf in snakemake.input.pmc_urls:
        jfs = jf.split("/")[-1]
        if not "test" in jfs:
            file_ym = jfs.split(".")[0]
            with open(jf, "r") as ji:
                url_data = json.load(ji)
            ym_data = url_data[file_ym]
            files = ym_data["files"]
            for pdf_name in files.keys():
                pmc_id = pdf_name.split(".")[-2]
                if not pmc_id in lut_data.keys():
                    missing.append(pmc_id)

    sickle = Sickle("https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi")

    for missing_id in missing:
        numerical_id = missing_id.lstrip("CMP")
        record = sickle.GetRecord(identifier=f"oai:pubmedcentral.nih.gov:{numerical_id}", metadataPrefix="oai_dc")
        if record is not None:
            dates = record.metadata["date"]
            _dates = []
            for date in dates:
                d = date.split("-")
                month = d[1] if len(d) > 1 else None
                day = d[2] if len(d) > 2 else None
                _dates.append(Date(d[0], month, day))

            lut_data = lut_data | {missing_id: get_earliest_date_str(_dates)}

    lut_data = patch_manually(lut_data)

    return lut_data

if __name__ == '__main__':
    patched_lut = assert_versions()
    with open(snakemake.output[0], "w") as out:
        json.dump(patched_lut, out)