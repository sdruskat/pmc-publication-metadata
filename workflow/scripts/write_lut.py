# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: MIT

import os
from datetime import datetime, MAXYEAR
import json
import logging
import calendar

from lxml import etree

log = logging.getLogger(__name__)
file_handler = logging.FileHandler(str(snakemake.log), mode="w")
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


def construct_lut(xml_dir: str) -> dict[str, str]:
    """
    Constructs a lookup table for publication dates for PMC publications.
    Takes a directory with PMC metadata files in JATS XML, retrieves the publication
    dates from the files, and writes the earliest publication date to a lookup table that is returned.

    :param xml_dir: A string representing a path to a directory containing PMC metadata XML files
    :return: A lookup table with the earliest recorded publication date for each PMC document identifier
    for which metadata is provided in an XML file in the input folder
    """
    in_dir = xml_dir
    _lut = {}
    for _file in os.listdir(in_dir):
        if _file.endswith(".xml"):
            # Parse XML
            parser = etree.XMLParser(encoding='utf-8', recover=True)
            xml_tree = etree.parse(os.path.join(in_dir, _file), parser=parser)
            root = xml_tree.getroot()

            pub_dates = root.findall("front/article-meta/pub-date")
            extracted_dates = []
            for pub_date in pub_dates:
                a = pub_date.find(".//year")
                year = a.text if a is not None and a.text is not None else None
                m = pub_date.find(".//month")
                month = m.text if m is not None and m.text is not None else None
                d = pub_date.find(".//day")
                day = str(d.text) if d is not None and d.text is not None else None
                extracted_dates.append(Date(year, month, day))
            earliest_date = get_earliest_date_str(extracted_dates)
            if earliest_date is not None:
                _lut[_file.rstrip("lmx.")] = earliest_date
            else:
                log.error(f"Could not determine an earliest publication date for {_file.rstrip('lmx.')}.")
    return _lut


if __name__ == '__main__':
    lut = construct_lut(str(snakemake.input))
    if lut:
        with open(str(snakemake.output), "w") as out:
            json.dump(lut, out)
