import os
from datetime import datetime
import json
import logging

from lxml import etree

log = logging.getLogger(__name__)
file_handler = logging.FileHandler(str(snakemake.log), mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)


def get_earliest_date(_full_dates: list[datetime], _month_dates: list[datetime], _year_dates: list[datetime]) -> datetime:
    """
    Determines the earliest date from a set of lists of datetime objects.
    The lists are ordered by the level of heuristic determination (or reliability) of dates.

    :param _full_dates: A list of datetime objects for which all parameters (year, month, day) were available
    :param _month_dates: A list of datetime objects for which 2/3 parameters (year, month) were available
    :param _year_dates: A list of datetime objects for which only year was available
    :return: the earliest date from the most reliable list, i.e., the list which contained the most complete
    data that could be determined without using fallback values.
    """
    if _full_dates:
        return min(_full_dates)
    elif _month_dates:
        return min(_month_dates)
    else:
        return min(_year_dates)


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
    lut = {}
    for _file in os.listdir(in_dir):
        if _file.endswith(".xml"):
            # Parse XML
            parser = etree.XMLParser(encoding='utf-8', recover=True)
            xml_tree = etree.parse(os.path.join(in_dir, _file), parser=parser)
            root = xml_tree.getroot()

            pub_dates = root.findall("front/article-meta/pub-date")
            full_dates = []
            month_dates = []
            year_dates = []
            for pub_date in pub_dates:
                # year = month = day = None
                # a = pub_date.find(".//year")
                a = pub_date.find(".//year")
                year = int(a.text) if a is not None and a.text is not None else None
                m = pub_date.find(".//month")
                month = int(m.text) if m is not None and m.text is not None else None
                d = pub_date.find(".//day")
                day = int(d.text) if d is not None and d.text is not None else None
                if day is not None:
                    full_dates.append(datetime(year, month, day))
                elif month is not None:
                    month_dates.append(datetime(year, month, 1))
                else:
                    year_dates.append(datetime(year, 1, 1))
                if day is not None:
                    full_dates.append(datetime(year, month, day))
                elif month is not None:
                    month_dates.append(datetime(year, month, 1))
                elif year is not None:
                    year_dates.append(datetime(year, 1, 1))
                else:
                    log.error(f"Could not find date data in {_file}: {etree.tostring(pub_date, pretty_print=True).decode()}.")
            earliest_date = get_earliest_date(full_dates, month_dates, year_dates)
            if earliest_date:
                lut[_file.rstrip("lmx.")] = earliest_date.strftime("%Y-%m-%d")
            else:
                log.error(f"No date found for {_file.rstrip('lmx.')}.")
    return lut


if __name__ == '__main__':
    lut = construct_lut(str(snakemake.input))
    if lut:
        with open(str(snakemake.output), "w") as out:
            json.dump(lut, out)
