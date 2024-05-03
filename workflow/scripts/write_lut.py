import os
from datetime import datetime
import json

from lxml import etree

def earliest_date(_full_dates: list[datetime], _month_dates: list[datetime], _year_dates: list[datetime]) -> datetime:
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
                year = month = day = None
                a = pub_date.find(".//year")
                if a and a.text:
                    year = int(a.text)
                m = pub_date.find(".//month")
                if m and m.text:
                    month = int(m.text)
                d = pub_date.find(".//day")
                if d and d.text:
                    day = int(d.text)
                if day is not None:
                    full_dates.append(datetime(year, month, day))
                elif month is not None:
                    month_dates.append(datetime(year, month, 1))
                else:
                    year_dates.append(datetime(year, 1, 1))
            lut[_file.rstrip("lmx.")] = earliest_date(full_dates, month_dates, year_dates).strftime("%Y-%m-%d")
    return lut


if __name__ == '__main__':
    with open(str(snakemake.output), "w") as out:
        json.dump(construct_lut(str(snakemake.input)), out)
