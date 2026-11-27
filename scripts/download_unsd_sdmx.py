# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
# SPDX-License-Identifier: AGPL-3.0-or-later

"""Download and cache country-specific UNSD energy statistics from SDMX."""

import argparse
import csv
import io
import logging
import os
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET

import country_converter as coco

logger = logging.getLogger(__name__)

BASE_URL = "https://data.un.org/legacy/ws/rest"
USER_AGENT = "PyPSA-Earth/UNSD-demand"

REQUIRED_COLUMNS = {
    "REF_AREA",
    "COMMODITY",
    "TRANSACTION",
    "TIME_PERIOD",
    "OBS_VALUE",
    "UNIT_MULT",
    "UNIT_MEASURE",
}


def country_to_unsd_area(country):
    """Convert an ISO 3166 alpha-2 code to its UN M49 numeric code."""
    # UNSD Energy Statistics uses specific geographic coverage codes
    # for some countries, differing from standard UN M49.
    unsd_area_overrides = {
        "FR": "251",  # France including Monaco
        "IT": "382",  # Italy including San Marino and Vatican City
    }

    country = country.upper()

    if country in unsd_area_overrides:
        return unsd_area_overrides[country]

    code = coco.convert(
        names=country.upper(),
        src="ISO2",
        to="UNcode",
        not_found="not found",
    )

    if code == "not found":
        raise ValueError(f"Unsupported ISO2 country code: {country}")

    return f"{int(code):03d}"


def validate_content(content, kind):
    """Reject incomplete or unexpected UNSD responses."""
    if not content:
        raise ValueError("Empty UNSD response")

    if kind == "structure":
        root = ET.fromstring(content)

        ns = {"s": ("http://www.sdmx.org/resources/" "sdmxml/schemas/v2_1/structure")}

        if root.find(".//s:Codelist[@id='CL_TRANSACTION_NRG']", ns) is None:
            raise ValueError("Missing UNSD transaction classification")

    else:
        first_line = content.decode("utf-8-sig").splitlines()[0]
        columns = set(next(csv.reader(io.StringIO(first_line))))

        missing = REQUIRED_COLUMNS - columns

        if missing:
            raise ValueError(f"Unexpected UNSD CSV columns: {sorted(missing)}")


def fetch(url, destination, kind, refresh=False):
    """Download atomically, reusing valid cached files."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not refresh:
        try:
            validate_content(destination.read_bytes(), kind)
            logger.info("Using cached UNSD file: %s", destination)
            return destination
        except (ValueError, ET.ParseError, UnicodeError):
            logger.warning("Invalid UNSD cache: %s", destination)

    headers = {"User-Agent": USER_AGENT}

    if kind == "data":
        headers["Accept"] = "text/csv"

    request = Request(url, headers=headers)

    last_error = None

    for attempt in range(3):
        try:
            with urlopen(request, timeout=180) as response:
                content = response.read()

            validate_content(content, kind)

            temporary = destination.with_suffix(destination.suffix + ".tmp")

            try:
                temporary.write_bytes(content)
                os.replace(temporary, destination)
            finally:
                temporary.unlink(missing_ok=True)

            logger.info(
                "Downloaded UNSD file: %s (%s bytes)",
                destination,
                len(content),
            )

            return destination

        except (HTTPError, URLError, TimeoutError, ValueError) as error:
            last_error = error

            if isinstance(error, HTTPError) and error.code in (400, 403, 404):
                break

            if attempt < 2:
                time.sleep(2 ** (attempt + 1))

    raise RuntimeError(f"Failed to retrieve UNSD data from {url}") from last_error


def download_structure(directory, refresh=False):
    """Download the official UNSD SDMX data structure."""
    directory = Path(directory)

    return fetch(
        (f"{BASE_URL}/datastructure/UNSD/DSD_ENERGY/1.2" "?references=children"),
        directory / "structure.xml",
        "structure",
        refresh,
    )


def download_country(code, directory, refresh=False):
    """Download UNSD energy observations for one UN M49 country code."""
    code = str(code)

    if len(code) != 3 or not code.isdigit():
        raise ValueError(f"Invalid UN M49 country code: {code}")

    directory = Path(directory)

    return fetch(
        f"{BASE_URL}/data/UNSD,DF_UNDATA_ENERGY,1.2/A.{code}..",
        directory / f"{code}.csv",
        "data",
        refresh,
    )


def download(countries, directory, refresh=False):
    """Download the structure and observations for selected countries."""
    structure = download_structure(directory, refresh)

    files = {}

    for country in sorted(set(countries)):
        files[country] = download_country(
            country_to_unsd_area(country),
            directory,
            refresh,
        )

    return structure, files


def main():
    parser = argparse.ArgumentParser()

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--countries", nargs="+")
    mode.add_argument("--country-code")
    mode.add_argument("--structure-only", action="store_true")

    parser.add_argument(
        "--output-dir",
        default="data/unsd_sdmx",
    )
    parser.add_argument("--refresh", action="store_true")

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if args.structure_only:
        download_structure(args.output_dir, args.refresh)

    elif args.country_code is not None:
        download_country(
            args.country_code,
            args.output_dir,
            args.refresh,
        )

    else:
        download(
            args.countries,
            args.output_dir,
            args.refresh,
        )


if __name__ == "__main__":
    main()
