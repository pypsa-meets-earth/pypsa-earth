# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Python interface to download OpenStreetMap data
Documented at https://github.com/pypsa-meets-earth/earth-osm

Relevant Settings
-----------------

*None*  # multiprocessing & infrastructure selection can be an option in future

Inputs
------

*None*

Outputs
-------

- ``data/osm/pbf``: Raw OpenStreetMap data as .pbf files per country
- ``data/osm/power``: Filtered power data as .json files per country
- ``data/osm/out``:  Prepared power data as .geojson and .csv files per country
- ``resources/osm/raw``: Prepared and per type (e.g. cable/lines) aggregated power data as .geojson and .csv files
"""
import inspect
import os
import shutil
from datetime import datetime
from pathlib import Path

from _helpers import BASE_DIR, configure_logging, create_logger, read_osm_config
from earth_osm import eo

logger = create_logger(__name__)


def country_list_to_geofk(country_list):
    """
    Convert the requested country list into geofk norm.

    Parameters
    ----------
    input : str
        Any two-letter country name or aggregation of countries given in the regions config file
        Country name duplications won't distort the result.
        Examples are:
        ["NG","ZA"], downloading osm data for Nigeria and South Africa
        ["SNGM"], downloading data for Senegal&Gambia shape
        ["NG","ZA","NG"], won't distort result.

    Returns
    -------
    full_codes_list : list
        Example ["NG","ZA"]
    """
    full_codes_list = [convert_iso_to_geofk(c_code) for c_code in set(country_list)]

    return full_codes_list


def convert_iso_to_geofk(
    iso_code, iso_coding=True, convert_dict=read_osm_config("iso_to_geofk_dict")
):
    """
    Function to convert the iso code name of a country into the corresponding
    geofabrik In Geofabrik, some countries are aggregated, thus if a single
    country is requested, then all the agglomeration shall be downloaded For
    example, Senegal (SN) and Gambia (GM) cannot be found alone in geofabrik,
    but they can be downloaded as a whole SNGM.

    The conversion directory, initialized to iso_to_geofk_dict is used to perform such conversion
    When a two-letter code country is found in convert_dict, and iso_coding is enabled,
    then that two-letter code is converted into the corresponding value of the dictionary

    Parameters
    ----------
    iso_code : str
        Two-code country code to be converted
    iso_coding : bool
        When true, the iso to geofk is performed
    convert_dict : dict
        Dictionary used to apply the conversion iso to geofk
        The keys correspond to the countries iso codes that need a different region to be downloaded
    """
    if iso_coding and iso_code in convert_dict:
        return convert_dict[iso_code]
    else:
        return iso_code


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("download_osm_data")
    configure_logging(snakemake)

    run = snakemake.config.get("run", {})
    RDIR = run["name"] + "/" if run.get("name") else ""
    store_path_resources = Path.joinpath(
        Path(BASE_DIR), "resources", RDIR, "osm", "raw"
    )
    store_path_data = Path.joinpath(Path(BASE_DIR), "data", "osm")
    country_list = country_list_to_geofk(snakemake.params.countries)
    custom_data = snakemake.config.get("custom_data", {}).get("osm_data", {})
    set_custom_data = custom_data.get("set", False)
    custom_data_path = custom_data.get("custom_path", "data/custom/osm")

    # Check for historical date configuration
    historical_config = snakemake.config.get("historical_osm_data", {})
    target_date = historical_config.get("osm_date", None)

    # Parse target_date if provided as string
    if target_date and isinstance(target_date, str):
        try:
            target_date = datetime.strptime(target_date, "%Y-%m-%d")
            logger.info(
                f"Using historical OSM data for date: {target_date.strftime('%Y-%m-%d')}"
            )
        except ValueError:
            logger.warning(
                f"Invalid date format '{target_date}', expected YYYY-MM-DD. Using latest data."
            )
            target_date = None
    elif target_date:
        logger.info(
            f"Using historical OSM data for date: {target_date.strftime('%Y-%m-%d')}"
        )

    # allow for custom data usage
    if set_custom_data and not os.path.exists(custom_data_path):
        raise FileNotFoundError(
            f"Custom OSM data path {custom_data_path} does not exist. Please provide a valid path."
        )
    elif set_custom_data and os.path.exists(custom_data_path):
        logger.info(
            "Custom OSM data usage is activated. Skipping the download of OSM data."
        )
        os.makedirs(store_path_data, exist_ok=True)
        shutil.copytree(custom_data_path, store_path_data, dirs_exist_ok=True)
    else:
        # Prepare save_osm_data arguments
        save_args = {
            "primary_name": "power",
            "region_list": country_list,
            "feature_list": ["substation", "line", "cable", "generator"],
            "update": False,
            "mp": True,
            "data_dir": store_path_data,
            "out_dir": store_path_resources,
            "out_format": ["csv", "geojson"],
            "out_aggregate": True,
            "progress_bar": snakemake.config["enable"]["progress_bar"],
        }

        # Add historical date support if available and target_date is provided
        if (
            target_date
            and "target_date" in inspect.signature(eo.save_osm_data).parameters
        ):
            save_args["target_date"] = target_date
            logger.info("Historical data download enabled")
        elif target_date:
            logger.warning(
                "Historical date requested but earth-osm version doesn't support target_date parameter. Downloading the latest OSM data."
            )

        eo.save_osm_data(**save_args)

    out_path = Path.joinpath(store_path_resources, "out")
    names = ["generator", "cable", "line", "substation"]
    out_formats = ["csv", "geojson"]
    new_files = os.listdir(out_path)  # list downloaded osm files

    # earth-osm (eo) only outputs files with content
    # If the file is empty, it is not created
    # This is a workaround to create empty files for the workflow

    # Rename and move osm files to the resources folder output
    for name in names:
        for f in out_formats:
            new_file_name = Path.joinpath(store_path_resources, f"all_raw_{name}s.{f}")
            old_files = list(Path(out_path).glob(f"*{name}.{f}"))
            # if file is missing, create empty file, otherwise rename it an move it
            if not old_files:
                with open(new_file_name, "w") as f:
                    pass
            else:
                logger.info(f"Move {old_files[0]} to {new_file_name}")
                shutil.move(old_files[0], new_file_name)
