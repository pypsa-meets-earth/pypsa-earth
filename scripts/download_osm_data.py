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

- ``data/osm/<date>/pbf``: Raw OpenStreetMap data as .pbf files per country
  - ``<date>`` is formatted as YYYYMM/ for historical data or "latest/" for current data
- ``data/osm/<date>/power``: Filtered power data as .json files per country
  - ``<date>`` is formatted as YYYYMM/ for historical data or "latest/" for current data
- ``data/osm/<date>/out``:  Prepared power data as .geojson and .csv files per country
  - ``<date>`` is formatted as YYYYMM/ for historical data or "latest/" for current data
- ``resources/osm/raw``: Prepared and per type (e.g. cable/lines) aggregated power data as .geojson and .csv files
"""
import inspect
import os
import shutil
from datetime import datetime
from pathlib import Path

from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
    read_osm_config,
    two_digits_2_name_country,
)
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
    country_list = country_list_to_geofk(snakemake.params.countries)

    # Get unified OSM data configuration
    osm_config = snakemake.config.get("osm_data", {})
    source = osm_config.get("source", "latest").lower()

    # Validate source option
    valid_sources = ["latest", "historical", "custom"]
    if source not in valid_sources:
        logger.warning(
            f"Invalid osm_data source '{source}'. Must be one of {valid_sources}. Defaulting to 'latest'."
        )
        source = "latest"

    # Extract configuration based on source
    custom_path_config = osm_config.get("custom_path", {})
    target_date = str(osm_config.get("target_date", None))

    # Handle custom_path - dict with pbf/power keys
    if isinstance(custom_path_config, dict):
        custom_pbf_path = custom_path_config.get("pbf", "data/custom/osm/pbf")
        custom_power_path = custom_path_config.get("power", "data/custom/osm/power")
    else:
        # Default paths
        custom_pbf_path = "data/custom/osm/pbf"
        custom_power_path = "data/custom/osm/power"

    # Create date-specific subdirectory for OSM data if historical date is provided
    osm_subdir = ""
    if source == "historical" and target_date:
        if isinstance(target_date, str):
            # Extract YYYYMM from date string
            osm_subdir = target_date.replace("-", "")[:6]  # Format: YYYYMM
        elif isinstance(target_date, datetime):
            osm_subdir = target_date.strftime("%Y%m")
    elif source == "custom":
        osm_subdir = "custom"
    else:
        osm_subdir = "latest"

    store_path_resources = Path.joinpath(
        Path(BASE_DIR), "resources", RDIR, "osm", "raw"
    )
    store_path_data = Path.joinpath(Path(BASE_DIR), "data", "osm", osm_subdir)

    # Parse target_date if provided as string (for historical source)
    if source == "historical" and target_date:
        if isinstance(target_date, str):
            try:
                target_date = datetime.strptime(target_date, "%Y-%m-%d")
                logger.info(
                    f"Historical OSM data mode: Using data from {target_date.strftime('%Y-%m-%d')}"
                )
                # Update osm_subdir after parsing
                osm_subdir = target_date.strftime("%Y%m")
                store_path_data = Path.joinpath(
                    Path(BASE_DIR), "data", "osm", osm_subdir
                )
            except ValueError:
                logger.warning(
                    f"Invalid date format '{target_date}', expected YYYY-MM-DD. Falling back to latest data download."
                )
                source = "latest"
                target_date = None
                osm_subdir = "latest"
                store_path_data = Path.joinpath(
                    Path(BASE_DIR), "data", "osm", osm_subdir
                )
        else:
            logger.info(
                f"Historical OSM data mode: Using data from {target_date.strftime('%Y-%m-%d')}"
            )
    elif source == "historical" and not target_date:
        logger.warning(
            "Historical source selected but no target_date provided. Falling back to latest data download."
        )
        source = "latest"

    # Log the data storage path
    logger.info(f"OSM data will be stored in: {store_path_data}")

    # Process OSM data based on selected source
    if source == "custom":
        # Custom data mode: copy user-provided OSM files then process them
        # PBF path is required, power path is optional

        # Check if pbf path exists (required)
        if not os.path.exists(custom_pbf_path):
            raise FileNotFoundError(
                f"Custom OSM pbf data path '{custom_pbf_path}' does not exist. "
                f"This path is required when using source: 'custom'. "
                f"Please provide valid pbf files or change osm_data.source to 'latest'."
            )

        # Validate that required pbf files exist for the specified countries
        # country_list contains ISO codes or special region codes
        # Earth-osm converts these to lowercase country names for pbf files
        missing_files = []
        for country_code in country_list:
            # Convert ISO code to full country name, then lowercase for geofabrik format
            # e.g., "BO" -> "Bolivia" -> "bolivia"
            # Special regions like "SN-GM" are handled by the helper function
            try:
                country_name = two_digits_2_name_country(country_code).lower()
            except:
                # If conversion fails, use the code itself lowercased
                country_name = country_code.lower()

            expected_file = f"{country_name}-latest.osm.pbf"
            expected_path = os.path.join(custom_pbf_path, expected_file)
            if not os.path.exists(expected_path):
                missing_files.append(expected_file)

        if missing_files:
            raise FileNotFoundError(
                f"Custom OSM data mode requires pbf files for all specified countries. "
                f"Missing files in '{custom_pbf_path}': {', '.join(missing_files)}. "
                f"Expected files should be named '<country>-latest.osm.pbf' (e.g., 'bolivia-latest.osm.pbf'). "
                f"Please provide the required pbf files or change osm_data.source to 'latest'."
            )

        # Create the data directories
        os.makedirs(store_path_data, exist_ok=True)
        pbf_target_path = Path.joinpath(store_path_data, "pbf")
        power_target_path = Path.joinpath(store_path_data, "power")

        # Copy pbf files (required)
        logger.info(
            f"Custom OSM data mode: Copying pbf files from '{custom_pbf_path}' to '{pbf_target_path}'"
        )
        shutil.copytree(custom_pbf_path, pbf_target_path, dirs_exist_ok=True)

        # Copy power files if path exists (optional)
        if os.path.exists(custom_power_path):
            logger.info(
                f"Custom OSM data mode: Copying power files from '{custom_power_path}' to '{power_target_path}'"
            )
            shutil.copytree(custom_power_path, power_target_path, dirs_exist_ok=True)
        else:
            logger.info(
                f"Custom power files path '{custom_power_path}' not found. "
                f"OSM data will be processed from pbf files only."
            )

    # Prepare save_osm_data arguments for all sources
    # earth-osm will skip downloading if pbf files already exist (custom mode)
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

    # Add historical date support if source is 'historical' and earth-osm supports it
    if source == "historical" and target_date:
        if "target_date" in inspect.signature(eo.save_osm_data).parameters:
            save_args["target_date"] = target_date
            logger.info(
                f"Downloading historical OSM data for {target_date.strftime('%Y-%m-%d')}"
            )
        else:
            logger.warning(
                f"Historical date requested but earth-osm version doesn't support target_date parameter. "
                f"Downloading the latest OSM data instead."
            )
    elif source == "latest":
        logger.info("Downloading latest OSM data")
    elif source == "custom":
        logger.info(
            "Processing custom OSM data (earth-osm will use existing pbf files)"
        )

    # Process OSM data (download or process existing pbf files)
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
