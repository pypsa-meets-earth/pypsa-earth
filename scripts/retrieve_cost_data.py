# -*- coding: utf-8 -*-
# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import datetime as dt
import os
import re
from zipfile import ZipFile

import geopandas as gpd
import pandas as pd
import yaml
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_country_list,
    create_logger,
    progress_retrieve,
)
from google_drive_downloader import GoogleDriveDownloader as gdd
from tqdm import tqdm

logger = create_logger(__name__)


def load_databundle_config(config):
    "Load databundle configurations from path file or dictionary"

    if type(config) is str:
        with open(config) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)["databundles"]
    elif type(config) is not dict:
        logger.error("Impossible to load the databundle configuration")

    # parse the "countries" list specified in the file before processing
    for bundle_name in config:
        config[bundle_name]["countries"] = create_country_list(
            config[bundle_name]["countries"], iso_coding=False
        )

    return config


if __name__ == "__main__":
    if "snakemake" not in globals():

        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_cost_data.py")

    # TODO Make logging compatible with progressbar (see PR #102, PyPSA-Eur)
    configure_logging(snakemake)

    rootpath = "."
    tutorial = snakemake.params.tutorial
    countries = snakemake.params.countries
    logger.info(f"Retrieving data for {len(countries)} countries.")

    config_enable = snakemake.config["enable"]
    config_bundles = load_databundle_config(snakemake.config["databundles"])
    disable_progress = not config_enable["progress_bar"]

    logger.warning(
        "DISCLAIMER LICENSES: the use of PyPSA-Earth is conditioned \n \
        to the acceptance of its multiple licenses.\n \
        The use of the code automatically implies that you accept all the licenses.\n \
        See our documentation for more information. \n \
        Link: https://pypsa-earth.readthedocs.io/en/latest/introduction.html#licence"
    )

    """
    TODO:
    retrieve_cost_data.py shall:
    1) get the list of countries from config
    2) get list of years for which the cost assumptions shall be considered
    3) loop through countries e years and fetch the corresponding technology-data file
    4) create resources/RDIR/costs subfolder. For the US create a subfolder with the costs fetched from technology-data. 
    
    --> For any other country copy the cost files under resources/RDIR/costs
    """
