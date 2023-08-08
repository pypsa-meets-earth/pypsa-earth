# -*- coding: utf-8 -*-


import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def build_ship_profile(export_volume, ship_opts):

    ship_capacity = ship_opts["ship_capacity"]
    travel_time = ship_opts["travel_time"]
    fill_time = ship_opts["fill_time"]
    unload_time = ship_opts["unload_time"]


    # Check profile

    return export_profile



if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake, sets_path_to_root

        snakemake = mock_snakemake(
            "build_ship_profile",
            h2export="9",
        )
        sets_path_to_root("pypsa-earth-sec")

    # Get parameters from config and wildcard
    ship_opts = snakemake.config["export"]["ship"]
    export_volume = eval(snakemake.wildcards.h2export)

    # Create export profile
    export_profile = build_ship_profile(export_volume, ship_opts)

    # Save export profile
    export_profile.to_csv(snakemake.output.ship_profile, header=False)

    logger.info("Ship profile successfully created")
