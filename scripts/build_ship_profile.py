# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

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

    landing = export_volume / ship_capacity  # fraction of max delivery
    pause_time = 8760 / landing - (fill_time + travel_time)
    full_cycle = fill_time + travel_time + unload_time + pause_time

    max_transport = ship_capacity * 8760 / (fill_time + travel_time + unload_time)
    print(f"The maximum transport capacity per ship is {max_transport:.2f} TWh/year")

    # throw error if max_transport < export_volume
    if max_transport < export_volume:
        ships = np.ceil(export_volume / max_transport)
        print(f"Number of ships needed to export {export_volume} TWh/year is {ships}")
        logger.info(
            "Not enough ship capacity to export all hydrogen in one ship. Extending the number of shipts to {}".format(
                ships
            )
        )

    # Set fill_time ->  1 and travel_time, unload_time, pause_time -> 0
    ship = pd.Series(
        [1.0] * fill_time + [0.0] * int(travel_time + unload_time + pause_time)
    )  # , index)
    ship.name = "profile"
    ship = pd.concat(
        [ship] * 1000, ignore_index=True
    )  # extend ship series to above 8760 hours

    # Add index, cut profile after length of snapshots
    snapshots = pd.date_range(freq="h", **snakemake.params.snapshots)
    ship = ship[: len(snapshots)]
    ship.index = snapshots

    # Scale ship profile to export_volume
    export_profile = ship / ship.sum() * export_volume * 1e6  # in MWh

    # Check profile
    if abs(export_profile.sum() / 1e6 - export_volume) > 0.001:
        raise ValueError(
            f"Sum of ship profile ({export_profile.sum()/1e6} TWh) does not match export demand ({export_volume} TWh)"
        )

    return export_profile


if __name__ == "__main__":
    if "snakemake" not in globals():

        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_ship_profile",
            h2export="120",
        )

    # Get parameters from config and wildcard
    ship_opts = snakemake.params.ship_opts
    export_volume = eval(snakemake.wildcards.h2export)

    # Create export profile
    if export_volume > 0:
        export_profile = build_ship_profile(export_volume, ship_opts)
    else:
        export_profile = pd.Series(
            0,
            index=pd.date_range(freq="h", **snakemake.params.snapshots),
            name="profile",
        )

    # Save export profile
    export_profile.to_csv(snakemake.output.ship_profile)  # , header=False)

    logger.info("Ship profile successfully created")
