# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import logging
import os
import shutil
from pathlib import Path

import pandas as pd
from _helpers import BASE_DIR, read_csv_nafix

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_transport_data_input")

    # configure_logging(snakemake)

    # run = snakemake.config.get("run", {})
    # RDIR = run["name"] + "/" if run.get("name") else ""
    # store_path_data = Path.joinpath(Path().cwd(), "data")
    # country_list = country_list_to_geofk(snakemake.config["countries"])'

    nbr_vehicles = read_csv_nafix(snakemake.input.n_vehicles_raw).copy()

    CO2_emissions = read_csv_nafix(snakemake.input.transport_emissions_raw).copy()

    if nbr_vehicles.empty or CO2_emissions.empty:
        # In case one of the urls is not working, we can use the hard-coded data
        src = BASE_DIR + "/data/temp_hard_coded/transport_data.csv"
        dest = snakemake.output.transport_data_input
        shutil.copy(src, dest)
    else:
        # Join the DataFrames by the 'country' column to prepare the tabular transport_data,
        # which will be saved as transport_data.csv in the resource folder and used
        # to prepare further (nodal) transport data in prepare_transport_data
        # and to scale the e-mob parameters in prepare_sector_network.
        transport = pd.merge(nbr_vehicles, CO2_emissions, on="country")
        transport = transport[["country", "number cars", "average fuel efficiency"]]

        missing = transport.index[transport["average fuel efficiency"].isna()]
        if not missing.empty:
            print(
                "Missing data on fuel efficiency from:\n",
                f"{list(transport.loc[missing].country)}.",
                "\nFilling gaps with averaged data.",
            )

            fill_value = transport["average fuel efficiency"].mean()
            transport.loc[missing, "average fuel efficiency"] = fill_value

        transport.loc[:, "average fuel efficiency"] = transport[
            "average fuel efficiency"
        ].round(3)

        transport.to_csv(
            snakemake.output.transport_data_input,
            sep=",",
            encoding="utf-8",
            header="true",
            index=False,
        )
