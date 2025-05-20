# -*- coding: utf-8 -*-
# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import pathlib

import requests
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
)

logger = create_logger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_cost_data")

    configure_logging(snakemake)
    model_country_list = list(snakemake.params.countries)
    costs_year = snakemake.params.costs["year"]
    url_technology_data = snakemake.params.costs["technology_data_url"]
    version_technology_data = snakemake.params.costs["version"]
    base_url_to_use = "/".join(
        (url_technology_data, version_technology_data, "outputs")
    )
    base_output_path = pathlib.Path(BASE_DIR, snakemake.output.output_path)
    cost_file_name = f"costs_{costs_year}.csv"
    base_output_path.mkdir(parents=True, exist_ok=True)

    if len(model_country_list) == 1 and "US" in model_country_list and snakemake.params.costs["technology_data_US"]:
        # take technology-data/outputs/US/costs_yyyy.csv
        output_path = pathlib.Path(base_output_path, cost_file_name)
        url_to_use = "/".join((base_url_to_use, "US", cost_file_name))
    else:
        # take technology-data/outputs/costs_yyyy.csv
        url_to_use = "/".join((base_url_to_use, cost_file_name))
        output_path = pathlib.Path(base_output_path, cost_file_name)

    response = requests.get(url_to_use)
    with open(output_path, mode="wb") as output_cost_file:
        output_cost_file.write(response.content)
    output_cost_file.close()
