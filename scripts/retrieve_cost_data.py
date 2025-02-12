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
    future_years_list = list(snakemake.params.costs["year"])
    url_technology_data = snakemake.params.costs["technology_data_url"]
    version_technology_data = snakemake.params.costs["version"]
    country_list_technology_data = list(
        snakemake.params.costs["technology_data_countries"]
    )
    base_url_to_use = "/".join(
        (url_technology_data, version_technology_data, "outputs")
    )
    base_output_path = pathlib.Path(BASE_DIR, snakemake.output.output_path)
    base_output_path.mkdir(parents=True, exist_ok=True)

    for country_val in model_country_list:
        for year_val in future_years_list:
            cost_file_name = f"costs_{year_val}.csv"
            if country_val in country_list_technology_data:
                base_output_path_with_country = pathlib.Path(
                    base_output_path, str(country_val)
                )
                base_output_path_with_country.mkdir(parents=True, exist_ok=True)
                output_path = pathlib.Path(
                    base_output_path_with_country, cost_file_name
                )
                url_to_use = "/".join(
                    (base_url_to_use, str(country_val), cost_file_name)
                )
            else:
                url_to_use = "/".join((base_url_to_use, cost_file_name))
                output_path = pathlib.Path(base_output_path, cost_file_name)
            response = requests.get(url_to_use)
            with open(output_path, mode="wb") as output_cost_file:
                output_cost_file.write(response.content)
            output_cost_file.close()
