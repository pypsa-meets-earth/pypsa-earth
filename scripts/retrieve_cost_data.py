# -*- coding: utf-8 -*-
# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import requests
import pathlib
from _helpers import (
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
    country_list_technology_data = list(snakemake.params.costs["technology_data_countries"])
    base_url_to_use = "/".join((url_technology_data, version_technology_data, "outputs"))
    base_output_path = snakemake.params.output_directory

    for country_val in model_country_list:
        for year_val in future_years_list:
            cost_file_name = f"costs_{year_val}.csv"
            if country_val in country_list_technology_data:
                output_path = pathlib.Path(base_output_path, str(country_val), cost_file_name)
                url_to_use = "/".join((base_url_to_use, str(country_val), cost_file_name))
            else:
                url_to_use = "/".join((base_url_to_use, cost_file_name))
                output_path = pathlib.Path(base_output_path, cost_file_name)
            response = requests.get(url_to_use)
            with open(output_path, mode="wb") as output_cost_file:
                output_cost_file.write(response.content)
            output_cost_file.close()

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
