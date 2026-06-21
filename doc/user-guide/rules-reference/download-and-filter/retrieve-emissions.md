<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Rule `retrieve_emissions`

The `retrieve_emissions` rule downloads the EDGAR v6.0 fossil CO2 emissions dataset used to derive automatic electricity-sector CO2 limits.

It retrieves the EDGAR archive from the Joint Research Centre open-data server, stores the downloaded ZIP file under `data/co2_emissions/`, and unpacks the Excel workbook consumed by `build_co2_emissions`.

This rule is active when `electricity.automatic_emission` is enabled.

## Outputs

- `data/co2_emissions/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip`
- `data/co2_emissions/v60_CO2_excl_short-cycle_org_C_1970_2018.xls`

## Downstream Rules

- [`build_co2_emissions`](../populate/build-co2-emissions.md) filters the retrieved workbook to public electricity and heat production emissions and writes the cleaned CSV used by `prepare_network`.
