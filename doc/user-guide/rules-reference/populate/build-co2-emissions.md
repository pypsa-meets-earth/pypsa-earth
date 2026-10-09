<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Rule `build_co2_emissions`

The `build_co2_emissions` rule prepares EDGAR fossil CO2 emissions for the automatic electricity-sector CO2 limit.

It reads the EDGAR Excel workbook retrieved by `retrieve_emissions`, detects the fossil CO2 emissions sheet, filters the data to `Public electricity and heat production`, and writes a cleaned CSV with country identifiers and yearly emission values.

This rule is active when `electricity.automatic_emission` is enabled.

## Inputs

- `data/co2_emissions/v60_CO2_excl_short-cycle_org_C_1970_2018.xls`

## Outputs

- `resources/{run}/co2_emissions_elec_and_heat.csv`

The output contains:

- `country_code_a3`: three-letter country code
- `country_code_a2`: two-letter country code
- `country_name`: country name
- `Y_YYYY`: yearly CO2 emissions in kt CO2

## Script Documentation

::: scripts.build_co2_emissions
    options:
        show_root_heading: false
        show_source: false
