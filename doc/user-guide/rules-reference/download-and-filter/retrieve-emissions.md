<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Rule `retrieve_emissions`

The `retrieve_emissions` rule downloads the EDGAR v6.0 CO2 emissions dataset used by the workflow to automatically derive electricity-sector CO2 limits.

EDGAR archive was created by the Joint Research Centre and is downloaded as zip file and placed under `data/co2_emissions/`. After that it is unpacked into the Excel workbook used as an input by `build_co2_emissions`.

This rule is active when `electricity.automatic_emission` is enabled.

## Outputs

- `data/co2_emissions/v60_GHG_CO2_excl_short-cycle_org_C_1970_2018.zip`
- `data/co2_emissions/v60_CO2_excl_short-cycle_org_C_1970_2018.xls`

## Downstream Rules

- [`build_co2_emissions`](../populate/build-co2-emissions.md) filters the retrieved CO2 dataset to electricity sector emissions and stores the cleaned dataset as a csv file used by `prepare_network`.
