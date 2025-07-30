# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Append cost data by readjusting to the currency of the main cost data (EUR or USD)
"""
import pandas as pd
from _helpers import create_logger, get_yearly_currency_exchange_average, mock_snakemake

logger = create_logger(__name__)

rename_subtech = {
    "Coal average": "coal",
    "Oil average": "oil",
    "Gas Turbine": "OCGT",
    "Gas Combined Cycle": "CCGT",
    "Nuclear average": "nuclear",
    "Large Hydro": "hydro",
    "Small Hydro": "ror",
    "Geothermal average": "geothermal",
    "Solar average": "solar",
    "Solar PV": "solar-utility",
    "Solar PV Rooftop": "solar-rooftop",
    # "Solar CSP": "central solar thermal", # Different unit and workflow
    "Wind Onshore": "onwind",
    "Wind Offshore": "offwind",
    "Bioenergy average": "biomass",
    "Biogas": "biogas",
    "Lithium Ion Batteries": "battery inverter",  # Only impacts the charge discharge capacity, not storage capacity
    "Pumped Hydro": "Pumped-Storage-Hydro-bicharger",  # Only impacts the charge discharge capacity, not storage capacity
}


def get_cost_if_exist(df, param, col):
    answer = df.loc[df.parameter == param, col]
    return answer.iloc[0] if not answer.empty else "None"


if __name__ == "__main__":
    if "snakemake" not in globals():
        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake(
            "append_cost_data",
            planning_horizons="2030",
        )

    # Get wildcard year
    year = snakemake.wildcards.year

    # Load declining factors
    declining = pd.read_csv(snakemake.input.declining_factor)
    declining = declining.set_index(["Technology", "Parameter"])[year].unstack()

    regional_config = snakemake.params.regional_factor
    if regional_config:
        if isinstance(regional_config, bool):
            regional = pd.read_csv(snakemake.input.regional_factor)
            regional_factor = float(regional.loc[regional.index[0], year])
        elif isinstance(regional_config, (float, int)) and regional_config > 0:
            regional_factor = regional_config
        else:
            raise ValueError(
                "regional_factor must be a boolean or a number larger than 0"
            )

        logger.info(
            f"Applying regional factor {regional_factor} to all investment costs"
        )
    else:
        logger.info("No regional factor is set")
        regional_factor = 1

    # Load cost summary
    app_costs = pd.read_csv(snakemake.input.app_costs)

    # Readjust USD to EUR
    app_costs[["capital_cost", "FOM", "VOM"]] *= get_yearly_currency_exchange_average(
        "USD", "EUR", 2020
    )

    # Include declining factors
    for t in declining.index:
        app_costs.loc[app_costs.technology == t, "capital_cost"] *= (
            1 - declining.loc[t, "Capex"]
        )
        app_costs.loc[app_costs.technology == t, "VOM"] *= 1 - declining.loc[t, "Opex"]

    # Load base cost data for lifetime and FOM
    base_costs = pd.read_csv(snakemake.input.costs)

    discount_rate = snakemake.params.discount_rate[0]
    source = "8th ASEAN Energy Outlook (D.15 D.17 D.18), https://aseanenergy.org/wp-content/uploads/2024/09/8th-ASEAN-Energy-Outlook.pdf"

    # Step 1 & 2: Convert capital cost to investment using lifetime from cost_20xx.csv
    rows = []
    for _, row in app_costs.iterrows():
        subtech = row["sub-technology"]
        tech = rename_subtech.get(subtech)
        if tech is None:
            continue
        techs = tech if isinstance(tech, list) else [tech]
        for t in techs:
            row_base = base_costs.loc[base_costs.technology == t]
            FOM = row["FOM"] / row["capital_cost"] * 100
            # Append investement, FOM and VOM costs
            rows.append(
                {
                    "technology": t,
                    "parameter": "investment",
                    "value": row["capital_cost"],  # Thousand USD/MW = USD/kW
                    "unit": "EUR/kW",
                    "source": source,
                    "currency_year": 2020,
                }
            )
            rows.append(
                {
                    "technology": t,
                    "parameter": "FOM",
                    "value": FOM,
                    "unit": "%/year",
                    "source": source,
                    "currency_year": 2020,
                }
            )
            rows.append(
                {
                    "technology": t,
                    "parameter": "VOM",
                    "value": row["VOM"],
                    "unit": "EUR/MWh",
                    "source": source,
                    "currency_year": 2020,
                }
            )

            old_investement = get_cost_if_exist(row_base, "investment", "value")
            old_investement_unit = get_cost_if_exist(row_base, "investment", "unit")
            old_investement *= 1e-3 if "/MW" in old_investement_unit else 1

            old_FOM = get_cost_if_exist(row_base, "FOM", "value")
            old_VOM = get_cost_if_exist(row_base, "VOM", "value")

            df_log = pd.DataFrame(
                {
                    "Before": [old_investement, old_FOM, old_VOM],
                    "After": [row["capital_cost"] * regional_factor, FOM, row["VOM"]],
                    "Unit": ["EUR/kW", "%/year", "EUR/MWh"],
                },
                index=["investment", "FOM", "VOM"],
            )

            logger.info(f"Updated costs for {t}:\n{df_log}\n")

    new_costs = pd.DataFrame(rows)

    # --- Combine with base_costs ---
    # Drop overlapping parameters for the technologies covered by new_costs
    mask_to_drop = base_costs.parameter.isin(
        new_costs.parameter
    ) & base_costs.technology.isin(new_costs.technology)

    base_costs_filtered = base_costs.loc[~mask_to_drop]

    # Concatenate and adjust investment by regional factor
    costs = pd.concat([base_costs_filtered, new_costs], ignore_index=True)
    costs.loc[costs.parameter == "investment", "value"] *= regional_factor

    # Group by to aggregate duplicates if any
    costs = (
        costs.groupby(["technology", "parameter", "value", "unit"]).sum().reset_index()
    )

    # Save the costs
    costs.to_csv(snakemake.output[0], index=False)
