# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Adds existing power and heat generation capacities for initial planning
horizon.
"""

import logging
import os
from types import SimpleNamespace

import country_converter as coco
import numpy as np
import pandas as pd
import powerplantmatching as pm
import pypsa
import xarray as xr
from _helpers import sanitize_carriers, sanitize_locations

# from _helpers import (
#     configure_logging,
#     set_scenario_config,
#     update_config_from_wildcards,
# )
# from add_electricity import sanitize_carriers
from prepare_sector_network import define_spatial, prepare_costs  # , cluster_heat_buses

logger = logging.getLogger(__name__)
cc = coco.CountryConverter()
idx = pd.IndexSlice
spatial = SimpleNamespace()


def add_build_year_to_new_assets(n: pypsa.Network, baseyear: int):
    """
    Assigns build_year to assets which have lifetimes and no build year, assuming they are
    built in the baseyear. This is necessary for the correct handling of existing assets.

    Parameters
    ----------
    n : pypsa.Network
    baseyear : int
        year in which optimized assets are built

    Returns
    -------
    None
    """
    # Set build_year for assets with lifetimes and no build year to the baseyear
    for c in n.iterate_components(["Link", "Generator", "Store"]):
        assets = c.df.index[(c.df.lifetime != np.inf) & (c.df.build_year == 0)]
        c.df.loc[assets, "build_year"] = baseyear

        # add -baseyear to name
        rename = pd.Series(c.df.index, c.df.index)
        rename[assets] += f"-{str(baseyear)}"
        c.df.rename(index=rename, inplace=True)

        # rename time-dependent
        selection = n.component_attrs[c.name].type.str.contains(
            "series"
        ) & n.component_attrs[c.name].status.str.contains("Input")
        for attr in n.component_attrs[c.name].index[selection]:
            c.pnl[attr] = c.pnl[attr].rename(columns=rename)


def add_heating_capacities_installed_before_baseyear(
    n,
    baseyear,
    grouping_years,
    ashp_cop,
    gshp_cop,
    time_dep_hp_cop,
    costs,
    default_lifetime,
):
    """
    Parameters
    ----------
    n : pypsa.Network
    baseyear : last year covered in the existing capacities database
    grouping_years : intervals to group existing capacities
        linear decommissioning of heating capacities from 2020 to 2045 is
        currently assumed heating capacities split between residential and
        services proportional to heating load in both 50% capacities
        in rural busess 50% in urban buses
    """
    logger.debug(f"Adding heating capacities installed before {baseyear}")

    existing_heating = pd.read_csv(
        snakemake.input.existing_heating_distribution, header=[0, 1], index_col=0
    )

    techs = existing_heating.columns.get_level_values(1).unique()

    for name in existing_heating.columns.get_level_values(0).unique():
        name_type = "central" if name == "urban central" else "decentral"

        nodes = pd.Index(n.buses.location[n.buses.index.str.contains(f"{name} heat")])

        if (name_type != "central") and options["electricity_distribution_grid"]:
            nodes_elec = nodes + " low voltage"
        else:
            nodes_elec = nodes

        heat_pump_type = "air" if "urban" in name else "ground"

        # Add heat pumps
        costs_name = f"decentral {heat_pump_type}-sourced heat pump"

        cop = {"air": ashp_cop, "ground": gshp_cop}

        if time_dep_hp_cop:
            efficiency = cop[heat_pump_type][nodes]
        else:
            efficiency = costs.at[costs_name, "efficiency"]

        for i, grouping_year in enumerate(grouping_years):
            if int(grouping_year) + default_lifetime <= int(baseyear):
                continue

            # installation is assumed to be linear for the past default_lifetime years
            ratio = (int(grouping_year) - int(grouping_years[i - 1])) / default_lifetime

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} {heat_pump_type} heat pump-{grouping_year}",
                bus0=nodes_elec,
                bus1=nodes + " " + name + " heat",
                carrier=f"{name} {heat_pump_type} heat pump",
                efficiency=efficiency,
                capital_cost=costs.at[costs_name, "efficiency"]
                * costs.at[costs_name, "fixed"],
                p_nom=existing_heating.loc[nodes, (name, f"{heat_pump_type} heat pump")]
                * ratio
                / costs.at[costs_name, "efficiency"],
                build_year=int(grouping_year),
                lifetime=costs.at[costs_name, "lifetime"],
            )

            # add resistive heater, gas boilers and oil boilers
            n.madd(
                "Link",
                nodes,
                suffix=f" {name} resistive heater-{grouping_year}",
                bus0=nodes_elec,
                bus1=nodes + " " + name + " heat",
                carrier=name + " resistive heater",
                efficiency=costs.at[f"{name_type} resistive heater", "efficiency"],
                capital_cost=(
                    costs.at[f"{name_type} resistive heater", "efficiency"]
                    * costs.at[f"{name_type} resistive heater", "fixed"]
                ),
                p_nom=(
                    existing_heating.loc[nodes, (name, "resistive heater")]
                    * ratio
                    / costs.at[f"{name_type} resistive heater", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} resistive heater", "lifetime"],
            )

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} gas boiler-{grouping_year}",
                bus0="EU gas" if "EU gas" in spatial.gas.nodes else nodes + " gas",
                bus1=nodes + " " + name + " heat",
                bus2="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[f"{name_type} gas boiler", "efficiency"],
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=(
                    costs.at[f"{name_type} gas boiler", "efficiency"]
                    * costs.at[f"{name_type} gas boiler", "fixed"]
                ),
                p_nom=(
                    existing_heating.loc[nodes, (name, "gas boiler")]
                    * ratio
                    / costs.at[f"{name_type} gas boiler", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} gas boiler", "lifetime"],
            )

            n.madd(
                "Link",
                nodes,
                suffix=f" {name} oil boiler-{grouping_year}",
                bus0=spatial.oil.nodes,
                bus1=nodes + " " + name + " heat",
                bus2="co2 atmosphere",
                carrier=name + " oil boiler",
                efficiency=costs.at["decentral oil boiler", "efficiency"],
                efficiency2=costs.at["oil", "CO2 intensity"],
                capital_cost=costs.at["decentral oil boiler", "efficiency"]
                * costs.at["decentral oil boiler", "fixed"],
                p_nom=(
                    existing_heating.loc[nodes, (name, "oil boiler")]
                    * ratio
                    / costs.at["decentral oil boiler", "efficiency"]
                ),
                build_year=int(grouping_year),
                lifetime=costs.at[f"{name_type} gas boiler", "lifetime"],
            )

            # delete links with p_nom=nan corresponding to extra nodes in country
            n.mremove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index and np.isnan(n.links.p_nom[index])
                ],
            )

            # delete links with capacities below threshold
            threshold = snakemake.params.existing_capacities["threshold_capacity"]
            n.mremove(
                "Link",
                [
                    index
                    for index in n.links.index.to_list()
                    if str(grouping_year) in index and n.links.p_nom[index] < threshold
                ],
            )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_existing_baseyear",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
            demand="AB",
            h2export="120",
        )

    # configure_logging(snakemake)
    # set_scenario_config(snakemake)

    # update_config_from_wildcards(snakemake.config, snakemake.wildcards)

    options = snakemake.params.sector

    baseyear = snakemake.params.baseyear

    n = pypsa.Network(snakemake.input.network)

    # define spatial resolution of carriers
    spatial = define_spatial(n.buses[n.buses.carrier == "AC"].index, options)
    add_build_year_to_new_assets(n, baseyear)

    Nyears = n.snapshot_weightings.generators.sum() / 8760.0
    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"],
        snakemake.params.costs["output_currency"],
        snakemake.params.costs["fill_values"],
        Nyears,
        snakemake.params.costs["default_exchange_rate"],
        snakemake.params.costs["future_exchange_rate_strategy"],
        snakemake.params.costs["custom_future_exchange_rate"],
    )

    grouping_years_heat = snakemake.params.existing_capacities["grouping_years_heat"]

    # TODO: not implemented in -sec yet
    # if options["enable"]["heat"]:
    #     time_dep_hp_cop = options["time_dep_hp_cop"]
    #     ashp_cop = (
    #         xr.open_dataarray(snakemake.input.cop_air_total)
    #         .to_pandas()
    #         .reindex(index=n.snapshots)
    #     )
    #     gshp_cop = (
    #         xr.open_dataarray(snakemake.input.cop_soil_total)
    #         .to_pandas()
    #         .reindex(index=n.snapshots)
    #     )
    #     default_lifetime = snakemake.params.existing_capacities[
    #         "default_heating_lifetime"
    #     ]
    #     add_heating_capacities_installed_before_baseyear(
    #         n,
    #         baseyear,
    #         grouping_years_heat,
    #         ashp_cop,
    #         gshp_cop,
    #         time_dep_hp_cop,
    #         costs,
    #         default_lifetime,
    #     )

    # if options.get("cluster_heat_buses", False):
    #     cluster_heat_buses(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))

    sanitize_carriers(n, snakemake.config)
    sanitize_locations(n)

    n.export_to_netcdf(snakemake.output[0])
