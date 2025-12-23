# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import os
from itertools import product

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from _helpers import mock_snakemake


def generate_periodic_profiles(dt_index, nodes, weekly_profile, localize=None):
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer
    time.
    """

    weekly_profile = pd.Series(weekly_profile, range(24 * 7))

    week_df = pd.DataFrame(index=dt_index, columns=nodes)

    for node in nodes:
        timezone = pytz.timezone(pytz.country_timezones[node[:2]][0])
        tz_dt_index = dt_index.tz_convert(timezone)
        week_df[node] = [24 * dt.weekday() + dt.hour for dt in tz_dt_index]
        week_df[node] = week_df[node].map(weekly_profile)

    week_df = week_df.tz_localize(localize)

    return week_df


def prepare_heat_data(n):
    ashp_cop = (
        xr.open_dataarray(snakemake.input.cop_air_total)
        .to_pandas()
        .reindex(index=n.snapshots)
    )
    gshp_cop = (
        xr.open_dataarray(snakemake.input.cop_soil_total)
        .to_pandas()
        .reindex(index=n.snapshots)
    )

    solar_thermal = (
        xr.open_dataarray(snakemake.input.solar_thermal_total)
        .to_pandas()
        .reindex(index=n.snapshots)
    )
    # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    solar_thermal = options["solar_cf_correction"] * solar_thermal / 1e3

    energy_totals = pd.read_csv(
        snakemake.input.energy_totals_name,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.0)
    nodal_energy_totals.index = pop_layout.index
    # # district heat share not weighted by population
    district_heat_share = nodal_energy_totals["district heat share"]  # .round(2)
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction, axis=0)

    # copy forward the daily average heat demand into each hour, so it can be multiplied by the intraday profile
    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand_total)
        .to_pandas()
        .reindex(index=n.snapshots, method="ffill")
    )

    intraday_profiles = pd.read_csv(
        snakemake.input.heat_profile, index_col=0
    )  # TODO GHALAT

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
    electric_heat_supply = {}
    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles[f"{sector} {use} weekend"])
        weekly_profile = weekday * 5 + weekend * 2
        intraday_year_profile = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile,
        )

        if use == "space":
            heat_demand_shape = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand_shape = intraday_year_profile

        heat_demand[f"{sector} {use}"] = (
            heat_demand_shape / heat_demand_shape.sum()
        ).multiply(
            nodal_energy_totals[f"total {sector} {use}"]
        ) * 1e6  # TODO v0.0.2
        electric_heat_supply[f"{sector} {use}"] = (
            heat_demand_shape / heat_demand_shape.sum()
        ).multiply(
            nodal_energy_totals[f"electricity {sector} {use}"]
        ) * 1e6  # TODO v0.0.2

    heat_demand = pd.concat(heat_demand, axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand #TODO v0.1
    # electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    # n.loads_t.p_set[electric_nodes] = (
    #     n.loads_t.p_set[electric_nodes]
    #     - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]
    # )

    return (
        nodal_energy_totals,
        heat_demand,
        ashp_cop,
        gshp_cop,
        solar_thermal,
        district_heat_share,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_heat_data",
            simpl="",
            clusters="4",
            planning_horizons=2030,
            demand="AB",
        )

    n = pypsa.Network(snakemake.input.network)

    # Get pop_layout
    pop_layout = pd.read_csv(
        snakemake.input.clustered_pop_layout,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )

    # Add options
    options = snakemake.config["sector"]

    # Get Nyears
    Nyears = n.snapshot_weightings.generators.sum() / 8760

    # Prepare transport data
    (
        nodal_energy_totals,
        heat_demand,
        ashp_cop,
        gshp_cop,
        solar_thermal,
        district_heat_share,
    ) = prepare_heat_data(n)

    # Save the generated output files to snakemake paths
    nodal_energy_totals.to_csv(snakemake.output.nodal_energy_totals)
    heat_demand.to_csv(snakemake.output.heat_demand)
    ashp_cop.to_csv(snakemake.output.ashp_cop)
    gshp_cop.to_csv(snakemake.output.gshp_cop)
    solar_thermal.to_csv(snakemake.output.solar_thermal)
    district_heat_share.to_csv(snakemake.output.district_heat_share)
