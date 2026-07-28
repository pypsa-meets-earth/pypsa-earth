# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Prepare heating-sector time series and summary tables used by the PyPSA-Earth workflow.

Relevant Settings
-----------------

None

Inputs
------

- ``networks/{RDIR}/elec_s{simpl}_{clusters}.nc``: Clustered PyPSA network
- ``resources/{SECDIR}/energy_totals_{demand}_{planning_horizons}.csv``: Energy totals for each sector
- ``resources/{SECDIR}/population_shares/pop_layout_elec_s{simpl}_{clusteres}_{planning_horizons}.csv``: Population shares indexed by node
- ``resources/{SECDIR}/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Hourly air temperature times series
- ``resources/{SECDIR}/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Ground/soil source heat pump COP time series aligned to the network snapshots.
- ``resources/{SECDIR}/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Air source heat pump COP time series aligned to the network snapshots.
- ``resources/{SECDIR}/demand/heat/solar_thermal_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Solar thermal irradiance or generation used to generate solar thermal time series per node.
- ``resources/{SECDIR}/demand/heat/heat_demand_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Daily average heat demand per node
- ``data/heat_load_profile_BDEW.csv``: Representative weekday/weekend heat load profiles by sector and use from BDEW

Outputs
-------

- ``resources/{SECDIR}/demand/heat/nodal_energy_heat_totals_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: per node annual energy totals CSV
- ``resources/{SECDIR}/demand/heat/heat_demand_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: hourly heat demand time series
- ``resources/{SECDIR}/demand/heat/ashp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: air source heat pump COP time series CSV
- ``resources/{SECDIR}/demand/heat/gshp_cop_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: ground source heat pump COP time series CSV
- ``resources/{SECDIR}/demand/heat/solar_thermal_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: Solar thermal demand time series CSV
- ``resources/{SECDIR}/demand/heat/district_heat_share_{demand}_s{simpl}_{clusters}_{planning_horizons}.csv``: per node district heat share CSV

Description
-----------
This module builds nodal heating energy totals, hourly heat demand time series, coefficient of performance (COP)
profiles for air-source and ground-source heat pumps, solar thermal time series and the district heating share per node.
"""
import os
from itertools import product

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from _helpers import mock_snakemake, read_csv_nafix


def generate_periodic_profiles(
    dt_index: pd.DatetimeIndex,
    nodes: list[str],
    weekly_profile: list,
    localize: str = None,
) -> pd.DataFrame:
    """
    Create per node hourly profiles from a 7-day (168 hour) template.

    The provided `weekly_profile` is a list or sequence of 168 hourly values
    (24 * 7). For each node the function converts the provided UTC-aware
    `dt_index` to the node's local timezone and maps each timestamp to the corresponding hour in the
    weekly profile. The resulting DataFrame is indexed by `dt_index` and has
    a column per node containing the hourly profile values localized for that node.

    Parameters
    ----------
    dt_index : pd.DatetimeIndex
        Time index for the target period (should be timezone-aware in UTC).
    nodes : list[str]
        Iterable of node identifiers where the first two characters are an
        ISO2 country code used to determine the node timezone.
    weekly_profile : list
        Iterable of 168 hourly values representing a typical week (24 * 7 hours).
    localize : str, optional
        Optional timezone string to localize the final DataFrame to

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by `dt_index` with one column per node containing
        the localized hourly profile values.
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


def prepare_heat_data(n: pypsa.Network) -> tuple:
    """Prepare heating sector inputs for a PyPSA network.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA `Network` object whose `snapshots` and other metadata are
        used to align time series and perform reindexing.

    Returns
    -------
    tuple
    A tuple with the following elements in order:
        nodal_energy_totals: pd.DataFrame
            per-node annual energy totals
        heat_demand: pd.DataFrame
            hourly heat demand time series
        ashp_cop: pd.DataFrame
            air-source heat pump COP time series
        gshp_cop: pd.DataFrame
            ground-source heat pump COP time series
        solar_thermal:pd.DataFrame
            solar thermal generation time series
        district_heat_share: float
            per-node fraction of heat served by district heat
    """
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
    solar_thermal = (
        options["solar_thermal_collector"]["cf_correction"] * solar_thermal / 1e3
    )

    energy_totals = read_csv_nafix(
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

    intraday_profiles = read_csv_nafix(
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
    pop_layout = read_csv_nafix(
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
