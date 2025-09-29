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
from _helpers import locate_bus, mock_snakemake

# TODO Add to Snakemake
CALIBR_DIR = "data/custom"
CALIBR_HEAT_FL = "mod_heating_calibr_clean.csv"
CALIBR_COOL_FL = "mod_cooling_calibr_clean.csv"

# TODO Add to the parameters into a configuration file
HEAT_DEMAND_TOTAL = 1
COOL_DEMAND_TOTAL = 1

SHARE_HEAT_RESID_DEMAND = 0.12
SHARE_WATER_RESID_DEMAND = 0.12

SHARE_HEAT_SERVICES_DEMAND = 0.32
SHARE_WATER_SERVICES_DEMAND = 0.05

SHARE_COOL_RESID_DEMAND = 0.19
SHARE_COOL_SERVICES_DEMAND = 0.11

SHARE_ELECTRICITY_RESID_SPACE = 0.10
SHARE_ELECTRICITY_SERVICES_SPACE = 0.13

SHARE_DISTRICT_HEAT = 0.1

CALIBRATE_LOAD = False
CALIBRATION_LEVEL = 1


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


def scale_demand(data_df, calibr_df, load_mode, geom_id, k=1):
    """
    Apply a linear transformation to the demand dataframe to match
    the demand values with the measured or modeled ones.

    Parameters
    ----------
    data_df : DataFrame
        Original demand
    calibr_df : DataFrame
        Calibration coefficients
    load_mode : string
        "heating" or "cooling" as is defined in the calibration data
    geom_id : string
        ID of the bus shapes
    """

    data_col = "scaling_" + load_mode

    if calibr_df is not None:
        scale_coeff = calibr_df[[data_col, geom_id]].set_index(geom_id)
        # some GADM sub-regions can be missed from the scaling factors calculations
        subregions_missed = data_df.columns.difference(scale_coeff.index)
        default_scaling = pd.DataFrame(
            data=[1.0] * len(subregions_missed),
            index=subregions_missed,
            columns=[data_col],
        )
        scale_coeff_clean = pd.concat([scale_coeff[data_col], default_scaling])
        data_df.mul(scale_coeff_clean.to_dict()[data_col], axis=1)
    # TODO Avoid using the global variables
    elif load_mode == "heating":
        k = HEAT_DEMAND_TOTAL / data_df.sum().sum()
        data_df = k * data_df
    elif load_mode == "cooling":
        k = COOL_DEMAND_TOTAL / data_df.sum().sum()
        data_df = k * data_df
    else:
        data_df = k * data_df
    return data_df


def prepare_heat_data(n, snapshots, countries):
    # heating
    ashp_cop = (
        xr.open_dataarray(snakemake.input.cop_air_total)
        .to_pandas()
        .reindex(index=snapshots)
    )
    gshp_cop = (
        xr.open_dataarray(snakemake.input.cop_soil_total)
        .to_pandas()
        .reindex(index=snapshots)
    )
    solar_thermal = (
        xr.open_dataarray(snakemake.input.solar_thermal_total)
        .to_pandas()
        .reindex(index=snapshots)
    )
    # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    solar_thermal = options["solar_cf_correction"] * solar_thermal / 1e3

    # cooling
    hp_cooling_total_cop = (
        xr.open_dataarray(snakemake.input.cop_hp_cooling_total)
        .to_pandas()
        .reindex(index=snapshots)
    )
    ac_cooling_total_cop = (
        xr.open_dataarray(snakemake.input.cop_ac_cooling_total)
        .to_pandas()
        .reindex(index=snapshots)
    )
    apft_abch_cooling_total_cop = (
        xr.open_dataarray(snakemake.input.capft_abch_cooling_total)
        .to_pandas()
        .reindex(index=snapshots)
    )

    # energy balance
    energy_totals = pd.read_csv(
        snakemake.input.energy_totals_name,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.0)
    nodal_energy_totals.index = pop_layout.index
    # TODO keeping the solution consistent
    nodal_energy_totals["district heat share"] = SHARE_DISTRICT_HEAT
    # district heat share not weighted by population
    district_heat_share = nodal_energy_totals["district heat share"]  # .round(2)
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction, axis=0)

    energy_col_names = nodal_energy_totals.columns

    # TODO reimplement using DRY approach
    # a proper structure of the inputs is essential
    residential_cols = energy_col_names[energy_col_names.str.contains("residential")]
    services_cols = energy_col_names[energy_col_names.str.contains("services")]

    energy_residential = nodal_energy_totals.loc[:, residential_cols].sum().sum()
    energy_services = nodal_energy_totals.loc[:, services_cols].sum().sum()
    energy_total = nodal_energy_totals.sum().sum()

    nodal_energy_totals[f"heat residential space"] = (
        SHARE_HEAT_RESID_DEMAND * energy_residential
    )
    nodal_energy_totals[f"heat residential water"] = (
        SHARE_WATER_RESID_DEMAND * energy_residential
    )
    nodal_energy_totals[f"heat services space"] = (
        SHARE_HEAT_SERVICES_DEMAND * energy_services
    )
    nodal_energy_totals[f"heat services water"] = (
        SHARE_WATER_SERVICES_DEMAND * energy_services
    )

    nodal_energy_totals[f"cool residential space"] = (
        SHARE_COOL_RESID_DEMAND * energy_residential
    )
    nodal_energy_totals[f"cool services space"] = (
        SHARE_COOL_SERVICES_DEMAND * energy_services
    )

    nodal_energy_totals[f"electricity residential space"] = (
        SHARE_ELECTRICITY_RESID_SPACE * energy_residential
    )
    nodal_energy_totals[f"electricity services space"] = (
        SHARE_ELECTRICITY_SERVICES_SPACE * energy_services
    )

    # heating/cooling demand profiles
    # copy forward the daily average heat demand into each hour, so it can be multiplied by the intraday profile
    daily_space_heat_demand = (
        xr.open_dataarray(snakemake.input.heat_demand_total)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )
    daily_space_cooling_demand = (
        xr.open_dataarray(snakemake.input.cooling_demand_total)
        .to_pandas()
        .reindex(index=snapshots, method="ffill")
    )

    intraday_profiles_heating = pd.read_csv(
        snakemake.input.heat_profile, index_col=0
    )  # TODO GHALAT
    intraday_profiles_cooling = pd.read_csv(
        snakemake.input.cooling_profile, index_col=0
    )  # TODO GHALAT

    if CALIBRATE_LOAD:
        calibr_heat_df = pd.read_csv(
            os.path.join(CALIBR_DIR, CALIBR_HEAT_FL), index_col=0
        )
        calibr_cool_df = pd.read_csv(
            os.path.join(CALIBR_DIR, CALIBR_COOL_FL), index_col=0
        )

        calibr_heat_buses_df = locate_bus(
            df=calibr_heat_df,
            countries=countries,
            gadm_level=1,
            path_to_gadm=snakemake.input.shapes_path,
            gadm_clustering=True,
            dropnull=True,
            col_out=None,
        )

        calibr_cool_buses_df = locate_bus(
            df=calibr_cool_df,
            countries=countries,
            gadm_level=1,
            path_to_gadm=snakemake.input.shapes_path,
            gadm_clustering=True,
            dropnull=True,
            col_out=None,
        )

    else:
        calibr_heat_buses_df = None
        calibr_cool_buses_df = None

    loads = ["heating", "cooling"]
    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
    cooling_demand = {}
    electric_heat_supply = {}

    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles_heating[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles_heating[f"{sector} {use} weekend"])
        weekly_profile_heating = weekday * 5 + weekend * 2

        intraday_year_profile_heating = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile_heating,
        )

        if use == "space":
            heating_demand_shape = (
                daily_space_heat_demand * intraday_year_profile_heating
            )
        else:
            heating_demand_shape = intraday_year_profile_heating

        heat_demand[f"{sector} {use}"] = (
            heating_demand_shape / heating_demand_shape.sum()
        ).multiply(
            nodal_energy_totals[f"total {sector} {use}"]
        )  # * 1e6  # TODO v0.0.2

        electric_heat_supply[f"{sector} {use}"] = (
            heating_demand_shape / heating_demand_shape.sum()
        ).multiply(
            nodal_energy_totals[f"electricity {sector} {use}"]
        )  # * 1e6  # TODO v0.0.2

        electric_heat_supply[f"total {use}"] = (
            heating_demand_shape / heating_demand_shape.sum()
        ).multiply(
            nodal_energy_totals[f"electricity residential {use}"]
            + nodal_energy_totals[f"electricity services {use}"]
        )  # * 1e6  # TODO v0.0.2

    heat_demand = pd.concat(heat_demand, axis=1)

    # TODO Account for the differences in a column name alternative/Voronoi clustering
    heat_demand = scale_demand(
        data_df=heat_demand,
        calibr_df=calibr_heat_buses_df,
        load_mode="heating",
        geom_id=CALIBRATION_LEVEL,
    )

    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand #TODO v0.1
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = (
        n.loads_t.p_set[electric_nodes]
        - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]
    )

    # TODO it's possible to account for weekday/weekend differences
    for use in ["space"]:
        day_cooling = list(intraday_profiles_cooling[f"cooling {use}"])
        weekly_profile_cooling = day_cooling * 7
        intraday_year_profile_cooling = generate_periodic_profiles(
            daily_space_cooling_demand.index.tz_localize("UTC"),
            nodes=daily_space_cooling_demand.columns,
            weekly_profile=weekly_profile_cooling,
        )

        if use == "space":
            cooling_demand_shape = (
                daily_space_cooling_demand * intraday_year_profile_cooling
            )
        else:
            cooling_demand_shape = intraday_year_profile_cooling

        cooling_demand[f"{use}"] = (
            cooling_demand_shape / cooling_demand_shape.sum()
        ).multiply(
            (
                nodal_energy_totals[f"cool residential space"]
                + nodal_energy_totals[f"cool services space"]
            )
        )  # * 1e6

    cooling_demand = pd.concat(cooling_demand, axis=1)

    # TODO Account for the differences in a column name alternative/Voronoi clustering
    cooling_demand = scale_demand(
        data_df=cooling_demand,
        calibr_df=calibr_cool_buses_df,
        load_mode="cooling",
        geom_id=CALIBRATION_LEVEL,
    )

    return (
        nodal_energy_totals,
        heat_demand,
        cooling_demand,
        ashp_cop,
        gshp_cop,
        solar_thermal,
        hp_cooling_total_cop,
        ac_cooling_total_cop,
        apft_abch_cooling_total_cop,
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
    country_list = snakemake.params.countries

    # Get Nyears
    Nyears = n.snapshot_weightings.generators.sum() / 8760

    (
        nodal_energy_totals,
        heat_demand,
        cooling_demand,
        ashp_cop,
        gshp_cop,
        solar_thermal,
        hp_cooling_total_cop,
        ac_cooling_total_cop,
        apft_abch_cooling_total_cop,
        district_heat_share,
    ) = prepare_heat_data(n, n.snapshots, country_list)

    # Save the generated output files to snakemake paths
    nodal_energy_totals.to_csv(snakemake.output.nodal_energy_totals)
    heat_demand.to_csv(snakemake.output.heat_demand)
    cooling_demand.to_csv(snakemake.output.cooling_demand)

    ashp_cop.to_csv(snakemake.output.ashp_cop)
    gshp_cop.to_csv(snakemake.output.gshp_cop)
    solar_thermal.to_csv(snakemake.output.solar_thermal)

    hp_cooling_total_cop.to_csv(snakemake.output.cop_hp_cooling_total)
    ac_cooling_total_cop.to_csv(snakemake.output.cop_ac_cooling_total)
    apft_abch_cooling_total_cop.to_csv(snakemake.output.capft_abch_cooling_total)

    district_heat_share.to_csv(snakemake.output.district_heat_share)
