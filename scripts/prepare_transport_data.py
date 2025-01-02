# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
import os

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr


def transport_degree_factor(
    temperature,
    deadband_lower=15,
    deadband_upper=20,
    lower_degree_factor=0.5,
    upper_degree_factor=1.6,
):
    """
    Work out how much energy demand in vehicles increases due to heating and
    cooling.

    There is a deadband where there is no increase. Degree factors are %
    increase in demand compared to no heating/cooling fuel consumption.
    Returns per unit increase in demand for each place and time
    """

    dd = temperature.copy()

    dd[(temperature > deadband_lower) & (temperature < deadband_upper)] = 0.0

    dT_lower = deadband_lower - temperature[temperature < deadband_lower]
    dd[temperature < deadband_lower] = lower_degree_factor / 100 * dT_lower

    dT_upper = temperature[temperature > deadband_upper] - deadband_upper
    dd[temperature > deadband_upper] = upper_degree_factor / 100 * dT_upper

    return dd


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


def prepare_transport_data(n):
    """
    Function to prepare the data required for the (land) transport sector.
    """

    energy_totals = pd.read_csv(
        snakemake.input.energy_totals_name,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )  # TODO change with real numbers

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.0)
    nodal_energy_totals.index = pop_layout.index
    # # district heat share not weighted by population
    # district_heat_share = nodal_energy_totals["district heat share"].round(2)
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction, axis=0)

    # Get overall demand curve for all vehicles

    traffic = pd.read_csv(
        snakemake.input.traffic_data_KFZ, skiprows=2, usecols=["count"]
    ).squeeze("columns")

    # Generate profiles
    transport_shape = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=traffic.values,
    )

    nodal_transport_shape = transport_shape / transport_shape.sum().sum()
    transport_shape = transport_shape / transport_shape.sum()

    transport_data = pd.read_csv(
        snakemake.input.transport_name, index_col=0, keep_default_na=False
    )

    nodal_transport_data = transport_data.reindex(pop_layout.ct, fill_value=0.0)

    nodal_transport_data.index = pop_layout.index
    nodal_transport_data["number cars"] = (
        pop_layout["fraction"] * nodal_transport_data["number cars"]
    )
    nodal_transport_data.loc[
        nodal_transport_data["average fuel efficiency"] == 0.0,
        "average fuel efficiency",
    ] = transport_data["average fuel efficiency"].mean()

    # electric motors are more efficient, so alter transport demand

    plug_to_wheels_eta = options.get("bev_plug_to_wheel_efficiency", 0.2)
    battery_to_wheels_eta = plug_to_wheels_eta * options.get(
        "bev_charge_efficiency", 0.9
    )

    efficiency_gain = (
        nodal_transport_data["average fuel efficiency"] / battery_to_wheels_eta
    )

    # get heating demand for correction to demand time series
    temperature = xr.open_dataarray(snakemake.input.temp_air_total).to_pandas()

    # correction factors for vehicle heating
    dd_ICE = transport_degree_factor(
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["ICE_lower_degree_factor"],
        options["ICE_upper_degree_factor"],
    )

    dd_EV = transport_degree_factor(
        temperature,
        options["transport_heating_deadband_lower"],
        options["transport_heating_deadband_upper"],
        options["EV_lower_degree_factor"],
        options["EV_upper_degree_factor"],
    )

    # divide out the heating/cooling demand from ICE totals
    # and multiply back in the heating/cooling demand for EVs
    ice_correction = (transport_shape * (1 + dd_ICE)).sum() / transport_shape.sum()

    if snakemake.config["custom_data"]["transport_demand"]:
        energy_totals_transport = nodal_energy_totals["total road"]

        transport = transport_shape.multiply(energy_totals_transport) * 1e6 * Nyears
    else:
        energy_totals_transport = (
            nodal_energy_totals["total road"]
            + nodal_energy_totals["total rail"]
            - nodal_energy_totals["electricity rail"]
        )
        transport = (
            (transport_shape.multiply(energy_totals_transport) * 1e6 * Nyears)
            .divide(efficiency_gain * ice_correction)
            .multiply(1 + dd_EV)
        )

    # derive plugged-in availability for PKW's (cars)

    traffic = pd.read_csv(
        snakemake.input.traffic_data_Pkw, skiprows=2, usecols=["count"]
    ).squeeze("columns")

    avail_max = options.get("bev_avail_max", 0.95)
    avail_mean = options.get("bev_avail_mean", 0.8)

    avail = avail_max - (avail_max - avail_mean) * (traffic - traffic.min()) / (
        traffic.mean() - traffic.min()
    )

    avail_profile = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=avail.values,
    )

    dsm_week = np.zeros((24 * 7,))

    dsm_week[(np.arange(0, 7, 1) * 24 + options["bev_dsm_restriction_time"])] = options[
        "bev_dsm_restriction_value"
    ]

    dsm_profile = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=dsm_week,
    )

    return (
        nodal_energy_totals,
        transport,
        avail_profile,
        dsm_profile,
        nodal_transport_data,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_transport_data",
            simpl="",
            clusters="4",
            planning_horizons="2030",
            demand="AB",
        )

    n = pypsa.Network(snakemake.input.network)

    # Get population layout
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
        transport,
        avail_profile,
        dsm_profile,
        nodal_transport_data,
    ) = prepare_transport_data(n)

    # Save the generated output files to snakemake paths

    # Transport demand per node per timestep
    transport.to_csv(snakemake.output.transport)

    # Available share of the battery to be used by the grid
    avail_profile.to_csv(snakemake.output.avail_profile)

    # Restrictions on state of charge of EVs
    dsm_profile.to_csv(snakemake.output.dsm_profile)

    # Nodal data on number of cars
    nodal_transport_data.to_csv(snakemake.output.nodal_transport_data)
