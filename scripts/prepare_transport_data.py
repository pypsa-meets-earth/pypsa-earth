# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepares the nodal energy demand, vehicle availability and demand-side
management profiles for the (land) transport sector.

Relevant Settings
-----------------

```yaml
custom_data:
    transport_demand:

sector:
    bev_plug_to_wheel_efficiency:
    bev_charge_efficiency:
    transport_heating_deadband_upper:
    transport_heating_deadband_lower:
    ICE_lower_degree_factor:
    ICE_upper_degree_factor:
    EV_lower_degree_factor:
    EV_upper_degree_factor:
    bev_avail_max:
    bev_avail_mean:
    bev_dsm_restriction_value:
    bev_dsm_restriction_time:
```

The remaining keys of the ``sector`` land transport block in
``config.default.yaml`` (``v2g``, ``bev_dsm``, ``bev_energy``,
``bev_availability``, ``transport_fuel_cell_efficiency`` and
``transport_internal_combustion_efficiency``) are not read here; they are
consumed by ``prepare_sector_network``.

Inputs
------

- ``networks/elec_s{simpl}_{clusters}.nc``: clustered base network.
- ``resources/energy_totals_{planning_horizons}.csv``: yearly energy totals per country, including road and rail transport demand.
- ``data/emobility/KFZ__count``: weekly traffic count profile for all motor vehicles, used to shape the overall transport demand time series.
- ``data/emobility/Pkw__count``: weekly traffic count profile for cars (Pkw), used to derive the BEV plugged-in availability profile.
- ``resources/transport_data.csv``: number of registered cars and average fuel efficiency per country.
- ``resources/population_shares/pop_layout_elec_s{simpl}_{clusters}_{planning_horizons}.csv``: clustered population layout per node.
- ``resources/temperatures/temp_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: nodal air temperature time series, used to correct demand for vehicle heating/cooling.

Outputs
-------

- ``resources/demand/transport_s{simpl}_{clusters}_{planning_horizons}.csv``: hourly land transport energy demand per node, in MWh.
- ``resources/pattern_profiles/avail_profile_s{simpl}_{clusters}_{planning_horizons}.csv``: share of the battery electric vehicle (BEV) fleet available to the grid, per node and timestep.
- ``resources/pattern_profiles/dsm_profile_s{simpl}_{clusters}_{planning_horizons}.csv``: restrictions on the state of charge of BEVs for demand-side management (DSM).
- ``resources/demand/nodal_transport_data_s{simpl}_{clusters}_{planning_horizons}.csv``: nodal data on number of cars and average fuel efficiency.

Description
-----------

The rule :mod:`prepare_transport_data` combines population layout, traffic
count profiles and country-level energy totals into nodal time series for
the land transport sector. It:

- shapes the yearly road (and rail) transport demand into an hourly profile using the weekly traffic count pattern,
- corrects the demand for the extra heating/cooling energy required by internal combustion engine (ICE) and battery electric vehicles (BEV) based on nodal air temperature,
- derives the share of the BEV fleet plugged in and available to the grid at each timestep, and
- derives the demand-side management (DSM) restriction profile that constrains when BEVs must be charged.
"""

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from _helpers import BASE_DIR, read_csv_nafix


def transport_degree_factor(
    temperature: pd.DataFrame,
    deadband_lower: float = 15,
    deadband_upper: float = 20,
    lower_degree_factor: float = 0.5,
    upper_degree_factor: float = 1.6,
) -> pd.DataFrame:
    """
    Work out how much energy demand in vehicles increases due to heating and
    cooling.

    There is a deadband where there is no increase. Degree factors are %
    increase in demand per degree Celsius outside the deadband, compared to
    no heating/cooling fuel consumption.

    Parameters
    ----------
    temperature : pd.DataFrame
        Air temperature per node (columns) and snapshot (index), in degree Celsius.
    deadband_lower : float
        Temperature below which heating demand starts to increase (default 15C).
    deadband_upper : float
        Temperature above which cooling demand starts to increase (default 20C).
    lower_degree_factor : float
        Percentage increase in demand per degree C below ``deadband_lower`` (default 0.5).
    upper_degree_factor : float
        Percentage increase in demand per degree C above ``deadband_upper`` (default 1.6).

    Returns
    -------
    pd.DataFrame
        Fractional increase in transport energy demand for each node (columns)
        and snapshot (index), relative to demand with no heating or cooling.
        The value is dimensionless: it already accounts for the temperature
        deviation from the deadband, so 0.05 means 5% additional demand rather
        than 5% per degree Celsius.
    """

    dd = temperature.copy()

    dd[(temperature > deadband_lower) & (temperature < deadband_upper)] = 0.0

    dT_lower = deadband_lower - temperature[temperature < deadband_lower]
    dd[temperature < deadband_lower] = lower_degree_factor / 100 * dT_lower

    dT_upper = temperature[temperature > deadband_upper] - deadband_upper
    dd[temperature > deadband_upper] = upper_degree_factor / 100 * dT_upper

    return dd


def generate_periodic_profiles(
    dt_index: pd.DatetimeIndex,
    nodes: pd.Index,
    weekly_profile: np.ndarray,
    localize: str | None = None,
) -> pd.DataFrame:
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer
    time.

    Each snapshot is converted to the local time of the node's country to
    decide which hour of ``weekly_profile`` applies, so the same snapshot can
    select different profile hours for nodes in different time zones.

    Parameters
    ----------
    dt_index : pd.DatetimeIndex
        Timezone-aware snapshots for which the profile should be generated.
    nodes : pd.Index
        Node names, whose first two characters are used as the country code to look up the timezone.
    weekly_profile : np.ndarray
        24*7 long array of hourly values for a typical week, indexed by ``24 * weekday + hour``.
    localize : str | None
        Timezone to attach to the index of the result. The default ``None``
        strips the timezone instead, preserving the wall-clock time of
        ``dt_index``. Local time is therefore used only to look up the profile
        hour, not for the index of the result.

    Returns
    -------
    pd.DataFrame
        Weekly profile per node (columns) and snapshot (index). With the
        default ``localize=None`` the index holds the same instants as
        ``dt_index`` with no timezone attached, so passing UTC snapshots
        returns UTC times rather than local ones.
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


def prepare_transport_data(
    n: pypsa.Network,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Function to prepare the data required for the (land) transport sector.

    Parameters
    ----------
    n : pypsa.Network
        Clustered network, used for its snapshots.

    Returns
    -------
    nodal_energy_totals : pd.DataFrame
        Energy totals per node, weighted by population share.
    transport : pd.DataFrame
        Nodal transport demand time series, corrected for vehicle heating/cooling.
    avail_profile : pd.DataFrame
        Share of the BEV fleet available to the grid, per node and timestep.
    dsm_profile : pd.DataFrame
        DSM restriction profile, per node and timestep.
    nodal_transport_data : pd.DataFrame
        Nodal data on number of cars and average fuel efficiency.
    """

    energy_totals = read_csv_nafix(
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

    traffic = read_csv_nafix(
        snakemake.input.traffic_data_KFZ, skiprows=2, usecols=["count"]
    ).squeeze("columns")

    # Generate profiles
    transport_shape = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=traffic.values,
    )

    snapshot_weights = n.snapshot_weightings.generators.reindex(transport_shape.index)

    if snapshot_weights.isna().any():
        raise ValueError(
            "Transport-profile snapshots do not match the network snapshot weights."
        )

    weighted_sum = transport_shape.mul(snapshot_weights, axis=0).sum(axis=0)

    if (weighted_sum <= 0).any():
        raise ValueError(
            "Transport profiles must have a positive weighted sum for every node."
        )

    transport_shape = transport_shape.div(weighted_sum, axis=1)

    transport_data = read_csv_nafix(
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
    weighted_transport_shape = transport_shape.mul(snapshot_weights, axis=0)

    ice_correction = (
        weighted_transport_shape.mul(1 + dd_ICE).sum() / weighted_transport_shape.sum()
    )

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

    traffic = read_csv_nafix(
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
        )

    n = pypsa.Network(snakemake.input.network)

    # Get population layout
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
