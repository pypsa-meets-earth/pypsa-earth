
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr

from helpers import mock_snakemake
from helpers import create_temperature_dummy
from helpers import create_transport_data_dummy
from helpers import create_energy_totals_dummy

def transport_degree_factor(
    temperature,
    deadband_lower=15,
    deadband_upper=20,
    lower_degree_factor=0.5,
    upper_degree_factor=1.6):
    """
    Work out how much energy demand in vehicles increases due to heating and cooling.
    There is a deadband where there is no increase.
    Degree factors are % increase in demand compared to no heating/cooling fuel consumption.
    Returns per unit increase in demand for each place and time
    """

    dd = temperature.copy()

    dd[(temperature > deadband_lower) & (temperature < deadband_upper)] = 0.

    dT_lower = deadband_lower - temperature[temperature < deadband_lower]
    dd[temperature < deadband_lower] = lower_degree_factor / 100 * dT_lower

    dT_upper = temperature[temperature > deadband_upper] - deadband_upper
    dd[temperature > deadband_upper] = upper_degree_factor / 100 * dT_upper

    return dd

def generate_periodic_profiles(dt_index, nodes, weekly_profile, localize=None):
    """
    Give a 24*7 long list of weekly hourly profiles, generate this for each
    country for the period dt_index, taking account of time zones and summer time.
    """

    weekly_profile = pd.Series(weekly_profile, range(24*7))

    week_df = pd.DataFrame(index=dt_index, columns=nodes)

    for node in nodes:
        timezone = pytz.timezone(pytz.country_timezones[node[:2]][0])
        tz_dt_index = dt_index.tz_convert(timezone)
        week_df[node] = [24 * dt.weekday() + dt.hour for dt in tz_dt_index]
        week_df[node] = week_df[node].map(weekly_profile)

    week_df = week_df.tz_localize(localize)

    return week_df


def prepare_transport_data(n):


    ##############
    #Heating
    ##############


    # ashp_cop = xr.open_dataarray(snakemake.input.cop_air_total).to_pandas().reindex(index=n.snapshots)
    # gshp_cop = xr.open_dataarray(snakemake.input.cop_soil_total).to_pandas().reindex(index=n.snapshots)

    # solar_thermal = xr.open_dataarray(snakemake.input.solar_thermal_total).to_pandas().reindex(index=n.snapshots)
    # # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    # solar_thermal = options['solar_cf_correction'] * solar_thermal / 1e3

    energy_totals = pd.read_csv(snakemake.input.energy_totals_name, index_col=0)

    # Create energy_totals data dummy for Morocco. TODO Remove once real data is available
    energy_totals = create_energy_totals_dummy(pop_layout, energy_totals)

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.)
    nodal_energy_totals.index = pop_layout.index
    # # district heat share not weighted by population
    # district_heat_share = nodal_energy_totals["district heat share"].round(2)
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction, axis=0)

    # # copy forward the daily average heat demand into each hour, so it can be multipled by the intraday profile
    # daily_space_heat_demand = xr.open_dataarray(snakemake.input.heat_demand_total).to_pandas().reindex(index=n.snapshots, method="ffill")

    # intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    # sectors = ["residential", "services"]
    # uses = ["water", "space"]

    # heat_demand = {}
    # electric_heat_supply = {}
    # for sector, use in product(sectors, uses):
    #     weekday = list(intraday_profiles[f"{sector} {use} weekday"])
    #     weekend = list(intraday_profiles[f"{sector} {use} weekend"])
    #     weekly_profile = weekday * 5 + weekend * 2
    #     intraday_year_profile = generate_periodic_profiles(
    #         daily_space_heat_demand.index.tz_localize("UTC"),
    #         nodes=daily_space_heat_demand.columns,
    #         weekly_profile=weekly_profile
    #     )

    #     if use == "space":
    #         heat_demand_shape = daily_space_heat_demand * intraday_year_profile
    #     else:
    #         heat_demand_shape = intraday_year_profile

    #     heat_demand[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(nodal_energy_totals[f"total {sector} {use}"]) * 1e6
    #     electric_heat_supply[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(nodal_energy_totals[f"electricity {sector} {use}"]) * 1e6

    # heat_demand = pd.concat(heat_demand, axis=1)
    # electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # # subtract from electricity load since heat demand already in heat_demand
    # electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    # n.loads_t.p_set[electric_nodes] = n.loads_t.p_set[electric_nodes] - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]

    ##############
    #Transport
    ##############

    ## Get overall demand curve for all vehicles

    traffic = pd.read_csv(snakemake.input.traffic_data_KFZ, skiprows=2, usecols=["count"], squeeze=True)

    #Generate profiles
    transport_shape = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=traffic.values
    )
    transport_shape = transport_shape / transport_shape.sum()

    transport_data = pd.read_csv(snakemake.input.transport_name, index_col=0)

    # Create transport data dummy for Morocco. TODO Remove once real data is available
    transport_data = create_transport_data_dummy(pop_layout, transport_data, cars = 4000000, average_fuel_efficiency = 0.7)

    nodal_transport_data = transport_data.loc[pop_layout.ct].fillna(0.)
    nodal_transport_data.index = pop_layout.index
    nodal_transport_data["number cars"] = pop_layout["fraction"] * nodal_transport_data["number cars"]
    nodal_transport_data.loc[nodal_transport_data["average fuel efficiency"] == 0., "average fuel efficiency"] = transport_data["average fuel efficiency"].mean()


    # electric motors are more efficient, so alter transport demand

    plug_to_wheels_eta = options.get("bev_plug_to_wheel_efficiency", 0.2)
    battery_to_wheels_eta = plug_to_wheels_eta * options.get("bev_charge_efficiency", 0.9)

    efficiency_gain = nodal_transport_data["average fuel efficiency"] / battery_to_wheels_eta

    #get heating demand for correction to demand time series
    temperature = xr.open_dataarray(snakemake.input.temp_air_total).to_pandas()

    # Create temperature data dummy for Morocco. TODO Remove once real data is available
    temperature = create_temperature_dummy(pop_layout, temperature)

    # correction factors for vehicle heating
    dd_ICE = transport_degree_factor(
        temperature,
        options['transport_heating_deadband_lower'],
        options['transport_heating_deadband_upper'],
        options['ICE_lower_degree_factor'],
        options['ICE_upper_degree_factor']
    )

    dd_EV = transport_degree_factor(
        temperature,
        options['transport_heating_deadband_lower'],
        options['transport_heating_deadband_upper'],
        options['EV_lower_degree_factor'],
        options['EV_upper_degree_factor']
    )

    # divide out the heating/cooling demand from ICE totals
    # and multiply back in the heating/cooling demand for EVs
    ice_correction = (transport_shape * (1 + dd_ICE)).sum() / transport_shape.sum()

    energy_totals_transport = nodal_energy_totals["total road"] + nodal_energy_totals["total rail"] - nodal_energy_totals["electricity rail"]

    transport = (transport_shape.multiply(energy_totals_transport) * 1e6 * Nyears).divide(efficiency_gain * ice_correction).multiply(1 + dd_EV)

    ## derive plugged-in availability for PKW's (cars)

    traffic = pd.read_csv(snakemake.input.traffic_data_Pkw, skiprows=2, usecols=["count"], squeeze=True)

    avail_max = options.get("bev_avail_max", 0.95)
    avail_mean = options.get("bev_avail_mean", 0.8)

    avail = avail_max - (avail_max - avail_mean) * (traffic - traffic.min()) / (traffic.mean() - traffic.min())

    avail_profile = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=avail.values
    )

    dsm_week = np.zeros((24*7,))

    dsm_week[(np.arange(0,7,1) * 24 + options['bev_dsm_restriction_time'])] = options['bev_dsm_restriction_value']

    dsm_profile = generate_periodic_profiles(
        dt_index=n.snapshots.tz_localize("UTC"),
        nodes=pop_layout.index,
        weekly_profile=dsm_week
    )


    return nodal_energy_totals, transport, avail_profile, dsm_profile, nodal_transport_data


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake("prepare_transport_data",
                                   simpl="",
                                   clusters="4")

    
    n = pypsa.Network(snakemake.input.network)

    # Get pop_layout (dummy)
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout_dummy, index_col=0)

    # Add options
    options = snakemake.config["sector"]

    # Get Nyears
    Nyears = n.snapshot_weightings.generators.sum() / 8760

    # Prepare transport data
    nodal_energy_totals, transport, avail_profile, dsm_profile, nodal_transport_data = prepare_transport_data(n)


    #Create create_temperature_dummy
    # temperature = xr.open_dataarray(snakemake.input.temp_air_total).to_pandas()
    # temperature_dummy = create_temperature_dummy(pop_layout, temperature)

    # Test temperature dummy


    # Test transport data dummy
    # transport_data = pd.read_csv(snakemake.input.transport_name, index_col=0)
    # transport_data_dummy = create_transport_data_dummy(pop_layout, transport_data)

    # Test energy_totals_dummy
    # energy_totals = pd.read_csv(snakemake.input.energy_totals_name, index_col=0)
    # energy_totals_dummy = create_energy_totals_dummy(pop_layout, energy_totals)


    print('successfull run')


