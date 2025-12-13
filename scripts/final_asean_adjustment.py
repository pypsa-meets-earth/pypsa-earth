# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Final adjustment to the ASEAN model
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)

carrier_to_keep = [
    # Power carriers (Renewable & Non-renewable)
    "nuclear",
    "CCGT",
    "OCGT",
    "geothermal",
    "hydro",
    "offwind-ac",
    "offwind-dc",
    "onwind",
    "solar",
    "solar rooftop",
    "biomass",
    "biomass EOP",
    "solid biomass",
    # "biogas", # Only used to reduce gas consumption
    # "biogas to gas", # Only used to reduce gas consumption
    "PHS",
    "ror",
    "DC",
    "AC",
    "B2B",
    # Electricity carriers
    "rail transport electricity",
    "industry electricity",
    "agriculture electricity",
    "service electricity",
    "electricity distribution grid",
    "low voltage",
    # Battery-related
    "battery",
    "battery charger",
    "battery discharger",
    "home battery",
    "home battery charger",
    "home battery discharger",
    # Hydrogen-related
    "H2",
    "H2 fuel cell",
    "H2 electrolysis",
    "H2 Electrolysis",
    "H2 Electrolysis",
    "H2 Store Tank",
    # Transport-related
    "Li ion",  # EV batteries
    "land transport EV",
    "BEV charger",
    "V2G",
    # CO₂ / Carbon-related
    "co2",
    "co2 stored",
    "co2 vent",
    # Fossil Fuel Technologies
    "gas",
    "oil",
    "coal",
    "lignite",
    "SMR CC",
    "helmeth",
    "SMR",
    "Fischer-Tropsch",
    "Sabatier",
]


elec_carrier = [
    "AC",
    "industry electricity",
    "agriculture electricity" "rail transport electricity",
    # 'land transport EV', # Electricity growth leads to nonlinearity
]


def extend_carrier_list(config_elec):

    config_carriers = [
        carrier
        for key in ["conventional_carriers", "renewable_carriers"]
        for carrier in config_elec[key]
    ]

    storage_carriers = [
        suffix
        for key in ["StorageUnit", "Store", "Link"]
        for carrier in config_elec["extendable_carriers"][key]
        for suffix in [carrier, carrier + " charger", carrier + " discharger"]
    ]
    return config_carriers + storage_carriers


def strip_network(n, carriers):
    """
    Strip network to core electricity related components

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network object to modify.
    carriers : list
        The carrier to be kept

    Returns
    -------
    pypsa.Network
    """

    m = n.copy()

    nodes_to_keep = m.buses[m.buses.carrier.isin(carriers)].index
    carriers_to_keep = m.carriers[m.carriers.index.isin(carriers)].index

    m.mremove("Bus", n.buses.index.symmetric_difference(nodes_to_keep))
    m.mremove("Carrier", n.carriers.index.symmetric_difference(carriers_to_keep))

    for c in m.iterate_components(
        ["Generator", "Link", "Line", "Store", "StorageUnit", "Load"]
    ):
        if c.name in ["Link", "Line"]:
            location_boolean = c.df.bus0.isin(nodes_to_keep) & c.df.bus1.isin(
                nodes_to_keep
            )
        else:
            location_boolean = c.df.bus.isin(nodes_to_keep)
        to_keep = c.df.index[location_boolean & c.df.carrier.isin(carrier_to_keep)]
        to_drop = c.df.index.symmetric_difference(to_keep)
        m.mremove(c.name, to_drop)

    m.links.loc[m.links.carrier.isin(["H2 Electrolysis", "H2 Fuel Cell"]), "bus2"] = ""

    logger.info("Strip network to core electricity related components")

    return m


def retrieve_population_forcast(file_path, target="TPopulation1Jan"):
    """
    Retrieve and process population forecast data by country and year.

    This function loads population data either from a local cached CSV file
    or downloads it from the UN World Population Prospects (WPP 2024) dataset if not available locally.
    It then groups the data by ISO2 country codes and time, returning a DataFrame of the selected population metric.
    """
    import os

    url = "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/WPP2024_Demographic_Indicators_Medium.csv.gz"

    if os.path.exists(file_path):
        data = pd.read_csv(file_path)
    else:
        data = pd.read_csv(url)
        data.to_csv(file_path)

    data = data.groupby(["ISO2_code", "Time"]).sum()[target].unstack("Time")

    return data


def include_electricity_growth(n, options, planning_horizons):
    """
    Adjust electricity demand for residential, industrial, agricultural, and rail transport sectors
    based on projected population growth and fixed per-capita electricity consumption.

    This function compensates for the lack of year-specific demand scaling in the build_demand_profiles scripts.

    Electricity demand per capita data are sourced from:
    - IEA: https://www.iea.org/data-and-statistics/data-tools/energy-statistics-data-browser
    - East Timor: https://www.ebsco.com/research-starters/power-and-energy/timor-leste-and-renewable-energy
    """
    pop_forecast_path = options["pop_forecast_path"]
    pop_forecast = retrieve_population_forcast(pop_forecast_path)
    pop_forecast = pop_forecast[int(planning_horizons)]

    elec_per_capita = pd.Series(options["elec_per_capita"])
    demand_forecast = pop_forecast.T[elec_per_capita.index] * 1e3 * elec_per_capita

    if options["total_elec_demand"]:
        total_elec_demand = pd.Series(options["total_elec_demand"])[int(planning_horizons)]
        demand_forecast = total_elec_demand * demand_forecast / demand_forecast.sum()

    loads = n.loads[n.loads.carrier.isin(elec_carrier)].copy()

    loads_dynamic_i = n.loads_t.p_set.columns.intersection(loads.index)
    loads_static_i = loads.index.difference(loads_dynamic_i)

    loads_t = n.loads_t.p_set[loads_dynamic_i]

    loads["p_set_static"] = loads["p_set"] * 8760.0
    loads["p_set_dynamic"] = n.snapshot_weightings.objective @ loads_t

    loads["p_set_sum"] = loads["p_set_static"].fillna(0) + loads[
        "p_set_dynamic"
    ].fillna(0)

    loads["country"] = loads.bus.map(n.buses.country)
    country_totals = loads.groupby("country")["p_set_sum"].transform("sum")

    loads["fraction"] = loads["p_set_sum"] / country_totals

    loads["p_set_sum_country_new"] = loads["fraction"] * loads["country"].map(
        demand_forecast
    )
    loads["p_set_new"] = (
        loads["p_set_sum_country_new"] / n.snapshot_weightings.objective.sum()
    )

    n.loads.loc[loads_static_i, "p_set"] = loads.loc[loads_static_i, "p_set_new"]
    n.loads_t.p_set[loads_dynamic_i] = (
        loads.loc[loads_dynamic_i, "p_set_new"] * loads_t / loads_t.mean()
    )

    logger.info(f"Activate: Adjust electricity growth for the year {planning_horizons}")


def redistribute_industrial_load(n, industry_share_data):
    """
    Redistribute total electricity load across industrial and non-industrial sectors
    based on fixed country-specific industry shares.

    Electricity demand by sector data are sourced from:
    - IEA: https://www.iea.org/data-and-statistics/data-tools/energy-statistics-data-browser
    """

    industry_share = pd.Series(industry_share_data)

    loads = n.loads[n.loads.carrier.isin(elec_carrier)].copy()

    loads_dynamic_i = n.loads_t.p_set.columns.intersection(loads.index)
    loads_t = n.loads_t.p_set[loads_dynamic_i]

    loads["p_set_static"] = loads["p_set"] * 8760.0
    loads["p_set_dynamic"] = n.snapshot_weightings.objective @ loads_t

    loads["p_set_sum"] = loads["p_set_static"].fillna(0) + loads[
        "p_set_dynamic"
    ].fillna(0)

    loads["country"] = loads.bus.map(n.buses.country)
    loads["country_totals"] = loads.groupby("country")["p_set_sum"].transform("sum")

    loads["sector_share"] = np.where(
        loads["carrier"] == "industry electricity",
        loads["country"].map(industry_share),
        1 - loads["country"].map(industry_share),
    )

    for carrier_list in [
        ["industry electricity"],
        ["AC", "agriculture electricity", "rail transport electricity"],
    ]:
        loads_c = loads[loads.carrier.isin(carrier_list)].copy()

        loads_c["fraction"] = loads_c["p_set_sum"] / loads_c.groupby("country")[
            "p_set_sum"
        ].transform("sum")
        loads_c["p_set_sum_country_new"] = (
            loads_c["fraction"] * loads["sector_share"] * loads_c["country_totals"]
        )
        loads_c["p_set_new"] = (
            loads_c["p_set_sum_country_new"] / n.snapshot_weightings.objective.sum()
        )

        loads_dynamic_i_c = n.loads_t.p_set.columns.intersection(loads_c.index)
        loads_static_i_c = loads_c.index.difference(loads_dynamic_i_c)

        loads_t_c = n.loads_t.p_set[loads_dynamic_i_c]

        n.loads.loc[loads_static_i_c, "p_set"] = loads_c.loc[
            loads_static_i_c, "p_set_new"
        ]
        n.loads_t.p_set[loads_dynamic_i_c] = (
            loads_c.loc[loads_dynamic_i_c, "p_set_new"] * loads_t_c / loads_t_c.mean()
        )

    logger.info("Activate: Redistribute industrial load")


def readjust_existing_interconnections(n):
    """
    Adjust existing cross-border transmission capacities in the PyPSA network
    to match AIMS (ASEAN Interconnection Masterplan Study) 2024 targets.

    See: https://aseanenergy.org/wp-content/uploads/2024/11/ASEAN-Power-Grid-Interconnections-Project-Profiles.pdf
    """
    AIMS_exist_data = [
        ("ID", "MY", 230),
        ("KH", "LA", 300),
        ("KH", "TH", 250),
        ("KH", "VN", 200),
        ("MY", "SG", 525),
        ("MY", "TH", 380),
        ("LA", "TH", 955),
        ("LA", "MM", 30),
        ("LA", "VN", 1),  # no grid-to-grid connections in 2024.
    ]

    AIMS_exist = pd.DataFrame(
        AIMS_exist_data, columns=["country0", "country1", "s_nom"]
    )

    lines = n.lines.copy()

    # filter in the relevant lines
    lines["country0"] = lines.bus0.map(n.buses.country)
    lines["country1"] = lines.bus1.map(n.buses.country)
    lines["existing"] = lines.build_year <= 2024

    lines = lines[
        (lines.country0 != lines.country1)
        & lines.country0.isin(AIMS_exist.country0)
        & lines.country1.isin(AIMS_exist.country1)
        & lines.existing
    ]  # TODO: future-prone features if AIMS interconnections is in OpenStreetMaps

    lines["fraction"] = lines["length"] / lines.groupby(["country0", "country1"])[
        "length"
    ].transform("sum")

    # reapply the transmission capacity by AIMS_exist data, distribute the line capacity by length
    for i in AIMS_exist.index:
        lines_i = lines[
            (lines.country0 == AIMS_exist.loc[i, "country0"])
            & (lines.country1 == AIMS_exist.loc[i, "country1"])
        ].index

        n.lines.loc[lines_i, "s_nom"] = (
            AIMS_exist.loc[i, "s_nom"] * lines.loc[lines_i, "fraction"]
        )
        n.lines.loc[lines_i, "s_nom_min"] = n.lines.loc[lines_i, "s_nom"]

    logger.info("Activate: Readjust existing interconnections")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "final_asean_adjustment",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate="0.071",
            demand="AB",
            h2export="120",
        )

    configure_logging(snakemake)

    config_elec = snakemake.params.electricity
    options = snakemake.params.final_adjustment
    planning_horizons = snakemake.wildcards.planning_horizons

    n = pypsa.Network(snakemake.input.network)

    if options["only_elec_network"]:
        carriers = carrier_to_keep + extend_carrier_list(config_elec)
        n = strip_network(n, carriers)

    include_electricity_growth(n, options, planning_horizons)

    if options["redistribute_industry"]:
        redistribute_industrial_load(n, options["industry_share"])

    if options["readjust_existing_interconnections"]:
        readjust_existing_interconnections(n)

    n.export_to_netcdf(snakemake.output[0])
