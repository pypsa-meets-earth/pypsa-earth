# -*- coding: utf-8 -*-
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from helpers import (
    create_dummy_data,
    create_network_topology,
    cycling_shift,
    locate_bus,
    mock_snakemake,
    override_component_attrs,
    prepare_costs,
    three_2_two_digits_country,
)
from prepare_transport_data import prepare_transport_data

spatial = SimpleNamespace()


def add_carrier_buses(n, carriers):
    """
    Add buses to connect e.g. coal, nuclear and oil plants
    """
    if isinstance(carriers, str):
        carriers = [carriers]

    for carrier in carriers:

        n.add("Carrier", carrier)

        n.add("Bus", "Africa " + carrier, location="Africa", carrier=carrier)

        # capital cost could be corrected to e.g. 0.2 EUR/kWh * annuity and O&M
        n.add(
            "Store",
            "Africa " + carrier + " Store",
            bus="Africa " + carrier,
            e_nom_extendable=True,
            e_cyclic=True,
            carrier=carrier,
        )

        n.add(
            "Generator",
            "Africa " + carrier,
            bus="Africa " + carrier,
            p_nom_extendable=True,
            carrier=carrier,
            marginal_cost=costs.at[carrier, "fuel"],
        )


def add_generation(n, costs):
    """Adds conventional generation as specified in config


    Args:
        n (network): PyPSA prenetwork
        costs (dataframe): _description_

    Returns:
        _type_: _description_
    """ """"""

    print("adding electricity generation")

    # Not required, because nodes are already defined in "nodes"
    # nodes = pop_layout.index

    fallback = {"OCGT": "gas"}
    conventionals = options.get("conventional_generation", fallback)

    add_carrier_buses(n, np.unique(list(conventionals.values())))

    for generator, carrier in conventionals.items():

        n.madd(
            "Link",
            nodes + " " + generator,
            bus0="Africa " + carrier,
            bus1=nodes,
            bus2="co2 atmosphere",
            marginal_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "VOM"],  # NB: VOM is per MWel
            # NB: fixed cost is per MWel
            capital_cost=costs.at[generator, "efficiency"]
            * costs.at[generator, "fixed"],
            p_nom_extendable=True,
            carrier=generator,
            efficiency=costs.at[generator, "efficiency"],
            efficiency2=costs.at[carrier, "CO2 intensity"],
            lifetime=costs.at[generator, "lifetime"],
        )


def add_oil(n, costs):
    """
    Function to add oil carrier and bus to network. If-Statements are required in
    case oil was already added from config ['sector']['conventional_generation']
    Oil is copper plated
    """
    # TODO function will not be necessary if conventionals are added using "add_carrier_buses()"
    # TODO before using add_carrier_buses: remove_elec_base_techs(n), otherwise carriers are added double

    if "oil" not in n.carriers.index:
        n.add("Carrier", "oil")

    if "Africa oil" not in n.buses.index:

        n.add("Bus", "Africa oil", location="Africa", carrier="oil")

    if "Africa oil Store" not in n.stores.index:

        # could correct to e.g. 0.001 EUR/kWh * annuity and O&M
        n.add(
            "Store",
            "Africa oil Store",
            bus="Africa oil",
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="oil",
        )

    if "Africa oil" not in n.generators.index:

        n.add(
            "Generator",
            "Africa oil",
            bus="Africa oil",
            p_nom_extendable=True,
            carrier="oil",
            marginal_cost=costs.at["oil", "fuel"],
        )


def H2_liquid_fossil_conversions(n, costs):
    """
    Function to add conversions between H2 and liquid fossil
    Carrier and bus is added in add_oil, which later on might be switched to add_generation
    """

    n.madd(
        "Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1="Africa oil",
        bus2=spatial.co2.nodes,
        carrier="Fischer-Tropsch",
        efficiency=costs.at["Fischer-Tropsch", "efficiency"],
        capital_cost=costs.at["Fischer-Tropsch", "fixed"],
        efficiency2=-costs.at["oil", "CO2 intensity"]
        * costs.at["Fischer-Tropsch", "efficiency"],
        p_nom_extendable=True,
        lifetime=costs.at["Fischer-Tropsch", "lifetime"],
    )


def add_hydrogen(n, costs):
    "function to add hydrogen as an energy carrier with its conversion technologies from and to AC"

    n.add("Carrier", "H2")

    n.madd(
        "Bus",
        nodes + " H2",
        location=nodes,
        carrier="H2",
        x=n.buses.loc[list(nodes)].x.values,
        y=n.buses.loc[list(nodes)].y.values,
    )

    n.madd(
        "Link",
        nodes + " H2 Electrolysis",
        bus1=nodes + " H2",
        bus0=nodes,
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "fixed"],
        lifetime=costs.at["electrolysis", "lifetime"],
    )

    n.madd(
        "Link",
        nodes + " H2 Fuel Cell",
        bus0=nodes + " H2",
        bus1=nodes,
        p_nom_extendable=True,
        carrier="H2 Fuel Cell",
        efficiency=costs.at["fuel cell", "efficiency"],
        # NB: fixed cost is per MWel
        capital_cost=costs.at["fuel cell", "fixed"]
        * costs.at["fuel cell", "efficiency"],
        lifetime=costs.at["fuel cell", "lifetime"],
    )

    cavern_nodes = pd.DataFrame()
    if options["hydrogen_underground_storage"]:

        h2_salt_cavern_potential = pd.read_csv(
            snakemake.input.h2_cavern, index_col=0, squeeze=True
        )
        h2_cavern_ct = h2_salt_cavern_potential[~h2_salt_cavern_potential.isna()]
        cavern_nodes = n.buses[n.buses.country.isin(h2_cavern_ct.index)]

        h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

        # assumptions: weight storage potential in a country by population
        # TODO: fix with real geographic potentials
        # convert TWh to MWh with 1e6
        h2_pot = h2_cavern_ct.loc[cavern_nodes.country]
        h2_pot.index = cavern_nodes.index
        # h2_pot = h2_pot * cavern_nodes.fraction * 1e6

        n.madd(
            "Store",
            cavern_nodes.index + " H2 Store",
            bus=cavern_nodes.index + " H2",
            e_nom_extendable=True,
            e_nom_max=h2_pot.values,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost,
        )

    # hydrogen stored overground (where not already underground)
    h2_capital_cost = costs.at["hydrogen storage tank incl. compressor", "fixed"]
    nodes_overground = cavern_nodes.index.symmetric_difference(nodes)

    n.madd(
        "Store",
        nodes_overground + " H2 Store",
        bus=nodes_overground + " H2",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="H2 Store",
        capital_cost=h2_capital_cost,
    )

    attrs = ["bus0", "bus1", "length"]
    h2_links = pd.DataFrame(columns=attrs)

    candidates = pd.concat(
        {"lines": n.lines[attrs], "links": n.links.loc[n.links.carrier == "DC", attrs]}
    )

    for candidate in candidates.index:
        buses = [candidates.at[candidate, "bus0"], candidates.at[candidate, "bus1"]]
        buses.sort()
        name = f"H2 pipeline {buses[0]} -> {buses[1]}"
        if name not in h2_links.index:
            h2_links.at[name, "bus0"] = buses[0]
            h2_links.at[name, "bus1"] = buses[1]
            h2_links.at[name, "length"] = candidates.at[candidate, "length"]

    # TODO Add efficiency losses
    n.madd(
        "Link",
        h2_links.index,
        bus0=h2_links.bus0.values + " H2",
        bus1=h2_links.bus1.values + " H2",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=h2_links.length.values,
        capital_cost=costs.at["H2 (g) pipeline", "fixed"] * h2_links.length.values,
        carrier="H2 pipeline",
        lifetime=costs.at["H2 (g) pipeline", "lifetime"],
    )


def add_co2(n, costs):
    "add carbon carrier, it's networks and storage units"
    spatial.nodes = nodes

    spatial.co2 = SimpleNamespace()

    if options["co2_network"]:
        spatial.co2.nodes = nodes + " co2 stored"
        spatial.co2.locations = nodes
        spatial.co2.vents = nodes + " co2 vent"
        spatial.co2.x = (n.buses.loc[list(nodes)].x.values,)
        spatial.co2.y = (n.buses.loc[list(nodes)].y.values,)
    else:
        spatial.co2.nodes = ["co2 stored"]
        spatial.co2.locations = ["Africa"]
        spatial.co2.vents = ["co2 vent"]
        spatial.co2.x = (0,)
        spatial.co2.y = 0

    spatial.co2.df = pd.DataFrame(vars(spatial.co2), index=nodes)

    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2", co2_emissions=-1.0)

    # this tracks CO2 in the atmosphere
    n.add(
        "Bus",
        "co2 atmosphere",
        location="Africa",  # TODO Ignoed by pypsa chck
        carrier="co2",
    )

    # can also be negative
    n.add(
        "Store",
        "co2 atmosphere",
        e_nom_extendable=True,
        e_min_pu=-1,
        carrier="co2",
        bus="co2 atmosphere",
    )

    # this tracks CO2 stored, e.g. underground
    n.madd(
        "Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored",
        x=spatial.co2.x[0],
        y=spatial.co2.y[0],
    )

    co2_stored_x = n.buses.filter(like="co2 stored", axis=0).loc[:, "x"]
    co2_stored_y = n.buses.loc[n.buses[n.buses.carrier == "co2 stored"].location].y

    n.buses[n.buses.carrier == "co2 stored"].x = co2_stored_x.values
    n.buses[n.buses.carrier == "co2 stored"].y = co2_stored_y.values

    n.madd(
        "Store",
        spatial.co2.nodes.str[:-2] + "age",
        e_nom_extendable=True,
        e_nom_max=np.inf,
        capital_cost=options["co2_sequestration_cost"],
        carrier="co2 stored",
        bus=spatial.co2.nodes,
    )

    n.madd(
        "Link",
        spatial.co2.vents,
        bus0=spatial.co2.nodes,
        bus1="co2 atmosphere",
        carrier="co2 vent",
        efficiency=1.0,
        p_nom_extendable=True,
    )

    # logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (
        (1 - co2_links.underwater_fraction)
        * costs.at["CO2 pipeline", "fixed"]
        * co2_links.length
    )
    cost_submarine = (
        co2_links.underwater_fraction
        * costs.at["CO2 submarine pipeline", "fixed"]
        * co2_links.length
    )
    capital_cost = cost_onshore + cost_submarine

    n.madd(
        "Link",
        co2_links.index,
        bus0=co2_links.bus0.values + " co2 stored",
        bus1=co2_links.bus1.values + " co2 stored",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=co2_links.length.values,
        capital_cost=capital_cost.values,
        carrier="CO2 pipeline",
        lifetime=costs.at["CO2 pipeline", "lifetime"],
    )

    n.madd(
        "Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        e_nom_max=np.inf,
        capital_cost=options["co2_sequestration_cost"],
        carrier="co2 stored",
        bus=spatial.co2.nodes,
    )

    # logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (
        (1 - co2_links.underwater_fraction)
        * costs.at["CO2 pipeline", "fixed"]
        * co2_links.length
    )
    cost_submarine = (
        co2_links.underwater_fraction
        * costs.at["CO2 submarine pipeline", "fixed"]
        * co2_links.length
    )
    capital_cost = cost_onshore + cost_submarine


def add_aviation(n, cost):
    # all_aviation = ["total international aviation", "total domestic aviation"]

    airports = pd.read_csv(snakemake.input.airports, index_col=None, squeeze=True)

    gadm_level = options["gadm_level"]

    airports["gadm_{}".format(gadm_level)] = airports[["x", "y", "country"]].apply(
        lambda airport: locate_bus(airport[["x", "y"]], airport["country"], gadm_level),
        axis=1,
    )

    # To change 3 country code to 2
    # airports["gadm_{}".format(gadm_level)] = airports["gadm_{}".format(gadm_level)].apply(
    # lambda cocode: three_2_two_digits_country(cocode[:3]) + " " + cocode[4:-2])

    airports = airports.set_index("gadm_{}".format(gadm_level))

    ind = pd.DataFrame(n.buses.index[n.buses.carrier == "AC"])

    ind = ind.set_index(n.buses.index[n.buses.carrier == "AC"])

    airports["p_set"] = airports["fraction"].apply(
        lambda frac: frac * 1e6 / 8760
    )  # TODO use real data here

    airports = pd.concat([airports, ind])

    airports = airports[~airports.index.duplicated(keep="first")]

    airports = airports.fillna(0)

    n.madd(
        "Load",
        nodes,
        suffix=" kerosene for aviation",
        bus="Africa oil",
        carrier="kerosene for aviation",
        p_set=airports["p_set"],
    )

    # co2_release = ["kerosene for aviation"]
    co2 = airports["p_set"].sum() * costs.at["oil", "CO2 intensity"] / 8760

    n.add(
        "Load",
        "oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2,
    )


def add_storage(n, costs):
    "function to add the different types of storage systems"
    n.add("Carrier", "battery")

    n.madd(
        "Bus",
        nodes + " battery",
        location=nodes,
        carrier="battery",
        x=n.buses.loc[list(nodes)].x.values,
        y=n.buses.loc[list(nodes)].y.values,
    )

    n.madd(
        "Store",
        nodes + " battery",
        bus=nodes + " battery",
        e_cyclic=True,
        e_nom_extendable=True,
        carrier="battery",
        capital_cost=costs.at["battery storage", "fixed"],
        lifetime=costs.at["battery storage", "lifetime"],
    )

    n.madd(
        "Link",
        nodes + " battery charger",
        bus0=nodes,
        bus1=nodes + " battery",
        carrier="battery charger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        capital_cost=costs.at["battery inverter", "fixed"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )

    n.madd(
        "Link",
        nodes + " battery discharger",
        bus0=nodes + " battery",
        bus1=nodes,
        carrier="battery discharger",
        efficiency=costs.at["battery inverter", "efficiency"] ** 0.5,
        marginal_cost=options["marginal_cost_storage"],
        p_nom_extendable=True,
        lifetime=costs.at["battery inverter", "lifetime"],
    )


def h2_hc_conversions(n, costs):
    "function to add the conversion technologies between H2 and hydrocarbons"
    if options["methanation"]:

        n.madd(
            "Link",
            spatial.nodes,
            suffix=" Sabatier",
            bus0=nodes + " H2",
            bus1="Africa gas",
            bus2=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="Sabatier",
            efficiency=costs.at["methanation", "efficiency"],
            efficiency2=-costs.at["methanation", "efficiency"]
            * costs.at["gas", "CO2 intensity"],
            # costs given per kW_gas
            capital_cost=costs.at["methanation", "fixed"]
            * costs.at["methanation", "efficiency"],
            lifetime=costs.at["methanation", "lifetime"],
        )

    if options["helmeth"]:

        n.madd(
            "Link",
            spatial.nodes,
            suffix=" helmeth",
            bus0=nodes,
            bus1="Africa gas",
            bus2=spatial.co2.nodes,
            carrier="helmeth",
            p_nom_extendable=True,
            efficiency=costs.at["helmeth", "efficiency"],
            efficiency2=-costs.at["helmeth", "efficiency"]
            * costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["helmeth", "fixed"],
            lifetime=costs.at["helmeth", "lifetime"],
        )

    if options["SMR"]:

        n.madd(
            "Link",
            spatial.nodes,
            suffix=" SMR CC",
            bus0="Africa gas",
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="SMR CC",
            efficiency=costs.at["SMR CC", "efficiency"],
            efficiency2=costs.at["gas", "CO2 intensity"] * (1 - options["cc_fraction"]),
            efficiency3=costs.at["gas", "CO2 intensity"] * options["cc_fraction"],
            capital_cost=costs.at["SMR CC", "fixed"],
            lifetime=costs.at["SMR CC", "lifetime"],
        )

        n.madd(
            "Link",
            nodes + " SMR",
            bus0="Africa gas",
            bus1=nodes + " H2",
            bus2="co2 atmosphere",
            p_nom_extendable=True,
            carrier="SMR",
            efficiency=costs.at["SMR", "efficiency"],
            efficiency2=costs.at["gas", "CO2 intensity"],
            capital_cost=costs.at["SMR", "fixed"],
            lifetime=costs.at["SMR", "lifetime"],
        )


def add_shipping(n, costs):

    gadm_level = options["gadm_level"]

    efficiency = (
        options["shipping_average_efficiency"] / costs.at["fuel cell", "efficiency"]
    )

    # check whether item depends on investment year
    shipping_hydrogen_share = get(options["shipping_hydrogen_share"], investment_year)

    ports = pd.read_csv(snakemake.input.ports, index_col=None, squeeze=True)

    ports["gadm_{}".format(gadm_level)] = ports[["x", "y", "country"]].apply(
        lambda port: locate_bus(port[["x", "y"]], port["country"], gadm_level), axis=1
    )

    ports = ports.set_index("gadm_{}".format(gadm_level))

    ind = pd.DataFrame(n.buses.index[n.buses.carrier == "AC"])

    ind = ind.set_index(n.buses.index[n.buses.carrier == "AC"])

    ports["p_set"] = ports["fraction"].apply(
        lambda frac: shipping_hydrogen_share * frac * 1e6 * efficiency / 8760
    )  # TODO use real data here

    ports = pd.concat([ports, ind])

    ports = ports[~ports.index.duplicated(keep="first")]

    # TODO Where is MAR.7_1 ????? Only 13 nodes

    ports = ports.fillna(0)

    if options["shipping_hydrogen_liquefaction"]:  # how to implement options?

        n.madd("Bus", nodes, suffix=" H2 liquid", carrier="H2 liquid", location=nodes)

        # link the H2 supply to liquified H2
        n.madd(
            "Link",
            nodes + " H2 liquefaction",
            bus0=nodes + " H2",
            bus1=nodes + " H2 liquid",
            carrier="H2 liquefaction",
            efficiency=costs.at["H2 liquefaction", "efficiency"],
            capital_cost=costs.at["H2 liquefaction", "fixed"],
            p_nom_extendable=True,
            lifetime=costs.at["H2 liquefaction", "lifetime"],
        )

        shipping_bus = nodes + " H2 liquid"
    else:
        shipping_bus = nodes + " H2"

    n.madd(
        "Load",
        nodes,
        suffix=" H2 for shipping",
        bus=shipping_bus,
        carrier="H2 for shipping",
        p_set=ports["p_set"],
    )

    if shipping_hydrogen_share < 1:

        shipping_oil_share = 1 - shipping_hydrogen_share

        ports["p_set"] = ports["fraction"].apply(
            lambda frac: shipping_oil_share * frac * 1e6 / 8760
        )

        n.madd(
            "Load",
            nodes,
            suffix=" shipping oil",
            bus="Africa oil",
            carrier="shipping oil",
            p_set=ports["p_set"],
        )

        co2 = (
            shipping_oil_share
            * ports["p_set"].sum()
            * 1e6
            / 8760
            * costs.at["oil", "CO2 intensity"]
        )

        n.add(
            "Load",
            "shipping oil emissions",
            bus="co2 atmosphere",
            carrier="shipping oil emissions",
            p_set=-co2,
        )

    if "Africa oil" not in n.buses.index:

        n.add("Bus", "Africa oil", location="Africa", carrier="oil")

    if "Africa oil Store" not in n.stores.index:

        # could correct to e.g. 0.001 EUR/kWh * annuity and O&M
        n.add(
            "Store",
            "Africa oil Store",
            bus="Africa oil",
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="oil",
        )

    if "Africa oil" not in n.generators.index:

        n.add(
            "Generator",
            "Africa oil",
            bus="Africa oil",
            p_nom_extendable=True,
            carrier="oil",
            marginal_cost=costs.at["oil", "fuel"],
        )


def add_industry(n, costs):

    #     print("adding industrial demand")

    #     # 1e6 to convert TWh to MWh
    #     industrial_demand = pd.read_csv(snakemake.input.industrial_demand, index_col=0) * 1e6
    industrial_demand = create_dummy_data(n, "industry", "")

    # TODO carrier Biomass

    # CARRIER = FOSSIL GAS

    demand_locations = pd.read_csv(
        snakemake.input.industry_demands, index_col=None, squeeze=True
    )

    gadm_level = options["gadm_level"]
    # carrier = "kerosene for aviation"

    demand_locations["gadm_{}".format(gadm_level)] = demand_locations[
        ["x", "y", "country"]
    ].apply(lambda loc: locate_bus(loc[["x", "y"]], loc["country"], gadm_level), axis=1)

    demand_locations = demand_locations.set_index("gadm_{}".format(gadm_level))

    ind = pd.DataFrame(n.buses.index[n.buses.carrier == "AC"])

    ind = ind.set_index(n.buses.index[n.buses.carrier == "AC"])

    # demand_locations["p_set_{}".format(carrier)] = demand_locations["fraction"].apply(
    #        lambda frac: frac * 1e6 / 8760)

    demand_locations = pd.concat([demand_locations, ind])

    demand_locations = demand_locations[
        ~demand_locations.index.duplicated(keep="first")
    ]

    demand_locations = demand_locations.fillna(0)

    n.add("Bus", "gas for industry", location="Africa", carrier="gas for industry")

    n.madd(
        "Load",
        nodes,
        suffix=" gas for industry",
        bus="gas for industry",
        carrier="gas for industry",
        p_set=demand_locations["gas fraction"].apply(
            lambda frac: frac * 1e6 / 8760
        ),  # TODO Multiply by gas demand and find true gas fraction
    )

    n.add(
        "Link",
        "gas for industry",
        bus0="Africa gas",
        bus1="gas for industry",
        bus2="co2 atmosphere",
        carrier="gas for industry",
        p_nom_extendable=True,
        efficiency=1.0,
        efficiency2=costs.at["gas", "CO2 intensity"],
    )

    n.madd(
        "Link",
        spatial.co2.locations,
        suffix=" gas for industry CC",
        bus0="Africa gas",
        bus1="gas for industry",
        bus2="co2 atmosphere",
        bus3=spatial.co2.nodes,
        carrier="gas for industry CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"]
        * costs.at["gas", "CO2 intensity"],
        efficiency=0.9,
        efficiency2=costs.at["gas", "CO2 intensity"]
        * (1 - costs.at["cement capture", "capture_rate"]),
        efficiency3=costs.at["gas", "CO2 intensity"]
        * costs.at["cement capture", "capture_rate"],
        lifetime=costs.at["cement capture", "lifetime"],
    )

    #################################################### CARRIER = HYDROGEN
    n.madd(
        "Load",
        nodes,
        suffix=" H2 for industry",
        bus=nodes + " H2",
        carrier="H2 for industry",
        p_set=demand_locations["H2 fraction"].apply(
            lambda frac: frac * 1e6 / 8760
        ),  # TODO Multiply by gas demand and find true gas fraction
    )

    # CARRIER = LIQUID HYDROCARBONS
    n.madd(
        "Load",
        nodes,
        suffix=" naphtha for industry",
        bus="Africa oil",
        carrier="naphtha for industry",
        p_set=demand_locations["naphta fraction"].apply(
            lambda frac: frac * 1e6 / 8760
        ),  # TODO Multiply by gas demand and find true gas fraction
    )

    #     #NB: CO2 gets released again to atmosphere when plastics decay or kerosene is burned
    #     #except for the process emissions when naphtha is used for petrochemicals, which can be captured with other industry process emissions
    #     #tco2 per hour
    # TODO kerosene for aviation should be added too but in the right func.
    co2_release = [" naphtha for industry"]
    # check land tranport

    co2 = (
        n.loads.loc[nodes + co2_release, "p_set"].sum()
        * costs.at["oil", "CO2 intensity"]
        - demand_locations["feedstock emissions fraction"].sum() / 8760
    )

    n.add(
        "Load",
        "industry oil emissions",
        bus="co2 atmosphere",
        carrier="industry oil emissions",
        p_set=-co2,
    )

    # co2 = airports["p_set"].sum() * costs.at["oil", 'CO2 intensity'] / 8760

    # n.add("Load",
    #    "oil emissions",
    #    bus="co2 atmosphere",
    #    carrier="oil emissions",
    #    p_set=-co2
    # )

    ########################################################### CARIER = HEAT
    # TODO simplify bus expression
    n.madd(
        "Load",
        nodes,
        suffix=" low-temperature heat for industry",
        bus=[
            node + " urban central heat"
            if node + " urban central heat" in n.buses.index
            else node + " services urban decentral heat"
            for node in nodes
        ],
        carrier="low-temperature heat for industry",
        p_set=industrial_demand.loc[nodes, "low-temperature heat"] / 8760,
    )

    ################################################## CARRIER = ELECTRICITY

    #     # remove today's industrial electricity demand by scaling down total electricity demand
    for ct in n.buses.country.dropna().unique():
        # TODO map onto n.bus.country
        # TODO make sure to check this one, should AC have carrier pf "electricity"?
        loads_i = n.loads.index[
            (n.loads.index.str[:2] == ct) & (n.loads.carrier == "electricity")
        ]
        if n.loads_t.p_set[loads_i].empty:
            continue
        factor = (
            1
            - industrial_demand.loc[loads_i, "current electricity"].sum()
            / n.loads_t.p_set[loads_i].sum().sum()
        )
        n.loads_t.p_set[loads_i] *= factor

    n.madd(
        "Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=demand_locations["electricity fraction"].apply(
            lambda frac: frac * 1e6 / 8760
        ),  # TODO Multiply by demand and find true fraction
    )

    n.add("Bus", "process emissions", location="Africa", carrier="process emissions")

    # this should be process emissions fossil+feedstock
    # then need load on atmosphere for feedstock emissions that are currently going to atmosphere via Link Fischer-Tropsch demand
    n.madd(
        "Load",
        nodes,
        suffix=" process emissions",
        bus="process emissions",
        carrier="process emissions",
        p_set=-(
            demand_locations["feedstock emissions fraction"]
            + demand_locations["process emissions fraction"]
        )
        / 8760,  # TODO multiply by total emissions
    )

    n.add(
        "Link",
        "process emissions",
        bus0="process emissions",
        bus1="co2 atmosphere",
        carrier="process emissions",
        p_nom_extendable=True,
        efficiency=1.0,
    )

    # assume enough local waste heat for CC
    n.madd(
        "Link",
        spatial.co2.locations,
        suffix=" process emissions CC",
        bus0="process emissions",
        bus1="co2 atmosphere",
        bus2=spatial.co2.nodes,
        carrier="process emissions CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"],
        efficiency=1 - costs.at["cement capture", "capture_rate"],
        efficiency2=costs.at["cement capture", "capture_rate"],
        lifetime=costs.at["cement capture", "lifetime"],
    )


def get(item, investment_year=None):
    """Check whether item depends on investment year"""
    if isinstance(item, dict):
        return item[investment_year]
    else:
        return item


"""
Missing data:
 - transport
 - aviation data
 - nodal_transport_data
 - cycling_shift
 - dsm_profile
 - avail_profile
"""


def add_land_transport(n, costs):
    """
    Function to add land transport to network
    """
    # TODO options?

    print("adding land transport")

    fuel_cell_share = get(options["land_transport_fuel_cell_share"], investment_year)
    electric_share = get(options["land_transport_electric_share"], investment_year)
    ice_share = 1 - fuel_cell_share - electric_share

    print("FCEV share", fuel_cell_share)
    print("EV share", electric_share)
    print("ICEV share", ice_share)

    assert ice_share >= 0, "Error, more FCEV and EV share than 1."

    # Nodes are already defined, remove it from here
    # nodes = pop_layout.index

    if electric_share > 0:

        n.add("Carrier", "Li ion")

        n.madd(
            "Bus",
            nodes,
            location=nodes,
            suffix=" EV battery",
            carrier="Li ion",
            x=n.buses.loc[list(nodes)].x.values,
            y=n.buses.loc[list(nodes)].y.values,
        )

        p_set = (
            electric_share
            * (
                transport[nodes]
                + cycling_shift(transport[nodes], 1)
                + cycling_shift(transport[nodes], 2)
            )
            / 3
        )

        n.madd(
            "Load",
            nodes,
            suffix=" land transport EV",
            bus=nodes + " EV battery",
            carrier="land transport EV",
            p_set=p_set,
        )

        p_nom = (
            nodal_transport_data["number cars"]
            * options.get("bev_charge_rate", 0.011)
            * electric_share
        )

        n.madd(
            "Link",
            nodes,
            suffix=" BEV charger",
            bus0=nodes,
            bus1=nodes + " EV battery",
            p_nom=p_nom,
            carrier="BEV charger",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
            # These were set non-zero to find LU infeasibility when availability = 0.25
            # p_nom_extendable=True,
            # p_nom_min=p_nom,
            # capital_cost=1e6,  #i.e. so high it only gets built where necessary
        )

    if electric_share > 0 and options["v2g"]:

        n.madd(
            "Link",
            nodes,
            suffix=" V2G",
            bus1=nodes,
            bus0=nodes + " EV battery",
            p_nom=p_nom,
            carrier="V2G",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
        )

    if electric_share > 0 and options["bev_dsm"]:

        e_nom = (
            nodal_transport_data["number cars"]
            * options.get("bev_energy", 0.05)
            * options["bev_availability"]
            * electric_share
        )

        n.madd(
            "Store",
            nodes,
            suffix=" battery storage",
            bus=nodes + " EV battery",
            carrier="battery storage",
            e_cyclic=True,
            e_nom=e_nom,
            e_max_pu=1,
            e_min_pu=dsm_profile[nodes],
        )

    if fuel_cell_share > 0:

        n.madd(
            "Load",
            nodes,
            suffix=" land transport fuel cell",
            bus=nodes + " H2",
            carrier="land transport fuel cell",
            p_set=fuel_cell_share
            / options["transport_fuel_cell_efficiency"]
            * transport[nodes],
        )

    if ice_share > 0:

        if "Africa oil" not in n.buses.index:
            n.add("Bus", "Africa oil", location="Africa", carrier="oil")

        ice_efficiency = options["transport_internal_combustion_efficiency"]

        n.madd(
            "Load",
            nodes,
            suffix=" land transport oil",
            bus="Africa oil",
            carrier="land transport oil",
            p_set=ice_share / ice_efficiency * transport[nodes],
        )

        co2 = (
            ice_share
            / ice_efficiency
            * transport[nodes].sum().sum()
            / 8760
            * costs.at["oil", "CO2 intensity"]
        )

        n.add(
            "Load",
            "land transport oil emissions",
            bus="co2 atmosphere",
            carrier="land transport oil emissions",
            p_set=-co2,
        )


def create_nodes_for_heat_sector():
    # TODO pop_layout

    # rural are areas with low heating density and individual heating
    # urban are areas with high heating density
    # urban can be split into district heating (central) and individual heating (decentral)

    ct_urban = pop_layout.urban.groupby(pop_layout.ct).sum()
    # distribution of urban population within a country
    pop_layout["urban_ct_fraction"] = pop_layout.urban / pop_layout.ct.map(ct_urban.get)

    sectors = ["residential", "services"]

    nodes = {}
    urban_fraction = pop_layout.urban / pop_layout[["rural", "urban"]].sum(axis=1)

    for sector in sectors:
        nodes[sector + " rural"] = pop_layout.index
        nodes[sector + " urban decentral"] = pop_layout.index

    # maximum potential of urban demand covered by district heating
    central_fraction = options["district_heating"]["potential"]
    # district heating share at each node
    dist_fraction_node = (
        district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    )
    nodes["urban central"] = dist_fraction_node.index
    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node], axis=1).max(axis=1)
    # difference of max potential and today's share of district heating
    diff = (urban_fraction * central_fraction) - dist_fraction_node
    progress = get(options["district_heating"]["progress"], investment_year)
    dist_fraction_node += diff * progress
    print(
        "The current district heating share compared to the maximum",
        f"possible is increased by a progress factor of\n{progress}",
        f"resulting in a district heating share of\n{dist_fraction_node}",
    )

    return nodes, dist_fraction_node, urban_fraction


def add_heat(n, costs):
    # TODO options?
    # TODO pop_layout?

    print("adding heat")

    sectors = ["residential", "services"]

    nodes, dist_fraction, urban_fraction = create_nodes_for_heat_sector()

    # NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get(options["reduce_space_heat_exogenously_factor"], investment_year)
        print(f"assumed space heat reduction of {dE*100} %")
        for sector in sectors:
            heat_demand[sector + " space"] = (1 - dE) * heat_demand[sector + " space"]

    heat_systems = [
        "residential rural",
        "services rural",
        "residential urban decentral",
        "services urban decentral",
        "urban central",
    ]

    for name in heat_systems:

        name_type = "central" if name == "urban central" else "decentral"

        n.add("Carrier", name + " heat")

        n.madd(
            "Bus",
            nodes[name] + " {} heat".format(name),
            location=nodes[name],
            carrier=name + " heat",
        )

        ## Add heat load

        for sector in sectors:
            # heat demand weighting
            if "rural" in name:
                factor = 1 - urban_fraction[nodes[name]]
            elif "urban central" in name:
                factor = dist_fraction[nodes[name]]
            elif "urban decentral" in name:
                factor = urban_fraction[nodes[name]] - dist_fraction[nodes[name]]
            else:
                raise NotImplementedError(
                    f" {name} not in " f"heat systems: {heat_systems}"
                )

            if sector in name:
                heat_load = (
                    heat_demand[[sector + " water", sector + " space"]]
                    .groupby(level=1, axis=1)
                    .sum()[nodes[name]]
                    .multiply(factor)
                )

        if name == "urban central":
            heat_load = (
                heat_demand.groupby(level=1, axis=1)
                .sum()[nodes[name]]
                .multiply(
                    factor * (1 + options["district_heating"]["district_heating_loss"])
                )
            )

        n.madd(
            "Load",
            nodes[name],
            suffix=f" {name} heat",
            bus=nodes[name] + f" {name} heat",
            carrier=name + " heat",
            p_set=heat_load,
        )

        ## Add heat pumps

        heat_pump_type = "air" if "urban" in name else "ground"

        costs_name = f"{name_type} {heat_pump_type}-sourced heat pump"
        cop = {"air": ashp_cop, "ground": gshp_cop}
        efficiency = (
            cop[heat_pump_type][nodes[name]]
            if options["time_dep_hp_cop"]
            else costs.at[costs_name, "efficiency"]
        )

        n.madd(
            "Link",
            nodes[name],
            suffix=f" {name} {heat_pump_type} heat pump",
            bus0=nodes[name],
            bus1=nodes[name] + f" {name} heat",
            carrier=f"{name} {heat_pump_type} heat pump",
            efficiency=efficiency,
            capital_cost=costs.at[costs_name, "efficiency"]
            * costs.at[costs_name, "fixed"],
            p_nom_extendable=True,
            lifetime=costs.at[costs_name, "lifetime"],
        )

        if options["tes"]:

            n.add("Carrier", name + " water tanks")

            n.madd(
                "Bus",
                nodes[name] + f" {name} water tanks",
                location=nodes[name],
                carrier=name + " water tanks",
            )

            n.madd(
                "Link",
                nodes[name] + f" {name} water tanks charger",
                bus0=nodes[name] + f" {name} heat",
                bus1=nodes[name] + f" {name} water tanks",
                efficiency=costs.at["water tank charger", "efficiency"],
                carrier=name + " water tanks charger",
                p_nom_extendable=True,
            )

            n.madd(
                "Link",
                nodes[name] + f" {name} water tanks discharger",
                bus0=nodes[name] + f" {name} water tanks",
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " water tanks discharger",
                efficiency=costs.at["water tank discharger", "efficiency"],
                p_nom_extendable=True,
            )

            if isinstance(options["tes_tau"], dict):
                tes_time_constant_days = options["tes_tau"][name_type]
            else:  # TODO add logger
                # logger.warning("Deprecated: a future version will require you to specify 'tes_tau' ",
                # "for 'decentral' and 'central' separately.")
                tes_time_constant_days = (
                    options["tes_tau"] if name_type == "decentral" else 180.0
                )

            # conversion from EUR/m^3 to EUR/MWh for 40 K diff and 1.17 kWh/m^3/K
            capital_cost = (
                costs.at[name_type + " water tank storage", "fixed"] / 0.00117 / 40
            )

            n.madd(
                "Store",
                nodes[name] + f" {name} water tanks",
                bus=nodes[name] + f" {name} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=name + " water tanks",
                standing_loss=1 - np.exp(-1 / 24 / tes_time_constant_days),
                capital_cost=capital_cost,
                lifetime=costs.at[name_type + " water tank storage", "lifetime"],
            )

        if options["boilers"]:

            key = f"{name_type} resistive heater"

            n.madd(
                "Link",
                nodes[name] + f" {name} resistive heater",
                bus0=nodes[name],
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " resistive heater",
                efficiency=costs.at[key, "efficiency"],
                capital_cost=costs.at[key, "efficiency"] * costs.at[key, "fixed"],
                p_nom_extendable=True,
                lifetime=costs.at[key, "lifetime"],
            )

            key = f"{name_type} gas boiler"

            n.madd(
                "Link",
                nodes[name] + f" {name} gas boiler",
                p_nom_extendable=True,
                bus0="EU gas",
                bus1=nodes[name] + f" {name} heat",
                bus2="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[key, "efficiency"],
                efficiency2=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at[key, "efficiency"] * costs.at[key, "fixed"],
                lifetime=costs.at[key, "lifetime"],
            )

        if options["solar_thermal"]:

            n.add("Carrier", name + " solar thermal")

            n.madd(
                "Generator",
                nodes[name],
                suffix=f" {name} solar thermal collector",
                bus=nodes[name] + f" {name} heat",
                carrier=name + " solar thermal",
                p_nom_extendable=True,
                capital_cost=costs.at[name_type + " solar thermal", "fixed"],
                p_max_pu=solar_thermal[nodes[name]],
                lifetime=costs.at[name_type + " solar thermal", "lifetime"],
            )

        if options["chp"] and name == "urban central":

            # add gas CHP; biomass CHP is added in biomass section
            n.madd(
                "Link",
                nodes[name] + " urban central gas CHP",
                bus0="EU gas",
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                carrier="urban central gas CHP",
                p_nom_extendable=True,
                capital_cost=costs.at["central gas CHP", "fixed"]
                * costs.at["central gas CHP", "efficiency"],
                marginal_cost=costs.at["central gas CHP", "VOM"],
                efficiency=costs.at["central gas CHP", "efficiency"],
                efficiency2=costs.at["central gas CHP", "efficiency"]
                / costs.at["central gas CHP", "c_b"],
                efficiency3=costs.at["gas", "CO2 intensity"],
                lifetime=costs.at["central gas CHP", "lifetime"],
            )

            n.madd(
                "Link",
                nodes[name] + " urban central gas CHP CC",
                bus0="EU gas",
                bus1=nodes[name],
                bus2=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                bus4=spatial.co2.df.loc[nodes[name], "nodes"].values,
                carrier="urban central gas CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at["central gas CHP", "fixed"]
                * costs.at["central gas CHP", "efficiency"]
                + costs.at["biomass CHP capture", "fixed"]
                * costs.at["gas", "CO2 intensity"],
                marginal_cost=costs.at["central gas CHP", "VOM"],
                efficiency=costs.at["central gas CHP", "efficiency"]
                - costs.at["gas", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "electricity-input"]
                    + costs.at["biomass CHP capture", "compression-electricity-input"]
                ),
                efficiency2=costs.at["central gas CHP", "efficiency"]
                / costs.at["central gas CHP", "c_b"]
                + costs.at["gas", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "heat-output"]
                    + costs.at["biomass CHP capture", "compression-heat-output"]
                    - costs.at["biomass CHP capture", "heat-input"]
                ),
                efficiency3=costs.at["gas", "CO2 intensity"]
                * (1 - costs.at["biomass CHP capture", "capture_rate"]),
                efficiency4=costs.at["gas", "CO2 intensity"]
                * costs.at["biomass CHP capture", "capture_rate"],
                lifetime=costs.at["central gas CHP", "lifetime"],
            )

        if options["chp"] and options["micro_chp"] and name != "urban central":

            n.madd(
                "Link",
                nodes[name] + f" {name} micro gas CHP",
                p_nom_extendable=True,
                bus0="EU gas",
                bus1=nodes[name],
                bus2=nodes[name] + f" {name} heat",
                bus3="co2 atmosphere",
                carrier=name + " micro gas CHP",
                efficiency=costs.at["micro CHP", "efficiency"],
                efficiency2=costs.at["micro CHP", "efficiency-heat"],
                efficiency3=costs.at["gas", "CO2 intensity"],
                capital_cost=costs.at["micro CHP", "fixed"],
                lifetime=costs.at["micro CHP", "lifetime"],
            )
    ###############################
    ###############################
    #######Heat Retrofitting#######
    ###############################

    # if options['retrofitting']['retro_endogen']:

    #     print("adding retrofitting endogenously")

    #     # resample heat demand temporal 'heat_demand_r' depending on in config
    #     # specified temporal resolution, to not overestimate retrofitting
    #     hours = list(filter(re.compile(r'^\d+h$', re.IGNORECASE).search, opts))
    #     if len(hours)==0:
    #         hours = [n.snapshots[1] - n.snapshots[0]]
    #     heat_demand_r =  heat_demand.resample(hours[0]).mean()

    #     # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
    #     # demand 'dE' [per unit of original heat demand] for each country and
    #     # different retrofitting strengths [additional insulation thickness in m]
    #     retro_data = pd.read_csv(snakemake.input.retro_cost_energy,
    #                              index_col=[0, 1], skipinitialspace=True,
    #                              header=[0, 1])
    #     # heated floor area [10^6 * m^2] per country
    #     floor_area = pd.read_csv(snakemake.input.floor_area, index_col=[0, 1])

    #     n.add("Carrier", "retrofitting")

    #     # share of space heat demand 'w_space' of total heat demand
    #     w_space = {}
    #     for sector in sectors:
    #         w_space[sector] = heat_demand_r[sector + " space"] / \
    #             (heat_demand_r[sector + " space"] + heat_demand_r[sector + " water"])
    #     w_space["tot"] = ((heat_demand_r["services space"] +
    #                        heat_demand_r["residential space"]) /
    #                        heat_demand_r.groupby(level=[1], axis=1).sum())

    #     for name in n.loads[n.loads.carrier.isin([x + " heat" for x in heat_systems])].index:

    #         node = n.buses.loc[name, "location"]
    #         ct = pop_layout.loc[node, "ct"]

    #         # weighting 'f' depending on the size of the population at the node
    #         f = urban_fraction[node] if "urban" in name else (1-urban_fraction[node])
    #         if f == 0:
    #             continue
    #         # get sector name ("residential"/"services"/or both "tot" for urban central)
    #         sec = [x if x in name else "tot" for x in sectors][0]

    #         # get floor aread at node and region (urban/rural) in m^2
    #         floor_area_node = ((pop_layout.loc[node].fraction
    #                               * floor_area.loc[ct, "value"] * 10**6).loc[sec] * f)
    #         # total heat demand at node [MWh]
    #         demand = (n.loads_t.p_set[name].resample(hours[0])
    #                   .mean())

    #         # space heat demand at node [MWh]
    #         space_heat_demand = demand * w_space[sec][node]
    #         # normed time profile of space heat demand 'space_pu' (values between 0-1),
    #         # p_max_pu/p_min_pu of retrofitting generators
    #         space_pu = (space_heat_demand / space_heat_demand.max()).to_frame(name=node)

    #         # minimum heat demand 'dE' after retrofitting in units of original heat demand (values between 0-1)
    #         dE = retro_data.loc[(ct, sec), ("dE")]
    #         # get addtional energy savings 'dE_diff' between the different retrofitting strengths/generators at one node
    #         dE_diff = abs(dE.diff()).fillna(1-dE.iloc[0])
    #         # convert costs Euro/m^2 -> Euro/MWh
    #         capital_cost =  retro_data.loc[(ct, sec), ("cost")] * floor_area_node / \
    #                         ((1 - dE) * space_heat_demand.max())
    #         # number of possible retrofitting measures 'strengths' (set in list at config.yaml 'l_strength')
    #         # given in additional insulation thickness [m]
    #         # for each measure, a retrofitting generator is added at the node
    #         strengths = retro_data.columns.levels[1]

    #         # check that ambitious retrofitting has higher costs per MWh than moderate retrofitting
    #         if (capital_cost.diff() < 0).sum():
    #             print(f"Warning: costs are not linear for {ct} {sec}")
    #             s = capital_cost[(capital_cost.diff() < 0)].index
    #             strengths = strengths.drop(s)

    #         # reindex normed time profile of space heat demand back to hourly resolution
    #         space_pu = space_pu.reindex(index=heat_demand.index).fillna(method="ffill")

    #         # add for each retrofitting strength a generator with heat generation profile following the profile of the heat demand
    #         for strength in strengths:
    #             n.madd('Generator',
    #                 [node],
    #                 suffix=' retrofitting ' + strength + " " + name[6::],
    #                 bus=name,
    #                 carrier="retrofitting",
    #                 p_nom_extendable=True,
    #                 p_nom_max=dE_diff[strength] * space_heat_demand.max(), # maximum energy savings for this renovation strength
    #                 p_max_pu=space_pu,
    #                 p_min_pu=space_pu,
    #                 country=ct,
    #                 capital_cost=capital_cost[strength] * options['retrofitting']['cost_factor']
    #             )


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake(
            "prepare_sector_network", simpl="", clusters="4", planning_horizons="2030"
        )
    # TODO add mock_snakemake func

    # TODO fetch from config

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    nodes = n.buses.index

    Nyears = n.snapshot_weightings.generators.sum() / 8760

    # TODO fetch investment year from config

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"]["USD2013_to_EUR2013"],
        snakemake.config["costs"]["discountrate"],
        Nyears,
        snakemake.config["costs"]["lifetime"],
    )
    # TODO logging

    options = snakemake.config["sector"]

    # Get the data required for land transport
    nodal_energy_totals = pd.read_csv(
        snakemake.input.nodal_energy_totals, index_col=0
    )  # where is nodal_energy_totals_4
    transport = pd.read_csv(snakemake.input.transport, index_col=0, parse_dates=True)
    avail_profile = pd.read_csv(
        snakemake.input.avail_profile, index_col=0, parse_dates=True
    )
    dsm_profile = pd.read_csv(
        snakemake.input.dsm_profile, index_col=0, parse_dates=True
    )
    nodal_transport_data = pd.read_csv(
        snakemake.input.nodal_transport_data, index_col=0
    )

    heat_demand = pd.read_csv(snakemake.input.heat_demand, index_col=0, header=[0, 1])
    gshp_cop = pd.read_csv(snakemake.input.gshp_cop, index_col=0)
    ashp_cop = pd.read_csv(snakemake.input.ashp_cop, index_col=0)
    solar_thermal = pd.read_csv(snakemake.input.solar_thermal, index_col=0)

    district_heat_share = pd.read_csv(
        snakemake.input.district_heat_share, index_col=0
    ).iloc[:, 0]
    add_co2(n, costs)  # TODO add costs

    # Add_generation() currently adds gas carrier/bus, as defined in config "conventional_generation"
    add_generation(n, costs)

    # Add_oil() adds oil carrier/bus.
    # TODO This might be transferred to add_generation, but before apply remove_elec_base_techs(n) from PyPSA-Eur-Sec
    add_oil(n, costs)

    add_hydrogen(n, costs)  # TODO add costs

    add_storage(n, costs)

    H2_liquid_fossil_conversions(n, costs)

    h2_hc_conversions(n, costs)

    add_industry(n, costs)

    add_shipping(n, costs)

    # Add_aviation runs with dummy data
    add_aviation(n, costs)

    # prepare_transport_data(n)

    add_land_transport(n, costs)
    # add_heat(n, costs)
    n.export_to_netcdf(snakemake.output[0])

    # n.lopf()
    # TODO define spatial (for biomass and co2)

    # TODO changes in case of myopic oversight

    # TODO add co2 tracking function

    # TODO add generation

    # TODO add storage  HERE THE H2 CARRIER IS ADDED IN PYPSA-EUR-SEC

    # TODO add options as in PyPSA-EUR-SEC
