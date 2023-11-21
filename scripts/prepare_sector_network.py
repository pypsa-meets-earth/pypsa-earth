# -*- coding: utf-8 -*-
import os
import re
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
    two_2_three_digits_country,
)
from prepare_transport_data import prepare_transport_data

spatial = SimpleNamespace()


def add_carrier_buses(n, carrier, nodes=None):
    """
    Add buses to connect e.g. coal, nuclear and oil plants
    """

    if nodes is None:
        nodes = vars(spatial)[carrier].nodes
    location = vars(spatial)[carrier].locations

    # skip if carrier already exists
    if carrier in n.carriers.index:
        return

    if not isinstance(nodes, pd.Index):
        nodes = pd.Index(nodes)

    n.add("Carrier", carrier, co2_emissions=costs.at[carrier, "CO2 intensity"])

    n.madd("Bus", nodes, location=location, carrier=carrier)

    # capital cost could be corrected to e.g. 0.2 EUR/kWh * annuity and O&M
    n.madd(
        "Store",
        nodes + " Store",
        bus=nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier=carrier,
    )

    n.madd(
        "Generator",
        nodes,
        bus=nodes,
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

    for generator, carrier in conventionals.items():
        add_carrier_buses(n, carrier)
        carrier_nodes = vars(spatial)[carrier].nodes
        n.madd(
            "Link",
            nodes + " " + generator,
            bus0=carrier_nodes,
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
    # spatial.gas = SimpleNamespace()

    spatial.oil = SimpleNamespace()

    if options["oil_network"]:
        spatial.oil.nodes = nodes + " oil"
        spatial.oil.locations = nodes
    else:
        spatial.oil.nodes = ["Africa oil"]
        spatial.oil.locations = ["Africa"]

    if "oil" not in n.carriers.index:
        n.add("Carrier", "oil")

    # Set the "co2_emissions" of the carrier "oil" to 0, because the emissions of oil usage taken from the spatial.oil.nodes are accounted seperately (directly linked to the co2 atmosphere bus). Setting the carrier to 0 here avoids double counting. Be aware to link oil emissions to the co2 atmosphere bus.
    n.carriers.loc["oil", "co2_emissions"] = 0
    print("co2_emissions of oil set to 0 for testing")  # TODO add logger.info

    n.madd(
        "Bus",
        spatial.oil.nodes,
        location=spatial.oil.locations,
        carrier="oil",
    )

    # if "Africa oil" not in n.buses.index:

    #     n.add("Bus", "Africa oil", location="Africa", carrier="oil")

    # if "Africa oil Store" not in n.stores.index:

    e_initial = (snakemake.config["fossil_reserves"]).get("oil", 0) * 1e6
    # could correct to e.g. 0.001 EUR/kWh * annuity and O&M
    n.madd(
        "Store",
        [oil_bus + " Store" for oil_bus in spatial.oil.nodes],
        bus=spatial.oil.nodes,
        e_nom_extendable=True,
        e_cyclic=False,
        carrier="oil",
        e_initial=e_initial,
        marginal_cost=costs.at["oil", "fuel"],
    )

    n.madd(
        "Generator",
        spatial.oil.nodes,
        bus=spatial.oil.nodes,
        p_nom_extendable=True,
        carrier="oil",
        marginal_cost=costs.at["oil", "fuel"],
    )


def add_gas(n, costs):
    spatial.gas = SimpleNamespace()

    if options["gas_network"]:
        spatial.gas.nodes = nodes + " gas"
        spatial.gas.locations = nodes
        spatial.gas.biogas = nodes + " biogas"
        spatial.gas.industry = nodes + " gas for industry"
        if snakemake.config["sector"]["cc"]:
            spatial.gas.industry_cc = nodes + " gas for industry CC"
        spatial.gas.biogas_to_gas = nodes + " biogas to gas"
    else:
        spatial.gas.nodes = ["Africa gas"]
        spatial.gas.locations = ["Africa"]
        spatial.gas.biogas = ["Africa biogas"]
        spatial.gas.industry = ["gas for industry"]
        if snakemake.config["sector"]["cc"]:
            spatial.gas.industry_cc = ["gas for industry CC"]
        spatial.gas.biogas_to_gas = ["Africa biogas to gas"]

    spatial.gas.df = pd.DataFrame(vars(spatial.gas), index=nodes)

    gas_nodes = vars(spatial)["gas"].nodes

    add_carrier_buses(n, "gas", gas_nodes)


def H2_liquid_fossil_conversions(n, costs):
    """
    Function to add conversions between H2 and liquid fossil
    Carrier and bus is added in add_oil, which later on might be switched to add_generation
    """

    n.madd(
        "Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1=spatial.oil.nodes,
        bus2=spatial.co2.nodes,
        carrier="Fischer-Tropsch",
        efficiency=costs.at["Fischer-Tropsch", "efficiency"],
        capital_cost=costs.at["Fischer-Tropsch", "fixed"]
        * costs.at[
            "Fischer-Tropsch", "efficiency"
        ],  # Use efficiency to convert from EUR/MW_FT/a to EUR/MW_H2/a
        efficiency2=-costs.at["oil", "CO2 intensity"]
        * costs.at["Fischer-Tropsch", "efficiency"],
        p_nom_extendable=True,
        p_min_pu=options.get("min_part_load_fischer_tropsch", 0),
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

    if snakemake.config["hydrogen_underground_storage"]:
        if snakemake.config["custom_data"]["h2_underground"]:
            custom_cavern = pd.read_csv("resources/custom_data/h2_underground.csv")
            # countries = n.buses.country.unique().to_list()
            countries = snakemake.config["countries"]
            custom_cavern = custom_cavern[custom_cavern.country.isin(countries)]
            cavern_nodes = n.buses[n.buses.country.isin(custom_cavern.country.index)]
            h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

            h2_pot = custom_cavern.loc[cavern_nodes.country]
            h2_pot.index = cavern_nodes.index

            n.madd(
                "Store",
                cavern_nodes.index + " H2 Store",
                bus=cavern_nodes.index + " H2",
                e_nom_extendable=False,
                e_nom_max=h2_pot.values,
                e_cyclic=True,
                carrier="H2 Store",
                capital_cost=h2_capital_cost,
            )

        else:
            h2_salt_cavern_potential = pd.read_csv(
                snakemake.input.h2_cavern, index_col=0
            ).squeeze()
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
    h2_capital_cost = costs.at[
        "hydrogen storage tank type 1 including compressor", "fixed"
    ]
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
    if snakemake.config["custom_data"]["gas_grid"]:
        h2_links = pd.read_csv(snakemake.input.pipelines)

        # Order buses to detect equal pairs for bidirectional pipelines
        buses_ordered = h2_links.apply(lambda p: sorted([p.bus0, p.bus1]), axis=1)

        # Appending string for carrier specification '_AC'
        h2_links["bus0"] = buses_ordered.str[0] + "_AC"
        h2_links["bus1"] = buses_ordered.str[1] + "_AC"

        # Conversion of GADM id to from 3 to 2-digit
        h2_links["bus0"] = (
            h2_links["bus0"]
            .str.split(".")
            .apply(lambda id: three_2_two_digits_country(id[0]) + "." + id[1])
        )
        h2_links["bus1"] = (
            h2_links["bus1"]
            .str.split(".")
            .apply(lambda id: three_2_two_digits_country(id[0]) + "." + id[1])
        )

        # Create index column
        h2_links["buses_idx"] = (
            "H2 pipeline " + h2_links["bus0"] + " -> " + h2_links["bus1"]
        )

        # Aggregate pipelines applying mean on length and sum on capacities
        h2_links = h2_links.groupby("buses_idx").agg(
            {"bus0": "first", "bus1": "first", "length": "mean", "capacity": "sum"}
        )
    else:
        attrs = ["bus0", "bus1", "length"]
        h2_links = pd.DataFrame(columns=attrs)

        candidates = pd.concat(
            {
                "lines": n.lines[attrs],
                "links": n.links.loc[n.links.carrier == "DC", attrs],
            }
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
    if snakemake.config["H2_network"]:
        if snakemake.config["H2_repurposed_network"]:
            n.madd(
                "Link",
                h2_links.index + " repurposed",
                bus0=h2_links.bus0.values + " H2",
                bus1=h2_links.bus1.values + " H2",
                p_min_pu=-1,
                p_nom_extendable=True,
                p_nom_max=h2_links.capacity.values
                * 0.8,  # https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf
                length=h2_links.length.values,
                capital_cost=costs.at["H2 (g) pipeline repurposed", "fixed"]
                * h2_links.length.values,
                carrier="H2 pipeline repurposed",
                lifetime=costs.at["H2 (g) pipeline repurposed", "lifetime"],
            )
            n.madd(
                "Link",
                h2_links.index,
                bus0=h2_links.bus0.values + " H2",
                bus1=h2_links.bus1.values + " H2",
                p_min_pu=-1,
                p_nom_extendable=True,
                length=h2_links.length.values,
                capital_cost=costs.at["H2 (g) pipeline", "fixed"]
                * h2_links.length.values,
                carrier="H2 pipeline",
                lifetime=costs.at["H2 (g) pipeline", "lifetime"],
            )
        else:
            n.madd(
                "Link",
                h2_links.index,
                bus0=h2_links.bus0.values + " H2",
                bus1=h2_links.bus1.values + " H2",
                p_min_pu=-1,
                p_nom_extendable=True,
                length=h2_links.length.values,
                capital_cost=costs.at["H2 (g) pipeline", "fixed"]
                * h2_links.length.values,
                carrier="H2 pipeline",
                lifetime=costs.at["H2 (g) pipeline", "lifetime"],
            )


def define_spatial(nodes):
    """
    Namespace for spatial

    Parameters
    ----------
    nodes : list-like
    """

    global spatial
    global options

    spatial.nodes = nodes

    # biomass

    spatial.biomass = SimpleNamespace()

    if options["biomass_transport"]:
        spatial.biomass.nodes = nodes + " solid biomass"
        spatial.biomass.locations = nodes
        spatial.biomass.industry = nodes + " solid biomass for industry"
        spatial.biomass.industry_cc = nodes + " solid biomass for industry CC"
    else:
        spatial.biomass.nodes = ["Africa solid biomass"]
        spatial.biomass.locations = ["Africa"]
        spatial.biomass.industry = ["solid biomass for industry"]
        spatial.biomass.industry_cc = ["solid biomass for industry CC"]

    spatial.biomass.df = pd.DataFrame(vars(spatial.biomass), index=nodes)

    # co2

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


def add_biomass(n, costs):
    print("adding biomass")

    # TODO get biomass potentials dataset and enable spatially resolved potentials
    biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)

    # if options["biomass_transport"]:
    #     biomass_potentials_spatial = biomass_potentials.rename(
    #         index=lambda x: x + " solid biomass"
    #     )
    # else:
    #     biomass_potentials_spatial = biomass_potentials.sum()
    biomass_potentials_spatial = biomass_potentials.sum() / len(biomass_potentials)
    print(
        "Solid Biomass potentials are not spatially resolved, the dummy value is the average solid Biomass of European countries"
    )

    n.add("Carrier", "biogas")
    n.add("Carrier", "solid biomass")

    n.madd(
        "Bus", spatial.gas.biogas, location=spatial.biomass.locations, carrier="biogas"
    )

    n.madd(
        "Bus",
        spatial.biomass.nodes,
        location=spatial.biomass.locations,
        carrier="solid biomass",
    )

    n.madd(
        "Store",
        spatial.gas.biogas,
        bus=spatial.gas.biogas,
        carrier="biogas",
        e_nom=biomass_potentials["biogas"].sum(),
        marginal_cost=costs.at["biogas", "fuel"],
        e_initial=biomass_potentials["biogas"].sum(),
    )

    n.madd(
        "Store",
        spatial.biomass.nodes,
        bus=spatial.biomass.nodes,
        carrier="solid biomass",
        e_nom=biomass_potentials_spatial["solid biomass"],
        marginal_cost=costs.at["solid biomass", "fuel"],
        e_initial=biomass_potentials_spatial["solid biomass"],
    )

    biomass_gen = "biomass EOP"
    n.madd(
        "Link",
        nodes + " biomass EOP",
        bus0=spatial.biomass.nodes,
        bus1=nodes,
        # bus2="co2 atmosphere",
        marginal_cost=costs.at[biomass_gen, "efficiency"]
        * costs.at[biomass_gen, "VOM"],  # NB: VOM is per MWel
        # NB: fixed cost is per MWel
        capital_cost=costs.at[biomass_gen, "efficiency"]
        * costs.at[biomass_gen, "fixed"],
        p_nom_extendable=True,
        carrier=biomass_gen,
        efficiency=costs.at[biomass_gen, "efficiency"],
        # efficiency2=costs.at["solid biomass", "CO2 intensity"],
        lifetime=costs.at[biomass_gen, "lifetime"],
    )
    n.madd(
        "Link",
        spatial.gas.biogas_to_gas,
        bus0=spatial.gas.biogas,
        bus1=spatial.gas.nodes,
        bus2="co2 atmosphere",
        carrier="biogas to gas",
        capital_cost=costs.loc["biogas upgrading", "fixed"],
        marginal_cost=costs.loc["biogas upgrading", "VOM"],
        efficiency2=-costs.at["gas", "CO2 intensity"],
        p_nom_extendable=True,
    )

    if options["biomass_transport"]:
        # TODO add biomass transport costs
        transport_costs = pd.read_csv(
            snakemake.input.biomass_transport_costs,
            index_col=0,
            keep_default_na=False,
        ).squeeze()

        # add biomass transport
        biomass_transport = create_network_topology(
            n, "biomass transport ", bidirectional=False
        )

        # costs
        bus0_costs = biomass_transport.bus0.apply(lambda x: transport_costs[x[:2]])
        bus1_costs = biomass_transport.bus1.apply(lambda x: transport_costs[x[:2]])
        biomass_transport["costs"] = pd.concat([bus0_costs, bus1_costs], axis=1).mean(
            axis=1
        )

        n.madd(
            "Link",
            biomass_transport.index,
            bus0=biomass_transport.bus0 + " solid biomass",
            bus1=biomass_transport.bus1 + " solid biomass",
            p_nom_extendable=True,
            length=biomass_transport.length.values,
            marginal_cost=biomass_transport.costs * biomass_transport.length.values,
            capital_cost=1,
            carrier="solid biomass transport",
        )

    # n.madd(
    #         "Link",
    #         urban_central + " urban central solid biomass CHP",
    #         bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
    #         bus1=urban_central,
    #         bus2=urban_central + " urban central heat",
    #         carrier="urban central solid biomass CHP",
    #         p_nom_extendable=True,
    #         capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
    #         marginal_cost=costs.at[key, "VOM"],
    #         efficiency=costs.at[key, "efficiency"],
    #         efficiency2=costs.at[key, "efficiency-heat"],
    #         lifetime=costs.at[key, "lifetime"],
    #     )

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty and options["chp"]:
        urban_central = urban_central.str[: -len(" urban central heat")]

        key = "central solid biomass CHP"

        n.madd(
            "Link",
            urban_central + " urban central solid biomass CHP",
            bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
            bus1=urban_central,
            bus2=urban_central + " urban central heat",
            carrier="urban central solid biomass CHP",
            p_nom_extendable=True,
            capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"],
            marginal_cost=costs.at[key, "VOM"],
            efficiency=costs.at[key, "efficiency"],
            efficiency2=costs.at[key, "efficiency-heat"],
            lifetime=costs.at[key, "lifetime"],
        )

        if snakemake.config["sector"]["cc"]:
            n.madd(
                "Link",
                urban_central + " urban central solid biomass CHP CC",
                bus0=spatial.biomass.df.loc[urban_central, "nodes"].values,
                bus1=urban_central,
                bus2=urban_central + " urban central heat",
                bus3="co2 atmosphere",
                bus4=spatial.co2.df.loc[urban_central, "nodes"].values,
                carrier="urban central solid biomass CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at[key, "fixed"] * costs.at[key, "efficiency"]
                + costs.at["biomass CHP capture", "fixed"]
                * costs.at["solid biomass", "CO2 intensity"],
                marginal_cost=costs.at[key, "VOM"],
                efficiency=costs.at[key, "efficiency"]
                - costs.at["solid biomass", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "electricity-input"]
                    + costs.at["biomass CHP capture", "compression-electricity-input"]
                ),
                efficiency2=costs.at[key, "efficiency-heat"]
                + costs.at["solid biomass", "CO2 intensity"]
                * (
                    costs.at["biomass CHP capture", "heat-output"]
                    + costs.at["biomass CHP capture", "compression-heat-output"]
                    - costs.at["biomass CHP capture", "heat-input"]
                ),
                efficiency3=-costs.at["solid biomass", "CO2 intensity"]
                * costs.at["biomass CHP capture", "capture_rate"],
                efficiency4=costs.at["solid biomass", "CO2 intensity"]
                * costs.at["biomass CHP capture", "capture_rate"],
                lifetime=costs.at[key, "lifetime"],
            )


def add_co2(n, costs):
    "add carbon carrier, it's networks and storage units"

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
        # x=spatial.co2.x[0],
        # y=spatial.co2.y[0],
    )
    """
    co2_stored_x = n.buses.filter(like="co2 stored", axis=0).loc[:, "x"]
    co2_stored_y = n.buses.loc[n.buses[n.buses.carrier == "co2 stored"].location].y

    n.buses[n.buses.carrier == "co2 stored"].x = co2_stored_x.values
    n.buses[n.buses.carrier == "co2 stored"].y = co2_stored_y.values
    """

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
    all_aviation = ["total international aviation", "total domestic aviation"]

    aviation_demand = (
        energy_totals.loc[countries, all_aviation].sum(axis=1).sum()  # * 1e6 / 8760
    )

    airports = pd.read_csv(snakemake.input.airports, keep_default_na=False)
    airports = airports[airports.country.isin(countries)]

    gadm_level = options["gadm_level"]

    airports["gadm_{}".format(gadm_level)] = airports[["x", "y", "country"]].apply(
        lambda airport: locate_bus(
            airport[["x", "y"]],
            airport["country"],
            gadm_level,
            snakemake.input.shapes_path,
            snakemake.config["clustering_options"]["alternative_clustering"],
        ),
        axis=1,
    )
    # To change 3 country code to 2
    # airports["gadm_{}".format(gadm_level)] = airports["gadm_{}".format(gadm_level)].apply(
    # lambda cocode: three_2_two_digits_country(cocode[:3]) + " " + cocode[4:-2])

    airports = airports.set_index("gadm_{}".format(gadm_level))

    ind = pd.DataFrame(n.buses.index[n.buses.carrier == "AC"])

    ind = ind.set_index(n.buses.index[n.buses.carrier == "AC"])
    airports["p_set"] = airports["fraction"].apply(
        lambda frac: frac * aviation_demand * 1e6 / 8760
    )

    airports = pd.concat([airports, ind])

    # airports = airports.fillna(0)

    airports = airports.groupby(airports.index).sum()
    n.madd(
        "Load",
        nodes,
        suffix=" kerosene for aviation",
        bus=spatial.oil.nodes,
        carrier="kerosene for aviation",
        p_set=airports["p_set"],
    )

    co2 = airports["p_set"].sum() * costs.at["oil", "CO2 intensity"]

    n.add(
        "Load",
        "aviation oil emissions",
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
            bus1=spatial.gas.nodes,
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
            bus1=spatial.gas.nodes,
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
            bus0=spatial.gas.nodes,
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
            bus0=spatial.gas.nodes,
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
    ports = pd.read_csv(
        snakemake.input.ports, index_col=None, keep_default_na=False
    ).squeeze()
    ports = ports[ports.country.isin(countries)]

    gadm_level = options["gadm_level"]

    all_navigation = ["total international navigation", "total domestic navigation"]

    navigation_demand = (
        energy_totals.loc[countries, all_navigation].sum(axis=1).sum()  # * 1e6 / 8760
    )

    efficiency = (
        options["shipping_average_efficiency"] / costs.at["fuel cell", "efficiency"]
    )

    # check whether item depends on investment year
    shipping_hydrogen_share = get(
        options["shipping_hydrogen_share"], demand_sc + "_" + str(investment_year)
    )

    ports["gadm_{}".format(gadm_level)] = ports[["x", "y", "country"]].apply(
        lambda port: locate_bus(
            port[["x", "y"]],
            port["country"],
            gadm_level,
            snakemake.input["shapes_path"],
            snakemake.config["clustering_options"]["alternative_clustering"],
        ),
        axis=1,
    )

    ports = ports.set_index("gadm_{}".format(gadm_level))

    ind = pd.DataFrame(n.buses.index[n.buses.carrier == "AC"])
    ind = ind.set_index(n.buses.index[n.buses.carrier == "AC"])

    ports["p_set"] = ports["fraction"].apply(
        lambda frac: shipping_hydrogen_share
        * frac
        * navigation_demand
        * efficiency
        * 1e6
        / 8760
        # TODO double check the use of efficiency
    )  # TODO use real data here

    ports = pd.concat([ports, ind]).drop("Bus", axis=1)

    # ports = ports.fillna(0.0)
    ports = ports.groupby(ports.index).sum()

    if options["shipping_hydrogen_liquefaction"]:
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
            lambda frac: shipping_oil_share * frac * navigation_demand * 1e6 / 8760
        )

        n.madd(
            "Load",
            nodes,
            suffix=" shipping oil",
            bus=spatial.oil.nodes,
            carrier="shipping oil",
            p_set=ports["p_set"],
        )

        co2 = (
            shipping_oil_share * ports["p_set"].sum() * costs.at["oil", "CO2 intensity"]
        )

        n.add(
            "Load",
            "shipping oil emissions",
            bus="co2 atmosphere",
            carrier="shipping oil emissions",
            p_set=-co2,
        )

    if "oil" not in n.buses.carrier.unique():
        n.madd("Bus", spatial.oil.nodes, location=spatial.oil.locations, carrier="oil")
    if "oil" not in n.stores.carrier.unique():
        # could correct to e.g. 0.001 EUR/kWh * annuity and O&M
        n.madd(
            "Store",
            [oil_bus + " Store" for oil_bus in spatial.oil.nodes],
            bus=spatial.oil.nodes,
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="oil",
        )

    if "oil" not in n.generators.carrier.unique():
        n.madd(
            "Generator",
            spatial.oil.nodes,
            bus=spatial.oil.nodes,
            p_nom_extendable=True,
            carrier="oil",
            marginal_cost=costs.at["oil", "fuel"],
        )


def add_industry(n, costs):
    print("adding industrial demand")
    # 1e6 to convert TWh to MWh

    # industrial_demand.reset_index(inplace=True)

    # Add carrier Biomass

    n.madd(
        "Bus",
        spatial.biomass.industry,
        location=spatial.biomass.locations,
        carrier="solid biomass for industry",
    )

    if options["biomass_transport"]:
        p_set = (
            industrial_demand.loc[spatial.biomass.locations, "solid biomass"].rename(
                index=lambda x: x + " solid biomass for industry"
            )
            / 8760
        )
    else:
        p_set = industrial_demand["solid biomass"].sum() / 8760

    n.madd(
        "Load",
        spatial.biomass.industry,
        bus=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_set=p_set,
    )

    n.madd(
        "Link",
        spatial.biomass.industry,
        bus0=spatial.biomass.nodes,
        bus1=spatial.biomass.industry,
        carrier="solid biomass for industry",
        p_nom_extendable=True,
        efficiency=1.0,
    )
    if snakemake.config["sector"]["cc"]:
        n.madd(
            "Link",
            spatial.biomass.industry_cc,
            bus0=spatial.biomass.nodes,
            bus1=spatial.biomass.industry,
            bus2="co2 atmosphere",
            bus3=spatial.co2.nodes,
            carrier="solid biomass for industry CC",
            p_nom_extendable=True,
            capital_cost=costs.at["cement capture", "fixed"]
            * costs.at["solid biomass", "CO2 intensity"],
            efficiency=0.9,  # TODO: make config option
            efficiency2=-costs.at["solid biomass", "CO2 intensity"]
            * costs.at["cement capture", "capture_rate"],
            efficiency3=costs.at["solid biomass", "CO2 intensity"]
            * costs.at["cement capture", "capture_rate"],
            lifetime=costs.at["cement capture", "lifetime"],
        )

    # CARRIER = FOSSIL GAS

    nodes = (
        pop_layout.index
    )  # TODO where to change country code? 2 letter country codes.

    # industrial_demand['TWh/a (MtCO2/a)'] = industrial_demand['TWh/a (MtCO2/a)'].apply(
    #     lambda cocode: two_2_three_digits_country(cocode[:2]) + "." + cocode[3:])

    # industrial_demand.set_index("TWh/a (MtCO2/a)", inplace=True)

    # n.add("Bus", "gas for industry", location="Africa", carrier="gas for industry")
    n.madd(
        "Bus",
        spatial.gas.industry,
        location=spatial.gas.locations,
        carrier="gas for industry",
    )

    gas_demand = industrial_demand.loc[nodes, "gas"] / 8760.0

    if options["gas_network"]:
        spatial_gas_demand = gas_demand.rename(index=lambda x: x + " gas for industry")
    else:
        spatial_gas_demand = gas_demand.sum()

    n.madd(
        "Load",
        spatial.gas.industry,
        bus=spatial.gas.industry,
        carrier="gas for industry",
        p_set=spatial_gas_demand,
    )

    n.madd(
        "Link",
        spatial.gas.industry,
        # bus0="Africa gas",
        bus0=spatial.gas.nodes,
        # bus1="gas for industry",
        bus1=spatial.gas.industry,
        bus2="co2 atmosphere",
        carrier="gas for industry",
        p_nom_extendable=True,
        efficiency=1.0,
        efficiency2=costs.at["gas", "CO2 intensity"],
    )
    if snakemake.config["sector"]["cc"]:
        n.madd(
            "Link",
            spatial.gas.industry_cc,
            # suffix=" gas for industry CC",
            # bus0="Africa gas",
            bus0=spatial.gas.nodes,
            bus1=spatial.gas.industry,
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
        p_set=industrial_demand["hydrogen"].apply(lambda frac: frac / 8760),
    )

    # CARRIER = LIQUID HYDROCARBONS
    n.madd(
        "Load",
        nodes,
        suffix=" naphtha for industry",
        bus=spatial.oil.nodes,
        carrier="naphtha for industry",
        p_set=industrial_demand["oil"].apply(lambda frac: frac / 8760),
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
        # - industrial_demand["process emission from feedstock"].sum()
        # / 8760
    )

    n.add(
        "Load",
        "industry oil emissions",
        bus="co2 atmosphere",
        carrier="industry oil emissions",
        p_set=-co2,
    )

    co2 = (
        industrial_demand["coal"].sum()
        * costs.at["coal", "CO2 intensity"]
        # - industrial_demand["process emission from feedstock"].sum()
        / 8760
    )

    n.add(
        "Load",
        "industry coal emissions",
        bus="co2 atmosphere",
        carrier="industry coal emissions",
        p_set=-co2,
    )

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
            (n.loads.index.str[:2] == ct) & (n.loads.carrier == "AC")
        ]
        if n.loads_t.p_set[loads_i].empty:
            continue

    # if not snakemake.config["custom_data"]["elec_demand"]:
    #     # if electricity demand is provided by pypsa-earth, the electricty used
    #     # in industry is included, and need to be removed from the default elec
    #     # demand here, and added as "industry electricity"
    #     factor = (
    #         1
    #         - industrial_demand.loc[loads_i, "current electricity"].sum()
    #         / n.loads_t.p_set[loads_i].sum().sum()
    #     )
    #     n.loads_t.p_set[loads_i] *= factor
    #     industrial_elec = industrial_demand["current electricity"].apply(
    #         lambda frac: frac / 8760
    #     )

    # else:
    industrial_elec = industrial_demand["electricity"].apply(lambda frac: frac / 8760)

    n.madd(
        "Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=industrial_elec,
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
            #    industrial_demand["process emission from feedstock"]+
            industrial_demand["process emissions"]
        )
        / 8760,
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
    if snakemake.config["sector"]["cc"]:
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

    fuel_cell_share = get(
        options["land_transport_fuel_cell_share"],
        demand_sc + "_" + str(investment_year),
    )
    electric_share = get(
        options["land_transport_electric_share"], demand_sc + "_" + str(investment_year)
    )
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
        if "oil" not in n.buses.carrier.unique():
            n.madd(
                "Bus", spatial.oil.nodes, location=spatial.oil.locations, carrier="oil"
            )
        ice_efficiency = options["transport_internal_combustion_efficiency"]

        n.madd(
            "Load",
            nodes,
            suffix=" land transport oil",
            bus=spatial.oil.nodes,
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
        district_heat_share["district heat share"]
        * pop_layout["urban_ct_fraction"]
        / pop_layout["fraction"]
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
        "resulting in a district heating share of",  # "\n{dist_fraction_node}", #TODO fix district heat share
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
                bus0=spatial.gas.nodes,
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
                bus0=spatial.gas.nodes,
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
            if snakemake.config["sector"]["cc"]:
                n.madd(
                    "Link",
                    nodes[name] + " urban central gas CHP CC",
                    # bus0="Africa gas",
                    bus0=spatial.gas.nodes,
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
                        + costs.at[
                            "biomass CHP capture", "compression-electricity-input"
                        ]
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
                # bus0="Africa gas",
                bus0=spatial.gas.nodes,
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


def average_every_nhours(n, offset):
    # logger.info(f'Resampling the network to {offset}')
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                if c.list_name == "stores" and k == "e_max_pu":
                    pnl[k] = df.resample(offset).min()
                elif c.list_name == "stores" and k == "e_min_pu":
                    pnl[k] = df.resample(offset).max()
                else:
                    pnl[k] = df.resample(offset).mean()

    return m


def add_dac(n, costs):
    heat_carriers = ["urban central heat", "services urban decentral heat"]
    heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]
    locations = n.buses.location[heat_buses]

    efficiency2 = -(
        costs.at["direct air capture", "electricity-input"]
        + costs.at["direct air capture", "compression-electricity-input"]
    )
    efficiency3 = -(
        costs.at["direct air capture", "heat-input"]
        - costs.at["direct air capture", "compression-heat-output"]
    )

    n.madd(
        "Link",
        heat_buses.str.replace(" heat", " DAC"),
        bus0="co2 atmosphere",
        bus1=spatial.co2.df.loc[locations, "nodes"].values,
        bus2=locations.values,
        bus3=heat_buses,
        carrier="DAC",
        capital_cost=costs.at["direct air capture", "fixed"],
        efficiency=1.0,
        efficiency2=efficiency2,
        efficiency3=efficiency3,
        p_nom_extendable=True,
        lifetime=costs.at["direct air capture", "lifetime"],
    )


def add_services(n, costs):
    # TODO make compatible with more counties
    profile_residential = n.loads_t.p_set[nodes] / n.loads_t.p_set[nodes].sum().sum()

    p_set_elec = (
        profile_residential
        * energy_totals.loc[countries, "services electricity"].sum()
        * 1e6
    )

    n.madd(
        "Load",
        nodes,
        suffix=" services electricity",
        bus=nodes,
        carrier="services electricity",
        p_set=p_set_elec,
    )
    p_set_biomass = (
        profile_residential
        * energy_totals.loc[countries, "services biomass"].sum()
        * 1e6
    )

    n.madd(
        "Load",
        nodes,
        suffix=" services biomass",
        bus=spatial.biomass.nodes,
        carrier="services biomass",
        p_set=p_set_biomass,
    )

    # co2 = (
    #     p_set_biomass.sum().sum() * costs.at["solid biomass", "CO2 intensity"]
    # ) / 8760

    # n.add(
    #     "Load",
    #     "services biomass emissions",
    #     bus="co2 atmosphere",
    #     carrier="biomass emissions",
    #     p_set=-co2,
    # )
    p_set_oil = (
        profile_residential * energy_totals.loc[countries, "services oil"].sum() * 1e6
    )

    n.madd(
        "Load",
        nodes,
        suffix=" services oil",
        bus=spatial.oil.nodes,
        carrier="services oil",
        p_set=p_set_oil,
    )

    co2 = (p_set_oil.sum().sum() * costs.at["oil", "CO2 intensity"]) / 8760

    n.add(
        "Load",
        "services oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2,
    )

    p_set_gas = (
        profile_residential * energy_totals.loc[countries, "services gas"].sum() * 1e6
    )

    n.madd(
        "Load",
        nodes,
        suffix=" services gas",
        bus=spatial.gas.nodes,
        carrier="services gas",
        p_set=p_set_gas,
    )

    co2 = (p_set_gas.sum().sum() * costs.at["gas", "CO2 intensity"]) / 8760

    n.add(
        "Load",
        "services gas emissions",
        bus="co2 atmosphere",
        carrier="gas emissions",
        p_set=-co2,
    )


def add_agriculture(n, costs):
    n.madd(
        "Load",
        nodes,
        suffix=" agriculture electricity",
        bus=nodes,
        carrier="agriculture electricity",
        p_set=nodal_energy_totals.loc[nodes, "agriculture electricity"] * 1e6 / 8760,
    )

    n.madd(
        "Load",
        nodes,
        suffix=" agriculture oil",
        bus=spatial.oil.nodes,
        carrier="agriculture oil",
        p_set=nodal_energy_totals.loc[nodes, "agriculture oil"] * 1e6 / 8760,
    )
    co2 = (
        nodal_energy_totals.loc[nodes, "agriculture oil"]
        * 1e6
        * costs.at["oil", "CO2 intensity"]
    )

    n.add(
        "Load",
        "agriculture oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2,
    )


def add_residential(n, costs):
    # need to adapt for many countries #TODO

    # if snakemake.config["custom_data"]["heat_demand"]:
    # heat_demand_index=n.loads_t.p.filter(like='residential').filter(like='heat').dropna(axis=1).index
    # oil_res_index=n.loads_t.p.filter(like='residential').filter(like='oil').dropna(axis=1).index

    heat_ind = (
        n.loads_t.p_set.filter(like="residential")
        .filter(like="heat")
        .dropna(axis=1)
        .columns
    )
    oil_ind = (
        n.loads_t.p_set.filter(like="residential")
        .filter(like="oil")
        .dropna(axis=1)
        .columns
    )
    bio_ind = (
        n.loads_t.p_set.filter(like="residential")
        .filter(like="biomass")
        .dropna(axis=1)
        .columns
    )

    gas_ind = (
        n.loads_t.p_set.filter(like="residential")
        .filter(like="gas")
        .dropna(axis=1)
        .columns
    )

    heat_shape = (
        n.loads_t.p_set.loc[:, heat_ind] / n.loads_t.p_set.loc[:, heat_ind].sum().sum()
    )

    heat_shape = heat_shape.groupby(
        lambda x: next((substring for substring in nodes if substring in x), x), axis=1
    ).sum()

    heat_oil_demand = (
        heat_shape * energy_totals.loc[countries[0], "residential heat oil"] * 1e6
    )

    heat_biomass_demand = (
        heat_shape * energy_totals.loc[countries[0], "residential heat biomass"] * 1e6
    )

    heat_gas_demand = (
        heat_shape * energy_totals.loc[countries[0], "residential heat gas"] * 1e6
    )

    n.loads_t.p_set.loc[:, heat_ind] = (
        heat_shape
        * (
            energy_totals.loc[countries, "total residential space"].sum()
            + energy_totals.loc[countries, "total residential water"].sum()
            - energy_totals.loc[countries[0], "residential heat biomass"]
            - energy_totals.loc[countries[0], "residential heat oil"]
            - energy_totals.loc[countries[0], "residential heat gas"]
        )
        * 1e6
    )

    print("#############################################")
    print("#############################################")
    print("#############################################")
    print(nodes)
    print("#############################################")
    print("#############################################")
    print("#############################################")
    print(heat_shape)
    print("#############################################")
    print("#############################################")
    print("#############################################")

    # TODO make compatible with more counties
    profile_residential = n.loads_t.p_set[nodes] / n.loads_t.p_set[nodes].sum().sum()

    p_set_oil = (
        profile_residential
        * energy_totals.loc[countries, "residential oil"].sum()
        * 1e6
    ) + heat_oil_demand.values

    p_set_biomass = (
        profile_residential
        * energy_totals.loc[countries, "residential biomass"].sum()
        * 1e6
    ) + heat_biomass_demand.values

    p_set_gas = (
        profile_residential
        * energy_totals.loc[countries, "residential gas"].sum()
        * 1e6
    ) + heat_gas_demand.values

    n.madd(
        "Load",
        nodes,
        suffix=" residential oil",
        bus=spatial.oil.nodes,
        carrier="residential oil",
        p_set=p_set_oil,
    )
    co2 = (p_set_oil.sum().sum() * costs.at["oil", "CO2 intensity"]) / 8760

    n.add(
        "Load",
        "residential oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2,
    )
    n.madd(
        "Load",
        nodes,
        suffix=" residential biomass",
        bus=spatial.biomass.nodes,
        carrier="residential biomass",
        p_set=p_set_biomass,
    )

    n.madd(
        "Load",
        nodes,
        suffix=" residential gas",
        bus=spatial.gas.nodes,
        carrier="residential gas",
        p_set=p_set_gas,
    )
    co2 = (p_set_gas.sum().sum() * costs.at["gas", "CO2 intensity"]) / 8760

    n.add(
        "Load",
        "residential gas emissions",
        bus="co2 atmosphere",
        carrier="gas emissions",
        p_set=-co2,
    )

    co2 = (p_set_oil.sum().sum() * costs.at["solid biomass", "CO2 intensity"]) / 8760

    n.add(
        "Load",
        "residential biomass emissions",
        bus="co2 atmosphere",
        carrier="biomass emissions",
        p_set=-co2,
    )

    for country in countries:
        rem_heat_demand = (
            energy_totals.loc[country, "total residential space"]
            + energy_totals.loc[country, "total residential water"]
            - energy_totals.loc[country, "residential heat biomass"]
            - energy_totals.loc[country, "residential heat oil"]
            - energy_totals.loc[country, "residential heat gas"]
        )

        heat_buses = (n.loads_t.p_set.filter(regex="heat").filter(like=country)).columns
        n.loads_t.p_set.loc[:, heat_buses] = (
            (
                n.loads_t.p_set.filter(like=country)[heat_buses]
                / n.loads_t.p_set.filter(like=country)[heat_buses].sum().sum()
            )
            * rem_heat_demand
            * 1e6
        )

        # if snakemake.config["custom_data"]["elec_demand"]:
    for country in countries:
        # indd=n.loads_t.p_set[n.loads_t.p_set.columns.str.contains(country)]

        buses = n.buses[(n.buses.carrier == "AC") & (n.buses.country == country)].index

        n.loads_t.p_set.loc[:, buses] = (
            (
                n.loads_t.p_set.filter(like=country)[buses]
                / n.loads_t.p_set.filter(like=country)[buses].sum().sum()
            )
            * energy_totals.loc[country, "electricity residential"]
            * 1e6
        )


# def add_co2limit(n, Nyears=1.0, limit=0.0):
#     print("Adding CO2 budget limit as per unit of 1990 levels of", limit)

#     countries = n.buses.country.dropna().unique()

#     sectors = emission_sectors_from_opts(opts)

#     # convert Mt to tCO2
#     co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

#     co2_limit = co2_totals.loc[countries, sectors].sum().sum()

#     co2_limit *= limit * Nyears

#     n.add(
#         "GlobalConstraint",
#         "CO2Limit",
#         carrier_attribute="co2_emissions",
#         sense="<=",
#         constant=co2_limit,
#     )


def add_custom_water_cost(n):
    for country in countries:
        water_costs = pd.read_csv(
            "resources/custom_data/{}_water_costs.csv".format(country),
            sep=",",
            index_col=0,
        )
        water_costs = water_costs.filter(like=country, axis=0).loc[nodes]
        electrolysis_links = n.links.filter(like=country, axis=0).filter(
            like="lectrolysis", axis=0
        )

        elec_index = n.links[
            (n.links.carrier == "H2 Electrolysis")
            & (n.links.bus0.str.contains(country))
        ].index
        n.links.loc[elec_index, "marginal_cost"] = water_costs.values
        # n.links.filter(like=country, axis=0).filter(like='lectrolysis', axis=0)["marginal_cost"] = water_costs.values
        # n.links.filter(like=country, axis=0).filter(like='lectrolysis', axis=0).apply(lambda x: water_costs[x.index], axis=0)
        # print(n.links.filter(like=country, axis=0).filter(like='lectrolysis', axis=0).marginal_cost)


def add_rail_transport(n, costs):
    p_set_elec = nodal_energy_totals.loc[nodes, "electricity rail"]
    p_set_oil = (nodal_energy_totals.loc[nodes, "total rail"]) - p_set_elec

    n.madd(
        "Load",
        nodes,
        suffix=" rail transport oil",
        bus=spatial.oil.nodes,
        carrier=" rail transport oil",
        p_set=p_set_oil * 1e6 / 8760,
    )

    n.madd(
        "Load",
        nodes,
        suffix=" rail transport electricity",
        bus=nodes,
        carrier=" rail transport electricity",
        p_set=p_set_elec * 1e6 / 8760,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake(
            "prepare_sector_network",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L",
            planning_horizons="2030",
            sopts="144H",
            discountrate="0.071",
            demand="AB",
        )

    # Load population layout
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    # Load all sector wildcards
    options = snakemake.config["sector"]

    # Load input network
    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    # Fetch the coutry list from the network
    # countries = list(n.buses.country.unique())
    countries = snakemake.config["countries"]
    # Locate all the AC buses
    nodes = n.buses[
        n.buses.carrier == "AC"
    ].index  # TODO if you take nodes from the index of buses of n it's more than pop_layout
    # clustering of regions must be double checked.. refer to regions onshore

    # Set carrier of AC loads
    n.loads.loc[nodes, "carrier"] = "AC"

    Nyears = n.snapshot_weightings.generators.sum() / 8760

    # Fetch wildcards
    investment_year = int(snakemake.wildcards.planning_horizons[-4:])
    demand_sc = snakemake.wildcards.demand  # loading the demand scenrario wildcard
    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    # Prepare the costs dataframe
    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"]["USD2013_to_EUR2013"],
        eval(snakemake.wildcards.discountrate),
        Nyears,
        snakemake.config["costs"]["lifetime"],
    )

    # Define spatial for biomass and co2. They require the same spatial definition
    define_spatial(pop_layout.index)

    # TODO logging

    nodal_energy_totals = pd.read_csv(
        snakemake.input.nodal_energy_totals,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )
    energy_totals = pd.read_csv(
        snakemake.input.energy_totals,
        index_col=0,
        keep_default_na=False,
        na_values=[""],
    )
    # Get the data required for land transport
    # TODO Leon, This contains transport demand, right? if so let's change it to transport_demand?
    transport = pd.read_csv(snakemake.input.transport, index_col=0, parse_dates=True)

    avail_profile = pd.read_csv(
        snakemake.input.avail_profile, index_col=0, parse_dates=True
    )
    dsm_profile = pd.read_csv(
        snakemake.input.dsm_profile, index_col=0, parse_dates=True
    )
    nodal_transport_data = pd.read_csv(  # TODO This only includes no. of cars, change name to something descriptive?
        snakemake.input.nodal_transport_data, index_col=0
    )

    # Load data required for the heat sector
    heat_demand = pd.read_csv(
        snakemake.input.heat_demand, index_col=0, header=[0, 1], parse_dates=True
    ).fillna(0)
    # Ground-sourced heatpump coefficient of performance
    gshp_cop = pd.read_csv(
        snakemake.input.gshp_cop, index_col=0, parse_dates=True
    )  # only needed with heat dep. hp cop allowed from config
    # TODO add option heat_dep_hp_cop to the config

    # Air-sourced heatpump coefficient of performance
    ashp_cop = pd.read_csv(
        snakemake.input.ashp_cop, index_col=0, parse_dates=True
    )  # only needed with heat dep. hp cop allowed from config

    # Solar thermal availability profiles
    solar_thermal = pd.read_csv(
        snakemake.input.solar_thermal, index_col=0, parse_dates=True
    )
    gshp_cop = pd.read_csv(snakemake.input.gshp_cop, index_col=0, parse_dates=True)

    # Share of district heating at each node
    district_heat_share = pd.read_csv(snakemake.input.district_heat_share, index_col=0)

    # Load data required for aviation and navigation
    # TODO follow the same structure as land transport and heat

    # Load industry demand data
    industrial_demand = pd.read_csv(
        snakemake.input.industrial_demand, index_col=0, header=0
    )  # * 1e6


    ##########################################################################
    ############## Functions adding different carrires and sectors ###########
    ##########################################################################

    add_co2(n, costs)  # TODO add costs

    # TODO This might be transferred to add_generation, but before apply remove_elec_base_techs(n) from PyPSA-Eur-Sec
    add_oil(n, costs)

    add_gas(n, costs)
    add_generation(n, costs)

    add_hydrogen(n, costs)  # TODO add costs

    add_storage(n, costs)

    H2_liquid_fossil_conversions(n, costs)

    h2_hc_conversions(n, costs)
    add_heat(n, costs)
    add_biomass(n, costs)

    add_industry(n, costs)

    add_shipping(n, costs)

    # Add_aviation runs with dummy data
    add_aviation(n, costs)

    # prepare_transport_data(n)

    add_land_transport(n, costs)

    # if snakemake.config["custom_data"]["transport_demand"]:
    add_rail_transport(n, costs)

    # if snakemake.config["custom_data"]["custom_sectors"]:
    add_agriculture(n, costs)
    add_residential(n, costs)
    add_services(n, costs)

    sopts = snakemake.wildcards.sopts.split("-")

    for o in sopts:
        m = re.match(r"^\d+h$", o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break

    # TODO add co2 limit here, if necessary
    # co2_limit_pu = eval(sopts[0][5:])
    # co2_limit = co2_limit_pu *
    # # Add co2 limit
    # co2_limit = 1e9
    # n.add(
    #     "GlobalConstraint",
    #     "CO2Limit",
    #     carrier_attribute="co2_emissions",
    #     sense="<=",
    #     constant=co2_limit,
    # )

    if options["dac"]:
        add_dac(n, costs)

    if snakemake.config["custom_data"]["water_costs"]:
        add_custom_water_cost(n)

    n.export_to_netcdf(snakemake.output[0])

    # TODO changes in case of myopic oversight
