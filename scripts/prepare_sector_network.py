import os
import pandas as pd
import numpy as np

from helpers import mock_snakemake, prepare_costs, create_network_topology  

def add_hydrogen(n, costs):
    "function to add hydrogen as an energy carrier with its conversion technologies from and to AC"

    n.add("Carrier", "H2")

    n.madd("Bus", nodes + " H2", location=nodes, carrier="H2")

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
        capital_cost=costs.at["fuel cell", "fixed"] *
        costs.at["fuel cell", "efficiency"],
        lifetime=costs.at["fuel cell", "lifetime"],
    )

    cavern_nodes = pd.DataFrame()
    if options['hydrogen_underground_storage']:
        
          h2_salt_cavern_potential = pd.read_csv(snakemake.input.h2_cavern, index_col=0, squeeze=True)
          h2_cavern_ct = h2_salt_cavern_potential[~h2_salt_cavern_potential.isna()]
          cavern_nodes = n.buses[n.buses.country.isin(h2_cavern_ct.index)]

          h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

          # assumptions: weight storage potential in a country by population
          # TODO: fix with real geographic potentials
          # convert TWh to MWh with 1e6
          h2_pot = h2_cavern_ct.loc[cavern_nodes.country]
          h2_pot.index = cavern_nodes.index
          # h2_pot = h2_pot * cavern_nodes.fraction * 1e6

          n.madd("Store",
            cavern_nodes.index + " H2 Store",
            bus=cavern_nodes.index + " H2",
            e_nom_extendable=True,
            e_nom_max=h2_pot.values,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost
        )

    # hydrogen stored overground (where not already underground)
    h2_capital_cost = costs.at["hydrogen storage tank incl. compressor", "fixed"]
    nodes_overground = cavern_nodes.index.symmetric_difference(nodes)

    n.madd("Store",
            nodes_overground + " H2 Store",
            bus=nodes_overground + " H2",
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost
            )

    attrs = ["bus0", "bus1", "length"]
    h2_links = pd.DataFrame(columns=attrs)

    candidates = pd.concat({"lines": n.lines[attrs],
                            "links": n.links.loc[n.links.carrier == "DC", attrs]})

    for candidate in candidates.index:
        buses = [candidates.at[candidate, "bus0"],
                  candidates.at[candidate, "bus1"]]
        buses.sort()
        name = f"H2 pipeline {buses[0]} -> {buses[1]}"
        if name not in h2_links.index:
            h2_links.at[name, "bus0"] = buses[0]
            h2_links.at[name, "bus1"] = buses[1]
            h2_links.at[name, "length"] = candidates.at[candidate, "length"]

    # TODO Add efficiency losses
    n.madd("Link",
            h2_links.index,
            bus0=h2_links.bus0.values + " H2",
            bus1=h2_links.bus1.values + " H2",
            p_min_pu=-1,
            p_nom_extendable=True,
            length=h2_links.length.values,
            capital_cost=costs.at['H2 (g) pipeline',
                                  'fixed'] * h2_links.length.values,
            carrier="H2 pipeline",
            lifetime=costs.at['H2 (g) pipeline', 'lifetime']
            )


def add_co2(n, costs):
    from types import SimpleNamespace
    spatial = SimpleNamespace()

    spatial.nodes = nodes

    spatial.co2 = SimpleNamespace()

    if options["co2_network"]:
        spatial.co2.nodes = nodes + " co2 stored"
        spatial.co2.locations = nodes
        spatial.co2.vents = nodes + " co2 vent"
    else:
        spatial.co2.nodes = ["co2 stored"]
        spatial.co2.locations = ["Africa"]
        spatial.co2.vents = ["co2 vent"]

    spatial.co2.df = pd.DataFrame(vars(spatial.co2), index=nodes)

    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2",
          co2_emissions=-1.)

    # this tracks CO2 in the atmosphere
    n.add("Bus",
        "co2 atmosphere",
        location="Africa",
        carrier="co2"
    )

    # can also be negative
    n.add("Store",
        "co2 atmosphere",
        e_nom_extendable=True,
        e_min_pu=-1,
        carrier="co2",
        bus="co2 atmosphere"
    )

    # this tracks CO2 stored, e.g. underground
    n.madd("Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored"
    )

    n.madd("Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        e_nom_max=np.inf,
        capital_cost=options['co2_sequestration_cost'],
        carrier="co2 stored",
        bus=spatial.co2.nodes
    )

   
    n.madd("Link",
        spatial.co2.vents,
        bus0=spatial.co2.nodes,
        bus1="co2 atmosphere",
        carrier="co2 vent",
        efficiency=1.,
        p_nom_extendable=
        True
    )

    #logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (1 - co2_links.underwater_fraction) * costs.at['CO2 pipeline', 'fixed'] * co2_links.length
    cost_submarine = co2_links.underwater_fraction * costs.at['CO2 submarine pipeline', 'fixed'] * co2_links.length
    capital_cost = cost_onshore + cost_submarine

    n.madd("Link",
        co2_links.index,
        bus0=co2_links.bus0.values + " co2 stored",
        bus1=co2_links.bus1.values + " co2 stored",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=co2_links.length.values,
        capital_cost=capital_cost.values,
        carrier="CO2 pipeline",
        lifetime=costs.at['CO2 pipeline', 'lifetime'])

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        # from helper import mock_snakemake #TODO remove func from here to helper script
        snakemake = mock_snakemake("prepare_sector_network",
                                   simpl="",
                                   clusters="4")
    # TODO add mock_snakemake func

    # TODO fetch from config

    n = pypsa.Network(snakemake.input.network)

    nodes = n.buses.index

    # costs = pd.read_csv( "{}/pypsa-earth-sec/data/costs.csv".format(os.path.dirname(os.getcwd())))

    Nyears = n.snapshot_weightings.generators.sum() / 8760

    costs = prepare_costs(
        snakemake.input.costs,
        snakemake.config["costs"]["USD2013_to_EUR2013"],
        snakemake.config["costs"]["discountrate"],
        Nyears,
        snakemake.config["costs"]["lifetime"],
    )
    # TODO logging

    # TODO fetch options from the config file

    options = {"co2_network": True,
               "co2_sequestration_potential": 200,  #MtCO2/a sequestration potential for Europe
                "co2_sequestration_cost": 10,
                "hydrogen_underground_storage": True,
                "h2_cavern": True}   #EUR/tCO2 for sequestration of CO2}
    
    add_hydrogen(n, costs)      #TODO add costs
    
    add_co2(n, costs)      #TODO add costs

    # TODO define spatial (for biomass and co2)

    # TODO changes in case of myopic oversight

    # TODO add co2 tracking function

    # TODO add generation

    # TODO add storage  HERE THE H2 CARRIER IS ADDED IN PYPSA-EUR-SEC

    # TODO add options as in PyPSA-EUR-SEC
