def add_hydrogen(n, costs):
    "function to add hydrogen as an energy carrier with its conversion technologies from and to AC"

    n.add('Carrier', 'H2')

    n.madd("Bus",
        nodes + " H2",
        location=nodes,
        carrier="H2"
    )    

    n.madd("Link",
        nodes + " H2 Electrolysis",
        bus1=nodes + " H2",
        bus0=nodes,
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=costs.at["electrolysis", "efficiency"],
        capital_cost=costs.at["electrolysis", "fixed"],
        lifetime=costs.at['electrolysis', 'lifetime']
    )

    n.madd("Link",
        nodes + " H2 Fuel Cell",
        bus0=nodes + " H2",
        bus1=nodes,
        p_nom_extendable=True,
        carrier ="H2 Fuel Cell",
        efficiency=costs.at["fuel cell", "efficiency"],
        capital_cost=costs.at["fuel cell", "fixed"] * costs.at["fuel cell", "efficiency"], #NB: fixed cost is per MWel
        lifetime=costs.at['fuel cell', 'lifetime']
    )

    cavern_nodes = pd.DataFrame()
    if options['hydrogen_underground_storage']:
         h2_salt_cavern_potential = pd.read_csv(snakemake.input.h2_cavern, index_col=0, squeeze=True)
         h2_cavern_ct = h2_salt_cavern_potential[~h2_salt_cavern_potential.isna()]
         cavern_nodes = pop_layout[pop_layout.ct.isin(h2_cavern_ct.index)]

         h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

         # assumptions: weight storage potential in a country by population
         # TODO: fix with real geographic potentials
         # convert TWh to MWh with 1e6
         h2_pot = h2_cavern_ct.loc[cavern_nodes.ct]
         h2_pot.index = cavern_nodes.index
         h2_pot = h2_pot * cavern_nodes.fraction * 1e6

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

           )
