# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Created on Sun May 30 18:11:07 2021.

@author: haz43975
"""


# -*- coding: utf-8 -*-
"""
Created on Tue May  4 10:22:36 2021

@author: haz43975
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

# %%

# base_path = os.path.dirname(os.path.realpath(__file__))
# dataset_paths = {'IGG': os.path.join(base_path, 'IGG', 'data'),
#                      'EMAP': os.path.join(base_path, 'EMAP', 'data')}

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_db",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
            demand="AB",
            h2export="120",
        )

    n0 = pypsa.Network(snakemake.input.network)

    tech_colors = snakemake.params.tech_colors


# %%
# def summary_h2(n, t):
t = 720

n = n0.copy()
# n = pypsa.Network("../results/MA_REALISTIC_2030/postnetworks/elec_s_195_ec_lc1.0_Co2L_3H_2030_0.071_AP_428export.nc")
# n = pypsa.Network("../results/MA_REALISTIC_2030_Q0_NoGreeness/postnetworks/elec_s_198_ec_lc1.0_Co2L_3H_2030_0.071_AP_0export.nc")
# n = pypsa.Network("../results/MA_REALISTIC_2030_Q0_oilnew_13/postnetworks/elec_s_213_ec_lc1.0_Co2L_720H_2030_0.071_AP_0export.nc")
summary_index = (n0.buses.loc[n0.buses.carrier == "AC"].index).sort_values()

nodes = n.buses.loc[n.buses.carrier == "AC"].index.tolist()
gens = n.generators_t.p.rename_axis(None, axis=1) * t  # /1e3
loads = n.loads_t.p.rename_axis(None, axis=1) * t  # /1e3
stores = n.stores_t.p.rename_axis(None, axis=1) * t  # /1e3
storage = n.storage_units_t.p.rename_axis(None, axis=1) * t  # /1e3

pipelines_h2 = n.links_t.p0.filter(like="H2 pipeline")
ac_lines = n.lines_t.p0.rename(columns=dict(n.lines.bus0 + " -> " + n.lines.bus1))

dc_lines = n.links_t.p0[
    n.links_t.p0.columns.intersection(n.links[n.links.carrier == "DC"].index.tolist())
].rename(
    columns=dict(
        n.links[n.links.carrier == "DC"].bus0
        + " -> "
        + n.links[n.links.carrier == "DC"].bus1
    )
)

summary_h2 = pd.DataFrame(index=n0.buses.loc[n0.buses.carrier == "AC"].index)

solar = (gens.filter(regex="solar$")).reset_index

summary_elec = pd.DataFrame(index=n0.buses.loc[n0.buses.carrier == "AC"].index)

db = pd.DataFrame(columns=["node_id", "carrier", "flow", "tech", "value"])

names = {"g": "Generator"}


def populate_db(tech_col, carrier, flow, tech, ngv=False):  # TODO Add scenario id
    global db
    dbf = tech_col.copy()
    # if tech != 'ac':
    #     dbf.name=dbf.name.str.replace(' '+tech, '')
    dbf = (
        dbf.stack()
        .reset_index(level=0)
        .rename(columns={"snapshot": "DateTime", 0: "value"})
        .reset_index()
        .rename(columns={"index": "node_id"})
    )
    dbf.node_id = dbf.node_id.str.replace(" " + tech, "")
    # dbf.columns = ['node_id', 'value']
    dbf["carrier"] = carrier
    dbf["flow"] = flow
    dbf["tech"] = tech
    if flow == "s":
        dbf["value"] = dbf["value"]
    else:
        if ngv == True:
            dbf["value"] = -1 * abs(dbf["value"])
        elif ngv == False:
            dbf["value"] = abs(dbf["value"])

    db = db.append(dbf)


def add_gen(tech, carrier, reg=False):
    if not reg:
        tech_col = gens.filter(like=tech)
    else:
        tech_col = gens.filter(regex=tech + "$")

    # tech_col.columns = tech_col.columns.str.replace(' '+tech, '')
    populate_db(tech_col, carrier, "g", tech, ngv=False)
    # summary_elec['{0}_g_{1}'.format(carrier, tech.replace(' ', '_'))] = tech_col.sum()


def add_load(tech, carrier, reg=False):
    global db
    if tech == "ac":
        ac_labels = loads.stack().reset_index(level=1).level_1
        ac_labels = ac_labels[ac_labels.str.len() < 11].unique()
        tech_col = loads.filter(ac_labels.tolist())
        # ac_labels = loads.reset_index()[loads.reset_index().name.str.len()<7].name.tolist() #TODO hard coded
        # tech_col = loads[ac_labels]

        # summary_elec['elec_l_ac'] = loads
    else:
        if reg == False:
            tech_col = loads.filter(regex=tech)  #
        else:
            tech_col = loads.filter(regex=tech + "$")  #

    populate_db(tech_col, carrier, "l", tech, ngv=True)

    # tech_col.index = tech_col.name.apply(lambda x: x.replace(' '+tech, ''))
    # summary_elec['{0}_l_{1}'.format(carrier, tech.replace(' ', '_'))] = -tech_col[0]


# Check node_id in db


def add_conv(tech, carrier, p, ngv, reg=False):
    global db
    if p == 0:
        links = n.links_t.p0.rename_axis(None, axis=1) * t  # /1e3
    elif p == 1:
        links = n.links_t.p1.rename_axis(None, axis=1) * t  # /1e3
    elif p == 2:
        links = n.links_t.p2.rename_axis(None, axis=1) * t  # /1e3
    elif p == 3:
        links = n.links_t.p3.rename_axis(None, axis=1) * t  # /1e3
    else:
        links = n.links_t.p4.rename_axis(None, axis=1) * t  # /1e3

    if tech == "battery charger" or tech == "battery discharger":
        drop_list = links.filter(like="home battery").columns.tolist()
        tech_col = links.drop(drop_list, axis=1).filter(like=tech)
    else:
        if not reg:
            tech_col = links.filter(like=tech)
        else:
            tech_col = links.filter(regex=tech + "$")
    populate_db(tech_col, carrier, "c", tech, ngv)


#    tech_col.index = tech_col.name.apply(lambda x: x.replace(' '+tech, ''))
#    summary_elec['{0}_c_{1}'.format(carrier, tech.replace(' ', '_'))] = -tech_col[0]


def add_store(tech, carrier, reg=False):
    global db
    if not reg:
        tech_col = stores.filter(like=tech)
    else:
        tech_col = stores.filter(regex=tech + "$")

    if tech == "co2 atmosphere" or tech == "co2 stored":
        tech_col *= -1
    populate_db(tech_col, carrier, "s", tech)


# tech_col.index = tech_col.name.apply(lambda x: x.replace(' '+tech, ''))
# summary_elec['{0}_s_{1}'.format(carrier, tech.replace(' ', '_'))] = tech_col[0]


def add_storage(
    tech, carrier, reg=False
):  # TODO commented out because there is no storage untis
    global db
    if not reg:
        tech_col = storage.filter(like=tech)
    else:
        tech_col = storage.filter(regex=tech + "$")
    populate_db(tech_col, carrier, "s", tech)
    # tech_col.index = tech_col.name.apply(lambda x: x.replace(' '+tech, ''))
    # summary_elec['{0}_s_{1}'.format(carrier, tech.replace(' ', '_'))] = tech_col[0]


def net_flow(co_code, tech, carrier, flow):
    global db
    if tech == "h2flow":
        tech_df = pipelines_h2
    elif tech == "acflow":
        tech_df = ac_lines
    elif tech == "dcflow":
        tech_df = dc_lines
    else:
        pass

    inflow = tech_df.filter(regex="{}$".format(co_code)).sum(axis=1) * t
    outflow = tech_df.filter(like="{} -> ".format(co_code)).sum(axis=1) * t

    dbf = pd.DataFrame()
    dbf["DateTime"] = inflow.index.copy()
    dbf["node_id"] = co_code
    dbf["carrier"] = carrier
    dbf["flow"] = flow
    dbf["tech"] = tech
    dbf["value"] = (inflow - outflow).reset_index(drop=True)  # /10**6
    db = db.append(dbf)
    return dbf


temp = pd.DataFrame(data=nodes)
temp[0].apply(net_flow, args=("h2flow", "h2", "t"))
temp[0].apply(net_flow, args=("acflow", "hv", "t"))
temp[0].apply(net_flow, args=("dcflow", "hv", "t"))


add_gen("solar", "hv", reg=True)
add_gen("onwind", "hv")
add_gen("offwind-ac", "hv")
add_gen("offwind-dc", "hv")
add_gen("ror", "hv")

add_conv("H2 export", "h2", 0, True)

add_conv("H2 Electrolysis", "hv", 0, True)
add_conv("H2 Fuel Cell", "hv", 1, False)
# add_conv("H2 Export", "hv", 1, False)
add_conv("DAC", "hv", 2, True)
add_conv("helmeth", "hv", 0, True)
add_conv("electricity distribution grid", "hv", 0, True)
add_conv("OCGT", "hv", 1, False)
add_conv("battery charger", "hv", 0, True, reg=True)
add_conv("battery discharger", "hv", 1, False, reg=True)
add_conv("urban central gas CHP", "hv", 1, False, reg=True)
add_conv("urban central gas CHP CC", "hv", 1, False)
add_conv("urban central solid biomass CHP", "hv", 1, False, reg=True)
add_conv("urban central solid biomass CHP CC", "hv", 1, False)

add_store("battery", "hv", reg=True)
add_store("battery storage", "hv")
add_store("home battery", "hv")

add_storage("PHS", "hv")  # TODO commented out because there is no storage untis
add_storage("hydro", "hv")


add_load("H2 for shipping", "h2")
add_load("H2 for industry", "h2")
add_load("land transport fuel cell", "h2")

add_conv("H2 Electrolysis", "h2", 1, False)
add_conv("H2 Fuel Cell", "h2", 0, True)
add_conv("Fischer-Tropsch", "h2", 0, True)
add_conv("Sabatier", "h2", 1, True)
add_conv("SMR", "h2", 1, False, reg=True)
add_conv("SMR CC", "h2", 1, False)

add_store("H2", "h2", reg=True)
add_store("H2 Store", "h2")


add_gen("solar rooftop", "hv")

add_conv("electricity distribution grid", "hv", 1, False)
add_conv("BEV charger", "hv", 0, True)
add_conv("V2G", "hv", 1, False)
add_conv("residential rural ground heat pump", "hv", 0, True)
add_conv("residential rural resistive heater", "hv", 0, True)
add_conv("services rural ground heat pump", "hv", 0, True)
add_conv("residential rural resistive heater", "hv", 0, True)
add_conv("urban central air heat pump", "hv", 0, True)
add_conv("urban central resistive heater", "hv", 0, True)
add_conv("home battery charger", "hv", 0, True)
add_conv("home battery discharger", "hv", 1, False)


add_load("ac", "hv")
add_load("industry electricity", "hv")

add_gen("residential rural solar thermal collector", "heat")
add_gen("services rural solar thermal collector", "heat")
add_gen("urban central solar thermal collector", "heat")

add_load("residential rural heat", "heat")
add_load("services rural heat", "heat")
add_load("urban central heat", "heat")
add_load("low-temperature heat for industry", "heat")

add_conv("residential rural water tanks charger", "heat", 0, True)
add_conv("services rural water tanks charger", "heat", 0, True)
add_conv("urban central water tanks charger", "heat", 0, True)
add_conv("residential rural ground heat pump", "heat", 1, False)
add_conv("residential rural water tanks discharger", "heat", 1, False)
add_conv("residential rural resistive heater", "heat", 1, False)
add_conv("residential rural gas boiler", "heat", 1, False)
add_conv("services rural ground heat pump", "heat", 1, False)
add_conv("services rural water tanks discharger", "heat", 1, False)
add_conv("services rural resistive heater", "heat", 1, False)
add_conv("services rural gas boiler", "heat", 1, False)
add_conv("urban central air heat pump", "heat", 1, False)
add_conv("urban central water tanks discharger", "heat", 1, False)
add_conv("urban central resistive heater", "heat", 1, False)
add_conv("urban central gas boiler", "heat", 1, False)
add_conv("H2 Fuel Cell", "heat", 2, False)
add_conv("urban central gas CHP", "heat", 2, False)
add_conv("urban central gas CHP CC", "heat", 2, False)
add_conv("urban central solid biomass CHP", "heat", 2, False)
add_conv("urban central solid biomass CHP CC", "heat", 2, False)
add_conv("DAC", "heat", 3, True)
add_conv("Fischer-Tropsch", "heat", 3, False)


add_conv("services urban decentral DAC", "co2", 0, False)
add_conv("urban central DAC", "co2", 0, False)

add_conv("process emissions", "co2", 1, True, True)
add_conv("process emissions CC", "co2", 1, True)
add_conv("co2 vent", "co2", 1, True)

add_conv("OCGT", "co2", 2, True)
add_conv("biogas to gas", "co2", 2, False)
add_conv("biomass EOP", "co2", 2, True)
add_conv("residential rural gas boiler", "co2", 2, True)
add_conv("services rural gas boiler", "co2", 2, True)
add_conv("residential urban decentral gas boiler", "co2", 2, True)
add_conv("services urban decentral gas boiler", "co2", 2, True)
add_conv("urban central gas boiler", "co2", 2, True)
add_conv("solid biomass for industry CC", "co2", 2, False)
add_conv("gas for industry", "co2", 2, True)
add_conv("gas for industry CC", "co2", 2, True)


add_load("industry oil emissions", "co2")
add_load("shipping oil emissions", "co2")
add_load("land transport oil emissions", "co2")
add_load("aviation oil emissions", "co2")
add_load("residential oil emissions", "co2")
add_load("residential biomass emissions", "co2")
add_load("services biomass emissions", "co2")

# add_store("co2 stored", "co2")
add_store("co2 atmosphere", "co2")

add_conv("Fischer-Tropsch", "oil", 1, False)

add_load("naphtha for industry", "oil")
add_load("residential oil", "oil", reg=True)
add_load("rail transport oil", "oil", reg=True)
add_load("agriculture oil", "oil", reg=True)
add_load("shipping oil", "oil", reg=True)
add_load("land transport oil", "oil", reg=True)  # mistakenly add oil emissions
add_load("kerosene for aviation", "oil")
add_store("oil Store", "oil", 1)
add_gen("oil", "oil")

# add_load("gas for industry", "gas")
add_conv("OCGT", "gas", 0, True)
add_conv("residential rural gas boiler", "gas", 0, True)
add_conv("services rural gas boiler", "gas", 0, True)
add_conv("residential urban decentral gas boiler", "gas", 0, True)
add_conv("services urban decentral gas boiler", "gas", 0, True)
add_conv("gas for industry", "gas", 0, True, True)
add_conv("gas for industry CC", "gas", 0, True)
add_conv("urban central gas boiler", "gas", 0, True)

add_conv("Sabatier", "gas", 1, False)
add_conv("helmeth", "gas", 1, False)

add_store("gas Store", "gas")
add_gen("gas", "gas")

add_conv("biogas to gas", "gas", 1, False)

add_conv("gas for industry", "gas", 1, True)
add_conv("gas for industry CC", "gas", 1, True)


add_conv("co2 vent", "co2 stored", 0, True)
# add_conv("CO2 pipeline", "co2 stored", 0, True)
add_conv("DAC", "co2 stored", 1, True)

add_conv("Fischer-Tropsch", "co2 stored", 2, False)
add_conv("Sabatier", "co2 stored", 2, False)
add_conv("helmeth", "co2 stored", 2, False)
add_conv("process emissions CC", "co2 stored", 2, True)

add_conv("solid biomass for industry CC", "co2 stored", 3, True)
add_conv("gas for industry CC", "co2 stored", 3, True)

# add_store("co2 stored")
add_store("co2 stored", "co2 stored")


# add
# REMEMBER TO ADD ALL DECENTRAL SHIT
# %%

# summary_elec['h2_t_pipeline'] = summary_elec.apply(lambda row: h2_net_flow(n0, row.name, 24), axis=1)

h2_flows = pd.DataFrame(index=pipelines_h2.index.copy(), columns=["node_id", "flow"])

# summary_elec['h2_balance'] = summary_elec.sum(axis=1)
# summary_elec=summary_elec.apply(lambda x: round(x, 2))

db.reset_index(drop=True, inplace=True)
# round(db).to_csv('db_fraction.csv')
round(db).to_csv(snakemake.output.db)
yearly_agg = round(db.groupby([db.node_id, db.carrier, db.flow, db.tech]).sum() / 1e3)


# yearly_agg.to_csv('summary_db.csv')
# yearly_agg.to_csv(snakemake.output.yr_agg)
# %%
def calc_energy_flow(carrier, node_id):
    agg = yearly_agg.reset_index()
    agg = agg[(agg.carrier == carrier)]
    agg.value = agg.value.apply(int)
    if node_id == "all":
        agg = agg.groupby("tech").sum().reset_index()
    else:
        agg = agg[agg.node_id.str.contains(node_id)].groupby("tech").sum().reset_index()
    return agg


def fetch_data_2(carrier, node_id):
    agg = db[db.DateTime == "2013-02-01"].drop("DateTime", axis=1)

    return agg[(agg.carrier == carrier) & (agg.node_id == node_id)]


def energy_pie(carrier, node_id, sign):
    sign_dict = {1: "generation", -1: "consumption"}
    agg = yearly_agg.reset_index()
    agg = agg[(agg.carrier == carrier) & (agg.value * sign > 0)]

    if node_id == "all":
        agg = agg[(agg.carrier == carrier) & (agg.value * sign > 0) & (agg.flow != "t")]
    else:
        agg = agg[agg.node_id.str.contains(node_id)]
    agg = agg.groupby("tech").sum().reset_index()
    agg["pct"] = round(agg["value"] / agg.value.sum(), 3)
    if agg.pct.sum() < 1:
        agg = agg.append(
            pd.DataFrame([["other", 0, 1 - agg.pct.sum()]], columns=agg.columns)
        )
    agg = agg[agg.pct > 0.009]

    fig1, ax1 = plt.subplots()  # figsize=(6, 4))
    ax1.pie(
        agg.pct,
        labels=agg.tech,
        autopct="%1.0f%%",
        colors=[tech_colors.get(key) for key in agg.tech.values.tolist()],
        explode=[0.05] * len(agg),
    )
    ax1.axis("equal")
    plt.title(
        "Yearly aggregate {0} of {1} at {2} node(s)\n".format(
            sign_dict[sign], carrier, node_id
        )
        + "Value = {} GWh".format(round(agg.value.sum(), 1)),
        # bbox={"facecolor": "0.8", "pad": 5},
    )
    plt.show()
    fig1.savefig(
        "Yearly_aggregate_{0}_of_{1}_at_{2}_node(s).png".format(
            sign_dict[sign], carrier, node_id
        ),
        dpi=100,
    )
