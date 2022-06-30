#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:50:42 2022

@author: user
"""
import os

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle, Ellipse


def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):

        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)

        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue

            names = ifind.index[ifind == i]

            c.df.loc[names, "location"] = names.str[:i]


def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    fig = ax.get_figure()

    def axes2pt():
        return np.diff(ax.transData.transform([(0, 0), (1, 1)]), axis=0)[0] * (
            72.0 / fig.dpi
        )

    ellipses = []
    if not dont_resize_actively:

        def update_width_height(event):
            dist = axes2pt()
            for e, radius in ellipses:
                e.width, e.height = 2.0 * radius * dist

        fig.canvas.mpl_connect("resize_event", update_width_height)
        ax.callbacks.connect("xlim_changed", update_width_height)
        ax.callbacks.connect("ylim_changed", update_width_height)

    def legend_circle_handler(
        legend, orig_handle, xdescent, ydescent, width, height, fontsize
    ):
        w, h = 2.0 * orig_handle.get_radius() * axes2pt()
        e = Ellipse(
            xy=(0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent),
            width=w,
            height=w,
        )
        ellipses.append((e, orig_handle.get_radius()))
        return e

    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}


def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0, 0), radius=(s / scale) ** 0.5, **kw) for s in sizes]


#############################################
# plot Hydrogen infrastructure map
#############################################
def plot_h2_infra(network):
    n = network.copy()

    # assign_location(n)

    bus_size_factor = 1e6
    linewidth_factor = 1e5
    # MW below which not drawn
    line_lower_threshold = 1e1
    bus_color = "m"
    link_color = "c"

    n.links.loc[:, "p_nom_opt"] = n.links.loc[:, "p_nom_opt"]

    # n.links.loc[n.links.carrier == "H2 Electrolysis"].p_nom_opt

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    elec = n.links.index[n.links.carrier == "H2 Electrolysis"]

    bus_sizes = (
        n.links.loc[elec, "p_nom_opt"].groupby(n.links.loc[elec, "bus0"]).sum()
        / bus_size_factor
    )

    # make a fake MultiIndex so that area is correct for legend
    bus_sizes.index = pd.MultiIndex.from_product([bus_sizes.index, ["electrolysis"]])

    n.links.drop(n.links.index[n.links.carrier != "H2 pipeline"], inplace=True)

    link_widths = n.links.p_nom_opt / linewidth_factor
    link_widths[n.links.p_nom_opt < line_lower_threshold] = 0.0

    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    print(link_widths.sort_values())

    print(n.links[["bus0", "bus1"]])

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

    fig.set_size_inches(10.5, 9)

    n.plot(
        bus_sizes=bus_sizes,
        bus_colors={"electrolysis": bus_color},
        link_colors=link_color,
        link_widths=link_widths,
        branch_components=["Link"],
        color_geomap={"ocean": "lightblue", "land": "oldlace"},
        ax=ax,
        boundaries=(-20, 0, 25, 40),
    )

    handles = make_legend_circles_for(
        [50000, 10000], scale=bus_size_factor, facecolor=bus_color
    )
    labels = ["{} GW".format(s) for s in (50, 10)]
    l2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.01, 1.01),
        labelspacing=0.8,
        framealpha=1.0,
        title="Electrolyzer capacity",
        handler_map=make_handler_map_to_scale_circles_as_in(ax),
    )
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in (50, 10):
        handles.append(
            plt.Line2D([0], [0], color=link_color, linewidth=s * 1e3 / linewidth_factor)
        )
        labels.append("{} GW".format(s))
    l1_1 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.32, 1.01),
        framealpha=1,
        labelspacing=0.8,
        handletextpad=1.5,
        title="H2 pipeline capacity",
    )
    ax.add_artist(l1_1)

    # fig.savefig(snakemake.output.hydrogen, bbox_inches='tight', transparent=True,
    fig.savefig(
        snakemake.output.map.replace("-costs-all", "-h2_network"), bbox_inches="tight"
    )


def plot_transmission_topology(network):

    n = network.copy()
    bus_size_factor = 1e5  # Def 1e5
    linewidth_factor = 2e4  # Def 1e4
    line_lower_threshold = 1e2  # MW below which not drawn. Def 1e3

    DC_lines = n.links[n.links.carrier == "DC"]

    n.links = n.links[n.links.carrier == "H2 pipeline"]
    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    n.lines.append(DC_lines[["bus0", "bus1"]])

    n.madd("Line", names=DC_lines.index, bus0=DC_lines.bus0, bus1=DC_lines.bus1)

    fig = plt.figure()
    fig.set_size_inches(10.5, 9)

    n.plot(
        branch_components=["Link", "Line"],
        boundaries=(-20, 0, 25, 40),
        color_geomap={"ocean": "lightblue", "land": "oldlace"},
        line_colors="darkblue",
        link_colors="turquoise",
        link_widths=5,
        bus_sizes=0.03,
        bus_colors="red",
        line_widths=1,
    )

    # Legend
    Elec_Circle = plt.Line2D(
        [0],
        [0],
        marker="o",
        color="darkblue",
        label="Clustered node",
        markerfacecolor="red",
        markersize=10,
    )
    elec_Line = plt.Line2D(
        [0],
        [0],
        marker="_",
        color="darkblue",
        label="Electricity Lines",
        markerfacecolor="w",
        markersize=16,
        lw=4,
    )

    H2_Line = plt.Line2D(
        [0],
        [0],
        marker="_",
        color="turquoise",
        label="H2 Pipeline",
        markerfacecolor="w",
        markersize=16,
        lw=4,
    )

    plt.legend(handles=[Elec_Circle, elec_Line, H2_Line], loc="upper left")

    fig.savefig(
        snakemake.output.map.replace("-costs-all", "-full_topology"),
        bbox_inches="tight",
    )


preferred_order = pd.Index(
    [
        "transmission lines",
        "hydroelectricity",
        "hydro reservoir",
        "run of river",
        "pumped hydro storage",
        "solid biomass",
        "biogas",
        "onshore wind",
        "offshore wind",
        "offshore wind (AC)",
        "offshore wind (DC)",
        "solar PV",
        "solar thermal",
        "solar",
        "building retrofitting",
        "ground heat pump",
        "air heat pump",
        "heat pump",
        "resistive heater",
        "power-to-heat",
        "gas-to-power/heat",
        "CHP",
        "OCGT",
        "gas boiler",
        "gas",
        "natural gas",
        "helmeth",
        "methanation",
        "hydrogen storage",
        "power-to-gas",
        "power-to-liquid",
        "battery storage",
        "hot water storage",
        "CO2 sequestration",
    ]
)


def rename_techs(label):

    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]

    rename_if_contains = [
        "CHP",
        "gas boiler",
        "biogas",
        "solar thermal",
        "air heat pump",
        "ground heat pump",
        "resistive heater",
        "Fischer-Tropsch",
    ]

    rename_if_contains_dict = {
        "water tanks": "hot water storage",
        "retrofitting": "building retrofitting",
        "H2": "hydrogen storage",
        "battery": "battery storage",
        "CCS": "CCS",
    }

    rename = {
        "solar": "solar PV",
        "Sabatier": "methanation",
        "offwind": "offshore wind",
        "offwind-ac": "offshore wind (AC)",
        "offwind-dc": "offshore wind (DC)",
        "onwind": "onshore wind",
        "ror": "hydroelectricity",
        "hydro": "hydroelectricity",
        "PHS": "hydroelectricity",
        "co2 Store": "DAC",
        "co2 stored": "CO2 sequestration",
        "AC": "transmission lines",
        "DC": "transmission lines",
        "B2B": "transmission lines",
    }

    for ptr in prefix_to_remove:
        if label[: len(ptr)] == ptr:
            label = label[len(ptr) :]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old, new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old, new in rename.items():
        if old == label:
            label = new
    return label


def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech or "resistive heater" in tech:
        return "power-to-heat"
    elif tech in ["methanation", "hydrogen storage", "helmeth"]:
        return "power-to-gas"
    elif tech in ["OCGT", "CHP", "gas boiler"]:
        return "gas-to-power/heat"
    elif "solar" in tech:
        return "solar"
    elif tech == "Fischer-Tropsch":
        return "power-to-liquid"
    elif "offshore wind" in tech:
        return "offshore wind"
    else:
        return tech


def plot_map(
    network,
    components=[
        "links",
        "generators",
        "stores",
    ],  # "storage_units"], #TODO uncomment after adding storage units
    bus_size_factor=1.7e10,
    transmission=False,
    geometry=True,
):

    n = network.copy()
    assign_location(n)
    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    costs = pd.DataFrame(index=n.buses.index)

    for comp in components:
        df_c = getattr(n, comp)
        df_c["nice_group"] = df_c.carrier.map(rename_techs_tyndp)

        attr = "e_nom_opt" if comp == "stores" else "p_nom_opt"

        costs_c = (
            (df_c.capital_cost * df_c[attr])
            .groupby([df_c.location, df_c.nice_group])
            .sum()
            .unstack()
            .fillna(0.0)
        )
        costs = pd.concat([costs, costs_c], axis=1)

        print(comp, costs)
    costs = costs.groupby(costs.columns, axis=1).sum()

    costs.drop(list(costs.columns[(costs == 0.0).all()]), axis=1, inplace=True)

    new_columns = preferred_order.intersection(costs.columns).append(
        costs.columns.difference(preferred_order)
    )
    costs = costs[new_columns]

    for item in new_columns:
        if item not in tech_colors:
            print("Warning!", item, "not in config/plotting/tech_colors")

    costs = costs.stack()  # .sort_index()

    n.links.drop(
        n.links.index[(n.links.carrier != "DC") & (n.links.carrier != "B2B")],
        inplace=True,
    )

    # drop non-bus
    to_drop = costs.index.levels[0].symmetric_difference(n.buses.index)
    if len(to_drop) != 0:
        print("dropping non-buses", list(to_drop))
        costs.drop(to_drop, level=0, inplace=True, axis=0)

    # make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)

    # PDF has minimum width, so set these to zero
    line_lower_threshold = 500.0
    line_upper_threshold = 1e4
    linewidth_factor = 2e3
    ac_color = "gray"
    dc_color = "m"

    # if snakemake.wildcards["lv"] == "1.0":         #TODO when we add wildcard lv
    # should be zero
    line_widths = n.lines.s_nom_opt - n.lines.s_nom
    link_widths = n.links.p_nom_opt - n.links.p_nom
    title = "Technologies"

    if transmission:
        line_widths = n.lines.s_nom_opt
        link_widths = n.links.p_nom_opt
        linewidth_factor = 2e3
        line_lower_threshold = 0.0
        title = "Technologies"
    else:
        line_widths = (
            n.lines.s_nom_opt - n.lines.s_nom_min
        )  # TODO when we add wildcard lv
        link_widths = n.links.p_nom_opt - n.links.p_nom_min
        title = "Transmission reinforcement"

        if transmission:
            line_widths = n.lines.s_nom_opt
            link_widths = n.links.p_nom_opt
            title = "Total transmission"

    line_widths.loc[line_widths < line_lower_threshold] = 0.0
    link_widths.loc[link_widths < line_lower_threshold] = 0.0

    line_widths.loc[line_widths > line_upper_threshold] = line_upper_threshold
    link_widths.loc[link_widths > line_upper_threshold] = line_upper_threshold

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    fig.set_size_inches(10.5, 9)

    n.plot(
        bus_sizes=costs / bus_size_factor,
        bus_colors=tech_colors,
        line_colors=ac_color,
        link_colors=dc_color,
        line_widths=line_widths / linewidth_factor,
        link_widths=link_widths / linewidth_factor,
        ax=ax,
        boundaries=(-20, 0, 25, 40),
        geomap="10m",
        color_geomap={"ocean": "lightblue", "land": "oldlace"},
    )

    handles = make_legend_circles_for(
        [5e9, 1e9], scale=bus_size_factor, facecolor="gray"
    )
    labels = ["{} b€/a".format(s) for s in (5, 1)]
    l2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.33, 1.005),
        labelspacing=1.0,
        framealpha=1.0,
        title="System cost",
        fontsize=12,
        handler_map=make_handler_map_to_scale_circles_as_in(ax),
    )
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in list(plot_labeles.keys()):
        handles.append(plt.Line2D([0], [0], color=tech_colors[s], linewidth=5))
        labels.append("{}".format(s))

    l1_1 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.001, 1.002),
        framealpha=1,
        labelspacing=0.4,
        handletextpad=1.5,
        fontsize=10,
    )

    ax.add_artist(l1_1)

    # import matplotlib.patches as mpatches
    # red_patch = mpatches.Patch(color='red', label='The red data')
    # plt.legend(handles=[red_patch])

    fig.savefig(snakemake.output.map, transparent=True, bbox_inches="tight")
    # fig.savefig('plot_map.pdf', transparent=True,
    #         bbox_inches="tight")#, dpi=300)


plot_labeles = {
    "onshore wind": "b",
    "offshore wind": "c",
    "hydroelectricity": "",
    "solar": "y",
    "power-to-gas": "#FF1493",
    "gas-to-power/heat": "orange",
    "power-to-heat": "",
    "power-to-liquid": "",
    "DAC": "",
    "electricity distribution grid": "",
}


nice_names = {
    # OCGT: "Gas",
    # OCGT marginal: "Gas (marginal)",
    "offwind": "offshore wind",
    "onwind": "onshore wind",
    "battery": "Battery storage",
    "lines": "Transmission lines",
    "AC line": "AC lines",
    "AC-AC": "DC lines",
    "ror": "Run of river",
}

nice_names_n = {
    "offwind": "offshore\nwind",
    "onwind": "onshore\nwind",
    # OCGT: "Gas"
    "H2": "Hydrogen\nstorage",
    # OCGT marginal: "Gas (marginal)"
    "lines": "transmission\nlines",
    "ror": "run of river",
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_network",
            simpl="",
            clusters="72",
            ll="c1",
            opts="Co2L-720H",
            planning_horizons="2030",
        )

    n = pypsa.Network(snakemake.input.network)

    tech_colors = snakemake.config["plotting"]["tech_colors"]
    plot_map(n, transmission=True)
    plot_transmission_topology(n)
    plot_h2_infra(n)
