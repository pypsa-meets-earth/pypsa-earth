# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Plots map with pie charts and cost box bar charts.

Relevant Settings
-----------------
Inputs
------
Outputs
-------
Description
-----------
"""

import os

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    aggregate_costs,
    aggregate_p,
    configure_logging,
    create_logger,
    load_network_for_plots,
)
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Circle, Ellipse

to_rgba = mpl.colors.colorConverter.to_rgba

logger = create_logger(__name__)


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


def set_plot_style():
    plt.style.use(
        [
            "classic",
            "seaborn-white",
            {
                "axes.grid": False,
                "grid.linestyle": "--",
                "grid.color": "0.6",
                "hatch.color": "white",
                "patch.linewidth": 0.5,
                "font.size": 12,
                "legend.fontsize": "medium",
                "lines.linewidth": 1.5,
                "pdf.fonttype": 42,
            },
        ]
    )


def plot_map(n, ax=None, attribute="p_nom", opts={}):
    if ax is None:
        ax = plt.gca()

    # DATA
    line_colors = {
        "cur": "purple",
        "exp": mpl.colors.rgb2hex(to_rgba("red", 0.7), True),
    }
    tech_colors = opts["tech_colors"]

    if attribute == "p_nom":
        # bus_sizes = n.generators_t.p.sum().loc[n.generators.carrier == "load"].groupby(n.generators.bus).sum()
        bus_sizes = pd.concat(
            (
                n.generators.query('carrier != "load"')
                .groupby(["bus", "carrier"])
                .p_nom_opt.sum(),
                n.storage_units.groupby(["bus", "carrier"]).p_nom_opt.sum(),
            )
        )
        line_widths_exp = n.lines.s_nom_opt
        line_widths_cur = n.lines.s_nom_min
        link_widths_exp = n.links.p_nom_opt
        link_widths_cur = n.links.p_nom_min
    else:
        logger.error("plotting of {} has not been implemented yet".format(attribute))

    line_colors_with_alpha = (line_widths_cur / n.lines.s_nom > 1e-3).map(
        {True: line_colors["cur"], False: to_rgba(line_colors["cur"], 0.0)}
    )
    link_colors_with_alpha = (link_widths_cur / n.links.p_nom > 1e-3).map(
        {True: line_colors["cur"], False: to_rgba(line_colors["cur"], 0.0)}
    )

    # FORMAT
    linewidth_factor = opts["map"][attribute]["linewidth_factor"]
    bus_size_factor = opts["map"][attribute]["bus_size_factor"]

    # PLOT
    n.plot(
        line_widths=line_widths_exp / linewidth_factor,
        link_widths=link_widths_exp / linewidth_factor,
        line_colors=line_colors["exp"],
        link_colors=line_colors["exp"],
        bus_sizes=bus_sizes / bus_size_factor,
        bus_colors=tech_colors,
        # boundaries=map_boundaries,
        color_geomap=True,
        geomap=True,
        ax=ax,
    )
    n.plot(
        line_widths=line_widths_cur / linewidth_factor,
        link_widths=link_widths_cur / linewidth_factor,
        line_colors=line_colors_with_alpha,
        link_colors=link_colors_with_alpha,
        bus_sizes=0,
        # boundaries=map_boundaries,
        color_geomap=True,
        geomap=False,
        ax=ax,
    )
    ax.set_aspect("equal")
    ax.axis("off")

    # Rasterize basemap
    # TODO : Check if this also works with cartopy
    for c in ax.collections[:2]:
        c.set_rasterized(True)

    # LEGEND
    handles = []
    labels = []

    for s in (10, 1):
        handles.append(
            plt.Line2D(
                [0], [0], color=line_colors["exp"], linewidth=s * 1e3 / linewidth_factor
            )
        )
        labels.append("{} GW".format(s))
    l1_1 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.24, 1.01),
        frameon=False,
        labelspacing=0.8,
        handletextpad=1.5,
        title="Transmission Exp./Exist.             ",
    )
    ax.add_artist(l1_1)

    handles = []
    labels = []
    for s in (10, 5):
        handles.append(
            plt.Line2D(
                [0], [0], color=line_colors["cur"], linewidth=s * 1e3 / linewidth_factor
            )
        )
        labels.append("/")
    l1_2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.26, 1.01),
        frameon=False,
        labelspacing=0.8,
        handletextpad=0.5,
        title=" ",
    )
    ax.add_artist(l1_2)

    handles = make_legend_circles_for(
        [10e3, 5e3, 1e3], scale=bus_size_factor, facecolor="w"
    )
    labels = ["{} GW".format(s) for s in (10, 5, 3)]
    l2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.01, 1.01),
        frameon=False,
        labelspacing=1.0,
        title="Generation",
        handler_map=make_handler_map_to_scale_circles_as_in(ax),
    )
    ax.add_artist(l2)

    techs = (bus_sizes.index.levels[1]).intersection(
        pd.Index(opts["vre_techs"] + opts["conv_techs"] + opts["storage_techs"])
    )
    handles = []
    labels = []
    for t in techs:
        handles.append(
            plt.Line2D(
                [0], [0], color=tech_colors[t], marker="o", markersize=8, linewidth=0
            )
        )
        labels.append(opts["nice_names"].get(t, t))
    l3 = ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.0),  # bbox_to_anchor=(0.72, -0.05),
        handletextpad=0.0,
        columnspacing=0.5,
        ncol=4,
        title="Technology",
    )

    return fig


def plot_total_energy_pie(n, ax=None):
    if ax is None:
        ax = plt.gca()

    ax.set_title("Energy per technology", fontdict=dict(fontsize="medium"))

    e_primary = aggregate_p(n).drop("load", errors="ignore").loc[lambda s: s > 0]

    patches, texts, autotexts = ax.pie(
        e_primary,
        startangle=90,
        labels=e_primary.rename(opts["nice_names"]).index,
        autopct="%.0f%%",
        shadow=False,
        colors=[opts["tech_colors"][tech] for tech in e_primary.index],
    )
    for t1, t2, i in zip(texts, autotexts, e_primary.index):
        if e_primary.at[i] < 0.04 * e_primary.sum():
            t1.remove()
            t2.remove()


def plot_total_cost_bar(n, ax=None):
    if ax is None:
        ax = plt.gca()

    total_load = (n.snapshot_weightings.generators * n.loads_t.p.sum(axis=1)).sum()
    tech_colors = opts["tech_colors"]

    def split_costs(n):
        costs = aggregate_costs(n).reset_index(level=0, drop=True)
        costs_ex = aggregate_costs(n, existing_only=True).reset_index(
            level=0, drop=True
        )
        return (
            costs["capital"].add(costs["marginal"], fill_value=0.0),
            costs_ex["capital"],
            costs["capital"] - costs_ex["capital"],
            costs["marginal"],
        )

    costs, costs_cap_ex, costs_cap_new, costs_marg = split_costs(n)

    costs_graph = pd.DataFrame(
        dict(a=costs.drop("load", errors="ignore")),
        index=[
            "AC-AC",
            "AC line",
            "onwind",
            "offwind-ac",
            "offwind-dc",
            "solar",
            "OCGT",
            "CCGT",
            "coal",
            "oil",
            "battery",
            "H2",
        ],
    ).dropna()
    bottom = np.array([0.0, 0.0])
    texts = []

    for i, ind in enumerate(costs_graph.index):
        data = np.asarray(costs_graph.loc[ind]) / total_load
        ax.bar([0.5], data, bottom=bottom, color=tech_colors[ind], width=0.7, zorder=-1)
        bottom_sub = bottom
        bottom = bottom + data

        if ind in opts["conv_techs"] + ["AC line"]:
            for c in [costs_cap_ex, costs_marg]:
                if ind in c:
                    data_sub = np.asarray([c.loc[ind]]) / total_load
                    ax.bar(
                        [0.5],
                        data_sub,
                        linewidth=0,
                        bottom=bottom_sub,
                        color=tech_colors[ind],
                        width=0.7,
                        zorder=-1,
                        alpha=0.8,
                    )
                    bottom_sub += data_sub

        if abs(data[-1]) < 5:
            continue

        text = ax.text(
            1.1, (bottom - 0.5 * data)[-1] - 3, opts["nice_names"].get(ind, ind)
        )
        texts.append(text)

    ax.set_ylabel("Average system cost [Eur/MWh]")
    ax.set_ylim([0, opts.get("costs_max", 80)])
    ax.set_xlim([0, 1])
    ax.set_xticklabels([])
    ax.grid(True, axis="y", color="k", linestyle="dotted")


#############################################
# plot Hydrogen infrastructure map
#############################################

# TODO function redundant with plot_h2_infra
# def plot_h2_infra(network):
#     n = network.copy()

#     # assign_location(n)

#     bus_size_factor = 1e5
#     linewidth_factor = 1e3
#     # MW below which not drawn
#     line_lower_threshold = 1e2
#     bus_color = "m"
#     link_color = "c"

#     n.links.loc[:, "p_nom_opt"] = n.links.loc[:, "p_nom_opt"]
#     # n.links.loc[n.links.carrier == "H2 Electrolysis"].p_nom_opt

#     # Drop non-electric buses so they don't clutter the plot
#     n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

#     elec = n.links.index[n.links.carrier == "SMR"]

#     bus_sizes = (
#         n.links.loc[elec, "p_nom_opt"].groupby(n.links.loc[elec, "bus0"]).sum()
#         / bus_size_factor
#     )

#     # make a fake MultiIndex so that area is correct for legend
#     bus_sizes.index = pd.MultiIndex.from_product([bus_sizes.index, ["SMR"]])

#     # n.links.drop(n.links.index[n.links.carrier != "H2 pipeline"], inplace=True)

#     # link_widths = n.links.p_nom_opt / linewidth_factor
#     # link_widths[n.links.p_nom_opt < line_lower_threshold] = 0.0

#     # n.links.bus0 = n.links.bus0.str.replace(" H2", "")
#     # n.links.bus1 = n.links.bus1.str.replace(" H2", "")

#     # print(link_widths.sort_values())

#     # print(n.links[["bus0", "bus1"]])

#     fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

#     fig.set_size_inches(10.5, 9)
#     bus_sizes.index = bus_sizes.index.set_levels(
#         bus_sizes.index.levels[0].str.replace(" gas", ""), level=0
#     )
#     n.plot(
#         bus_sizes=bus_sizes,
#         bus_colors={"SMR": "darkolivegreen"},
#         # link_colors=link_color,
#         # link_widths=link_widths,
#         branch_components=["Link"],
#         color_geomap={"ocean": "lightblue", "land": "oldlace"},
#         ax=ax,
#         boundaries=(-20, 0, 25, 40),
#     )

#     handles = make_legend_circles_for(
#         [5000, 1000], scale=bus_size_factor, facecolor="darkolivegreen"
#     )
#     labels = ["{} GW".format(s) for s in (5, 1)]
#     l2 = ax.legend(
#         handles,
#         labels,
#         loc="upper left",
#         bbox_to_anchor=(0.01, 1.01),
#         labelspacing=0.8,
#         framealpha=1.0,
#         title="SMR capacity",
#         handler_map=make_handler_map_to_scale_circles_as_in(ax),
#     )
#     ax.add_artist(l2)

#     handles = []
#     labels = []

#     for s in (5, 1):
#         handles.append(
#             plt.Line2D([0], [0], color=link_color, linewidth=s * 1e3 / linewidth_factor)
#         )
#         labels.append("{} GW".format(s))
#     l1_1 = ax.legend(
#         handles,
#         labels,
#         loc="upper left",
#         bbox_to_anchor=(0.32, 1.01),
#         framealpha=1,
#         labelspacing=0.8,
#         handletextpad=1.5,
#         title="H2 pipeline capacity",
#     )
#     ax.add_artist(l1_1)

#     # fig.savefig(snakemake.output.hydrogen, bbox_inches='tight', transparent=True,
#     fig.savefig(
#         snakemake.output.map.replace("-costs-all", "-h2_network"), bbox_inches="tight"
#     )


def plot_h2_infra(network):
    n = network.copy()

    # assign_location(n)

    bus_size_factor = 1e5
    linewidth_factor = 4e2
    # MW below which not drawn
    line_lower_threshold = 1e2
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
        # boundaries=(-20, 0, 25, 40),
    )

    handles = make_legend_circles_for(
        [5000, 1000], scale=bus_size_factor, facecolor=bus_color
    )
    labels = ["{} GW".format(s) for s in (5, 1)]
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

    for s in (5, 1):
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


def plot_smr(network):
    n = network.copy()

    # assign_location(n)

    bus_size_factor = 1e5
    linewidth_factor = 1e3
    # MW below which not drawn
    line_lower_threshold = 1e2
    bus_color = "m"
    link_color = "c"

    n.links.loc[:, "p_nom_opt"] = n.links.loc[:, "p_nom_opt"]
    # n.links.loc[n.links.carrier == "H2 Electrolysis"].p_nom_opt

    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    elec = n.links.index[n.links.carrier == "SMR"]

    bus_sizes = (
        n.links.loc[elec, "p_nom_opt"].groupby(n.links.loc[elec, "bus0"]).sum()
        / bus_size_factor
    )

    # make a fake MultiIndex so that area is correct for legend
    bus_sizes.index = pd.MultiIndex.from_product([bus_sizes.index, ["SMR"]])

    # n.links.drop(n.links.index[n.links.carrier != "H2 pipeline"], inplace=True)

    # link_widths = n.links.p_nom_opt / linewidth_factor
    # link_widths[n.links.p_nom_opt < line_lower_threshold] = 0.0

    # n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    # n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    # print(link_widths.sort_values())

    # print(n.links[["bus0", "bus1"]])

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

    fig.set_size_inches(10.5, 9)
    bus_sizes.index = bus_sizes.index.set_levels(
        bus_sizes.index.levels[0].str.replace(" gas", ""), level=0
    )
    n.plot(
        bus_sizes=bus_sizes,
        bus_colors={"SMR": "darkolivegreen"},
        # link_colors=link_color,
        # link_widths=link_widths,
        branch_components=["Link"],
        color_geomap={"ocean": "lightblue", "land": "oldlace"},
        ax=ax,
        # boundaries=(-20, 0, 25, 40),
    )

    handles = make_legend_circles_for(
        [5000, 1000], scale=bus_size_factor, facecolor="darkolivegreen"
    )
    labels = ["{} GW".format(s) for s in (5, 1)]
    l2 = ax.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.01, 1.01),
        labelspacing=0.8,
        framealpha=1.0,
        title="SMR capacity",
        handler_map=make_handler_map_to_scale_circles_as_in(ax),
    )
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in (5, 1):
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
    fig.savefig(snakemake.output.map.replace("-costs-all", "-SMR"), bbox_inches="tight")


def plot_transmission_topology(network):
    n = network.copy()
    bus_size_factor = 1e5  # Def 1e5
    linewidth_factor = 2e4  # Def 1e4
    line_lower_threshold = 1e2  # MW below which not drawn. Def 1e3

    DC_lines = n.links[n.links.carrier == "DC"]

    n.links = n.links[n.links.carrier == "H2 pipeline"]
    n.links.bus0 = n.links.bus0.str.replace(" H2", "")
    n.links.bus1 = n.links.bus1.str.replace(" H2", "")

    n.lines = pd.concat([n.lines, DC_lines[["bus0", "bus1"]]])

    n.madd("Line", names=DC_lines.index, bus0=DC_lines.bus0, bus1=DC_lines.bus1)

    fig = plt.figure()
    fig.set_size_inches(10.5, 9)

    n.plot(
        branch_components=["Link", "Line"],
        # boundaries=(-20, 0, 25, 40),
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
        label="Existing Power Lines",
        markerfacecolor="w",
        markersize=16,
        lw=4,
    )

    H2_Line = plt.Line2D(
        [0],
        [0],
        marker="_",
        color="turquoise",
        label="Allowed H2 Pipeline Routes",
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


def plot_sector_map(
    network,
    components=[
        "links",
        "generators",
        "stores",
    ],  # "storage_units"], #TODO uncomment after adding storage units
    bus_size_factor=2e10,
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
        line_widths = n.lines.s_nom_opt - n.lines.s_nom_min
        line_widths = (
            n.lines.s_nom_opt - n.lines.s_nom_opt
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
        # boundaries=(-20, 0, 25, 40),
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
    fig.savefig(
        snakemake.output.map.replace("pdf", "png"),
        transparent=True,
        bbox_inches="tight",
    )
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
    "OCGT": "Gas",
    "OCGT marginal": "Gas (marginal)",
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
    "OCGT": "Gas",
    "H2": "Hydrogen\nstorage",
    "OCGT marginal": "Gas (marginal)",
    "lines": "transmission\nlines",
    "ror": "run of river",
}


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os

        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_network",
            network="elec",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            attr="p_nom",
            ext="pdf",
        )

    configure_logging(snakemake)

    if snakemake.rule == "plot_network":

        # load africa shape to identify borders of the image
        africa_shape = gpd.read_file(snakemake.input.africa_shape)["geometry"].iloc[0]

        set_plot_style()

        opts = snakemake.params.plotting
        map_figsize = opts["map"]["figsize"]
        map_boundaries = opts["map"]["boundaries"]

        if len(map_boundaries) != 4:
            map_boundaries = africa_shape.boundary.bounds

        n = load_network_for_plots(
            snakemake.input.network,
            snakemake.input.tech_costs,
            snakemake.params.costs,
            snakemake.params.electricity,
        )

        scenario_opts = snakemake.wildcards.opts.split("-")

        fig, ax = plt.subplots(
            figsize=map_figsize, subplot_kw={"projection": ccrs.PlateCarree()}
        )
        plot_map(n, ax, snakemake.wildcards.attr, opts)

        fig.savefig(snakemake.output.only_map, dpi=150, bbox_inches="tight")

        ax1 = fig.add_axes([-0.115, 0.625, 0.2, 0.2])
        plot_total_energy_pie(n, ax1)

        ax2 = fig.add_axes([-0.075, 0.1, 0.1, 0.45])
        plot_total_cost_bar(n, ax2)

        ll = snakemake.wildcards.ll
        ll_type = ll[0]
        ll_factor = ll[1:]
        lbl = dict(c="line cost", v="line volume")[ll_type]
        amnt = (
            "{ll} x today's".format(ll=ll_factor) if ll_factor != "opt" else "optimal"
        )
        fig.suptitle(
            "Expansion to {amount} {label} at {clusters} clusters".format(
                amount=amnt, label=lbl, clusters=snakemake.wildcards.clusters
            )
        )

        fig.savefig(snakemake.output.ext, transparent=True, bbox_inches="tight")

    if snakemake.rule == "plot_sector_network":

        n = pypsa.Network(snakemake.input.network)

        tech_colors = snakemake.config["plotting"]["tech_colors"]
        plot_sector_map(n, transmission=False)
        plot_transmission_topology(n)
        if snakemake.config["sector"]["SMR"]:
            plot_smr(n)
        plot_h2_infra(n)
