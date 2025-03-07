# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepare gas network.
"""

import logging

logger = logging.getLogger(__name__)

import geopandas as gpd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
from _helpers import (
    three_2_two_digits_country,
)
from matplotlib.lines import Line2D
from shapely.geometry import LineString
from shapely import wkt


def plot_gas_network(pipelines, country_borders, bus_regions_onshore, gas_network_data):
    df = pipelines.copy()
    df = gpd.overlay(df, country_borders, how="intersection")

    if gas_network_data == "IGGIELGN":
        df = df.rename({"p_nom": "capacity [MW]"}, axis=1)

    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(12, 7)
    bus_regions_onshore.to_crs(epsg=3857).plot(
        ax=ax, color="white", edgecolor="darkgrey", linewidth=0.5
    )
    df.loc[(df.amount_states_passed > 1)].to_crs(epsg=3857).plot(
        ax=ax,
        column="capacity [MW]",
        linewidth=2.5,
        # linewidth=df['capacity [MW]'],
        # alpha=0.8,
        categorical=False,
        cmap="viridis_r",
        # legend=True,
        # legend_kwds={'label':'Pipeline capacity [MW]'},
    )

    df.loc[(df.amount_states_passed <= 1)].to_crs(epsg=3857).plot(
        ax=ax,
        column="capacity [MW]",
        linewidth=2.5,
        # linewidth=df['capacity [MW]'],
        alpha=0.5,
        categorical=False,
        # color='darkgrey',
        ls="dotted",
    )

    # # Create custom legend handles for line types
    # line_types = [ 'solid', 'dashed', 'dotted'] # solid
    # legend_handles = [Line2D([0], [0], color='black', linestyle=line_type) for line_type in line_types]

    # Define line types and labels
    line_types = ["solid", "dotted"]
    line_labels = ["Operating", "Not considered \n(within-state)"]

    # Create custom legend handles for line types
    legend_handles = [
        Line2D([0], [0], color="black", linestyle=line_type, label=line_label)
        for line_type, line_label in zip(line_types, line_labels)
    ]

    # Add the line type legend
    ax.legend(
        handles=legend_handles,
        title="Status",
        borderpad=1,
        title_fontproperties={"weight": "bold"},
        fontsize=11,
        loc=1,
    )

    # # create the colorbar
    norm = colors.Normalize(
        vmin=df["capacity [MW]"].min(), vmax=df["capacity [MW]"].max()
    )
    cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis_r")
    # fig.colorbar(cbar, ax=ax).set_label('Capacity [MW]')

    # add colorbar
    ax_cbar = fig.colorbar(cbar, ax=ax, location="left", shrink=0.8, pad=0.01)
    # add label for the colorbar
    ax_cbar.set_label("Natural gas pipeline capacity [MW]", fontsize=15)

    ax.set_axis_off()
    fig.savefig(snakemake.output.existing_gas_pipelines, dpi=300, bbox_inches="tight")

def plot_clustered_gas_network(pipelines, bus_regions_onshore):
    # Create a new GeoDataFrame with centroids
    centroids = bus_regions_onshore.copy()
    centroids["geometry"] = centroids["geometry"].centroid
    centroids["gadm_id"] = centroids["gadm_id"].apply(
        lambda id: three_2_two_digits_country(id[:3]) + id[3:]
    )
    gdf1 = pd.merge(
        pipelines, centroids, left_on=["bus0"], right_on=["gadm_id"], how="left"
    )
    gdf1.rename(columns={"geometry": "geometry_bus0"}, inplace=True)
    pipe_links = pd.merge(
        gdf1, centroids, left_on=["bus1"], right_on=["gadm_id"], how="left"
    )
    pipe_links.rename(columns={"geometry": "geometry_bus1"}, inplace=True)

    # Create LineString geometries from the points
    pipe_links["geometry"] = pipe_links.apply(
        lambda row: LineString([row["geometry_bus0"], row["geometry_bus1"]]), axis=1
    )

    clustered = gpd.GeoDataFrame(pipe_links, geometry=pipe_links["geometry"])

    # Optional: Set the coordinate reference system (CRS) if needed
    clustered.crs = "EPSG:3857"  # For example, WGS84

    # plot pipelines
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(12, 7)
    bus_regions_onshore.to_crs(epsg=3857).plot(
        ax=ax, color="white", edgecolor="darkgrey", linewidth=0.5
    )
    clustered.to_crs(epsg=3857).plot(
        ax=ax,
        column="capacity",
        linewidth=1.5,
        categorical=False,
        cmap="plasma",
    )

    centroids.to_crs(epsg=3857).plot(
        ax=ax,
        color="red",
        markersize=35,
        alpha=0.5,
    )

    norm = colors.Normalize(
        vmin=pipelines["capacity"].min(), vmax=pipelines["capacity"].max()
    )
    cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma")

    # add colorbar
    ax_cbar = fig.colorbar(cbar, ax=ax, location="left", shrink=0.8, pad=0.01)
    # add label for the colorbar
    ax_cbar.set_label("Natural gas pipeline capacity [MW]", fontsize=15)

    ax.set_axis_off()
    fig.savefig(snakemake.output.clustered_gas_pipelines, dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_gas_network",
            simpl="",
            clusters="4",
        )
    # configure_logging(snakemake)

    gas_network_data = snakemake.params.gas_network_data

    bus_regions_onshore = gpd.read_file(snakemake.input.bus_regions_onshore)
    country_borders = gpd.read_file(snakemake.input.country_borders)

    pipelines = pd.read_csv(snakemake.input.preclustered_gas_network)
    pipelines['geometry'] = pipelines['geometry'].apply(wkt.loads)
    pipelines = gpd.GeoDataFrame(pipelines, geometry = "geometry", crs="epsg:3857")

    clustered_pipelines = pd.read_csv(snakemake.input.clustered_gas_network)

    plot_gas_network(pipelines, country_borders, bus_regions_onshore, gas_network_data)
    plot_clustered_gas_network(clustered_pipelines, bus_regions_onshore)