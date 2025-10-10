# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Plot water network for Desalination .
"""

import logging

logger = logging.getLogger(__name__)

import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely import wkt
import numpy as np
from pathlib import Path
from _helpers import (
    BASE_DIR,
    content_retrieve,
    progress_retrieve,
    two_2_three_digits_country,
)



def plot_water_networks (regions_onshore_aqueduct, regions_onshore_aqueduct_desalination, clustered_water_network, shorelines_natura):
    """
    Plots the water network visualization including altitude changes along pipelines.

    Parameters:
    - regions_onshore_aqueduct (GeoDataFrame): GeoDataFrame containing onshore aqueduct regions with polygon geometries.
    - regions_onshore_aqueduct_desalination (GeoDataFrame): GeoDataFrame containing desalination regions with centroid geometries.
    - clustered_water_network (GeoDataFrame): GeoDataFrame representing the clustered water network with altitude change data.
    - shorelines_natura (GeoDataFrame): GeoDataFrame containing shoreline geometries.

    The function creates a plot with the following elements:
    - Clipped shorelines in blue.
    - Original aqueduct polygons with specified colors and transparency.
    - Desalination centroids in red.
    - Clustered water network lines colored by total positive altitude change.

    The plot is saved to a file specified by `snakemake.output.water_network`.
    """

    # Plot the lines with altitude changes
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    shorelines_natura.to_crs(epsg=4326).plot(ax=ax, color='blue', label='Clipped Shoreline')
    regions_onshore_aqueduct.to_crs(epsg=4326).plot(ax=ax, edgecolor='black', color=regions_onshore_aqueduct['color'], alpha=0.4, label="Original Polygons")

#-----------
    # regions_onshore_aqueduct.to_crs(epsg=4326).plot(
    #     ax=ax, 
    #     edgecolor='black', 
    #     color=regions_onshore_aqueduct['color'], 
    #     alpha=0.4, 
    #     label="Original Polygons"
    # )

    # # Add labels for each polygon using the 'name' column
    # for idx, row in regions_onshore_aqueduct.iterrows():
    #     ax.annotate(
    #         text=row['name'], 
    #         xy=(row.geometry.centroid.x, row.geometry.centroid.y),
    #         xytext=(-9, -9), 
    #         textcoords="offset points",
    #         fontsize=8,
    #         color='black',
    #         bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white", alpha=0.7),
    #         arrowprops=dict(arrowstyle="->", color='black', lw=0.5)
    #     )
#----------
    regions_onshore_aqueduct_desalination.plot(ax=ax, color='red', marker='o', label='Centroids')

    clustered_water_network.to_crs(epsg=4326).plot(
        column='total_positive_altitude_change',
        cmap='viridis',
        linewidth=2,
        ax=ax,
        legend=True
    )
    plt.legend(loc='upper left')
    plt.title("Total positive altitude change along water pipelines", fontsize=16)
    
    # Save the figure to a file
    plt.savefig(snakemake.output.water_network, dpi=300, bbox_inches='tight')  # Save with high resolution and tight layout
    # plt.show()



def plot_altitude_profile(df, resolution=1000):
    """
    Makes a simple plot showing:
      (1) elevation along line
      (2) running head and per-step added head (bars)
    """
    k = 1.5          # scaling factor for distance matching prepare_water_network
    inv = 1.0 / k    # scaling factor for effective_increment_m and head matching prepare_water_network
    factors = {
        "distance_m": k,
        # "elevation_m": inv,
        "head_clamped_m": inv,
        "effective_increment_m": inv,
    }

    # scale only the columns that exist (avoids KeyError)
    cols = [c for c in factors if c in df.columns]
    df.loc[:, cols] = df.loc[:, cols].mul(pd.Series({c: factors[c] for c in cols}), axis=1)

    # Running head + added head per step
    plt.figure(figsize=(10, 4.2))
    plt.plot(df["distance_m"], df["head_clamped_m"], label="Running head (clamped)")
    plt.plot(df["distance_m"], df["elevation_m"], label="Elevation", linestyle=":")
    # plt.axhline(df["max_altitude_cap"].unique(), linestyle="--", label="Head cap")
    plt.axhline(df["min_required_head"].unique(), linestyle=":", label="Min required (end-start)")

    # bar of what was added each step; align at segment ends
    if "effective_increment_m" in df:
        w = resolution * 0.9
        plt.bar(df["distance_m"].iloc[1:], df["effective_increment_m"].iloc[1:].fillna(0.0),
                width=w, align="edge", alpha=0.4, label="Added this step")
        
        # --- show sum of effective_increment_m on the plot ---
        total_added = float(np.nansum(df["effective_increment_m"].to_numpy()))
        plt.gca().text(
            0.99, 0.05, f"Σ added head = {total_added:.1f} m",
            transform=plt.gca().transAxes, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7, linewidth=0)
        )
        
    plt.ylabel("Head [m]")
    plt.xlabel("Distance along line [m]")
    plt.grid(True, alpha=0.3)
    plt.title(df["line_name"].unique())
    plt.legend()
    plt.tight_layout()

    # Save the figure to a file
    path_to_save = Path(os.path.join(BASE_DIR, "results/plots/desalination/{}_profile.png".format(df["line_name"].unique()[0])))
    plt.savefig(path_to_save, dpi=300, bbox_inches='tight')  # Save with high resolution and tight layout
    # plt.show()



if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_water_network",
            simpl="",
            clusters="31",
            ext="png",
        )


    # Load input data
    regions_onshore_aqueduct = gpd.read_file(snakemake.input.regions_onshore_aqueduct)
    regions_onshore_aqueduct_desalination = gpd.read_file(snakemake.input.regions_onshore_aqueduct_desalination)
    clustered_water_network = gpd.read_file(snakemake.input.clustered_water_network)
    shorelines_natura = gpd.read_file(snakemake.input.shorelines_natura)
    water_pipes_profiles = pd.read_csv(snakemake.input.water_pipes_profiles, index_col=0)


    # Convert the 'centroid' column from WKT to geometry
    regions_onshore_aqueduct_desalination['centroid'] = regions_onshore_aqueduct_desalination['centroid'].apply(wkt.loads)
    regions_onshore_aqueduct_desalination = regions_onshore_aqueduct_desalination.set_geometry('centroid', crs="EPSG:4326")



    # Plot the water networks
    logger.info("Plotting water networks...")
    plot_water_networks (regions_onshore_aqueduct, regions_onshore_aqueduct_desalination, clustered_water_network, shorelines_natura)

    # Plot altitude profiles for each water pipe
    logger.info("Plotting altitude profiles for water pipes...")
    for line in water_pipes_profiles.line_name.unique():
        to_plot = water_pipes_profiles[water_pipes_profiles.line_name == line]
        plot_altitude_profile(to_plot, resolution=1000)