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




if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "plot_water_network",
            simpl="",
            clusters="16",
            ext="png",
        )


    # Load input data
    regions_onshore_aqueduct = gpd.read_file(snakemake.input.regions_onshore_aqueduct)
    regions_onshore_aqueduct_desalination = gpd.read_file(snakemake.input.regions_onshore_aqueduct_desalination)
    clustered_water_network = gpd.read_file(snakemake.input.clustered_water_network)
    shorelines_natura = gpd.read_file(snakemake.input.shorelines_natura)


    # Convert the 'centroid' column from WKT to geometry
    regions_onshore_aqueduct_desalination['centroid'] = regions_onshore_aqueduct_desalination['centroid'].apply(wkt.loads)
    regions_onshore_aqueduct_desalination = regions_onshore_aqueduct_desalination.set_geometry('centroid', crs="EPSG:4326")




    plot_water_networks (regions_onshore_aqueduct, regions_onshore_aqueduct_desalination, clustered_water_network, shorelines_natura)
