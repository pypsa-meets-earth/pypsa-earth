# -*- coding: utf-8 -*-
"""Build EGS potentials for different depths and temperatures. Upper limits available
heat content by overlaying subsurface potential with geospatial demand data."""


import logging

logger = logging.getLogger(__name__)

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from _helpers import configure_logging
from rasterio.transform import xy
from shapely.geometry import Point
from tqdm import tqdm


name_transformer = lambda x: Path(str(x).split('/')[-1]).name.replace('.tif', '')

def tif_to_gdf(tif_files, name_transformer=name_transformer):
    """
    Convert a list of .tif files into a GeoDataFrame where each pixel is represented as a square polygon.
    
    Parameters:
    -----------
    tif_files : list
        List of paths to .tif files
    
    Returns:
    --------
    geopandas.GeoDataFrame
        GeoDataFrame with pixel values and geometries
    """
    
    # Dictionary to store values for each geometry
    geometry_data = {}
    
    for tif_file in tqdm(tif_files, desc="Joining tif files to GeoDataFrame"):
        with rasterio.open(tif_file) as src:
            data = src.read(1)  # Read the first band
            transform = src.transform
            
            # Get the file name without path and without .tif extension
            source_name = name_transformer(tif_file)
            
            # Iterate through each pixel
            for row in range(data.shape[0]):
                for col in range(data.shape[1]):
                    value = data[row, col]
                    
                    # Skip nodata values
                    if value == src.nodata or (src.nodata is None and value == 0):
                        continue
                    
                    # Get the pixel's coordinates
                    x1, y1 = transform * (col, row)
                    x2, y2 = transform * (col + 1, row + 1)
                    
                    # Create a box geometry for the pixel
                    geometry = box(min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))
                    
                    # Convert geometry to WKT for dictionary key
                    geom_wkt = geometry.wkt
                    
                    # Initialize entry if this geometry hasn't been seen before
                    if geom_wkt not in geometry_data:
                        geometry_data[geom_wkt] = {'geometry': geometry}
                    
                    # Add the value for this source
                    geometry_data[geom_wkt][source_name] = value
    
    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(geometry_data, orient='index')
    
    # Extract geometry column
    geometries = df.pop('geometry')
    
    # Create the GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry=geometries, crs=src.crs)
    gdf.index = range(len(gdf))
    
    return gdf


def myround(x):
    if x % 10 == 0:
        return int(x)
    elif x % 10 < 5:
        return int(x - x % 10)
    else:
        return int(x + 10 - x % 10)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_egs_potentials",
            simpl="",
            clusters=50,
        )

    configure_logging(snakemake)

    regions = gpd.read_file(snakemake.input.shapes).set_index("name")

    capex_gdf = (
        get_raster_file(snakemake.input.egs_capex)
        .set_index(["x", "y"])
        .rename(columns={"value": "capex"})
    )
    gen_gdf = (
        get_raster_file(snakemake.input.egs_gen)
        .set_index(["x", "y"])
        .rename(columns={"value": "gen"})
    )
    opex_gdf = (
        get_raster_file(snakemake.input.egs_opex)
        .set_index(["x", "y"])
        .rename(columns={"value": "opex"})
    )

    gdf = capex_gdf[["capex", "geometry"]].join(gen_gdf[["gen"]], how="inner")
    gdf = gdf.join(opex_gdf[["opex"]], how="inner")

    gdf["plant_capacity"] = gdf["gen"].div(8760)
    gdf["capex($/kWel)"] = gdf["capex"].div(gdf["plant_capacity"]).mul(1e3)
    gdf["plant_capacity(MWel)"] = gdf["plant_capacity"]
    gdf["opex[$/kWh]"] = gdf["opex"].mul(1e-2)

    gdf.rename(
        columns={
            "capex($/kWel)": "capex[$/kW]",
            "plant_capacity(MWel)": "available_capacity[MW]",
        },
        inplace=True,
    )

    gdf = gdf[["capex[$/kW]", "opex[$/kWh]", "available_capacity[MW]", "geometry"]]

    gdf = gdf.replace([np.inf, -np.inf], np.nan)

    nodal_egs_potentials = pd.DataFrame(
        np.nan,
        index=regions.index,
        columns=["capex[$/kW]", "opex[$/kWh]", "p_nom_max[MW]"],
    )

    config = snakemake.params["enhanced_geothermal"]

    regional_potentials = []

    for name, geom in tqdm(
        regions.geometry.items(),
        desc="Matching EGS potentials to network regions",
        ascii=True,
    ):

        ss = gdf.loc[gdf.geometry.within(geom)]
        if ss.empty:
            continue

        ss = (
            ss[["capex[$/kW]", "opex[$/kWh]", "available_capacity[MW]"]]
            .reset_index(drop=True)
            .sort_values(by="capex[$/kW]")
        )

        ss["agg_available_capacity[MW]"] = ss["available_capacity[MW]"].cumsum()

        bins = pd.Series(
            np.linspace(
                # ss["capex[$/kW]"].min(),
                0,
                ss["capex[$/kW]"].max(),
                config["max_levels"] + 1,
            )
        )

        labels = bins.rolling(2).mean().dropna().tolist()

        if ss.empty:
            continue

        if len(ss) > 1:
            ss["level"] = pd.cut(
                ss["capex[$/kW]"], bins=bins, labels=labels, duplicates="drop"
            )
        else:
            ss["level"] = ss["capex[$/kW]"].iloc[0]

        ss = ss.dropna()

        ss = ss.groupby("level", observed=False)[
            ["available_capacity[MW]", "opex[$/kWh]"]
        ].agg({"available_capacity[MW]": "sum", "opex[$/kWh]": "mean"})
        ss.index = list(map(myround, ss.index))
        ss.index = pd.MultiIndex.from_product(
            [[name], ss.index], names=["network_region", "capex[$/kW]"]
        )

        regional_potentials.append(ss)

    regional_potentials = pd.concat(regional_potentials).dropna()
    dr = snakemake.params["costs"]["fill_values"]["discount rate"]

    def to_capital_cost(index, lt=25):
        return (index[0], index[1] * dr / (1 - (1 + dr) ** (-lt)))

    regional_potentials.index = regional_potentials.index.map(to_capital_cost)
    regional_potentials.index.names = ["network_region", "capital_cost[$/kW]"]

    regional_potentials.to_csv(snakemake.output.egs_potentials)
