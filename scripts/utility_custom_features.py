# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: Contributors to PyPSA-Earth
# SPDX-FileCopyrightText: Open Energy Transition gGmbH
#
# SPDX-License-Identifier: AGPL-3.0-or-later


import logging

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from pypsa.io import import_components_from_dataframe, import_series_from_dataframe
from rasterio.features import rasterize
from rasterio.mask import mask
from rasterio.transform import from_bounds
from shapely import wkt

logger = logging.getLogger(__name__)


def annual_gwh_to_average_mw(energy_gwh, hours_per_year=8760):
    """Convert annual energy in GWh to average power in MW."""
    return energy_gwh * 1000 / hours_per_year


def load_interconnector_data(countries_path, links_path, substations_path):
    """Load interconnector input data from CSV files."""
    return (
        pd.read_csv(countries_path),
        pd.read_csv(links_path),
        pd.read_csv(substations_path),
    )


def find_nearest_bus(n, lat, lon, distance_crs="EPSG:20935"):
    """Return the nearest Zambian bus to a given latitude and longitude."""
    buses = n.buses[n.buses["country"] == "ZM"].copy()
    buses = gpd.GeoDataFrame(
        buses,
        geometry=gpd.points_from_xy(buses["x"], buses["y"]),
        crs="EPSG:4326",
    ).to_crs(distance_crs)
    target_point = (
        gpd.GeoSeries.from_xy([lon], [lat], crs="EPSG:4326")
        .to_crs(distance_crs)
        .iloc[0]
    )
    distances = buses.geometry.distance(target_point)
    return distances.idxmin()


def add_foreign_buses(n, power_pool_countries):
    """Add neighbouring-country buses to the network."""
    for _, row in power_pool_countries.iterrows():
        country = row["country"]
        if country not in n.buses.index:
            n.add("Bus", country, x=row["lon"], y=row["lat"], carrier="AC")
            n.buses.loc[country, "country"] = country
    return n


def add_cross_border_links(n, power_pool_links, substation_dict, distance_crs):
    """Add cross-border links to the network."""
    for _, row in power_pool_links.iterrows():
        name = row["name"]
        if row["from_country"] == "ZM":
            lat, lon = substation_dict[name]
            bus0 = find_nearest_bus(n, lat, lon, distance_crs)
        else:
            bus0 = row["from_country"]
        if row["to_country"] == "ZM":
            lat, lon = substation_dict[name]
            bus1 = find_nearest_bus(n, lat, lon, distance_crs)
        else:
            bus1 = row["to_country"]

        if name not in n.links.index:
            n.add(
                "Link",
                name,
                bus0=bus0,
                bus1=bus1,
                carrier="AC",
                p_nom=row["capacity_mw"],
                efficiency=1.0,
                p_min_pu=-1.0,
            )
    return n


def add_trade_components(n, power_pool_countries, hours_per_year=8760):
    """Add import loads and export generators for neighbouring countries."""
    for _, row in power_pool_countries.iterrows():
        country = row["country"]
        if country not in n.buses.index:
            continue

        load_name = f"import_{country}"
        gen_name = f"export_{country}"

        if load_name not in n.loads.index:
            n.add(
                "Load",
                load_name,
                bus=country,
                carrier="import",
                p_set=annual_gwh_to_average_mw(row["demand_gwh"], hours_per_year),
            )

        if gen_name not in n.generators.index:
            n.add(
                "Generator",
                gen_name,
                bus=country,
                carrier="export",
                p_nom=annual_gwh_to_average_mw(row["generation_gwh"], hours_per_year),
                marginal_cost=row["marginal_cost"],
            )
    return n


def add_interconnectors(
    n,
    power_pool_countries,
    power_pool_links,
    substations,
    distance_crs,
    hours_per_year=8760,
):
    """Add foreign buses, interconnectors, and trade components to the network."""
    substation_dict = {
        row["name"]: (row["lat"], row["lon"]) for _, row in substations.iterrows()
    }

    n = add_foreign_buses(n, power_pool_countries)
    n = add_cross_border_links(n, power_pool_links, substation_dict, distance_crs)
    n = add_trade_components(n, power_pool_countries, hours_per_year)

    return n


def load_custom_line_types(line_types: str) -> pd.DataFrame:
    """Load and format custom transmission line types for a PyPSA network."""
    custom_line_types = pd.read_csv(line_types)

    custom_line_types = custom_line_types.rename(
        columns={
            "PyPSA Type Name": "name",
            "f_nom (Hz)": "f_nom",
            "r_per_length (Ω/km)": "r_per_length",
            "x_per_length (Ω/km)": "x_per_length",
            "c_per_length (nF/km)": "c_per_length",
            "i_nom (kA)": "i_nom",
            "cross_section (mm²)": "cross_section",
        }
    )
    custom_line_types = custom_line_types.set_index("name")
    return custom_line_types


def add_custom_line_types(n, custom_line_types):
    """merge custom line_types into the pypsa network"""
    n.line_types = pd.concat([n.line_types, custom_line_types], axis=0)
    n.line_types = n.line_types[~n.line_types.index.duplicated(keep="last")]
    return n


def map_buses_from_coords(
    n, df, buses_df=None, distance_crs="EPSG:20935", geo_crs="EPSG:4326"
):
    """Find the nearest bus for each row in df using lat/lon.

    Parameters
    ----------
    n : pypsa.Network
        The network whose buses are used for geo-matching.
    df : pd.DataFrame
        Power plant table. Must contain ``lat`` and ``lon`` columns.
    buses_df : pd.DataFrame or None, optional
        Candidate buses for geo-matching. Defaults to all buses in ``n``.
    distance_crs : str, optional
        Geographic CRS used to measure distances.
    geo_crs : str, optional
        Geographic CRS used to construct GeoDataFrames before reprojection.
    """
    candidates = n.buses if buses_df is None else buses_df
    bus_points = gpd.GeoDataFrame(
        index=candidates.index,
        geometry=gpd.points_from_xy(candidates["x"], candidates["y"]),
        crs=geo_crs,
    ).to_crs(distance_crs)
    plant_points = gpd.GeoDataFrame(
        index=df.index,
        geometry=gpd.points_from_xy(df["lon"], df["lat"]),
        crs=geo_crs,
    ).to_crs(distance_crs)

    return plant_points.geometry.apply(
        lambda plant: bus_points.geometry.distance(plant).idxmin()
    )


def disaggregate_plants(
    n, df, name_fallback="plant", buses_df=None, geo_crs="EPSG:4326"
):
    """Rename df rows to real plant names and assign each to its nearest bus.

    Parameters
    ----------
    n : pypsa.Network
        The network whose buses are used for geo-matching.
    df : pd.DataFrame
        Power plant table. Must contain ``lat`` and ``lon`` columns.
        If a ``name`` column is present, its values are used as the plant name
        prefix in the new index; otherwise ``name_fallback`` is used.
    name_fallback : str, optional
        Prefix used in the new index when a plant has no name (default: "plant").
        E.g. ``name_fallback="hydro"`` produces index entries like ``"hydro-3"``.
    buses_df : pd.DataFrame or None, optional
        Candidate buses for geo-matching. Defaults to all buses in ``n``.
    geo_crs : str, optional
        Geographic CRS used to construct GeoDataFrames before reprojection
        (default: ``"EPSG:4326"``). Should match ``config["crs"]["geo_crs"]``.
    """
    new_names = []
    for i in df.index:
        if "name" in df.columns:
            plant_name = df.loc[i, "name"]
            if plant_name is not None and str(plant_name).strip() != "":
                new_name = str(plant_name) + "-" + str(i)
            else:
                new_name = name_fallback + "-" + str(i)
        else:
            new_name = name_fallback + "-" + str(i)
        new_names.append(new_name)
    df.index = new_names
    df["bus"] = map_buses_from_coords(n, df, buses_df=buses_df, geo_crs=geo_crs)
    return df


def save_excluded_components(n, component, busmap, exclude_carriers):
    """
    Save components we want to protect from aggregation.

    Works for any component type: "Generator", "Load", "StorageUnit", etc.
    Pulls out matching rows, remaps their buses, and saves their time-series.

    Parameters
    ----------
    n : pypsa.Network
        The network containing the components to save.
    component : str
        Component type to filter, e.g. ``"Generator"`` or ``"StorageUnit"``.
    busmap : pd.Series
        Maps old bus names to new bus names (as produced by clustering).
    exclude_carriers : list of str
        Carriers whose components should be saved and excluded from aggregation.

    Returns
    -------
    saved : pd.DataFrame
        Static component table for the excluded components, with buses remapped.
    saved_timeseries : dict of str -> pd.DataFrame
        Time-varying attributes (e.g. ``p_max_pu``) for the excluded components.
        Empty dict if none exist.
    """
    if not exclude_carriers:
        return pd.DataFrame(), {}
    component_table = n.df(component)
    if "carrier" not in component_table.columns:
        return pd.DataFrame(), {}
    is_excluded = component_table.carrier.isin(exclude_carriers)
    if not is_excluded.any():
        return pd.DataFrame(), {}
    saved = component_table[is_excluded].copy()
    saved["bus"] = saved["bus"].replace(busmap)
    list_name = n.components[component]["list_name"]
    ts_container = getattr(n, list_name + "_t")
    saved_timeseries = {}
    for attr in ts_container:
        ts_data = getattr(ts_container, attr)
        if ts_data.empty:
            continue
        cols = ts_data.columns.intersection(saved.index)
        if not cols.empty:
            saved_timeseries[attr] = ts_data[cols]
    return saved, saved_timeseries


def restore_excluded_components(n, component, saved_components, saved_timeseries):
    """
    Put protected components back into the network after aggregation
    """
    if saved_components.empty:
        return
    import_components_from_dataframe(n, saved_components, component)
    for attr, data in saved_timeseries.items():
        if not data.empty:
            import_series_from_dataframe(n, data, component, attr)


def busmap_keeps_topology(busmap):
    """
    Check if every bus maps to itself.
    """
    return (busmap.index == busmap.values).all()


def add_mining_data(df_gadm, mining_raster_path):
    """
    Add mining electricity demand (GWh/year) to each province.

    Uses a raster with intensity (MWh/km²/year), clips it to each province,
    converts to total demand, and stores results in df_gadm["mining"].
    """
    df_gadm["mining"] = 0.0
    with rasterio.open(mining_raster_path) as src:
        pixel_size_deg = abs(src.transform.a)
        for idx, row in df_gadm.iterrows():
            geom = row.geometry
            if geom is None or geom.is_empty:
                continue
            # Clip raster to province
            clipped, _ = mask(src, [geom], all_touched=True, nodata=0.0)
            intensity = clipped[0]  # MWh/km²/year
            if intensity.max() == 0:
                continue
            # Approximate pixel area (km²)
            lat = geom.centroid.y
            pixel_area_km2 = (pixel_size_deg * 111.0) ** 2 * np.cos(np.radians(lat))
            # Total demand (GWh/year)
            df_gadm.loc[idx, "mining"] = (intensity * pixel_area_km2).sum() / 1000.0
    return df_gadm


def load_mining_data(
    provincial_demand_path, mining_polygons_path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load mining raster input data from CSV files

    Returns:
    --------
    tuple
        Tuple of dataframes containing mining demand and mining polygons
    """
    return (pd.read_csv(provincial_demand_path), pd.read_csv(mining_polygons_path))


def build_mining_raster(
    provincial_demand,
    mining_polygons,
    output_path,
    resolution=1000,
    geo_crs="EPSG:4326",
    area_crs="ESRI:54009",
):
    """
    Create a mining demand raster for Zambia.

    Arguments
    ---------
    provincial_demand: pd.DataFrame
        columns [province, mining_demand_gwh]
    mining_polygons: pd.DataFrame
        columns [province, area_km2, geometry_wkt].
        Each polygon carries a native "province" field, so demand is assigned by
        direct attribute lookup — no spatial intersection is required.
    output_path: str
        path for the output GeoTIFF file
    resolution: int
        raster resolution in units of area_crs (default: 1000 m)
    geo_crs: str
        CRS of the input WKT geometries (default: EPSG:4326)
    area_crs: str
        CRS used for area calculations and raster output (default: ESRI:54009)

    Returns
    -------
    str
        Path to the saved mining demand raster file

    Notes
    -----
    Output values are in MWh/km²/year.
    """

    demand = provincial_demand.set_index("province")
    mines = mining_polygons
    # Total mining area per province
    mines["province_area_km2"] = mines.groupby("province")["area_km2"].transform("sum")
    # Demand intensity
    mines["demand_mwh_per_km2"] = (
        mines["province"].map(demand["mining_demand_gwh"])
        * 1000
        / mines["province_area_km2"]
    )
    # Convert WKT to geometries, then reproject to area_crs for rasterization
    gdf = gpd.GeoDataFrame(
        mines,
        geometry=mines["geometry_wkt"].apply(wkt.loads),
        crs=geo_crs,
    ).to_crs(area_crs)
    # Derive bounding box from reprojected data
    x_min, y_min, x_max, y_max = gdf.total_bounds
    width = round((x_max - x_min) / resolution)
    height = round((y_max - y_min) / resolution)
    transform = from_bounds(x_min, y_min, x_max, y_max, width, height)
    # Rasterize mining demand
    shapes = zip(gdf.geometry, gdf["demand_mwh_per_km2"])

    raster = rasterize(
        shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype="float32",
    )
    # Save raster
    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs=area_crs,
        transform=transform,
        nodata=0,
    ) as dst:
        dst.write(raster, 1)

    return output_path
