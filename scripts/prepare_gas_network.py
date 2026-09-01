# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepare gas network.


Relevant Settings
-----------------

```yaml

    sector:
        gas:
            spatial_gas:
            network:
            network_data:
            network_data_GGIT_status:

    clustering:
        alternative_clustering:

    custom_data:
        gas_network:

```

Inputs
------
- ``resources/{RDIR}/bus_regions/regions_onshore_elec_s{simpl}_{clusters}.geojson``: path to onshore region polygons used to assign pipelines to regions.

Outputs
-------
- ``resources/{SECDIR}/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv``: CSV containing aggregated pipeline capacities and lengths between onshore regions.

Description
-----------

Utilities to download, load and cluster global natural gas pipeline datasets
for use with PyPSA-Earth. This module supports two source datasets:
`GGIT` and `IGGIELGN`. It provides helpers to normalise pipeline diameters and
capacities, clean geometry, assign pipelines to onshore regions, and aggregate
inter-state pipeline capacities.
"""

import logging

logger = logging.getLogger(__name__)

import os
import zipfile
from pathlib import Path
from typing import Any, List, Union

import fiona
import geopandas as gpd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
from _helpers import (
    BASE_DIR,
    content_retrieve,
    progress_retrieve,
    three_2_two_digits_country,
    two_2_three_digits_country,
)
from build_shapes import gadm
from matplotlib.lines import Line2D
from pyproj import CRS
from pypsa.geo import haversine_pts
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import unary_union, linemerge
from shapely.validation import make_valid

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_gas_network",
            simpl="",
            clusters="10",
            H2_retrofit_capacity_per_CH4= "10"
        )

    # configure_logging(snakemake)

    # run = snakemake.config.get("run", {})
    # RDIR = run["name"] + "/" if run.get("name") else ""
    # store_path_data = Path.joinpath(Path().cwd(), "data")
    # country_list = country_list_to_geofk(snakemake.config["countries"])'


def download_IGGIELGN_gas_network() -> None:
    """
    Downloads a global dataset for gas networks as .xlsx.

    The following xlsx file was downloaded from the webpage
    https://globalenergymonitor.org/projects/global-gas-infrastructure-tracker/
    The dataset contains 3144 pipelines.
    """

    url = "https://zenodo.org/record/4767098/files/IGGIELGN.zip"

    # Save locations
    zip_fn = Path(os.path.join(BASE_DIR, "IGGIELGN.zip"))
    to_fn = Path(os.path.join(BASE_DIR, "data/gas_network/scigrid-gas"))

    logger.info(f"Downloading databundle from '{url}'.")
    progress_retrieve(url, zip_fn)

    logger.info(f"Extracting databundle.")
    zipfile.ZipFile(zip_fn).extractall(to_fn)

    zip_fn.unlink()

    logger.info(f"Gas infrastructure data available in '{to_fn}'.")


def download_GGIT_gas_network() -> pd.DataFrame:
    """
    Downloads a global dataset for gas networks as .xlsx.

    The following xlsx file was downloaded from the webpage
    https://globalenergymonitor.org/projects/global-gas-infrastructure-tracker/
    The dataset contains 3144 pipelines.
    """
    url = "https://github.com/pypsa-meets-earth/temporary_storage/raw/refs/heads/main/datasets/GEM-GGIT-Gas-Pipelines-December-2022.xlsx"
    
    GGIT_gas_pipeline = pd.read_excel(
        content_retrieve(url),
        index_col=0,
        sheet_name="Gas Pipelines 2022-12-16",
        header=0,
    )
    

    return GGIT_gas_pipeline


def diameter_to_capacity(pipe_diameter_mm: int) -> int:
    """
    Calculate pipe capacity in MW based on diameter in mm.

    20 inch (500 mm)  50 bar -> 1.5   GW CH4 pipe capacity (LHV) 24 inch
    (600 mm)  50 bar -> 5     GW CH4 pipe capacity (LHV) 36 inch (900
    mm)  50 bar -> 11.25 GW CH4 pipe capacity (LHV) 48 inch (1200 mm) 80
    bar -> 21.7  GW CH4 pipe capacity (LHV)

    Based on p.15 of
    https://gasforclimate2050.eu/wp-content/uploads/2020/07/2020_European-Hydrogen-Backbone_Report.pdf

    Parameters
    ----------
    pipe_diameter_mm: int
        Diameter of gas pipeline in mm

    Returns
    -------
    int
        Pipe capacity in MW
    """
    # slopes definitions
    m0 = (1500 - 0) / (500 - 0)
    m1 = (5000 - 1500) / (600 - 500)
    m2 = (11250 - 5000) / (900 - 600)
    m3 = (21700 - 11250) / (1200 - 900)

    # intercept
    a0 = 0
    a1 = -16000
    a2 = -7500
    a3 = -20100

    if pipe_diameter_mm < 500:
        return a0 + m0 * pipe_diameter_mm
    elif pipe_diameter_mm < 600:
        return a1 + m1 * pipe_diameter_mm
    elif pipe_diameter_mm < 900:
        return a2 + m2 * pipe_diameter_mm
    else:
        return a3 + m3 * pipe_diameter_mm


def inch_to_mm(len_inch: float) -> float:
    """
    Convert a length from inches to millimetres.

    Parameters
    ----------
    len_inch : float
        Length in inches.

    Returns
    -------
    float
        Length in millimetres.
    """
    return len_inch / 0.0393701


def bcm_to_MW(cap_bcm: float) -> float:
    """
    Convert volumetric capacity in bcm/year to power in MW.

    Parameters
    ----------
    cap_bcm : float
        Capacity in billion cubic metres per year.

    Returns
    -------
    float
        Equivalent average power in MW.
    """
    return cap_bcm * 9769444.44 / 8760


def correct_Diameter_col(value: Any) -> float:
    """
    Parse and average compound pipeline diameter values.

    Handles diameter strings containing commas, slashes or dashes by splitting
    the string into numeric parts and returning the mean.
    """
    value = str(value)
    # Check if the value contains a comma
    if "," in value:
        # Split the value by comma and convert each part to a float
        diameter_values = [float(val) for val in value.split(",")]
        # Calculate the mean of the values
        return sum(diameter_values) / len(diameter_values)
    elif "/" in value:
        # Split the value by slash and convert each part to a float
        diameter_values = [float(val) for val in value.split("/")]
        # Calculate the mean of the values
        return sum(diameter_values) / len(diameter_values)
    elif "-" in value:
        # Split the value by slash and convert each part to a float
        diameter_values = [float(val) for val in value.split("-")]
        # Calculate the mean of the values
        return sum(diameter_values) / len(diameter_values)
    else:
        # Return the original value for rows without a comma or slash
        return float(value)


def prepare_GGIT_data(GGIT_gas_pipeline: pd.DataFrame) -> gpd.GeoDataFrame:
    """
    Clean and normalise the GGIT pipeline dataset.

    Parameters
    ----------
    GGIT_gas_pipeline : pandas.DataFrame
        Raw GGIT pipeline table read from the Excel source.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame containing valid geometries, corrected diameter values,
        and capacities expressed in MW.
    """

    df = GGIT_gas_pipeline.copy().reset_index()

    # Drop rows containing "--" in the 'WKTFormat' column
    df = df[df["WKTFormat"] != "--"]

    # Keep pipelines that are as below
    df = df[df["Status"].isin(snakemake.params.gas_config["network_data_GGIT_status"])]

    # Convert the WKT column to a GeoDataFrame
    df = gpd.GeoDataFrame(
        df, geometry=gpd.GeoSeries.from_wkt(df["WKTFormat"], on_invalid="warn")
    )

    df = df[df.geometry.is_valid & ~df.geometry.is_empty]

    # Set the CRS to EPSG:4326
    df.crs = CRS.from_epsg(4326)

    # Convert CRS to EPSG:3857 so we can measure distances
    df = df.to_crs(epsg=3857)

    # Convert and correct diameter column to be in mm
    df.loc[df["DiameterUnits"] == "mm", "diameter_mm"] = df.loc[
        df["DiameterUnits"] == "mm", "Diameter"
    ].apply(correct_Diameter_col)
    df.loc[df["DiameterUnits"] == "in", "diameter_mm"] = (
        df.loc[df["DiameterUnits"] == "in", "Diameter"]
        .apply(correct_Diameter_col)
        .apply(
            lambda d: inch_to_mm(float(d))
        )  #         .apply(lambda ds: pd.Series(ds).apply(lambda d: inch_to_mm(float(d))))
    )

    # Convert Bcm/y to MW
    df["CapacityBcm/y"] = pd.to_numeric(df["CapacityBcm/y"], errors="coerce")
    df["capacity [MW]"] = df["CapacityBcm/y"].apply(lambda d: bcm_to_MW(d))

    # Get capacity from diameter for rows where no capacity is given
    df.loc[df["CapacityBcm/y"] == "--", "capacity [MW]"] = df.loc[
        df["CapacityBcm/y"] == "--", "diameter_mm"
    ].apply(lambda d: diameter_to_capacity(int(d)))
    df["diameter_mm"] = pd.to_numeric(
        df["diameter_mm"], errors="coerce", downcast="integer"
    )
    df.loc[pd.isna(df["CapacityBcm/y"]), "capacity [MW]"] = df.loc[
        pd.isna(df["CapacityBcm/y"]), "diameter_mm"
    ].apply(lambda d: diameter_to_capacity(d))

    return df


def load_IGGIELGN_data(fn: Path) -> gpd.GeoDataFrame:
    """
    Load and flatten the IGGIELGN gas pipeline dataset.

    Parameters
    ----------
    fn : pathlib.Path
        Path to the IGGIELGN GeoJSON/GeoPackage file.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame with flattened parameter columns and dropped raw
        metadata fields.
    """
    df = gpd.read_file(fn)
    param = df.param.apply(pd.Series)
    method = df.method.apply(pd.Series)[["diameter_mm", "max_cap_M_m3_per_d"]]
    method.columns = method.columns + "_method"
    df = pd.concat([df, param, method], axis=1)
    to_drop = ["param", "uncertainty", "method", "tags"]
    to_drop = df.columns.intersection(to_drop)
    df.drop(to_drop, axis=1, inplace=True)
    return df


def prepare_IGGIELGN_data(
    df: gpd.GeoDataFrame,
    length_factor: float = 1.5,
    correction_threshold_length: float = 4,
    correction_threshold_p_nom: float = 8,
    bidirectional_below: float = 10,
) -> gpd.GeoDataFrame:  # Taken from pypsa-eur and adapted
    """
    Process IGGIELGN pipeline data and infer missing attributes.

    Parameters
    ----------
    df : geopandas.GeoDataFrame
        Raw IGGIELGN pipeline GeoDataFrame.
    length_factor : float, optional
        Multiplier applied to haversine-based distances when correcting
        reported line lengths.
    correction_threshold_length : float, optional
        Threshold ratio for when the reported length is replaced by the
        haversine-based length.
    correction_threshold_p_nom : float, optional
        Threshold ratio for when reported capacity is corrected using
        diameter-based capacity estimates.
    bidirectional_below : float, optional
        Line length (km) below which pipelines are assumed bidirectional.

    Returns
    -------
    geopandas.GeoDataFrame
        Cleaned pipeline GeoDataFrame with normalized capacity, length, and
        bidirectionality attributes.
    """
    # extract start and end from LineString
    df["point0"] = df.geometry.apply(lambda x: Point(x.coords[0]))
    df["point1"] = df.geometry.apply(lambda x: Point(x.coords[-1]))

    conversion_factor = 437.5  # MCM/day to MWh/h
    df["p_nom"] = df.max_cap_M_m3_per_d * conversion_factor

    # for inferred diameters, assume 500 mm rather than 900 mm (more conservative)
    df.loc[df.diameter_mm_method != "raw", "diameter_mm"] = 500.0

    keep = [
        "name",
        "diameter_mm",
        "is_H_gas",
        "is_bothDirection",
        "length_km",
        "p_nom",
        "max_pressure_bar",
        "start_year",
        "point0",
        "point1",
        "geometry",
    ]
    to_rename = {
        "is_bothDirection": "bidirectional",
        "is_H_gas": "H_gas",
        "start_year": "build_year",
        "length_km": "length",
    }
    df = df[keep].rename(columns=to_rename)

    df.bidirectional = df.bidirectional.astype(bool)
    df.H_gas = df.H_gas.astype(bool)

    # short lines below 10 km are assumed to be bidirectional
    short_lines = df["length"] < bidirectional_below
    df.loc[short_lines, "bidirectional"] = True

    # correct all capacities that deviate correction_threshold factor
    # to diameter-based capacities, unless they are NordStream pipelines
    # also all capacities below 0.5 GW are now diameter-based capacities
    df["p_nom_diameter"] = df.diameter_mm.apply(diameter_to_capacity)
    ratio = df.p_nom / df.p_nom_diameter
    not_nordstream = df.max_pressure_bar < 220
    df.p_nom.update(
        df.p_nom_diameter.where(
            (df.p_nom <= 500)
            | ((ratio > correction_threshold_p_nom) & not_nordstream)
            | ((ratio < 1 / correction_threshold_p_nom) & not_nordstream)
        )
    )

    # lines which have way too discrepant line lengths
    # get assigned haversine length * length factor
    df["length_haversine"] = df.apply(
        lambda p: length_factor
        * haversine_pts([p.point0.x, p.point0.y], [p.point1.x, p.point1.y]),
        axis=1,
    )
    ratio = df.eval("length / length_haversine")
    df["length"].update(
        df.length_haversine.where(
            (df["length"] < 20)
            | (ratio > correction_threshold_length)
            | (ratio < 1 / correction_threshold_length)
        )
    )

    # Convert CRS to EPSG:3857 so we can measure distances
    df = df.to_crs(epsg=3857)

    return df


def load_bus_regions(onshore_path, offshore_path=None):
    """
    Load onshore and optional offshore bus regions.

    The function reads onshore bus regions, converts them to EPSG:3857, standardizes
    the region identifier column to ``gadm_id``, and marks them as onshore. If an
    offshore regions file is provided, it is converted to the same CRS, standardized
    to the same column layout, and marked as offshore. The onshore and offshore
    regions are then combined and dissolved into model and country border
    geometries.

    Parameters
    ----------
    onshore_path : str or path-like
        Path to the onshore bus regions file. The file must contain ``name``,
        ``country`` and ``geometry`` columns.
    offshore_path : str or path-like, optional
        Path to the offshore bus regions file. If provided, the file must contain
        either a ``name`` or ``gadm_id`` column and a geometry column. If no
        ``country`` column is present, it is filled with ``None``.

    Returns
    -------
    tuple
        Tuple containing ``bus_regions_onshore``, ``bus_regions_offshore``,
        ``bus_regions_all``, ``model_borders`` and ``country_borders``.

        ``bus_regions_onshore`` and ``bus_regions_offshore`` contain standardized
        bus regions with ``gadm_id``, ``country``, ``geometry`` and
        ``is_offshore`` columns. ``bus_regions_all`` contains the combined non-empty
        onshore and offshore regions. ``model_borders`` contains the dissolved
        geometry of all model regions, while ``country_borders`` contains the
        dissolved geometry of the onshore regions.

    Raises
    ------
    KeyError
        If offshore regions are provided without a ``name`` or ``gadm_id`` column.
    """
    
    bus_regions_onshore = gpd.read_file(onshore_path).to_crs(epsg=3857)
    bus_regions_onshore = (
        bus_regions_onshore.rename(columns={"name": "gadm_id"})
        .loc[:, ["gadm_id", "country", "geometry"]]
        .copy()
    )
    bus_regions_onshore["is_offshore"] = False

    if offshore_path is not None:
        bus_regions_offshore = gpd.read_file(offshore_path).to_crs(bus_regions_onshore.crs)

        if "name" in bus_regions_offshore.columns and "gadm_id" not in bus_regions_offshore.columns:
            bus_regions_offshore = bus_regions_offshore.rename(columns={"name": "gadm_id"})
        elif "gadm_id" not in bus_regions_offshore.columns:
            raise KeyError("Offshore regions need a 'name' or 'gadm_id' column")

        if "country" not in bus_regions_offshore.columns:
            bus_regions_offshore["country"] = None

        bus_regions_offshore = (
            bus_regions_offshore.loc[:, ["gadm_id", "country", "geometry"]].copy()
        )
        bus_regions_offshore["is_offshore"] = True
    else:
        bus_regions_offshore = gpd.GeoDataFrame(
            columns=["gadm_id", "country", "geometry", "is_offshore"],
            geometry="geometry",
            crs=bus_regions_onshore.crs,
        )

    bus_regions_all = gpd.GeoDataFrame(
        pd.concat([bus_regions_onshore, bus_regions_offshore], ignore_index=True),
        geometry="geometry",
        crs=bus_regions_onshore.crs,
    )
    bus_regions_all = bus_regions_all[
        bus_regions_all.geometry.notna() & ~bus_regions_all.geometry.is_empty
    ].copy()

    # Important: CRS of borders = CRS of the underlying geometries
    model_borders = gpd.GeoDataFrame(
        geometry=[unary_union(bus_regions_all.geometry)],
        crs=bus_regions_all.crs,
    )
    country_borders = gpd.GeoDataFrame(
        geometry=[unary_union(bus_regions_onshore.geometry)],
        crs=bus_regions_onshore.crs,
    )

    return bus_regions_onshore, bus_regions_offshore, bus_regions_all, model_borders, country_borders



def get_states_in_order(
    pipeline: Union[LineString, MultiLineString],
    bus_regions_onshore: gpd.GeoDataFrame,
) -> list[str]:
    """
    Determine the ordered onshore regions traversed by a pipeline.

    Parameters
    ----------
    pipeline : shapely.geometry.LineString or MultiLineString
        Pipeline geometry to sample along its route.
    bus_regions_onshore : geopandas.GeoDataFrame
        Onshore region geometries with a `gadm_id` column.

    Returns
    -------
    list[str]
        Ordered list of `gadm_id` values for the regions intersected by the
        pipeline geometry.
    """
    states_p: List[str] = []

    if pipeline.geom_type == "LineString":
        # Interpolate points along the LineString with a given step size (e.g., 5)
        step_size = 10000
        interpolated_points = [
            pipeline.interpolate(i) for i in range(0, int(pipeline.length), step_size)
        ]
        interpolated_points.append(
            pipeline.interpolate(pipeline.length)
        )  # Add the last point

    elif pipeline.geom_type == "MultiLineString":
        interpolated_points = []
        # Iterate over each LineString within the MultiLineString
        for line in pipeline.geoms:
            # Interpolate points along each LineString with a given step size (e.g., 5)
            step_size = 10000
            interpolated_points_line = [
                line.interpolate(i) for i in range(0, int(line.length), step_size)
            ]
            interpolated_points_line.append(
                line.interpolate(line.length)
            )  # Add the last point
            interpolated_points.extend(interpolated_points_line)
    
    else:
        return states_p  # Return an empty list if the geometry type is not supported    
    
    # Check each interpolated point against the state geometries
    for point in interpolated_points:
        for _, state_row in bus_regions_onshore.iterrows():
            if state_row.geometry.contains(point):
                gadm_id = state_row["gadm_id"]
                if gadm_id not in states_p:
                    states_p.append(gadm_id)
                break  # Stop checking other states once a match is found

    return states_p


def parse_states(
    pipelines: gpd.GeoDataFrame, bus_regions_onshore: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Parse which onshore regions each pipeline traverses.

    Parameters
    ----------
    pipelines : geopandas.GeoDataFrame
        Pipeline geometries with a `geometry` column.
    bus_regions_onshore : geopandas.GeoDataFrame
        Onshore regions with a `gadm_id` column.

    Returns
    -------
    geopandas.GeoDataFrame
        Input pipeline GeoDataFrame augmented with `states_passed`,
        `amount_states_passed`, and `nodes` columns.
    """
    # Parse the states of the points which are connected by the pipeline geometry object
    pipelines["nodes"] = None
    pipelines["states_passed"] = None
    pipelines["amount_states_passed"] = None

    for pipeline, row in pipelines.iterrows():
        states_p = get_states_in_order(row.geometry, bus_regions_onshore)
        # states_p = pd.unique(states_p)
        row["states_passed"] = states_p
        row["amount_states_passed"] = len(states_p)
        row["nodes"] = list(zip(states_p[0::1], states_p[1::1]))
        pipelines.loc[pipeline] = row
    print(
        "The maximum number of states which are passed by one single pipeline amounts to {}.".format(
            pipelines.states_passed.apply(lambda n: len(n)).max()
        )
    )
    return pipelines


def _countries_from_row(row):
    """
    Extract country names associated with a GGIT pipeline record.

    The function collects country names from the ``StartCountry`` and
    ``EndCountry`` columns and additionally parses comma-separated entries from the
    ``Countries`` column. Missing values and empty country names are ignored.

    Parameters
    ----------
    row : pandas.Series
        GGIT pipeline record containing country metadata.

    Returns
    -------
    set of str
        Unique country names associated with the pipeline record.
    """
    countries = set()

    for col in ["StartCountry", "EndCountry"]:
        val = row.get(col)
        if pd.notna(val):
            countries.add(str(val).strip())

    val = row.get("Countries")
    if pd.notna(val):
        countries.update(c.strip() for c in str(val).split(",") if c.strip())

    return countries


def _iter_lines(geom):
    """
    Iterate over valid line geometries contained in a geometry object.

    The function yields a single ``LineString`` directly, or merges and yields the
    individual line parts of a ``MultiLineString``. Empty or missing geometries are
    ignored.

    Parameters
    ----------
    geom : shapely geometry
        Geometry to iterate over. Supported geometry types are ``LineString`` and
        ``MultiLineString``.

    Yields
    ------
    shapely.geometry.LineString
        Valid line geometries extracted from the input geometry.
    """
    if geom is None or geom.is_empty:
        return

    if geom.geom_type == "LineString":
        yield geom

    elif geom.geom_type == "MultiLineString":
        merged = linemerge(geom)

        if merged.geom_type == "LineString":
            yield merged
        elif merged.geom_type == "MultiLineString":
            for part in merged.geoms:
                if part is not None and not part.is_empty:
                    yield part


def _orient_from_outside(line, polygon):
    """
    Orient a line from outside to inside a polygon.

    The function checks whether exactly one endpoint of the line lies inside or on
    the polygon boundary. If the line starts outside and ends inside, it is returned
    unchanged. If it starts inside and ends outside, the coordinate order is
    reversed. Lines with both endpoints inside, both endpoints outside, or fewer
    than two coordinates are ignored.

    Parameters
    ----------
    line : shapely.geometry.LineString
        Line geometry to orient.
    polygon : shapely geometry
        Polygon or multipolygon used to determine whether the line endpoints are
        inside the target area.

    Returns
    -------
    shapely.geometry.LineString or None
        Line oriented from outside to inside the polygon, or ``None`` if no clear
        entry direction can be determined.
    """
    coords = list(line.coords)
    if len(coords) < 2:
        return None

    start = Point(coords[0])
    end = Point(coords[-1])

    start_in = polygon.covers(start)
    end_in = polygon.covers(end)

    # Only one clear case of entry: exactly one end outside, the other inside
    if (not start_in) and end_in:
        return line

    if start_in and (not end_in):
        return LineString(coords[::-1])

    return None

def _first_inside_point(line, polygon, eps_m=1000.0):
    """
    Return a point slightly inside a polygon along an oriented line.

    The line is expected to be oriented from outside to inside the polygon. The
    function intersects the line with the polygon, selects the first line segment
    inside the polygon along the original line direction, and interpolates a point
    slightly after the entry point. This avoids selecting a point exactly on the
    polygon boundary.

    Parameters
    ----------
    line : shapely.geometry.LineString
        Pipeline geometry oriented from outside to inside the polygon.
    polygon : shapely geometry
        Polygon or multipolygon used to determine the inside section of the line.
    eps_m : float, default 1000.0
        Maximum inward offset distance along the line. The distance uses the units
        of the input geometries.

    Returns
    -------
    shapely.geometry.Point or None
        Point located slightly inside the polygon along the line, or ``None`` if
        the line does not enter the polygon with a valid line segment.
    """
    inside = line.intersection(polygon)
    if inside.is_empty:
        return None

    if inside.geom_type == "LineString":
        segments = [inside]
    elif inside.geom_type == "MultiLineString":
        segments = [seg for seg in inside.geoms if not seg.is_empty]
    else:
        # only point contact or similar
        return None

    first_seg = None
    first_proj = None

    for seg in segments:
        p0 = Point(seg.coords[0])
        p1 = Point(seg.coords[-1])

        proj = min(line.project(p0), line.project(p1))
        if first_proj is None or proj < first_proj:
            first_proj = proj
            first_seg = seg

    if first_seg is None:
        return None

    # move slightly inward so that the point lies securely within a region
    step = min(eps_m, max(first_seg.length * 0.01, 1.0))
    dist = min(first_proj + step, line.length)

    return line.interpolate(dist)


def build_h2_pipeline_import_nodes_from_GGIT(
    pipelines,
    bus_regions_all,
    bus_regions_onshore,
    h2_capacity_factor,
    target_country_name,
    entry_scope="country",   # "country" = first onshore Bus, "model" = first on/offshore Bus
):
    """
    Build hydrogen pipeline import nodes from GGIT pipeline data.

    The function identifies international pipeline projects that enter the target
    country or model region, assigns each valid pipeline to the first bus region at
    the selected entry boundary, converts the reported pipeline capacity to hydrogen
    import capacity, and aggregates the resulting capacities by model node.

    Parameters
    ----------
    pipelines : geopandas.GeoDataFrame
        GGIT pipeline geometries. The dataframe must contain a geometry column,
        pipeline capacity in ``capacity [MW]``, and country information readable by
        ``_countries_from_row``. The optional column ``ProjectID`` is retained as
        project metadata.
    bus_regions_all : geopandas.GeoDataFrame
        Onshore and offshore model bus regions. The dataframe must contain
        ``geometry`` and ``gadm_id`` columns. Optional columns include ``country``
        and ``is_offshore``.
    bus_regions_onshore : geopandas.GeoDataFrame
        Onshore model bus regions used to define the country entry polygon. The
        dataframe must contain a geometry column and should use the same region
        identifiers as ``bus_regions_all``.
    h2_capacity_factor : float
        Factor used to convert the reported pipeline capacity in MW to hydrogen
        import capacity.
    target_country_name : str
        Name of the target country, for example ``"Japan"`` or ``"Germany"``.
    entry_scope : {"country", "model"}, default "country"
        Boundary used to determine the first entry point of a pipeline. ``"country"``
        assigns the pipeline to the first onshore bus region entered within the
        target country. ``"model"`` assigns the pipeline to the first bus region
        entered within the full model domain, including offshore regions.

    Returns
    -------
    pandas.DataFrame
        Aggregated hydrogen pipeline import nodes with the columns ``bus``,
        ``p_nom``, ``source``, ``is_offshore``, ``country`` and ``project_id``.
        ``project_id`` contains a list of GGIT project IDs contributing to each
        aggregated node.

    Raises
    ------
    ValueError
        If ``entry_scope`` is neither ``"country"`` nor ``"model"``.
    """
    if pipelines.crs != bus_regions_all.crs:
        pipelines = pipelines.to_crs(bus_regions_all.crs)

    if bus_regions_onshore.crs != bus_regions_all.crs:
        bus_regions_onshore = bus_regions_onshore.to_crs(bus_regions_all.crs)

    model_geom = unary_union(bus_regions_all.geometry)
    country_geom = unary_union(bus_regions_onshore.geometry)

    if entry_scope == "model":
        entry_polygon = model_geom
        search_regions = bus_regions_all.copy()
    elif entry_scope == "country":
        entry_polygon = country_geom
        search_regions = bus_regions_onshore.copy()
    else:
        raise ValueError("entry_scope must be either 'country' or 'model'")

    search_regions = search_regions[
        search_regions["gadm_id"].notna()
        & search_regions.geometry.notna()
        & ~search_regions.geometry.is_empty
    ].copy()

    sindex = search_regions.sindex
    records = []

    for _, row in pipelines.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        # 1) Consider only international lines
        countries = _countries_from_row(row)
        if target_country_name not in countries:
            continue
        if len(countries) < 2:
            continue

        # 2) The capacity must be known
        capacity = row.get("capacity [MW]")
        if capacity is None or pd.isna(capacity):
            continue

        matched_region = None

        for part in _iter_lines(geom):
            line = _orient_from_outside(part, entry_polygon)
            if line is None:
                continue

            if not line.intersects(entry_polygon):
                continue

            sample_point = _first_inside_point(line, entry_polygon, eps_m=1000.0)
            if sample_point is None:
                continue

            # Search by region using the spatial index
            cand_idx = list(sindex.intersection(sample_point.bounds))
            if not cand_idx:
                continue

            candidates = search_regions.iloc[cand_idx]
            candidates = candidates[
                candidates.geometry.apply(lambda g: g.covers(sample_point))
            ]

            if candidates.empty:
                continue

            # When creating a model entry, optionally select "offshore" if the item is not yet "onshore"
            if entry_scope == "model" and "is_offshore" in candidates.columns:
                if not country_geom.covers(sample_point):
                    candidates = candidates.sort_values("is_offshore", ascending=False)
                else:
                    candidates = candidates.sort_values("is_offshore", ascending=True)

            matched_region = candidates.iloc[0]
            break

        if matched_region is None:
            continue

        records.append(
            {
                "bus": matched_region["gadm_id"],
                "p_nom": float(capacity) * h2_capacity_factor,
                "source": "pipeline",
                "is_offshore": bool(matched_region.get("is_offshore", False)),
                "country": matched_region.get("country", None),
                "project_id": row.get("ProjectID", None),
            }
        )

    if not records:
        return pd.DataFrame(
            columns=["bus", "p_nom", "source", "is_offshore", "country", "project_id"]
        )

    df = pd.DataFrame(records)

    df = (
        df.groupby(
            ["bus", "source", "is_offshore", "country"],
            as_index=False,
            dropna=False,
        )
        .agg(
            {
                "p_nom": "sum",
                "project_id": lambda x: [v for v in x if pd.notna(v)],
            }
        )
    )

    return df



def cluster_gas_network(pipelines, bus_regions_onshore, length_factor):
    # drop innerstatal pipelines
    pipelines_interstate = pipelines.drop(
        pipelines.loc[pipelines.amount_states_passed < 2].index
    )

    # Convert CRS to EPSG:3857 so we can measure distances
    pipelines_interstate = pipelines_interstate.to_crs(epsg=3857)  # 3857

    # Perform overlay operation to split lines by polygons
    pipelines_interstate = gpd.overlay(
        pipelines_interstate, bus_regions_onshore, how="intersection"
    )

    column_set = ["ProjectID", "nodes", "gadm_id", "capacity [MW]"]

    if snakemake.params.gas_config["network_data"] == "IGGIELGN":
        pipelines_per_state = (
            pipelines_interstate.rename(
                {"p_nom": "capacity [MW]", "name": "ProjectID"}, axis=1
            )
            .loc[:, column_set]
            .reset_index(drop=True)
        )
    elif snakemake.params.gas_config["network_data"] == "GGIT":
        pipelines_per_state = pipelines_interstate.loc[:, column_set].reset_index(
            drop=True
        )

    # Explode the column containing lists of tuples
    df_exploded = pipelines_per_state.explode("nodes").reset_index(drop=True)

    # Create new columns for the tuples
    df_exploded.insert(0, "bus1", pd.DataFrame(df_exploded["nodes"].tolist())[1])
    df_exploded.insert(0, "bus0", pd.DataFrame(df_exploded["nodes"].tolist())[0])

    # Drop the original column
    df_exploded.drop("nodes", axis=1, inplace=True)

    # Reset the index if needed
    df_exploded.reset_index(drop=True, inplace=True)

    # Custom function to check if value in column 'gadm_id' exists in either column 'bus0' or column 'bus1'
    def check_existence(row):
        return row["gadm_id"] in [row["bus0"], row["bus1"]]

    # Apply the custom function to each row and keep only the rows that satisfy the condition
    df_filtered = df_exploded[df_exploded.apply(check_existence, axis=1)]
    df_grouped = df_filtered.groupby(["bus0", "bus1", "ProjectID"], as_index=False).agg(
        {
            "capacity [MW]": "first",
        }
    )

    # Rename columns to match pypsa-earth-sec format
    df_grouped = df_grouped.rename({"capacity [MW]": "capacity"}, axis=1).loc[
        :, ["bus0", "bus1", "capacity"]
    ]
    # df_exploded = df_exploded.loc[:, ['bus0', 'bus1', 'length']] # 'capacity'

    # Group by buses to get average length and sum of capacites of all pipelines between any two states on the route.
    grouped = df_grouped.groupby(["bus0", "bus1"], as_index=False).agg(
        {"capacity": "sum"}
    )
    states1 = bus_regions_onshore.copy()
    states1 = states1.set_index("gadm_id")

    # Create center points for each polygon and store them in a new column 'center_point'
    states1["center_point"] = (
        states1["geometry"].to_crs(3857).centroid.to_crs(4326)
    )  # ----> If haversine_pts method  for length calc is used
    # states1['center_point'] = states1['geometry'].centroid

    # Create an empty DataFrame to store distances
    distance_data = []

    # Iterate over all combinations of polygons
    for i in range(len(states1)):
        for j in range(len(states1)):
            if i != j:
                polygon1 = states1.iloc[i]
                polygon2 = states1.iloc[j]

                # Calculate Haversine distance
                distance = haversine_pts(
                    [
                        Point(polygon1["center_point"].coords[0]).x,
                        Point(polygon1["center_point"].coords[-1]).y,
                    ],
                    [
                        Point(polygon2["center_point"].coords[0]).x,
                        Point(polygon2["center_point"].coords[-1]).y,
                    ],
                )  # ----> If haversine_pts method  for length calc is used

                # Store the distance along with polygon IDs or other relevant information
                polygon_id1 = states1.index[i]
                polygon_id2 = states1.index[j]
                distance_data.append([polygon_id1, polygon_id2, distance])

    # Create a DataFrame from the distance data
    distance_df = pd.DataFrame(distance_data, columns=["bus0", "bus1", "distance"])

    merged_df = pd.merge(grouped, distance_df, on=["bus0", "bus1"], how="left")

    length_factor = 1.25

    merged_df["length"] = merged_df["distance"] * length_factor

    merged_df = merged_df.drop("distance", axis=1)

    merged_df["GWKm"] = (merged_df["capacity"] / 1000) * merged_df["length"]

    return merged_df


# TODO: Move it to a separate plotting rule!
# def plot_gas_network(pipelines, country_borders, bus_regions_onshore):
#     df = pipelines.copy()
#     df = gpd.overlay(df, country_borders, how="intersection")

#     if snakemake.params.gas_config["network_data"] == "IGGIELGN":
#         df = df.rename({"p_nom": "capacity [MW]"}, axis=1)

#     fig, ax = plt.subplots(1, 1)
#     fig.set_size_inches(12, 7)
#     bus_regions_onshore.to_crs(epsg=3857).plot(
#         ax=ax, color="white", edgecolor="darkgrey", linewidth=0.5
#     )
#     df.loc[(df.amount_states_passed > 1)].to_crs(epsg=3857).plot(
#         ax=ax,
#         column="capacity [MW]",
#         linewidth=2.5,
#         # linewidth=df['capacity [MW]'],
#         # alpha=0.8,
#         categorical=False,
#         cmap="viridis_r",
#         # legend=True,
#         # legend_kwds={'label':'Pipeline capacity [MW]'},
#     )

#     df.loc[(df.amount_states_passed <= 1)].to_crs(epsg=3857).plot(
#         ax=ax,
#         column="capacity [MW]",
#         linewidth=2.5,
#         # linewidth=df['capacity [MW]'],
#         alpha=0.5,
#         categorical=False,
#         # color='darkgrey',
#         ls="dotted",
#     )

#     # # Create custom legend handles for line types
#     # line_types = [ 'solid', 'dashed', 'dotted'] # solid
#     # legend_handles = [Line2D([0], [0], color='black', linestyle=line_type) for line_type in line_types]

#     # Define line types and labels
#     line_types = ["solid", "dotted"]
#     line_labels = ["Operating", "Not considered \n(within-state)"]

#     # Create custom legend handles for line types
#     legend_handles = [
#         Line2D([0], [0], color="black", linestyle=line_type, label=line_label)
#         for line_type, line_label in zip(line_types, line_labels)
#     ]

#     # Add the line type legend
#     ax.legend(
#         handles=legend_handles,
#         title="Status",
#         borderpad=1,
#         title_fontproperties={"weight": "bold"},
#         fontsize=11,
#         loc=1,
#     )

#     # # create the colorbar
#     norm = colors.Normalize(
#         vmin=df["capacity [MW]"].min(), vmax=df["capacity [MW]"].max()
#     )
#     cbar = plt.cm.ScalarMappable(norm=norm, cmap="viridis_r")
#     # fig.colorbar(cbar, ax=ax).set_label('Capacity [MW]')

#     # add colorbar
#     ax_cbar = fig.colorbar(cbar, ax=ax, location="left", shrink=0.8, pad=0.01)
#     # add label for the colorbar
#     ax_cbar.set_label("Natural gas pipeline capacity [MW]", fontsize=15)

#     ax.set_axis_off()
#     fig.savefig(snakemake.output.gas_network_fig_1, dpi=300, bbox_inches="tight")

# TODO: Move it to a separate plotting rule!
# def plot_clustered_gas_network(pipelines, bus_regions_onshore):
#     # Create a new GeoDataFrame with centroids
#     centroids = bus_regions_onshore.copy()
#     centroids["geometry"] = centroids["geometry"].centroid
#     centroids["gadm_id"] = centroids["gadm_id"].apply(
#         lambda id: three_2_two_digits_country(id[:3]) + id[3:]
#     )
#     gdf1 = pd.merge(
#         pipelines, centroids, left_on=["bus0"], right_on=["gadm_id"], how="left"
#     )
#     gdf1.rename(columns={"geometry": "geometry_bus0"}, inplace=True)
#     pipe_links = pd.merge(
#         gdf1, centroids, left_on=["bus1"], right_on=["gadm_id"], how="left"
#     )
#     pipe_links.rename(columns={"geometry": "geometry_bus1"}, inplace=True)

#     # Create LineString geometries from the points
#     pipe_links["geometry"] = pipe_links.apply(
#         lambda row: LineString([row["geometry_bus0"], row["geometry_bus1"]]), axis=1
#     )

#     clustered = gpd.GeoDataFrame(pipe_links, geometry=pipe_links["geometry"])

#     # Optional: Set the coordinate reference system (CRS) if needed
#     clustered.crs = "EPSG:3857"  # For example, WGS84

#     # plot pipelines
#     fig, ax = plt.subplots(1, 1)
#     fig.set_size_inches(12, 7)
#     bus_regions_onshore.to_crs(epsg=3857).plot(
#         ax=ax, color="white", edgecolor="darkgrey", linewidth=0.5
#     )
#     clustered.to_crs(epsg=3857).plot(
#         ax=ax,
#         column="capacity",
#         linewidth=1.5,
#         categorical=False,
#         cmap="plasma",
#     )

#     centroids.to_crs(epsg=3857).plot(
#         ax=ax,
#         color="red",
#         markersize=35,
#         alpha=0.5,
#     )

#     norm = colors.Normalize(
#         vmin=pipelines["capacity"].min(), vmax=pipelines["capacity"].max()
#     )
#     cbar = plt.cm.ScalarMappable(norm=norm, cmap="plasma")

#     # add colorbar
#     ax_cbar = fig.colorbar(cbar, ax=ax, location="left", shrink=0.8, pad=0.01)
#     # add label for the colorbar
#     ax_cbar.set_label("Natural gas pipeline capacity [MW]", fontsize=15)

#     ax.set_axis_off()
#     fig.savefig(snakemake.output.gas_network_fig_2, dpi=300, bbox_inches="tight")


if not snakemake.params.custom_gas_network:
    if snakemake.params.gas_config["network_data"] == "GGIT":
        pipelines = download_GGIT_gas_network()
        pipelines = prepare_GGIT_data(pipelines)

    elif snakemake.params.gas_config["network_data"] == "IGGIELGN":
        download_IGGIELGN_gas_network()

        gas_network = os.path.join(
            BASE_DIR, "data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson"
        )

        pipelines = load_IGGIELGN_data(gas_network)
        pipelines = prepare_IGGIELGN_data(pipelines)

    (
        bus_regions_onshore,
        bus_regions_offshore,
        bus_regions_all,
        model_borders,
        country_borders,
    ) = load_bus_regions(
        snakemake.input.regions_onshore,
        snakemake.input.regions_offshore,
    )

    bus_regions_onshore.geometry = bus_regions_onshore.geometry.apply(make_valid)
    bus_regions_offshore.geometry = bus_regions_offshore.geometry.apply(make_valid)
    bus_regions_all.geometry = bus_regions_all.geometry.apply(make_valid)
    model_borders.geometry = model_borders.geometry.apply(make_valid)
    country_borders.geometry = country_borders.geometry.apply(make_valid)

    if snakemake.params.gas_config["network_data"] == "GGIT":
        
        h2_pipeline_import_nodes = build_h2_pipeline_import_nodes_from_GGIT(
        pipelines=pipelines,
        bus_regions_all=bus_regions_all,
        bus_regions_onshore=bus_regions_onshore,
        h2_capacity_factor=snakemake.params.hydrogen_config["H2_retrofit_capacity_per_CH4"],
        target_country_name=bus_regions_onshore["country"].dropna().iloc[0],
        entry_scope="country",
    )
    else:
        h2_pipeline_import_nodes = pd.DataFrame(
            columns=[
                "bus",
                "p_nom",
                "source",
                "is_offshore",
                "country",
                "project_id",
            ]
        )

    h2_pipeline_import_nodes.to_csv(
        snakemake.output.h2_pipeline_import_nodes,
        index=False,
    )

    pipelines = parse_states(pipelines, bus_regions_onshore)

    if len(pipelines.loc[pipelines.amount_states_passed >= 2]) > 0:
        pipelines = cluster_gas_network(
            pipelines, bus_regions_onshore, length_factor=1.25
        )

        pipelines.to_csv(snakemake.output.clustered_gas_network, index=False)

        average_length = pipelines["length"].mean()
        print("average_length = ", average_length)

        total_system_capacity = pipelines["GWKm"].sum()
        print("total_system_capacity = ", total_system_capacity)

    else:
        bus_regions_onshore["country"] = bus_regions_onshore["gadm_id"]
        print(
            "The following countries have no existing Natural Gas network between the chosen bus regions:\n"
            + ", ".join(bus_regions_onshore.country.unique().tolist())
        )

        pipelines = pd.DataFrame(
            {"bus0": [], "bus1": [], "capacity": [], "length": [], "GWKm": []}
        )
        pipelines.to_csv(snakemake.output.clustered_gas_network, index=False)
        
    
