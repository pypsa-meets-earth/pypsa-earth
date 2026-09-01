# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Build hydrogen import nodes from pipeline- and port-based import sources.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    imports:
        sources:
        limit:

.. seealso::

    Documentation of the configuration file ``config.yaml``.

Inputs
------

- ``h2_pipeline_import_nodes``: CSV file with fixed pipeline-based hydrogen
  import nodes and nominal import capacities.
- ``export_ports``: CSV file with prepared port locations, country codes and
  capacity allocation fractions.
- ``regions_onshore``: GeoJSON file with onshore bus regions used to assign
  port import locations to model nodes.

Outputs
-------

- ``h2_import_nodes``: CSV file with aggregated maximum hydrogen import
  capacities per model node, split into pipeline-based, port-based and total
  import capacity.

Description
-----------

The script builds model nodes for hydrogen imports from the import sources
selected in the configuration. Pipeline-based import capacities are read as
fixed capacities from the pipeline import table. Port-based import capacities
are derived from the remaining annual import budget after subtracting the fixed
pipeline import capacity, distributed across valid ports according to their
capacity fractions, and mapped to onshore bus regions.

The resulting table contains one row per import node with the maximum nominal
pipeline, port and total hydrogen import capacity.
"""

import logging
from os import path
from pathlib import Path

import fiona
import geopandas as gpd
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import pandas as pd
from pypsa.geo import haversine_pts
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
from shapely.validation import make_valid


logger = logging.getLogger(__name__)


def load_pipeline_import_nodes(path):
    """
    Load and validate pipeline-based hydrogen import nodes.

    Parameters
    ----------
    path : str or path-like
        Path to the CSV file containing pipeline-based hydrogen import nodes. The
        file must contain the columns ``bus``, ``country`` and ``p_nom``.

    Returns
    -------
    pandas.DataFrame
        Cleaned pipeline import table with valid bus names, country codes and
        non-negative nominal import capacities.

    Raises
    ------
    ValueError
        If required columns are missing or negative ``p_nom`` values are found.
    """
    df = pd.read_csv(
        path,
        index_col=None,
        keep_default_na=False,
    ).copy()

    required_columns = {"bus", "country", "p_nom"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in pipeline import file {path}: "
            f"{sorted(missing_columns)}"
        )

    df["bus"] = df["bus"].astype(str).str.strip()
    df["country"] = df["country"].astype(str).str.strip()
    df["p_nom"] = pd.to_numeric(df["p_nom"], errors="coerce")

    df = df.dropna(subset=["bus", "country", "p_nom"])
    df = df[df["bus"] != ""]
    df = df[df["country"] != ""]

    if (df["p_nom"] < 0).any():
        invalid = df.loc[df["p_nom"] < 0, ["bus", "country", "p_nom"]]
        raise ValueError(
            "Negative p_nom values found in pipeline import file:\n"
            f"{invalid.to_string(index=False)}"
        )

    return df
    


def load_ports(path):
    """
    Load and validate prepared hydrogen import ports.

    Parameters
    ----------
    path : str or path-like
        Path to the prepared ports CSV file. The file must contain the columns
        ``country``, ``x``, ``y`` and ``fraction``.

    Returns
    -------
    pandas.DataFrame
        Cleaned port table filtered to the model countries defined in the
        configuration.

    Raises
    ------
    ValueError
        If required columns are missing or negative ``fraction`` values are found.
    """
    ports = pd.read_csv(
            path,
            index_col=None,
            keep_default_na=False,
        ).copy()

    required_columns = {"country", "x", "y", "fraction"}
    missing_columns = required_columns.difference(ports.columns)
    if missing_columns:
        raise ValueError(
            f"Missing required columns in ports file {path}: "
            f"{sorted(missing_columns)}"
        )

    ports["country"] = ports["country"].astype(str).str.strip()
    ports["x"] = pd.to_numeric(ports["x"], errors="coerce")
    ports["y"] = pd.to_numeric(ports["y"], errors="coerce")
    ports["fraction"] = pd.to_numeric(ports["fraction"], errors="coerce")

    ports = ports.dropna(subset=["country", "x", "y", "fraction"])
    ports = ports[ports["country"] != ""]

    if (ports["fraction"] < 0).any():
        invalid = ports.loc[ports["fraction"] < 0, ["country", "x", "y", "fraction"]]
        raise ValueError(
            "Negative fraction values found in ports file:\n"
            f"{invalid.to_string(index=False)}"
        )

    model_countries = set(snakemake.config["countries"])
    ports = ports[ports["country"].isin(model_countries)].copy()

    return ports


def get_remaining_port_budget(pipeline_p_nom_max_total):
    """
    Calculate the remaining port-based hydrogen import capacity budget.

    The total annual hydrogen import limit is read from
    ``snakemake.params.imports_config['limit']`` and converted from TWh/a to MW.
    The fixed pipeline-based import capacity is subtracted from this limit. If the
    pipeline capacity already reaches or exceeds the import limit, the remaining
    port budget is set to zero.

    Parameters
    ----------
    pipeline_p_nom_max_total : float
        Total fixed pipeline-based hydrogen import capacity in MW.

    Returns
    -------
    float
        Remaining maximum nominal capacity in MW available for port-based hydrogen
        imports.

    Raises
    ------
    ValueError
        If the import limit or the pipeline import capacity is negative.
    """
    import_limit = snakemake.params.h2import
    import_limit = float(import_limit) * 1e6 / 8760  # Convert from TWh/a to MW

    if import_limit is None:
        raise ValueError(
            "Missing config entry for imports limit: "
            "snakemake.params.imports_config['limit']"
        )

    pipeline_p_nom_max_total = float(pipeline_p_nom_max_total)

    if import_limit < 0:
        raise ValueError("Import limit must be non-negative.")

    if pipeline_p_nom_max_total < 0:
        raise ValueError("pipeline_p_nom_max_total must be non-negative.")

    return max(import_limit - pipeline_p_nom_max_total, 0.0)


def assign_port_p_nom(ports, pipeline_p_nom_max_total):
    """
    Assign maximum port-based hydrogen import capacity to ports.

    The remaining port import budget is calculated after subtracting fixed
    pipeline-based import capacity from the total import limit. This remaining
    budget is distributed across ports according to their normalized ``fraction``
    values.

    Parameters
    ----------
    ports : pandas.DataFrame
        Port table with at least the column ``fraction``.
    pipeline_p_nom_max_total : float
        Total fixed pipeline-based hydrogen import capacity in MW.

    Returns
    -------
    pandas.DataFrame
        Copy of the port table with an additional ``p_nom_max`` column containing
        the assigned maximum nominal port import capacity in MW.

    Raises
    ------
    ValueError
        If the ``fraction`` column is missing, contains negative values, or sums to
        zero after cleaning.
    """
    ports = ports.copy()

    if "fraction" not in ports.columns:
        raise ValueError("Ports dataframe must contain column 'fraction'.")

    if ports.empty:
        ports["p_nom_max"] = pd.Series(dtype=float)
        return ports

    ports["fraction"] = pd.to_numeric(ports["fraction"], errors="coerce")
    ports = ports.dropna(subset=["fraction"])

    if (ports["fraction"] < 0).any():
        invalid = ports.loc[ports["fraction"] < 0, ["fraction"]]
        raise ValueError(
            "Negative values found in port 'fraction' column:\n"
            f"{invalid.to_string(index=False)}"
        )

    port_budget = get_remaining_port_budget(pipeline_p_nom_max_total)

    fraction_sum = ports["fraction"].sum()
    if fraction_sum <= 0:
        raise ValueError("Sum of port fractions must be greater than zero.")

    ports["fraction"] = ports["fraction"] / fraction_sum
    ports["p_nom_max"] = ports["fraction"] * port_budget

    return ports


def map_import_locations_to_nodes(df, regions_or_buses):
    """
    Map hydrogen import locations to onshore model bus regions.

    Import locations are first assigned to intersecting onshore bus regions. For
    locations that do not intersect any region, the nearest onshore bus region is
    used as a fallback. The result is restricted to one bus assignment per original
    import location.

    Parameters
    ----------
    df : pandas.DataFrame
        Import locations with at least the columns ``x``, ``y`` and ``country``.
    regions_or_buses : str or path-like
        Path to the onshore bus regions GeoJSON file. The file must contain the
        columns ``name`` and ``geometry``.

    Returns
    -------
    pandas.DataFrame
        Import locations with an added ``bus`` column identifying the assigned
        model node.

    Raises
    ------
    ValueError
        If required columns are missing from the import locations or the bus
        regions file.
    """
    df = df.reset_index(drop=True).copy()
    df["port_id"] = df.index

    required_columns = {"x", "y", "country"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            "Missing required columns for mapping import locations: "
            f"{sorted(missing_columns)}"
        )

    gdf_regions = gpd.read_file(regions_or_buses).copy()

    required_region_columns = {"name", "geometry"}
    missing_region_columns = required_region_columns.difference(gdf_regions.columns)
    if missing_region_columns:
        raise ValueError(
            "Missing required columns in bus regions file: "
            f"{sorted(missing_region_columns)}"
        )

    gdf_regions = gdf_regions.rename(columns={"name": "bus"})

    if gdf_regions.crs is None:
        gdf_regions = gdf_regions.set_crs("EPSG:4326")
    else:
        gdf_regions = gdf_regions.to_crs("EPSG:4326")

    gdf_points = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["x"], df["y"]),
        crs="EPSG:4326",
    )

    # 1) first try direct match
    strict_matches = gpd.sjoin(
        gdf_points,
        gdf_regions[["bus", "geometry"]],
        how="inner",
        predicate="intersects",
        rsuffix="region",
    )
    strict_matches["distance_to_region"] = 0.0

    matched_port_ids = set(strict_matches["port_id"].tolist())

    # 2) fallback: nearest only for still-unmatched ports
    unmatched_points = gdf_points.loc[~gdf_points["port_id"].isin(matched_port_ids)].copy()

    if not unmatched_points.empty:
        nearest_matches = gpd.sjoin_nearest(
            unmatched_points,
            gdf_regions[["bus", "geometry"]],
            how="left",
            rsuffix="region",
            distance_col="distance_to_region",
        )
        gdf_joined = pd.concat([strict_matches, nearest_matches], ignore_index=True)
    else:
        gdf_joined = strict_matches.copy()

    gdf_joined = gdf_joined.dropna(subset=["bus"]).copy()

    # IMPORTANT: enforce exactly one mapping per original port
    gdf_joined = gdf_joined.sort_values(
        ["port_id", "distance_to_region", "bus"]
    )
    gdf_joined = gdf_joined.drop_duplicates(subset=["port_id"], keep="first")

    return pd.DataFrame(
        gdf_joined.drop(
            columns=["geometry", "index_region", "distance_to_region"],
            errors="ignore",
        )
    )
 
 
def aggregate_h2_import_by_node(pipelines, ports, pipeline_p_nom_max_total, tol=1e-6):
    """
    Aggregate pipeline- and port-based hydrogen import capacities by model node.

    Pipeline import capacities and assigned port import capacities are grouped by
    ``bus`` and combined into one table. After aggregation, port import capacities
    are capped to the remaining port budget as a final safety check.

    Parameters
    ----------
    pipelines : pandas.DataFrame
        Pipeline import nodes with at least the columns ``bus`` and ``p_nom``.
    ports : pandas.DataFrame
        Mapped port import nodes with at least the columns ``bus`` and
        ``p_nom_max``.
    pipeline_p_nom_max_total : float
        Total fixed pipeline-based hydrogen import capacity in MW.
    tol : float, default 1e-6
        Numerical tolerance used when checking whether aggregated port capacity
        exceeds the remaining port budget.

    Returns
    -------
    pandas.DataFrame
        Aggregated hydrogen import capacities per model node with the columns
        ``bus``, ``p_nom_max_pipeline``, ``p_nom_max_port`` and
        ``p_nom_max_total``.

    Raises
    ------
    ValueError
        If required columns are missing from the pipeline or port tables.
    """
    pipeline_required = {"bus", "p_nom"}
    port_required = {"bus", "p_nom_max"}

    missing_pipeline = pipeline_required.difference(pipelines.columns)
    if missing_pipeline:
        raise ValueError(
            "Missing required columns in pipelines dataframe: "
            f"{sorted(missing_pipeline)}"
        )

    missing_ports = port_required.difference(ports.columns)
    if missing_ports:
        raise ValueError(
            "Missing required columns in ports dataframe: "
            f"{sorted(missing_ports)}"
        )

    pipeline_agg = (
        pipelines.copy()
        .groupby("bus", as_index=False)["p_nom"]
        .sum()
        .rename(columns={"p_nom": "p_nom_max_pipeline"})
    )

    port_agg = (
        ports.copy()
        .groupby("bus", as_index=False)["p_nom_max"]
        .sum()
        .rename(columns={"p_nom_max": "p_nom_max_port"})
    )

    # final safety cap after node aggregation
    port_budget = get_remaining_port_budget(pipeline_p_nom_max_total)

    if not port_agg.empty:
        total_port = port_agg["p_nom_max_port"].sum()

        if port_budget <= 0.0:
            port_agg["p_nom_max_port"] = 0.0
        elif total_port > port_budget + tol:
            scale = port_budget / total_port
            logger.warning(
                "Aggregated port import capacity %.3f exceeds remaining port "
                "budget %.3f. Scaling by factor %.6f.",
                total_port,
                port_budget,
                scale,
            )
            port_agg["p_nom_max_port"] = port_agg["p_nom_max_port"] * scale

    df = pipeline_agg.merge(
        port_agg,
        on="bus",
        how="outer",
    ).fillna(0.0)

    df["p_nom_max_total"] = df["p_nom_max_pipeline"] + df["p_nom_max_port"]
    df = df.sort_values("bus").reset_index(drop=True)

    return df

    
if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_h2_import_locations",
            clusters="10",
        )

    sources = set(snakemake.params.imports_config.get("sources", []))

    if "pipelines" in sources:
        pipelines = load_pipeline_import_nodes(snakemake.input.h2_pipeline_import_nodes)
    else:
        pipelines = pd.DataFrame(columns=["bus", "country", "p_nom"])

    pipeline_p_nom_max_total = pd.to_numeric(
        pipelines.get("p_nom", pd.Series(dtype=float)),
        errors="coerce",
    ).fillna(0.0).sum()

    if "ports" in sources:
        ports = load_ports(snakemake.input.export_ports)

        # first map to bus nodes
        ports = map_import_locations_to_nodes(
            ports,
            snakemake.input.regions_onshore,
        )

        # then distribute remaining budget only across valid mapped ports
        ports = assign_port_p_nom(ports, pipeline_p_nom_max_total)
    else:
        ports = pd.DataFrame(columns=["bus", "country", "fraction", "p_nom_max"])

    h2_import_nodes = aggregate_h2_import_by_node(
        pipelines,
        ports,
        pipeline_p_nom_max_total,
    )

    h2_import_nodes.to_csv(
        snakemake.output.h2_import_nodes,
        index=False,
    )