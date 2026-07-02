# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: Open Energy Transition
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import logging
import re
import sys
import tempfile
import zipfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pypsa
import requests
from scipy.spatial import cKDTree

PYPSA_EARTH_SCRIPTS = Path.cwd() / "scripts"
if not PYPSA_EARTH_SCRIPTS.exists():
    PYPSA_EARTH_SCRIPTS = (
        Path(__file__).resolve().parents[1] / "pypsa-earth" / "scripts"
    )
sys.path.insert(0, str(PYPSA_EARTH_SCRIPTS))

from _helpers import mock_snakemake

logger = logging.getLogger(__name__)

POA_SHAPES_URL = (
    "https://www.abs.gov.au/statistics/standards/"
    "australian-statistical-geography-standard-asgs-edition-3/"
    "jul2021-jun2026/access-and-downloads/digital-boundary-files/"
    "POA_2021_AUST_GDA2020_SHP.zip"
)

POA_SHAPES_DIR = Path("data/shapes/POA_2021_AUST_GDA2020_SHP")
POA_SHP_FILE = POA_SHAPES_DIR / "POA_2021_AUST_GDA2020.shp"


def ensure_poa_shapefile() -> Path:
    """
    Ensure the ABS POA shapefile is available locally.

    Returns
    -------
    Path
        Path to the POA shapefile.
    """
    if POA_SHP_FILE.exists():
        logger.info("Using existing POA shapefile at %s", POA_SHP_FILE)
        return POA_SHP_FILE

    logger.info("POA shapefile not found. Downloading from ABS.")
    POA_SHAPES_DIR.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        zip_path = tmpdir / "POA_2021_AUST_GDA2020_SHP.zip"

        response = requests.get(POA_SHAPES_URL, timeout=120)
        response.raise_for_status()
        zip_path.write_bytes(response.content)

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(POA_SHAPES_DIR.parent)

    if not POA_SHP_FILE.exists():
        raise FileNotFoundError(
            f"Expected shapefile not found after download: {POA_SHP_FILE}"
        )

    logger.info("Downloaded and extracted POA shapefile to %s", POA_SHAPES_DIR)
    return POA_SHP_FILE


def parse_month_column(col: str):
    """
    Parse CER monthly column names like:
    'Jan 2011 - Rated Power Output In kW'

    Returns
    -------
    pandas.Timestamp | None
        First day of the month if parsing succeeds, otherwise None.
    """
    match = re.match(r"^([A-Z][a-z]{2})\s+(\d{4})\s+-\s+Rated Power Output In kW$", col)
    if not match:
        return None
    month_abbr, year = match.groups()
    return pd.to_datetime(f"01 {month_abbr} {year}", format="%d %b %Y")


def detect_postcode_column(df: pd.DataFrame) -> str:
    """
    Detect postcode column in CER dataset.
    """
    candidates = [
        "Small Unit Postcode",
        "Small Unit Installation Postcode",
        "Small Unit Installation Post Code",
        "postcode",
        "Postcode",
        "POSTCODE",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    raise KeyError(
        "Could not find a postcode column in CER dataset. "
        f"Available columns: {list(df.columns)}"
    )


def detect_capacity_column(df: pd.DataFrame, base_year: int) -> str | None:
    """
    Detect a direct cumulative capacity column for the requested base year.

    The CER SGU monthly columns are monthly additions, not cumulative values,
    so they must not be used as direct cumulative capacity columns.
    """
    candidates = [
        f"Total Rated Power Output In kW",
    ]

    for col in candidates:
        if col in df.columns:
            return col

    return None


def build_cumulative_capacity_by_postcode(
    cer_path: str | Path,
    base_year: int,
) -> pd.DataFrame:
    """
    Build cumulative rooftop solar capacity by postcode up to Dec of base_year.

    Parameters
    ----------
    cer_path : str or Path
        Path to CER rooftop solar CSV.
    base_year : int
        Base year.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['postcode', 'capacity_kw'].
    """
    logger.info("Reading CER dataset from %s", cer_path)
    df = pd.read_csv(cer_path, thousands=",")

    postcode_col = detect_postcode_column(df)
    historic_col = "Historic Total Rated Power Output In kW (2001 - 2010)"
    direct_capacity_col = detect_capacity_column(df, base_year)

    if historic_col not in df.columns and direct_capacity_col is None:
        raise KeyError(
            f"Neither '{historic_col}' nor a direct Dec {base_year} capacity column "
            "was found in CER dataset."
        )

    logger.info("Using postcode column: %s", postcode_col)

    df = df.copy()
    df[postcode_col] = (
        df[postcode_col].astype(str).str.extract(r"(\d+)")[0].str.zfill(4)
    )

    if direct_capacity_col is not None:
        logger.info("Using direct cumulative capacity column: %s", direct_capacity_col)

        df[direct_capacity_col] = pd.to_numeric(
            df[direct_capacity_col], errors="coerce"
        ).fillna(0.0)

        out = (
            df[[postcode_col, direct_capacity_col]]
            .rename(
                columns={postcode_col: "postcode", direct_capacity_col: "capacity_kw"}
            )
            .groupby("postcode", as_index=False)["capacity_kw"]
            .sum()
        )

        logger.info(
            "Built cumulative rooftop capacity by postcode from direct Dec %s column: %.3f GW total",
            base_year,
            out["capacity_kw"].sum() / 1e6,
        )
        return out

    logger.info(
        "Direct Dec %s column not found. Falling back to historic + monthly additions.",
        base_year,
    )

    df[historic_col] = pd.to_numeric(df[historic_col], errors="coerce").fillna(0.0)

    monthly_cols = []
    monthly_dates = {}

    for col in df.columns:
        dt = parse_month_column(col)
        if dt is not None:
            monthly_cols.append(col)
            monthly_dates[col] = dt

    cutoff = pd.Timestamp(year=base_year, month=12, day=1)
    selected_cols = [col for col in monthly_cols if monthly_dates[col] <= cutoff]

    if not selected_cols:
        raise ValueError(f"No CER monthly columns found up to Dec {base_year}.")

    logger.info(
        "Using CER monthly additions from Jan 2011 to Dec %s (%d monthly columns)",
        base_year,
        len(selected_cols),
    )

    for col in selected_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    df["capacity_kw"] = df[historic_col] + df[selected_cols].sum(axis=1)

    out = (
        df[[postcode_col, "capacity_kw"]]
        .rename(columns={postcode_col: "postcode"})
        .groupby("postcode", as_index=False)["capacity_kw"]
        .sum()
    )

    logger.info(
        "Built cumulative rooftop capacity by postcode up to Dec %s: %.3f GW total",
        base_year,
        out["capacity_kw"].sum() / 1e6,
    )

    return out


def load_postcode_centroids(shapefile_path: str | Path) -> pd.DataFrame:
    """
    Load ABS POA shapefile and compute postcode centroids.

    Parameters
    ----------
    shapefile_path : str or Path
        Path to POA shapefile (.shp).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['postcode', 'lon', 'lat'].
    """
    logger.info("Reading POA shapefile from %s", shapefile_path)
    poa = gpd.read_file(shapefile_path)

    postcode_col = "POA_CODE21"
    if postcode_col not in poa.columns:
        raise KeyError(f"Column '{postcode_col}' not found in POA shapefile.")

    centroids = gpd.GeoSeries(poa.geometry.centroid, crs=poa.crs).to_crs(epsg=4326)

    out = pd.DataFrame(
        {
            "postcode": poa[postcode_col].astype(str).str.zfill(4),
            "lon": centroids.x,
            "lat": centroids.y,
        }
    ).dropna(subset=["postcode", "lon", "lat"])

    logger.info("Loaded %d postcode centroids from POA shapefile", len(out))
    return out


def map_postcodes_to_nearest_buses(
    rooftop_by_postcode: pd.DataFrame,
    postcode_centroids: pd.DataFrame,
    network_path: str | Path,
) -> pd.DataFrame:
    """
    Map postcode rooftop capacities to nearest AC buses in the network.

    Parameters
    ----------
    rooftop_by_postcode : pd.DataFrame
        DataFrame with ['postcode', 'capacity_kw'].
    postcode_centroids : pd.DataFrame
        DataFrame with ['postcode', 'lon', 'lat'].
    network_path : str or Path
        Path to the PyPSA network file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['node', 'p_nom'] in MW.
    """
    logger.info("Reading network from %s", network_path)
    n = pypsa.Network(network_path)

    ac_buses = n.buses[n.buses.carrier == "AC"].copy()
    if ac_buses.empty:
        raise ValueError("No AC buses found in the network.")

    ac_buses = ac_buses.dropna(subset=["x", "y"])
    if ac_buses.empty:
        raise ValueError("AC buses do not have valid coordinates.")

    merged = rooftop_by_postcode.merge(postcode_centroids, on="postcode", how="inner")
    merged = merged.dropna(subset=["lon", "lat", "capacity_kw"])

    if merged.empty:
        raise ValueError(
            "No postcode capacities could be matched to postcode centroids."
        )

    logger.info(
        "Matched %d postcodes with centroid geometries (%.3f GW total)",
        len(merged),
        merged["capacity_kw"].sum() / 1e6,
    )

    bus_names = ac_buses.index.to_numpy()
    bus_xy = ac_buses[["x", "y"]].to_numpy()
    tree = cKDTree(bus_xy)

    postcode_xy = merged[["lon", "lat"]].to_numpy()
    _, indices = tree.query(postcode_xy, k=1)

    merged["node"] = bus_names[indices]
    merged["p_nom"] = merged["capacity_kw"] / 1000.0  # kW -> MW

    out = merged.groupby("node", as_index=False)["p_nom"].sum()

    logger.info(
        "Mapped rooftop solar to %d AC buses (%.3f GW total)",
        len(out),
        out["p_nom"].sum() / 1000.0,
    )

    return out.sort_values("node").reset_index(drop=True)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_rooftop_solar_existing",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2025",
        )

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s:%(name)s:%(message)s",
    )

    base_year = int(str(snakemake.wildcards.planning_horizons)[-4:])

    rooftop_by_postcode = build_cumulative_capacity_by_postcode(
        cer_path=snakemake.input.cer,
        base_year=base_year,
    )

    poa_shapefile = ensure_poa_shapefile()

    postcode_centroids = load_postcode_centroids(
        shapefile_path=poa_shapefile,
    )

    rooftop_by_node = map_postcodes_to_nearest_buses(
        rooftop_by_postcode=rooftop_by_postcode,
        postcode_centroids=postcode_centroids,
        network_path=snakemake.input.network,
    )

    output_path = Path(snakemake.output[0])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rooftop_by_node.to_csv(output_path, index=False)

    logger.info("Saved rooftop solar existing capacities to %s", output_path)
