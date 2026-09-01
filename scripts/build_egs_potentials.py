# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Build enhanced geothermal system potentials and optional capacity-factor profiles.

Relevant Settings
-----------------
.. code:: yaml

    countries:
    snapshots:

.. seealso::

    Documentation of the configuration file ``config.yaml``.

Inputs
------
- ``egs_input``: CSV file with point-based enhanced geothermal system data from
  Franzmann, Heinrichs and Stolten (2025), "Global geothermal electricity
  potentials: A technical, economic, and thermal renewability assessment",
  Renewable Energy, 250, 123199, and the associated Mendeley Data dataset
  "Global Technical Energy Potentials of Enhanced Geothermal Systems".
  https://www.sciencedirect.com/science/article/pii/S0960148125008614
  The file includes coordinates, country identifiers, power output estimates
  and LCOE values for different EGS reservoir models.
- ``regions``: GeoJSON file with PyPSA-Earth network regions used to aggregate
  point-based EGS potentials to model nodes.
- ``temp_air_total``: Optional NetCDF file with hourly air temperatures per
  network region used to derive temperature-dependent EGS capacity factors.

Outputs
-------
- ``egs_potentials``: CSV file with maximum installable EGS capacity and CAPEX
  per network region.
- ``egs_capacity_factors``: Optional CSV file with hourly EGS capacity-factor
  profiles per network region.

Description
-----------
The script prepares enhanced geothermal system potentials for PyPSA-Earth from
the global point-based EGS potential dataset published by Franzmann, Heinrichs
and Stolten (2025). The underlying study assesses global geothermal electricity
potentials from a technical, economic and thermal renewability perspective and
compares three reservoir models: Gringarten, Volume Method and Sustainable
Approach.

For use in PyPSA-Earth, this script filters the input data to the configured
model countries, uses the Gringarten power and LCOE estimates, converts power
output from MW to GW, derives a CAPEX-like cost value from LCOE assumptions,
and assigns each valid EGS point to a PyPSA-Earth network region.

Point-level potentials are aggregated by network region. The resulting regional
potential table contains the maximum installable EGS capacity and CAPEX used by
the model. If air temperature input data and a corresponding output path are
available, the script additionally builds temperature-dependent EGS
capacity-factor profiles.
"""
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from _helpers import (
    configure_logging,
    create_logger,
    mock_snakemake,
    two_2_three_digits_country,
)
from shapely.geometry import Polygon

logger = create_logger(__name__)


CRF = 0.09  # calculated with 0.08 interest rate, assumtion from Franzmann et al. (https://doi.org/10.1016/j.renene.2025.123199) for calculation of the LCOE's.
# To calculate the CAPEX based on the paper, this exact value must be used.


def read_network_regions(network_regions_file):
    """
    Read PyPSA-Earth network regions and index them by region name.

    Parameters
    ----------
    network_regions_file : str or path-like
        Path to the network regions file. The file must contain a ``name`` column
        and a geometry column.

    Returns
    -------
    geopandas.GeoDataFrame
        Network regions indexed by ``name`` and assigned to CRS EPSG:4326.

    Raises
    ------
    ValueError
        If the network regions file does not contain a ``name`` column.
    """
    regions = gpd.read_file(network_regions_file)

    if "name" not in regions.columns:
        raise ValueError(
            f"Network regions file '{network_regions_file}' does not contain a 'name' column."
        )

    regions = regions.set_index("name", drop=True).set_crs(
        epsg=4326, allow_override=True
    )
    return regions


def prepare_egs_data(egs_file, countries, network_regions_file):
    """
    Prepare region-aggregated enhanced geothermal system potentials.

    The function reads point-based EGS input data, filters it to the configured
    countries, keeps the Gringarten power and LCOE estimates, converts the power
    output from MW to GW, derives CAPEX from LCOE and annual generation, and assigns
    each point to a PyPSA-Earth network region. Potentials are then aggregated by
    region.

    Parameters
    ----------
    egs_file : str or path-like
        Path to the EGS input CSV file. The file must contain coordinates,
        ``gid1`` country or region identifiers, and the required EGS power and LCOE
        columns.
    countries : str or list-like
        Model countries to keep. Two-letter country codes are converted to
        three-letter country codes before filtering.
    network_regions_file : str or path-like
        Path to the PyPSA-Earth network regions file used to assign EGS points to
        model regions.

    Returns
    -------
    pandas.DataFrame
        Region-indexed dataframe with mean ``CAPEX``, summed ``PowerSust``, number
        of assigned EGS points in ``n_points``, and ``p_nom_max`` equal to
        ``PowerSust``.

    Raises
    ------
    ValueError
        If required input columns are missing, no rows remain after country
        filtering, no valid positive-output EGS rows remain, or no EGS points can be
        assigned to the provided network regions.
    """
    df = pd.read_csv(egs_file)

    required_columns = {
        "lon",
        "lat",
        "gid1",
        "Pout_Gringarten_MW",
        "LCOE_Gringarten_Eurct_per_kWh",
        "Pout_VolumeMethod_MW",
        "LCOE_VolumeMethod_Eurct_per_kWh",
        "Pout_Sustainable_MW",
        "LCOE_Sustainable_Eurct_per_kWh",
    }
    missing_columns = sorted(required_columns - set(df.columns))
    if missing_columns:
        raise ValueError(f"Missing required input columns: {missing_columns}")

    countries = two_2_three_digits_country(countries)
    if isinstance(countries, str):
        countries = [countries]

    countries = [str(c).upper() for c in countries]
    df = df.drop(
        columns=[
            "Pout_VolumeMethod_MW",
            "LCOE_VolumeMethod_Eurct_per_kWh",
            "Pout_Sustainable_MW",
            "LCOE_Sustainable_Eurct_per_kWh",
        ]
    )

    df["gid1"] = df["gid1"].astype(str).str.strip().str.upper()
    df["country_code"] = df["gid1"].str.extract(r"^([A-Z]{3})", expand=False)
    df = df.loc[df["country_code"].isin(countries)].copy()

    if df.empty:
        raise ValueError(
            f"No EGS rows left after filtering by configured countries: {countries}"
        )

    df = df.rename(
        columns={
            "lon": "Lon",
            "lat": "Lat",
            "Pout_Gringarten_MW": "PowerSust_MW",
            "LCOE_Gringarten_Eurct_per_kWh": "LCOE_Eur_per_kWh",
        }
    )

    df = df.dropna(subset=["Lon", "Lat", "PowerSust_MW", "LCOE_Eur_per_kWh"])
    df = df.loc[df["PowerSust_MW"] > 0.0].copy()

    if df.empty:
        raise ValueError(
            "No valid EGS rows left after dropping invalid or zero-output rows."
        )

    # EUR/kWh => EUR/MWh
    df["LCOE_EUR_per_MWh"] = df["LCOE_Eur_per_kWh"] * 1000.0

    # Energy produced per year in MWh
    df["Leistung_MWh"] = df["PowerSust_MW"] * 8760.0

    # Annualized CAPEX-like quantity
    # Assumption: 2 % of CAPEX as OPEX, -> CRF + 0.02
    capex = df["LCOE_EUR_per_MWh"] * df["Leistung_MWh"] / (+0.02)

    # PowerSust in GW
    df["PowerSust"] = df["PowerSust_MW"] / 1000.0  # MW -> GW

    # CAPEX in EUR/GW
    df["CAPEX"] = capex / df["PowerSust"]

    point_gdf = gpd.GeoDataFrame(
        df[["gid1", "country_code", "Lon", "Lat", "CAPEX", "PowerSust"]],
        geometry=gpd.points_from_xy(df["Lon"], df["Lat"]),
        crs="EPSG:4326",
    )

    regions = read_network_regions(network_regions_file).reset_index()[
        ["name", "geometry"]
    ]

    # Assign each point directly to a network region
    joined = gpd.sjoin(
        point_gdf,
        regions,
        how="left",
        predicate="intersects",
    ).rename(columns={"name": "region"})

    missing_regions = int(joined["region"].isna().sum())
    if missing_regions > 0:
        logger.warning(
            f"{missing_regions} EGS points could not be assigned to any network region and will be dropped."
        )

    joined = joined.dropna(subset=["region"]).copy()

    if joined.empty:
        raise ValueError(
            "No EGS points could be assigned to the provided network regions."
        )

    # ------------------------------------------------------------
    # Power-weighted CAPEX by region
    # ------------------------------------------------------------
    def weighted_average_capex(group):
        weights = group["PowerSust"].fillna(0.0)

        if weights.sum() <= 0:
            return group["CAPEX"].mean()

        return np.average(
            group["CAPEX"],
            weights=weights,
        )

    # Aggregate standard sizes first
    region_agg = (
        joined.groupby("region", as_index=True)
        .agg(
            PowerSust=("PowerSust", "sum"),
            n_points=("region", "size"),
        )
        .sort_index()
    )

    # Then add CAPEX on a performance-weighted basis
    weighted_capex = (
        joined.groupby("region").apply(weighted_average_capex).rename("CAPEX")
    )

    region_agg = region_agg.join(weighted_capex)

    region_agg.index.name = "region"
    region_agg["p_nom_max"] = region_agg["PowerSust"]

    return region_agg


def get_capacity_factors(network_regions_file, air_temperatures_file, snapshots=None):
    """
    Build temperature-dependent EGS capacity-factor profiles.

    Enhanced geothermal system performance is adjusted based on deviations of
    regional air temperature from the regional mean temperature. Lower-than-average
    air temperatures increase the capacity factor, while higher-than-average
    temperatures decrease it. The adjustment curve follows the PyPSA-Eur
    temperature-based EGS capacity-factor approach.

    Parameters
    ----------
    network_regions_file : str or path-like
        Path to the network regions file. The file must contain a ``name`` column
        matching the regional names in the air temperature dataset.
    air_temperatures_file : str or path-like
        Path to a NetCDF file containing air temperatures by network region and
        time. The dataset must provide a ``temperature`` variable and a ``name``
        dimension or coordinate matching the network region names.
    snapshots : pandas.DatetimeIndex, optional
        Snapshots for which capacity factors should be returned. If omitted, the
        time index from the air temperature dataset is used.

    Returns
    -------
    pandas.DataFrame
        Capacity-factor time series with snapshots as index and network region
        names as columns. The index name is ``snapshot``.
    """

    # these values are taken from the paper's
    # Supplementary Figure 20 from https://zenodo.org/records/7093330
    # and relate deviations of the ambient temperature from the year-average
    # ambient temperature to EGS capacity factors.
    delta_t = [-15, -10, -5, 0, 5, 10, 15, 20]
    cf = [1.17, 1.13, 1.07, 1.0, 0.925, 0.84, 0.75, 0.65]

    x = np.linspace(-15, 20, 200)
    y = np.interp(x, delta_t, cf)

    upper_x = np.linspace(20, 25, 50)
    m_upper = (y[-1] - y[-2]) / (x[-1] - x[-2])
    upper_y = upper_x * m_upper - x[-1] * m_upper + y[-1]

    lower_x = np.linspace(-20, -15, 50)
    m_lower = (y[1] - y[0]) / (x[1] - x[0])
    lower_y = lower_x * m_lower - x[0] * m_lower + y[0]

    x = np.hstack((lower_x, x, upper_x))
    y = np.hstack((lower_y, y, upper_y))

    network_regions = gpd.read_file(network_regions_file).set_crs(
        epsg=4326, allow_override=True
    )
    region_names = network_regions["name"]

    air_temp = xr.open_dataset(air_temperatures_file)
    if snapshots is None:
        snapshots = pd.DatetimeIndex(air_temp.indexes["time"])

    capacity_factors = pd.DataFrame(index=snapshots)

    for bus in region_names:
        temp = air_temp.sel(name=bus).to_dataframe()["temperature"]
        temp = temp.reindex(snapshots)
        capacity_factors[bus] = np.interp((temp - temp.mean()).values, x, y)

    capacity_factors.index.name = "snapshot"
    return capacity_factors


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_egs_potentials",
        )

    configure_logging(snakemake)

    egs_data = prepare_egs_data(
        egs_file=snakemake.input.egs_input,
        countries=snakemake.params.countries,
        network_regions_file=snakemake.input.regions,
    )

    egs_data["p_nom_max"] = egs_data["PowerSust"]

    egs_data[["p_nom_max", "CAPEX"]].to_csv(snakemake.output.egs_potentials)

    if hasattr(snakemake.input, "temp_air_total") and hasattr(
        snakemake.output, "egs_capacity_factors"
    ):
        snapshots = pd.date_range(
            freq="h",
            **snakemake.params.snapshots,
        )
        capacity_factors = get_capacity_factors(
            snakemake.input.regions,
            snakemake.input.temp_air_total,
            snapshots,
        )
        capacity_factors.to_csv(snakemake.output.egs_capacity_factors)
