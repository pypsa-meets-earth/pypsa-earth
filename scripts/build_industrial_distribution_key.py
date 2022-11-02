# -*- coding: utf-8 -*-
"""Build industrial distribution keys from hotmaps database."""

import uuid
from distutils.version import StrictVersion
from itertools import product

import geopandas as gpd
import pandas as pd
from helpers import locate_bus, three_2_two_digits_country
from shapely.geometry import Point

gpd_version = StrictVersion(gpd.__version__)


def map_industry_to_buses(df):
    """
    Load hotmaps database of industrial sites and map onto bus regions.
    Build industrial demand... Change name and add other functions.
    Function similar to aviation/shipping. Use functions to disaggregate.
    Only cement not steel - proof of concept.
    Change hotmaps to more descriptive name, etc.
    """
    df = df[df.country.isin(snakemake.config["countries"])]
    df["gadm_{}".format(gadm_level)] = df[["x", "y", "country"]].apply(
        lambda site: locate_bus(
            site[["x", "y"]].astype("float"),
            site["country"],
            gadm_level,
            shapes_path,
            gadm_clustering,
        ),
        axis=1,
    )

    return df.set_index("gadm_" + str(snakemake.config["sector"]["gadm_level"]))



def build_nodal_distribution_key(
    industrial_database, regions
):  # returns percentage of co2 emissions
    """Build nodal distribution keys for each sector."""

    countries = regions["name"].str[:2].unique()

    keys = pd.DataFrame(index=regions.name, columns=technology, dtype=float)

    pop = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    pop["country"] = pop.index.str[:2]
    keys["population"] = pop["total"].values / pop["total"].sum()

    for tech, country in product(technology, countries):

        regions_ct = regions.name[regions.name.str.contains(country)]

        facilities = industrial_database.query(
            "country == @country and technology == @tech"
        )
        # TODO adapt for facilities with production values not emissions
        if not facilities.empty:
            indicator = facilities["capacity"]
            if indicator.sum() == 0:
                key = pd.Series(1 / len(facilities), facilities.index)
            else:
                # TODO BEWARE: this is a strong assumption
                # indicator = indicator.fillna(0)
                key = indicator / indicator.sum()
            key = (
                key.groupby(facilities.index).sum().reindex(regions_ct, fill_value=0.0)
            )
        else:
            key = keys.loc[regions_ct, "population"]

        keys.loc[regions_ct, tech] = key

    return keys


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_industrial_distribution_key",
            simpl="",
            clusters=4,
            demand="NZ",
            planning_horizons=2030,
        )

    options = snakemake.config["sector"]
    gadm_level = options["gadm_level"]

    regions = gpd.read_file(snakemake.input.regions_onshore)

    if regions["name"][0][
        :3
    ].isalpha():  # TODO clean later by changing all codes to 2 letters
        regions["name"] = regions["name"].apply(
            lambda name: three_2_two_digits_country(name[:3]) + name[3:]
        )

    geo_locs = pd.read_csv(
        snakemake.input.industrial_database, sep=",", header=0  # , index_col=0
    )

    gadm_clustering = snakemake.config["clustering_options"]["alternative_clustering"]

    geo_locs = geo_locs[geo_locs.quality != "nonexistent"]

    technology = geo_locs.technology.unique()

    shapes_path = snakemake.input.shapes_path

    industrial_database = map_industry_to_buses(
        geo_locs[geo_locs.quality != "unavailable"]
    )

    keys = build_nodal_distribution_key(industrial_database, regions)

    keys.to_csv(snakemake.output.industrial_distribution_key)
