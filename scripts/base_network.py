# SPDX-FileCopyrightText: : 2021 The PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
# coding: utf-8
"""
Creates the network topology from a OpenStreetMap

Relevant Settings
-----------------

.. code:: yaml

    snapshots:

    countries:

    electricity:
        voltages:

    lines:
        types:
        s_max_pu:
        under_construction:

    links:
        p_max_pu:
        under_construction:
        include_tyndp:

    transformers:
        x:
        s_nom:
        type:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`snapshots_cf`, :ref:`toplevel_cf`, :ref:`electricity_cf`, :ref:`load_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`transformers_cf`

Inputs
------



Outputs
-------

- ``networks/base.nc``

    .. image:: ../img/base.png
        :scale: 33 %

Description
-----------

"""
import logging
import os

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import scipy as sp
import shapely.prepared
import shapely.wkt
import yaml
from _helpers import _read_csv_nafix
from _helpers import configure_logging
from download_osm_data import create_country_list
from scipy.sparse import csgraph
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.ops import unary_union

logger = logging.getLogger(__name__)


def _get_oid(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"oid"=>"(\\d+)"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _get_country(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"country"=>"([A-Z]{2})"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _find_closest_links(links, new_links, distance_upper_bound=1.5):
    treecoords = np.asarray([
        np.asarray(shapely.wkt.loads(s))[[0, -1]].flatten()
        for s in links.geometry
    ])
    querycoords = np.vstack([
        new_links[["x1", "y1", "x2", "y2"]],
        new_links[["x2", "y2", "x1", "y1"]]
    ])
    tree = sp.spatial.KDTree(treecoords)
    dist, ind = tree.query(querycoords,
                           distance_upper_bound=distance_upper_bound)
    found_b = ind < len(links)
    found_i = np.arange(len(new_links) * 2)[found_b] % len(new_links)

    return (pd.DataFrame(
        dict(D=dist[found_b], i=links.index[ind[found_b] % len(links)]),
        index=new_links.index[found_i],
    ).sort_values(by="D")
        [lambda ds: ~ds.index.duplicated(keep="first")].sort_index()["i"])


def _load_buses_from_osm():
    buses = (_read_csv_nafix(
        snakemake.input.osm_buses).set_index("bus_id").drop(
            ["station_id"], axis=1).rename(columns=dict(voltage="v_nom")))

    buses = buses.loc[:, ~buses.columns.str.contains("^Unnamed")]
    buses["v_nom"] /= 1e3
    buses["carrier"] = buses.pop("dc").map({True: "DC", False: "AC"})
    buses["under_construction"] = buses["under_construction"].fillna(
        False).astype(bool)
    buses["x"] = buses["lon"]
    buses["y"] = buses["lat"]
    # TODO: Drop NAN maybe somewhere else?
    buses = buses.dropna(axis="index", subset=["x", "y", "country"])

    # Rebase all voltages to three levels
    buses = _rebase_voltage_to_config(buses)

    logger.info("Removing buses with voltages {}".format(
        pd.Index(buses.v_nom.unique()).dropna().difference(
            snakemake.config["electricity"]["voltages"])))

    return buses


def _load_lines_from_osm(buses):
    lines = (_read_csv_nafix(
        snakemake.input.osm_lines,
        dtype=dict(
            line_id="str",
            bus0="str",
            bus1="str",
            underground="bool",
            under_construction="bool",
        ),
    ).set_index("line_id").rename(
        columns=dict(voltage="v_nom", circuits="num_parallel")))

    lines["length"] /= 1e3  # m to km conversion
    lines["v_nom"] /= 1e3  # V to kV conversion
    lines = lines.loc[:, ~lines.columns.str.contains(
        "^Unnamed")]  # remove unnamed col
    lines = _rebase_voltage_to_config(lines)  # rebase voltage to config inputs
    # lines = _remove_dangling_branches(lines, buses)  # TODO: add dangling branch removal?

    return lines


def _set_electrical_parameters_lines(lines):
    v_noms = snakemake.config["electricity"]["voltages"]
    linetypes = snakemake.config["lines"]["types"]

    for v_nom in v_noms:
        lines.loc[lines["v_nom"] == v_nom, "type"] = linetypes[v_nom]

    lines["s_max_pu"] = snakemake.config["lines"]["s_max_pu"]

    return lines


def _set_lines_s_nom_from_linetypes(n):
    # Info: n.line_types is a lineregister from pypsa/pandapowers
    n.lines["s_nom"] = (np.sqrt(3) * n.lines["type"].map(n.line_types.i_nom) *
                        n.lines["v_nom"] * n.lines.num_parallel)


def _remove_dangling_branches(branches, buses):
    return pd.DataFrame(branches.loc[branches.bus0.isin(buses.index)
                                     & branches.bus1.isin(buses.index)])


def _set_countries_and_substations(n):

    buses = n.buses

    countries = create_country_list(snakemake.config["countries"])
    country_shapes = (gpd.read_file(snakemake.input.country_shapes).set_index(
        "name")["geometry"].set_crs(4326))
    offshore_shapes = unary_union(
        gpd.read_file(
            snakemake.input.offshore_shapes)["geometry"].set_crs(4326))

    bus_locations = buses
    bus_locations = gpd.GeoDataFrame(
        bus_locations,
        geometry=gpd.points_from_xy(bus_locations.x, bus_locations.y),
        crs=4326,
    )
    # Check if bus is in shape
    offshore_b = bus_locations.within(offshore_shapes)

    # Assumption that HV-bus qualifies as potential offshore bus. Offshore bus is empty otherwise.
    offshore_hvb = (
        buses["v_nom"] >=
        snakemake.config["base_network"]["min_voltage_substation_offshore"] /
        1000)
    # Compares two lists & makes list value true if at least one is true
    buses["substation_off"] = offshore_b | offshore_hvb

    # Busses without country tag are removed OR get a country tag if close to country
    c_nan_b = buses.country.isnull()
    if c_nan_b.sum() > 0:
        c_tag = _get_country(buses.loc[c_nan_b])
        c_tag.loc[~c_tag.isin(countries)] = np.nan
        n.buses.loc[c_nan_b, "country"] = c_tag

        c_tag_nan_b = n.buses.country.isnull()

        # Nearest country in path length defines country of still homeless buses
        # Work-around until commit 705119 lands in pypsa release
        # pypsa-africa comment: Important to connect 'homeless' offshore assets
        # Otherwise
        n.transformers["length"] = 0.0
        graph = n.graph(weight="length")
        n.transformers.drop("length", axis=1, inplace=True)

        for b in n.buses.index[c_tag_nan_b]:
            df = (pd.DataFrame(
                dict(pathlength=nx.single_source_dijkstra_path_length(
                    graph, b, cutoff=200))).join(n.buses.country).dropna())
            assert (
                not df.empty
            ), "No buses with defined country within 200km of bus `{}`".format(
                b)
            n.buses.at[b, "country"] = df.loc[df.pathlength.idxmin(),
                                              "country"]

        logger.warning(
            "{} buses are not in any country or offshore shape,"
            " {} have been assigned from the tag of the entsoe map,"
            " the rest from the next bus in terms of pathlength.".format(
                c_nan_b.sum(),
                c_nan_b.sum() - c_tag_nan_b.sum()))

    return buses


def _rebase_voltage_to_config(component):
    """
    Rebase the voltage of components to the config.yaml input

    Components such as line and buses have voltage levels between
    110 kV up to around 850 kV. PyPSA-Africa uses 3 voltages as config input.
    This function rebases all inputs to the lower, middle and upper voltage
    bound.

    Parameters
    ----------
    component : dataframe
    """
    v_min = (snakemake.config["base_network"]["min_voltage_rebase_voltage"] /
             1000)  # min. filtered value in dataset
    v_low = snakemake.config["electricity"]["voltages"][0]
    v_mid = snakemake.config["electricity"]["voltages"][1]
    v_up = snakemake.config["electricity"]["voltages"][2]
    v_low_mid = (v_mid - v_low) / 2 + v_low  # between low and mid voltage
    v_mid_up = (v_up - v_mid) / 2 + v_mid  # between mid and upper voltage
    component.loc[(v_min <= component["v_nom"]) &
                  (component["v_nom"] < v_low_mid), "v_nom"] = v_low
    component.loc[(v_low_mid <= component["v_nom"]) &
                  (component["v_nom"] < v_mid_up), "v_nom"] = v_mid
    component.loc[v_mid_up <= component["v_nom"], "v_nom"] = v_up

    return component


def base_network():
    buses = _load_buses_from_osm().reset_index()
    lines = _load_lines_from_osm(buses)
    lines = _set_electrical_parameters_lines(lines)

    n = pypsa.Network()
    n.name = "PyPSA-Eur"

    n.set_snapshots(pd.date_range(freq="h", **snakemake.config["snapshots"]))
    n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

    n.import_components_from_dataframe(buses, "Bus")
    n.import_components_from_dataframe(lines, "Line")

    _set_lines_s_nom_from_linetypes(n)

    _set_countries_and_substations(n)

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake("base_network")
    configure_logging(snakemake)

    n = base_network()
    n.buses = pd.DataFrame(n.buses.drop(columns="geometry"))
    n.export_to_netcdf(snakemake.output[0])
