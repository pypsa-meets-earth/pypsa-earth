# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Creates the network topology from a OpenStreetMap.

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
        p_nom_max:
        under_construction:

    transformers:
        x:
        s_nom:
        type:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`snapshots_cf`, :ref:`toplevel_cf`, :ref:`electricity_cf`, :ref:`load_options_cf`,
    :ref:`lines_cf`, :ref:`links_cf`, :ref:`transformers_cf`

Inputs
------



Outputs
-------

- ``networks/base.nc``

    .. image:: /img/base.png
        :width: 33 %

Description
-----------
"""
import os

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import pypsa
import scipy as sp
import shapely.prepared
import shapely.wkt
from _helpers import configure_logging, create_logger, read_csv_nafix
from shapely.ops import unary_union

logger = create_logger(__name__)


def _get_oid(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"oid"=>"(\\d+)"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def get_country(df):
    if "tags" in df.columns:
        return df.tags.str.extract('"country"=>"([A-Z]{2})"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _find_closest_links(links, new_links, distance_upper_bound=1.5):
    treecoords = np.asarray(
        [np.asarray(shapely.wkt.loads(s))[[0, -1]].flatten() for s in links.geometry]
    )
    querycoords = np.vstack(
        [new_links[["x1", "y1", "x2", "y2"]], new_links[["x2", "y2", "x1", "y1"]]]
    )
    tree = sp.spatial.KDTree(treecoords)
    dist, ind = tree.query(querycoords, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(links)
    found_i = np.arange(len(new_links) * 2)[found_b] % len(new_links)

    return (
        pd.DataFrame(
            dict(D=dist[found_b], i=links.index[ind[found_b] % len(links)]),
            index=new_links.index[found_i],
        )
        .sort_values(by="D")[lambda ds: ~ds.index.duplicated(keep="first")]
        .sort_index()["i"]
    )


def _load_buses_from_osm(fp_buses):
    buses = (
        read_csv_nafix(fp_buses, dtype=dict(bus_id="str", voltage="float"))
        .set_index("bus_id")
        .drop(["station_id"], axis=1)
        .rename(columns=dict(voltage="v_nom"))
    )

    buses = buses.loc[:, ~buses.columns.str.contains("^Unnamed")]
    buses["v_nom"] /= 1e3
    buses["carrier"] = buses.pop("dc").map({True: "DC", False: "AC"})
    buses["under_construction"] = buses["under_construction"].fillna(False).astype(bool)
    buses["x"] = buses["lon"]
    buses["y"] = buses["lat"]
    # TODO: Drop NAN maybe somewhere else?
    buses = buses.dropna(axis="index", subset=["x", "y", "country"])

    return buses


def add_underwater_links(n, fp_offshore_shapes):
    if not hasattr(n.links, "geometry"):
        n.links["underwater_fraction"] = 0.0
    else:
        offshore_shape = gpd.read_file(fp_offshore_shapes).unary_union
        if offshore_shape is None or offshore_shape.is_empty:
            n.links["underwater_fraction"] = 0.0
        else:
            links = gpd.GeoSeries(n.links.geometry.dropna().map(shapely.wkt.loads))
            n.links["underwater_fraction"] = (
                links.intersection(offshore_shape).length / links.length
            )


def _set_dc_underwater_fraction(lines_or_links, fp_offshore_shapes):
    # HVDC part always has some links as converters
    # excluding probably purely DC networks which are currently somewhat exotic
    if lines_or_links.empty:
        return

    if lines_or_links.loc[lines_or_links.carrier == "DC"].empty:
        # Add "underwater_fraction" both to lines and links
        lines_or_links["underwater_fraction"] = 0.0
        return

    if not hasattr(lines_or_links, "geometry"):
        lines_or_links["underwater_fraction"] = 0.0
    else:
        offshore_shape = gpd.read_file(fp_offshore_shapes).unary_union
        if offshore_shape is None or offshore_shape.is_empty:
            lines_or_links["underwater_fraction"] = 0.0
        else:
            branches = gpd.GeoSeries(
                lines_or_links.geometry.dropna().map(shapely.wkt.loads)
            )
            # fix to avoid NaN for links during augmentation
            if branches.empty:
                lines_or_links["underwater_fraction"] = 0
            else:
                lines_or_links["underwater_fraction"] = (
                    # TODO Check assumption that all underwater lines are DC
                    branches.intersection(offshore_shape).length
                    / branches.length
                )


def _load_lines_from_osm(fp_osm_lines):
    lines = (
        read_csv_nafix(
            fp_osm_lines,
            dtype=dict(
                line_id="str",
                bus0="str",
                bus1="str",
                underground="bool",
                under_construction="bool",
                voltage="float",
                circuits="float",
            ),
        )
        .set_index("line_id")
        .rename(columns=dict(voltage="v_nom", circuits="num_parallel"))
    )

    lines["length"] /= 1e3  # m to km conversion
    lines["v_nom"] /= 1e3  # V to kV conversion
    lines = lines.loc[:, ~lines.columns.str.contains("^Unnamed")]  # remove unnamed col
    # lines = _remove_dangling_branches(lines, buses)  # TODO: add dangling branch removal?

    return lines


# TODO Seems to be not needed anymore
def _load_links_from_osm(fp_osm_converters, base_network_config, voltages_config):
    # the links file can be empty
    if os.path.getsize(fp_osm_converters) == 0:
        links = pd.DataFrame()
        return links

    links = (
        read_csv_nafix(
            fp_osm_converters,
            dtype=dict(
                line_id="str",
                bus0="str",
                bus1="str",
                underground="bool",
                under_construction="bool",
            ),
        )
        .set_index("line_id")
        .rename(columns=dict(voltage="v_nom", circuits="num_parallel"))
    )

    links["length"] /= 1e3  # m to km conversion
    links["v_nom"] /= 1e3  # V to kV conversion
    links = links.loc[:, ~links.columns.str.contains("^Unnamed")]  # remove unnamed col
    # links = _remove_dangling_branches(links, buses)  # TODO: add dangling branch removal?

    return links


def _load_converters_from_osm(fp_osm_converters, buses):
    # the links file can be empty
    if os.path.getsize(fp_osm_converters) == 0:
        converters = pd.DataFrame()
        return converters

    converters = read_csv_nafix(
        fp_osm_converters,
        dtype=dict(converter_id="str", bus0="str", bus1="str"),
    ).set_index("converter_id")

    # converters = _remove_dangling_branches(converters, buses)

    converters["carrier"] = "B2B"
    converters["dc"] = True

    return converters


def _load_transformers_from_osm(fp_osm_transformers, buses):
    transformers = (
        read_csv_nafix(
            fp_osm_transformers,
            dtype=dict(transformer_id="str", bus0="str", bus1="str"),
        )
        .rename(columns=dict(line_id="transformer_id"))
        .set_index("transformer_id")
    )
    # transformers = _remove_dangling_branches(transformers, buses)  # TODO: add dangling branch removal?

    return transformers


def _load_linetypes_from_csv(path):
    """
    Load voltage-to-line-type mappings from a CSV file.

    The CSV must use voltages as index and countries as columns. The
    ``default`` column reproduces the current PyPSA-Earth mapping.

    Parameters
    ----------
    path : str
        Path to the CSV file containing the line type mapping.

    Returns
    -------
    pandas.DataFrame
        Line type mapping with nominal voltages in kV as index and countries as columns.
    """
    linetypes = pd.read_csv(path, index_col=0)
    linetypes.index = linetypes.index.astype(float)

    if "default" not in linetypes.columns:
        raise ValueError(f"Missing 'default' column in line type mapping file: {path}")

    return linetypes


def _add_custom_line_types(n, path):
    """
    Add custom line-type definitions to the PyPSA line-type register.

    Parameters
    ----------
    n : pypsa.Network
        Network whose line-type register is extended.
    path : str
        CSV file containing PyPSA-compatible line-type parameters.
    """
    line_types = pd.read_csv(path).set_index("type")

    required_columns = [
        "f_nom",
        "r_per_length",
        "x_per_length",
        "c_per_length",
        "i_nom",
    ]

    missing_columns = set(required_columns) - set(line_types.columns)
    if missing_columns:
        raise ValueError(
            f"Custom line type file {path} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    duplicated_types = line_types.index[line_types.index.duplicated()].unique()
    if len(duplicated_types):
        raise ValueError(
            f"Duplicate custom line types in {path}: " f"{duplicated_types.tolist()}"
        )

    existing_types = line_types.index.intersection(n.line_types.index)
    if not existing_types.empty:
        raise ValueError(
            "Custom line types already exist in the PyPSA register: "
            f"{existing_types.tolist()}"
        )

    line_types = line_types[required_columns]
    n.line_types = pd.concat([n.line_types, line_types])

    logger.info("Added %s custom line types from %s.", len(line_types), path)


def _get_linetype_by_voltage_and_country(v_nom, country, linetypes):
    """
    Return the closest available line type for a voltage and country.

    Country-specific values override the default mapping at the corresponding
    voltage. Missing country-specific values fall back to the default mapping.
    """
    mapping = linetypes["default"]

    if country in linetypes.columns:
        mapping = linetypes[country].combine_first(mapping)

    mapping = mapping.dropna()

    if mapping.empty:
        raise ValueError(
            f"No line type mapping found for voltage {v_nom} kV "
            f"and country '{country}'."
        )

    voltage = min(mapping.index, key=lambda value: abs(value - v_nom))
    return mapping.at[voltage]


def _set_electrical_parameters_ac_lines(lines_config, buses, lines, line_type_mapping):
    if lines.empty:
        lines["type"] = []
        return lines

    linetypes = _load_linetypes_from_csv(line_type_mapping)

    lines["carrier"] = "AC"
    lines["dc"] = False
    lines["country"] = lines["bus0"].map(buses["country"])

    lines.loc[:, "type"] = lines.apply(
        lambda x: _get_linetype_by_voltage_and_country(x.v_nom, x.country, linetypes),
        axis=1,
    )

    lines["s_max_pu"] = lines_config["s_max_pu"]

    return lines


def _set_electrical_parameters_dc_lines(lines_config, buses, lines, line_type_mapping):
    if lines.empty:
        lines["type"] = []
        return lines

    linetypes = _load_linetypes_from_csv(line_type_mapping)

    lines["carrier"] = "DC"
    lines["dc"] = True
    lines["country"] = lines["bus0"].map(buses["country"])

    lines.loc[:, "type"] = lines.apply(
        lambda x: _get_linetype_by_voltage_and_country(x.v_nom, x.country, linetypes),
        axis=1,
    )

    lines["s_max_pu"] = lines_config["s_max_pu"]

    return lines


def _set_electrical_parameters_links(links_config, links):
    if links.empty:
        return links

    p_max_pu = links_config.get("p_max_pu", 1.0)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = -p_max_pu

    links["carrier"] = "DC"
    links["dc"] = True

    return links


def _set_electrical_parameters_transformers(transformers_config, transformers):
    config = transformers_config

    # Add transformer parameters
    transformers["x"] = config.get("x", 0.1)
    transformers["s_nom"] = config.get("s_nom", 2000)
    transformers["type"] = config.get("type", "")

    return transformers


def _set_electrical_parameters_converters(links_config, converters):
    p_max_pu = links_config.get("p_max_pu", 1.0)
    converters["p_max_pu"] = p_max_pu
    converters["p_min_pu"] = -p_max_pu

    converters["p_nom"] = 2000  # [MW]?

    # Converters are combined with links
    converters["under_construction"] = False
    converters["underground"] = False

    return converters


def _set_lines_s_nom_from_linetypes(n):
    # Info: n.line_types is a lineregister from pypsa/pandapowers
    n.lines["s_nom"] = (
        np.sqrt(3)
        * n.lines["type"].map(n.line_types.i_nom)
        * n.lines.eval("v_nom * num_parallel")
    )
    # Re-define s_nom for DC lines
    n.lines.loc[n.lines["carrier"] == "DC", "s_nom"] = n.lines["type"].map(
        n.line_types.i_nom
    ) * n.lines.eval("v_nom * num_parallel")


def _remove_dangling_branches(branches, buses):
    return pd.DataFrame(
        branches.loc[branches.bus0.isin(buses.index) & branches.bus1.isin(buses.index)]
    )


def _set_countries_and_substations(inputs, base_network_config, countries_config, n):
    countries = countries_config
    country_shapes = gpd.read_file(inputs.country_shapes).set_index("name")["geometry"]

    offshore_shapes = unary_union(gpd.read_file(inputs.offshore_shapes)["geometry"])

    buses = n.buses
    bus_locations = buses
    bus_locations = gpd.GeoDataFrame(
        bus_locations,
        geometry=gpd.points_from_xy(bus_locations.x, bus_locations.y),
        crs=country_shapes.crs,  # the workflow sets the the same crs for buses and shapes
    )
    # Check if bus is in shape
    offshore_b = bus_locations.within(offshore_shapes)

    # Assumption that HV-bus qualifies as potential offshore bus. Offshore bus is empty otherwise.
    offshore_hvb = (
        buses["v_nom"] >= base_network_config["min_voltage_substation_offshore"] / 1000
    )
    # Compares two lists & makes list value true if at least one is true
    buses["substation_off"] = offshore_b | offshore_hvb

    # Buses without country tag are removed OR get a country tag if close to country
    c_nan_b = buses.country.isnull()
    if c_nan_b.sum() > 0:
        c_tag = get_country(buses.loc[c_nan_b])
        c_tag.loc[~c_tag.isin(countries)] = np.nan
        n.buses.loc[c_nan_b, "country"] = c_tag

        c_tag_nan_b = n.buses.country.isnull()

        # Nearest country in path length defines country of still homeless buses
        # Work-around until commit 705119 lands in pypsa release
        # pypsa-earth comment: Important to connect 'homeless' offshore assets
        # Otherwise
        n.transformers["length"] = 0.0
        graph = n.graph(weight="length")
        n.transformers.drop("length", axis=1, inplace=True)

        for b in n.buses.index[c_tag_nan_b]:
            df = (
                pd.DataFrame(
                    dict(
                        pathlength=nx.single_source_dijkstra_path_length(
                            graph, b, cutoff=200
                        )
                    )
                )
                .join(n.buses.country)
                .dropna()
            )
            assert (
                not df.empty
            ), "No buses with defined country within 200km of bus `{}`".format(b)
            n.buses.at[b, "country"] = df.loc[df.pathlength.idxmin(), "country"]

        logger.warning(
            "{} buses are not in any country or offshore shape,"
            " {} have been assigned from the tag of the entsoe map,"
            " the rest from the next bus in terms of pathlength.".format(
                c_nan_b.sum(), c_nan_b.sum() - c_tag_nan_b.sum()
            )
        )

    return buses


def base_network(
    inputs,
    base_network_config,
    countries_config,
    hvdc_as_lines_config,
    lines_config,
    links_config,
    snapshots_config,
    transformers_config,
    voltages_config,
):
    buses = _load_buses_from_osm(inputs.osm_buses).reset_index(drop=True)
    lines = _load_lines_from_osm(inputs.osm_lines).reset_index(drop=True)
    transformers = _load_transformers_from_osm(inputs.osm_transformers, buses)
    converters = _load_converters_from_osm(inputs.osm_converters, buses)

    lines_ac = lines[~lines.dc].copy()
    lines_dc = lines[lines.dc].copy()

    lines_ac = _set_electrical_parameters_ac_lines(
        lines_config,
        buses,
        lines_ac,
        inputs.line_type_mapping_ac,
    )

    lines_dc = _set_electrical_parameters_dc_lines(
        lines_config,
        buses,
        lines_dc,
        inputs.line_type_mapping_dc,
    )

    transformers = _set_electrical_parameters_transformers(
        transformers_config, transformers
    )
    converters = _set_electrical_parameters_converters(links_config, converters)

    n = pypsa.Network()
    n.name = "PyPSA-Earth"

    _add_custom_line_types(n, inputs.custom_line_types)

    n.set_snapshots(pd.date_range(freq="h", **snapshots_config))
    n.snapshot_weightings[:] *= 8760.0 / n.snapshot_weightings.sum()

    n.import_components_from_dataframe(buses, "Bus")

    if hvdc_as_lines_config:
        lines = pd.concat([lines_ac, lines_dc])
        n.import_components_from_dataframe(lines, "Line")
    else:
        lines_dc = _set_electrical_parameters_links(links_config, lines_dc)
        # parse line information into p_nom required for converters
        lines_dc["p_nom"] = lines_dc.apply(
            lambda x: x["v_nom"] * n.line_types.i_nom[x["type"]],
            axis=1,
            result_type="reduce",
        )
        n.import_components_from_dataframe(lines_ac, "Line")
        n.import_components_from_dataframe(lines_dc, "Link")

    n.import_components_from_dataframe(transformers, "Transformer")
    n.import_components_from_dataframe(converters, "Link")

    # greenfield capacity expansion is represented with null capacity using num_parallel==0
    n.lines["num_parallel"] = n.lines["num_parallel"].where(
        ~n.lines["under_construction"], 0.0
    )
    n.lines.drop(columns="under_construction", inplace=True, errors="ignore")

    _set_lines_s_nom_from_linetypes(n)

    _set_countries_and_substations(inputs, base_network_config, countries_config, n)

    _set_dc_underwater_fraction(n.lines, inputs.offshore_shapes)
    _set_dc_underwater_fraction(n.links, inputs.offshore_shapes)

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("base_network")

    configure_logging(snakemake)

    inputs = snakemake.input

    # Snakemake imports:
    base_network_config = snakemake.params.base_network
    countries = snakemake.params.countries
    hvdc_as_lines = snakemake.params.hvdc_as_lines
    lines = snakemake.params.lines
    links = snakemake.params.links
    snapshots = snakemake.params.snapshots
    transformers = snakemake.params.transformers
    voltages = snakemake.params.voltages

    n = base_network(
        inputs,
        base_network_config,
        countries,
        hvdc_as_lines,
        lines,
        links,
        snapshots,
        transformers,
        voltages,
    )

    n.buses = pd.DataFrame(n.buses.drop(columns="geometry"))
    n.meta = snakemake.config
    n.export_to_netcdf(snakemake.output[0])
