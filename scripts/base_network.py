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


def _get_oid(df: pd.DataFrame) -> pd.Series:
    """Extract the OpenStreetMap object ID (oid) from the tags column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame that may contain a ``tags`` column with OSM tag strings.

    Returns
    -------
    pd.Series
        Series of oid strings aligned to *df*'s index, or NaN where absent.
    """
    if "tags" in df.columns:
        return df.tags.str.extract('"oid"=>"(\\d+)"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def get_country(df: pd.DataFrame) -> pd.Series:
    """Extract the two-letter country code from the tags column.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame that may contain a ``tags`` column with OSM tag strings.

    Returns
    -------
    pd.Series
        Series of ISO 3166-1 alpha-2 country codes, or NaN where absent.
    """
    if "tags" in df.columns:
        return df.tags.str.extract('"country"=>"([A-Z]{2})"', expand=False)
    else:
        return pd.Series(np.nan, df.index)


def _find_closest_links(
    links: pd.DataFrame,
    new_links: pd.DataFrame,
    distance_upper_bound: float = 1.5,
) -> pd.Series:
    """Match each entry in *new_links* to the closest existing link via a KDTree.

    Endpoints are encoded as (x1, y1, x2, y2) coordinate tuples and queried in
    both forward and reverse direction so that link orientation does not matter.

    Parameters
    ----------
    links : pd.DataFrame
        Existing links with a WKT ``geometry`` column.
    new_links : pd.DataFrame
        Candidate links with columns ``x1``, ``y1``, ``x2``, ``y2``.
    distance_upper_bound : float, optional
        Maximum allowed endpoint distance (degrees) for a match, by default 1.5.

    Returns
    -------
    pd.Series
        Mapping from *new_links* index to the matched index in *links*.
    """
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


def _load_buses_from_osm(fp_buses: str) -> pd.DataFrame:
    """Load and normalise bus data from an OSM-derived CSV file.

    Converts voltage from V to kV, maps the ``dc`` boolean to a ``carrier``
    string (``"DC"`` / ``"AC"``), and drops rows with missing coordinates or
    country.

    Parameters
    ----------
    fp_buses : str
        Path to the CSV file produced by the OSM bus extraction step.

    Returns
    -------
    pd.DataFrame
        Cleaned bus table indexed by ``bus_id``.
    """
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


def add_underwater_links(n: "pypsa.Network", fp_offshore_shapes: str) -> None:
    """Compute and attach the underwater fraction for each HVDC link.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``links`` component is updated in-place.
    fp_offshore_shapes : str
        Path to the offshore shapes file (GeoPackage or shapefile).
    """
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


def _set_dc_underwater_fraction(
    lines_or_links: pd.DataFrame, fp_offshore_shapes: str
) -> None:
    """Set the ``underwater_fraction`` column on DC lines or links in-place.

    Skips computation when the DataFrame is empty, contains no DC carriers, or
    has no geometry column.  The fraction is the ratio of each branch's length
    that intersects the offshore polygon.

    Parameters
    ----------
    lines_or_links : pd.DataFrame
        Network component table (lines or links) modified in-place.
    fp_offshore_shapes : str
        Path to the offshore shapes file used to determine underwater sections.
    """
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


def _load_lines_from_osm(fp_osm_lines: str) -> pd.DataFrame:
    """Load and normalise AC/DC line data from an OSM-derived CSV file.

    Converts length from metres to kilometres and voltage from volts to
    kilovolts.  Unnamed columns produced by the CSV exporter are dropped.

    Parameters
    ----------
    fp_osm_lines : str
        Path to the CSV file produced by the OSM line extraction step.

    Returns
    -------
    pd.DataFrame
        Cleaned line table indexed by ``line_id``.
    """
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
def _load_links_from_osm(
    fp_osm_converters: str,
    base_network_config: dict,
    voltages_config: list,
) -> pd.DataFrame:
    """Load HVDC link data from an OSM-derived CSV file.

    Returns an empty DataFrame when the file is empty.  Converts length and
    voltage units and drops unnamed columns.

    .. deprecated::
        This function appears to be superseded by
        :func:`_load_converters_from_osm` and may be removed in a future
        release.

    Parameters
    ----------
    fp_osm_converters : str
        Path to the OSM converter/link CSV file.
    base_network_config : dict
        Base network configuration block from ``config.yaml``.
    voltages_config : list
        List of nominal voltages to consider.

    Returns
    -------
    pd.DataFrame
        Link table indexed by ``line_id``, or empty DataFrame if file is empty.
    """
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


def _load_converters_from_osm(
    fp_osm_converters: str, buses: pd.DataFrame
) -> pd.DataFrame:
    """Load back-to-back (B2B) converter data from an OSM-derived CSV file.

    Returns an empty DataFrame when the file is empty.  Converters are marked
    with ``carrier="B2B"`` and ``dc=True``.

    Parameters
    ----------
    fp_osm_converters : str
        Path to the OSM converter CSV file.
    buses : pd.DataFrame
        Bus table used for future dangling-branch removal (currently unused).

    Returns
    -------
    pd.DataFrame
        Converter table indexed by ``converter_id``, or empty DataFrame.
    """
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


def _load_transformers_from_osm(
    fp_osm_transformers: str, buses: pd.DataFrame
) -> pd.DataFrame:
    """Load transformer data from an OSM-derived CSV file.

    Parameters
    ----------
    fp_osm_transformers : str
        Path to the OSM transformer CSV file.
    buses : pd.DataFrame
        Bus table used for future dangling-branch removal (currently unused).

    Returns
    -------
    pd.DataFrame
        Transformer table indexed by ``transformer_id``.
    """
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


def _get_linetypes_config(line_types, voltages):
    """
    Return the dictionary of linetypes for selected voltages. The dictionary is
    a subset of the dictionary line_types, whose keys match the selected
    voltages.

    Parameters
    ----------
    line_types : dict
        Dictionary of linetypes: keys are nominal voltages and values are linetypes.
    voltages : list
        List of selected voltages.

    Returns
    -------
        Dictionary of linetypes for selected voltages.
    """
    # get voltages value that are not available in the line types
    vnoms_diff = set(voltages).symmetric_difference(set(line_types.keys()))
    if vnoms_diff:
        logger.warning(
            f"Voltages {vnoms_diff} not in the {line_types} or {voltages} list."
        )
    return {k: v for k, v in line_types.items() if k in voltages}


def _get_linetype_by_voltage(v_nom, d_linetypes):
    """
    Return the linetype of a specific line based on its voltage v_nom.

    Parameters
    ----------
    v_nom : float
        The voltage of the line.
    d_linetypes : dict
        Dictionary of linetypes: keys are nominal voltages and values are linetypes.

    Returns
    -------
        The linetype of the line whose nominal voltage is closest to the line voltage.
    """
    v_nom_min, line_type_min = min(
        d_linetypes.items(),
        key=lambda x: abs(x[0] - v_nom),
    )
    return line_type_min


def _set_electrical_parameters_lines(
    lines_config: dict, voltages: list, lines: pd.DataFrame
) -> pd.DataFrame:
    """Assign AC electrical parameters (carrier, type, s_max_pu) to lines.

    Parameters
    ----------
    lines_config : dict
        ``lines`` section of ``config.yaml``, must contain ``ac_types`` and
        ``s_max_pu``.
    voltages : list
        Nominal voltages to include.
    lines : pd.DataFrame
        AC line table modified in-place and returned.

    Returns
    -------
    pd.DataFrame
        Lines with ``carrier``, ``dc``, ``type``, and ``s_max_pu`` columns set.
    """
    if lines.empty:
        lines["type"] = []
        return lines

    linetypes = _get_linetypes_config(lines_config["ac_types"], voltages)

    lines["carrier"] = "AC"
    lines["dc"] = False

    lines.loc[:, "type"] = lines.v_nom.apply(
        lambda x: _get_linetype_by_voltage(x, linetypes)
    )

    lines["s_max_pu"] = lines_config["s_max_pu"]

    return lines


def _set_electrical_parameters_dc_lines(
    lines_config: dict, voltages: list, lines: pd.DataFrame
) -> pd.DataFrame:
    """Assign DC electrical parameters (carrier, type, s_max_pu) to lines.

    Parameters
    ----------
    lines_config : dict
        ``lines`` section of ``config.yaml``, must contain ``dc_types`` and
        ``s_max_pu``.
    voltages : list
        Nominal voltages to include.
    lines : pd.DataFrame
        DC line table modified in-place and returned.

    Returns
    -------
    pd.DataFrame
        Lines with ``carrier``, ``dc``, ``type``, and ``s_max_pu`` columns set.
    """
    if lines.empty:
        lines["type"] = []
        return lines

    linetypes = _get_linetypes_config(lines_config["dc_types"], voltages)

    lines["carrier"] = "DC"
    lines["dc"] = True
    lines.loc[:, "type"] = lines.v_nom.apply(
        lambda x: _get_linetype_by_voltage(x, linetypes)
    )

    lines["s_max_pu"] = lines_config["s_max_pu"]

    return lines


def _set_electrical_parameters_links(
    links_config: dict, links: pd.DataFrame
) -> pd.DataFrame:
    """Assign DC link parameters (p_max_pu, p_min_pu, carrier) to HVDC links.

    Parameters
    ----------
    links_config : dict
        ``links`` section of ``config.yaml``; ``p_max_pu`` defaults to 1.0.
    links : pd.DataFrame
        HVDC link table modified in-place and returned.

    Returns
    -------
    pd.DataFrame
        Links with ``p_max_pu``, ``p_min_pu``, ``carrier``, and ``dc`` set.
    """
    if links.empty:
        return links

    p_max_pu = links_config.get("p_max_pu", 1.0)
    links["p_max_pu"] = p_max_pu
    links["p_min_pu"] = -p_max_pu

    links["carrier"] = "DC"
    links["dc"] = True

    return links


def _set_electrical_parameters_transformers(
    transformers_config: dict, transformers: pd.DataFrame
) -> pd.DataFrame:
    """Assign electrical parameters (x, s_nom, type) to transformers.

    Parameters
    ----------
    transformers_config : dict
        ``transformers`` section of ``config.yaml``; keys ``x``, ``s_nom``,
        and ``type`` default to 0.1, 2000, and ``""`` respectively.
    transformers : pd.DataFrame
        Transformer table modified in-place and returned.

    Returns
    -------
    pd.DataFrame
        Transformers with ``x``, ``s_nom``, and ``type`` columns set.
    """
    config = transformers_config

    # Add transformer parameters
    transformers["x"] = config.get("x", 0.1)
    transformers["s_nom"] = config.get("s_nom", 2000)
    transformers["type"] = config.get("type", "")

    return transformers


def _set_electrical_parameters_converters(
    links_config: dict, converters: pd.DataFrame
) -> pd.DataFrame:
    """Assign parameters to back-to-back converters before merging with links.

    Sets ``p_max_pu``, ``p_min_pu``, a fixed ``p_nom`` of 2000 MW, and marks
    converters as neither under construction nor underground.

    Parameters
    ----------
    links_config : dict
        ``links`` section of ``config.yaml``; ``p_max_pu`` defaults to 1.0.
    converters : pd.DataFrame
        Converter table modified in-place and returned.

    Returns
    -------
    pd.DataFrame
        Converters ready for import into the PyPSA network as links.
    """
    p_max_pu = links_config.get("p_max_pu", 1.0)
    converters["p_max_pu"] = p_max_pu
    converters["p_min_pu"] = -p_max_pu

    converters["p_nom"] = 2000  # [MW]?

    # Converters are combined with links
    converters["under_construction"] = False
    converters["underground"] = False

    return converters


def _set_lines_s_nom_from_linetypes(n: "pypsa.Network") -> None:
    """Compute s_nom for all lines from their linetype current rating.

    Uses the PyPSA line type register (``n.line_types.i_nom``) to derive the
    apparent power rating.  DC lines use a single-phase formula; AC lines use
    the three-phase formula (√3 · i_nom · v_nom · num_parallel).

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``lines`` component is updated in-place.
    """
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


def _remove_dangling_branches(
    branches: pd.DataFrame, buses: pd.DataFrame
) -> pd.DataFrame:
    """Drop branches whose bus0 or bus1 is not present in *buses*.

    Parameters
    ----------
    branches : pd.DataFrame
        Line, link, transformer, or converter table with ``bus0`` and ``bus1``
        columns.
    buses : pd.DataFrame
        Bus table whose index defines the set of valid bus IDs.

    Returns
    -------
    pd.DataFrame
        Subset of *branches* where both endpoints exist in *buses*.
    """
    return pd.DataFrame(
        branches.loc[branches.bus0.isin(buses.index) & branches.bus1.isin(buses.index)]
    )


def _set_countries_and_substations(
    inputs: object,
    base_network_config: dict,
    countries_config: list,
    n: "pypsa.Network",
) -> pd.DataFrame:
    """Assign country tags and offshore substation flags to buses.

    Buses without a country tag are assigned one from the OSM tag if available,
    then by nearest-neighbour path length in the network graph.  Buses above
    the configured offshore voltage threshold that lie within offshore shapes
    are flagged as ``substation_off``.

    Parameters
    ----------
    inputs : snakemake.io.Namedlist
        Snakemake input object providing ``country_shapes`` and
        ``offshore_shapes`` file paths.
    base_network_config : dict
        Base network config block; must contain
        ``min_voltage_substation_offshore``.
    countries_config : list
        List of valid ISO 3166-1 alpha-2 country codes.
    n : pypsa.Network
        Network whose ``buses`` component is updated in-place.

    Returns
    -------
    pd.DataFrame
        Updated bus DataFrame (same object as ``n.buses``).
    """
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
    inputs: object,
    base_network_config: dict,
    countries_config: list,
    hvdc_as_lines_config: bool,
    lines_config: dict,
    links_config: dict,
    snapshots_config: dict,
    transformers_config: dict,
    voltages_config: list,
) -> "pypsa.Network":
    """Build the base PyPSA network from OSM-derived component tables.

    Loads buses, lines, transformers, and converters from CSV files, assigns
    electrical parameters, imports all components into a new
    :class:`pypsa.Network`, and returns it ready for further processing.

    When *hvdc_as_lines_config* is ``True``, DC lines are imported as PyPSA
    ``Line`` components; otherwise they are modelled as ``Link`` components.

    Parameters
    ----------
    inputs : snakemake.io.Namedlist
        Snakemake input paths: ``osm_buses``, ``osm_lines``,
        ``osm_transformers``, ``osm_converters``, ``country_shapes``,
        ``offshore_shapes``.
    base_network_config : dict
        ``base_network`` section of ``config.yaml``.
    countries_config : list
        List of ISO 3166-1 alpha-2 country codes to include.
    hvdc_as_lines_config : bool
        If ``True``, treat HVDC lines as PyPSA Line components.
    lines_config : dict
        ``lines`` section of ``config.yaml``.
    links_config : dict
        ``links`` section of ``config.yaml``.
    snapshots_config : dict
        ``snapshots`` section passed to :meth:`pypsa.Network.set_snapshots`.
    transformers_config : dict
        ``transformers`` section of ``config.yaml``.
    voltages_config : list
        Nominal voltages (kV) to retain.

    Returns
    -------
    pypsa.Network
        Fully assembled base network with buses, lines, links, transformers,
        converters, country assignments, and underwater fractions.
    """
    buses = _load_buses_from_osm(inputs.osm_buses).reset_index(drop=True)
    lines = _load_lines_from_osm(inputs.osm_lines).reset_index(drop=True)
    transformers = _load_transformers_from_osm(inputs.osm_transformers, buses)
    converters = _load_converters_from_osm(inputs.osm_converters, buses)

    lines_ac = lines[~lines.dc].copy()
    lines_dc = lines[lines.dc].copy()
    lines_ac = _set_electrical_parameters_lines(lines_config, voltages_config, lines_ac)

    lines_dc = _set_electrical_parameters_dc_lines(
        lines_config, voltages_config, lines_dc
    )

    transformers = _set_electrical_parameters_transformers(
        transformers_config, transformers
    )
    converters = _set_electrical_parameters_converters(links_config, converters)

    n = pypsa.Network()
    n.name = "PyPSA-Earth"

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
