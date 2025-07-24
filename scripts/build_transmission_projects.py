# SPDX-FileCopyrightText: PyPSA-ASEAN, PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


"""
Gets the transmission projects defined in the config file, concatenates and
deduplicates them. Projects are later included in :mod:`add_electricity.py`.

Inputs
------

- ``networks/base_network.nc``:  Base network topology for the electricity grid. This is processed in :mod:`base_network.py`.
- ``data/transmission_projects/"project_name"/``: Takes the transmission projects from the subfolder of data/transmission_projects. The subfolder name is the project name.
- ``offshore_shapes.geojson``: Shapefile containing the offshore regions. Used to determine if a new bus should be added for a new line or link.
- ``country_shapes.geojson``: Shapefile containing the shape of countries. Used to determine if a project is within the considered countries.

Outputs
-------

- ``transmission_projects/new_lines.csv``: New project lines to be added to the network. This includes new lines and upgraded lines.
- ``transmission_projects/new_links.csv``: New project links to be added to the network. This includes new links and upgraded links.
- ``transmission_projects/adjust_lines.csv``: For lines which are upgraded, the decommissioning year of the existing line is adjusted to the build year of the upgraded line.
- ``transmission_projects/adjust_links.csv``: For links which are upgraded, the decommissioning year of the existing link is adjusted to the build year of the upgraded link.
- ``transmission_projects/new_buses.csv``: For some links, we have to add new buses (e.g. North Sea Wind Power Hub).
"""

from pathlib import Path

import geopandas as gpd
from matplotlib import lines
import numpy as np
import pandas as pd
import pypsa
import shapely
from pypsa.descriptors import nominal_attrs
from scipy import spatial
from shapely.geometry import LineString, Point, Polygon, MultiPolygon
from _helpers import configure_logging, create_logger, read_csv_nafix

logger = create_logger(__name__)


def add_new_buses(n, new_ports):
    # Add new buses for the ports which do not have an existing bus close by. If there are multiple ports at the same location, only one bus is added.
    duplicated = new_ports.duplicated(subset=["x", "y"], keep="first")
    to_add = new_ports[~duplicated]
    added_buses = n.madd(
        "Bus",
        names=to_add.index,
        x=to_add.x,
        y=to_add.y,
        v_nom=380,
        under_construction=True,
        symbol="substation",
        substation_off=True,
        substation_lv=True,
        carrier="AC",
    )
    new_buses = n.buses.loc[added_buses].copy().dropna(axis=1, how="all")
    new_ports.loc[to_add.index, "neighbor"] = added_buses
    new_ports["neighbor"] = new_ports.groupby(["x", "y"])["neighbor"].transform("first")
    return new_ports, new_buses


def find_country_for_bus(bus, shapes):
    """
    Find the country of a bus based on its coordinates and the provided
    shapefile.

    Shapefile must contain a column "country" with the country names.
    """
    point = Point(bus.x, bus.y)
    country = shapes.loc[shapes.contains(point), "country"]
    return country.values[0]


def connect_new_lines(
    lines,
    n,
    new_buses_df, status,
    offshore_shapes=None,
    country_shapes=None,
    distance_upper_bound=np.inf,
    bus_carrier="AC"
):
    """
    Find the closest existing bus to the port of each line.

    If closest bus is further away than distance_upper_bound and is
    inside an offshore or onshore region, a new bus is created and the line is
    connected to it.
    """
    bus_carrier = np.atleast_1d(bus_carrier)
    buses = n.buses.query("carrier in @bus_carrier").copy()
    bus_tree = spatial.KDTree(buses[["x", "y"]])

    for port in [0, 1]:
        lines_port = lines["geometry"].apply(
            lambda x: pd.Series(
                get_bus_coords_from_port(x, port=port), index=["x", "y"]
            )
        )
        distances, indices = bus_tree.query(lines_port)
        lines_port["neighbor"] = buses.iloc[indices].index
        lines_port["match_distance"] = distances < distance_upper_bound

        # --- OFFSHORE NEW BUSES ---
        if not lines_port.match_distance.all() and offshore_shapes is not None:
            potential_new_buses = lines_port[~lines_port.match_distance]
            is_offshore = potential_new_buses.apply(
                lambda x: offshore_shapes.unary_union.contains(Point(x.x, x.y)), axis=1
            )
            new_buses = potential_new_buses[is_offshore]

            if not new_buses.empty:
                new_port, new_buses = add_new_buses(n, new_buses)

                new_buses["country"] = new_buses.apply(
                    lambda bus: find_country_for_bus(bus, offshore_shapes), axis=1
                )

                lines_port.loc[new_port.index, "match_distance"] = True
                lines_port.loc[new_port.index, "neighbor"] = new_port["neighbor"]
                new_buses_df = pd.concat([new_buses_df, new_buses])

        # --- ONSHORE NEW BUSES ---
        if not lines_port.match_distance.all() and country_shapes is not None:
            potential_new_buses = lines_port[~lines_port.match_distance]

            is_onshore = potential_new_buses.apply(
                lambda x: country_shapes.unary_union.contains(Point(x.x, x.y)),
                axis=1
            )
            new_buses = potential_new_buses[is_onshore]

            if not new_buses.empty:
                new_port, new_buses = add_new_buses(n, new_buses)

                new_buses["country"] = new_buses.apply(
                    lambda bus: find_country_for_bus(bus, country_shapes), axis=1
                )

                new_buses["tag_substation"] = "transmission"
                new_buses["tag_area"] = 0.0

                lines_port.loc[new_port.index, "match_distance"] = True
                lines_port.loc[new_port.index, "neighbor"] = new_port["neighbor"]

                new_buses_df = pd.concat([new_buses_df, new_buses])


        # Assign bus connections for this port
        lines.loc[lines_port.index, f"bus{port}"] = lines_port["neighbor"]

    lines["under_construction"] = lines["project_status"] != "existing"


    return lines, new_buses_df



def get_branch_coords_from_geometry(linestring, reversed=False):
    """
    Reduces a linestring to its start and end points. Used to simplify the
    linestring which can have more than two points.

    Parameters
    ----------
    linestring: Shapely linestring
    reversed (bool, optional): If True, returns the end and start points instead of the start and end points.
                               Defaults to False.

    Returns
    -------
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1] if not reversed else [-1, 0]
    start_end_coords = coords[ind]
    return start_end_coords.flatten()


def get_branch_coords_from_buses(line):
    """
    Gets line string for branch component in an pypsa network.

    Parameters
    ----------
    linestring: shapely linestring
    reversed (bool, optional): If True, returns the end and start points instead of the start and end points.
                               Defaults to False.

    Returns
    -------
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    start_coords = n.buses.loc[line.bus0, ["x", "y"]].values
    end_coords = n.buses.loc[line.bus1, ["x", "y"]].values
    return np.array([start_coords, end_coords]).flatten()


def get_bus_coords_from_port(linestring, port=0):
    """
    Extracts the coordinates of a specified port from a given linestring.

    Parameters
    ----------
    linestring: The shapely linestring.
    port (int): The index of the port to extract coordinates from. Default is 0.

    Returns
    -------
    tuple: The coordinates of the specified port as a tuple (x, y).
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1]
    coords = coords[ind]
    coords = coords[port]
    return coords


def find_closest_lines(lines, new_lines, distance_upper_bound=0.1, type="new"):
    """
    Find the closest lines in the existing set of lines to a set of new lines.

    Parameters
    ----------
    lines (pandas.DataFrame): DataFrame of the existing lines.
    new_lines (pandas.DataFrame): DataFrame with column geometry containing the new lines.
    distance_upper_bound (float, optional): Maximum distance to consider a line as a match. Defaults to 0.1 which corresponds to approximately 15 km.

    Returns
    -------
    pandas.Series: Series containing with index the new lines and values providing closest existing line.
    """

    # get coordinates of start and end points of all lines, for new lines we need to check both directions
    treelines = lines.apply(get_branch_coords_from_buses, axis=1)
    querylines = pd.concat(
        [
            new_lines["geometry"].apply(get_branch_coords_from_geometry),
            new_lines["geometry"].apply(get_branch_coords_from_geometry, reversed=True),
        ]
    )
    treelines = np.vstack(treelines)
    querylines = np.vstack(querylines)
    tree = spatial.KDTree(treelines)
    dist, ind = tree.query(querylines, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(lines)
    # since the new lines are checked in both directions, we need to find the correct index of the new line
    found_i = np.arange(len(querylines))[found_b] % len(new_lines)
    # create a DataFrame with the distances, new line and its closest existing line
    line_map = pd.DataFrame(
        dict(D=dist[found_b], existing_line=lines.index[ind[found_b] % len(lines)]),
        index=new_lines.index[found_i].rename("new_lines"),
    )
    if type == "new":
        if len(found_i) != 0:
            # compare if attribute of new line and existing line is similar
            attr = "p_nom" if "p_nom" in lines else "v_nom"
            # potential duplicates
            duplicated = line_map["existing_line"]
            # only if lines are similar in terms of p_nom or v_nom they are kept as duplicates
            to_keep = is_similar(
                new_lines.loc[duplicated.index, attr],
                duplicated.map(lines[attr]),
                percentage=10,
            )
            line_map = line_map[to_keep]
            if not line_map.empty:
                logger.warning(
                    "Found new lines similar to existing lines:\n"
                    + str(line_map["existing_line"].to_dict())
                    + "\n Lines are assumed to be duplicated and will be ignored."
                )
    elif type == "upgraded":
        if len(found_i) < len(new_lines):
            not_found = new_lines.index.difference(line_map.index)
            logger.warning(
                "Could not find upgraded lines close enough to existing lines:\n"
                + str(not_found.to_list())
                + "\n Lines will be ignored."
            )
    # only keep the closer line of the new line pair (since lines are checked in both directions)
    line_map = line_map.sort_values(by="D")[
        lambda ds: ~ds.index.duplicated(keep="first")
    ].sort_index()["existing_line"]
    return line_map


def adjust_decommissioning(upgraded_lines, line_map):
    """
    Adjust the decommissioning year of the existing lines to the built year of
    the upgraded lines.
    """
    to_update = pd.DataFrame(index=line_map)
    to_update["build_year"] = (
        1990  # dummy build_year to make existing lines decommissioned when upgraded lines are built
    )
    to_update["lifetime"] = (
        upgraded_lines.rename(line_map)["build_year"] - 1990
    )  # set lifetime to the difference between build year of upgraded line and existing line
    return to_update


def get_upgraded_lines(branch_component, n, upgraded_lines, line_map):
    """
    Get upgraded lines or links by merging info of existing and upgraded components.
    """

    # get first the information of the existing lines/links to be upgraded
    df_existing = n.df(branch_component).loc[line_map].copy()

    # get columns of upgraded lines/links which are not in existing lines
    new_columns = upgraded_lines.columns.difference(df_existing.columns)

    # rename upgraded lines to match existing lines
    upgraded_lines = upgraded_lines.rename(line_map)

    # set the same index names to be able to merge
    upgraded_lines.index.name = df_existing.index.name
    # Ensure type is preserved if missing from upgraded lines
    
    columns_to_update = upgraded_lines.columns.difference(["type"])
    df_existing.update(upgraded_lines[columns_to_update])


    # add columns which were new in upgraded_lines
    df_existing = pd.concat([df_existing, upgraded_lines[new_columns]], axis=1)
    if "type" not in df_existing.columns:
        df_existing["type"] = n.df(branch_component).loc[line_map, "type"].values
    else:
        # If type column exists but has nulls for some upgrades, fill from original
        mask_null_type = df_existing["type"].isna()
        if mask_null_type.any():
            df_existing.loc[mask_null_type, "type"] = (
                n.df(branch_component).loc[line_map[mask_null_type], "type"].values
            )

    # only consider columns of original upgraded lines/links and bus0 and bus1
    df_existing = df_existing.loc[:, ["bus0", "bus1", *upgraded_lines.columns]]

    # decide capacity attribute:
    capacity_attr = "s_nom" if branch_component == "Line" else "p_nom"

    for idx, row in df_existing.iterrows():
        original_idx = idx.replace("_upgraded", "")
        if original_idx in n.df(branch_component).index:
            cap_existing = n.df(branch_component).at[original_idx, capacity_attr]

            # defaults
            scaling_factor = 1.0

            if branch_component == "Line":
                v_nom_existing = n.df(branch_component).at[original_idx, "v_nom"]
                num_parallel_existing = n.df(branch_component).at[original_idx, "num_parallel"]

                # ensure no NaN or zero
                if pd.isna(v_nom_existing) or v_nom_existing == 0:
                    v_nom_existing = row.get("v_nom", 1.0)
                if pd.isna(num_parallel_existing) or num_parallel_existing == 0:
                    num_parallel_existing = 1.0

                scaling_factor = (
                    (row["v_nom"] / v_nom_existing)
                    * (row["num_parallel"] / num_parallel_existing)
                )

            elif branch_component == "Link":
                # only num_parallel applies (no voltage concept for links)
                num_parallel_existing = n.df(branch_component).at[original_idx, "num_parallel"] \
                    if "num_parallel" in n.df(branch_component).columns else 1.0
                if pd.isna(num_parallel_existing) or num_parallel_existing == 0:
                    num_parallel_existing = 1.0

                scaling_factor = row.get("num_parallel", 1.0) / num_parallel_existing

            new_capacity = cap_existing * scaling_factor
            new_capacity = min(new_capacity, 1000.0)

            df_existing.at[idx, capacity_attr] = round(new_capacity, 2)

        else:
            logger.warning(
                f"Original {branch_component} {original_idx} not found. "
                f"Leaving {capacity_attr} unchanged for upgraded {branch_component} {idx}."
            )

    # Ensure capacity attribute exists
    if capacity_attr not in df_existing.columns:
        df_existing[capacity_attr] = cap_existing

    if capacity_attr not in df_existing.columns:
        df_existing[capacity_attr] = cap_existing

    df_existing.index = df_existing.index.astype(str) + "_upgraded"

    return df_existing


    return df_existing



def get_project_files(path, skip=[]):
    path = Path(path)
    lines = {}
    files = [
        p
        for p in path.iterdir()
        if p.is_file()
        and p.suffix == ".csv"
        and not any(substring in p.name for substring in skip)
    ]
    if not files:
        logger.warning(f"No projects found for {path.parent.name}")
        return lines
    for file in files:
        df = pd.read_csv(file, index_col=0)
        df["geometry"] = df.apply(
            lambda x: LineString([[x.x0, x.y0], [x.x1, x.y1]]), axis=1
        )
        df.drop(columns=["x0", "y0", "x1", "y1"], inplace=True)
        lines[file.stem] = df
    return lines

def remove_projects_outside_countries(lines, region_shape):
    """
    Remove projects which are not in the considered countries.
    """
    region_shape_prepped = shapely.prepared.prep(region_shape)
    is_within_covered_countries = lines["geometry"].apply(
        lambda x: region_shape_prepped.contains(x)
    )

    if not is_within_covered_countries.all():
        logger.warning(
            "Project lines outside of the covered area (skipping): "
            + ", ".join(str(i) for i in lines.loc[~is_within_covered_countries].index)
        )

    lines = lines.loc[is_within_covered_countries]
    return lines


def is_similar(ds1, ds2, percentage=10):
    """
    Check if values in series ds2 are within a specified percentage of series
    ds1.

    Returns:
    - A boolean series where True indicates ds2 values are within the percentage range of ds2.
    """
    lower_bound = ds1 * (1 - percentage / 100)
    upper_bound = ds1 * (1 + percentage / 100)
    return np.logical_and(ds2 >= lower_bound, ds2 <= upper_bound)


def set_underwater_fraction(new_links, offshore_shapes):
    new_links_gds = gpd.GeoSeries(new_links["geometry"])
    new_links["underwater_fraction"] = (
        new_links_gds.intersection(offshore_shapes.union_all()).length
        / new_links_gds.length
    ).round(2)


def add_projects(
    n,
    new_lines_df,
    new_links_df,
    adjust_lines_df,
    adjust_links_df,
    new_buses_df,
    region_shape,
    country_shapes,
    offshore_shapes,
    path,
    plan,
    status=["confirmed", "under construction"],
    skip=[],
):
    lines_dict = get_project_files(path, skip=skip)
    for key, lines in lines_dict.items():
        logger.info(f"Processing {key.replace('_', ' ')} projects from {plan}.")
        lines = remove_projects_outside_countries(lines, region_shape)
        if isinstance(status, dict):
            status = status[plan]
        lines = lines.loc[lines.project_status.isin(status)]
        if lines.empty:
            continue
        if key == "new_lines":
            new_lines, new_buses_df = connect_new_lines(
                lines, 
                n, 
                new_buses_df, 
                status, 
                bus_carrier="AC", 
                country_shapes=country_shapes
            )
            duplicate_lines = find_closest_lines(
                n.lines, new_lines, distance_upper_bound=0.10, type="new"
            )
            new_lines = new_lines.drop(duplicate_lines.index, errors="ignore")
            new_lines_df = pd.concat([new_lines_df, new_lines])
            # add new lines to network to be able to find added duplicates
            #n.add("Line", new_lines.index, **new_lines)
            new_lines_df["dc"] = 0
            new_lines_df["underwater_fraction"] = 0.0 #only relevant for dc 
            n.madd("Line", new_lines.index, **new_lines.to_dict(orient="list"))
        elif key == "new_links":
            new_links, new_buses_df = connect_new_lines(
                lines,
                n,
                new_buses_df,
                status,
                offshore_shapes=offshore_shapes,
                distance_upper_bound=0.4,
                bus_carrier=["AC", "DC"],
                country_shapes=country_shapes
            )
            duplicate_links = find_closest_lines(
                n.links, new_links, distance_upper_bound=0.10, type="new"
            )
            new_links = new_links.drop(duplicate_links.index, errors="ignore")
            set_underwater_fraction(new_links, offshore_shapes)
            new_links_df = pd.concat([new_links_df, new_links])

            # Ensure bus columns are strings
            new_links["bus0"] = new_links["bus0"].astype(str)
            new_links["bus1"] = new_links["bus1"].astype(str)

            # CHECK FOR MISSING BUSES
            link_buses = pd.Series(
                pd.concat([new_links["bus0"], new_links["bus1"]])
                .unique(),
                name="bus"
            )
            missing_buses = link_buses[~link_buses.isin(n.buses.index)]

            if not missing_buses.empty:
                new_bus_rows = []
                for idx, row in new_links.iterrows():
                    for port, bus in zip([0, 1], [row["bus0"], row["bus1"]]):
                        if bus in missing_buses.values:
                            coords = get_bus_coords_from_port(row["geometry"], port)
                            new_bus_rows.append({
                                "bus": bus,
                                "x": coords[0],
                                "y": coords[1],
                                "v_nom": 380,
                                "under_construction": True,
                                "symbol": "substation",
                                "substation_off": True,
                                "substation_lv": True,
                                "carrier": "AC" if row["carrier"] == "AC" else "DC"
                            })

                if new_bus_rows:
                    new_buses = pd.DataFrame(new_bus_rows).set_index("bus")
                    new_buses.index = new_buses.index.map(str)
                    new_buses = new_buses[~new_buses.index.duplicated(keep="first")]
                    new_buses_df = pd.concat([new_buses_df, new_buses])
                    n.madd("Bus", new_buses.index, **new_buses.to_dict(orient="list"))
                    logger.info(f"Added missing buses to new_buses_df and network:\n{new_buses.index.tolist()}")

            # Add new links to the network
            n.madd("Link", new_links.index, **new_links.to_dict(orient="list"))


        elif key == "upgraded_lines":
            line_map = find_closest_lines(
                n.lines, lines, distance_upper_bound=0.30, type="upgraded"
            )
            upgraded_lines = lines.loc[line_map.index]
            lines_to_adjust = adjust_decommissioning(upgraded_lines, line_map)
            adjust_lines_df = pd.concat([adjust_lines_df, lines_to_adjust])
            upgraded_lines = get_upgraded_lines("Line", n, upgraded_lines, line_map)
            upgraded_lines["dc"] = 0  # ensure upgraded AC lines have dc=0
            upgraded_lines["underwater_fraction"] = 0.0 #only relevant for dc
            if "underground" in upgraded_lines.columns:
                upgraded_lines["underground"] = upgraded_lines["underground"].astype("boolean")
            upgraded_lines["under_construction"] = upgraded_lines["project_status"] != "existing"

            new_lines_df = pd.concat([new_lines_df, upgraded_lines])

        elif key == "upgraded_links":
            line_map = find_closest_lines(
                n.links.query("carrier=='DC'"),
                lines,
                distance_upper_bound=0.30,
                type="upgraded",
            )
            upgraded_links = lines.loc[line_map.index]
            links_to_adjust = adjust_decommissioning(upgraded_links, line_map)
            adjust_links_df = pd.concat([adjust_links_df, links_to_adjust])
            upgraded_links = get_upgraded_lines("Link", n, upgraded_links, line_map)
            new_links_df = pd.concat([new_links_df, upgraded_links])
            set_underwater_fraction(new_links_df, offshore_shapes)
        else:
            logger.warning(f"Unknown project type {key}")
            continue
    return new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df


def fill_length_from_geometry(line, line_factor=1.2):
    if not pd.isna(line.length):
        return line.length
    length = gpd.GeoSeries(line["geometry"], crs=4326).to_crs(3035).length.values[0]
    length = length / 1000 * line_factor
    return round(length, 1)

def remove_holes(geom):
    if geom.geom_type == 'Polygon':
        return Polygon(geom.exterior)
    elif geom.geom_type == 'MultiPolygon':
        return MultiPolygon([Polygon(p.exterior) for p in geom.geoms])
    else:
        return geom  # In case it's something else


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_transmission_projects")
    configure_logging(snakemake)

    line_factor = snakemake.params.line_factor
    s_max_pu = snakemake.params.s_max_pu

    n = pypsa.Network(snakemake.input.base_network)

    new_lines_df = pd.DataFrame()
    new_links_df = pd.DataFrame()
    adjust_lines_df = pd.DataFrame()
    adjust_links_df = pd.DataFrame()
    new_buses_df = pd.DataFrame()

    region_shape = remove_holes(gpd.read_file(snakemake.input.region_shape).geometry[0])
    country_shapes = gpd.read_file(snakemake.input.country_shapes).rename(columns={"name": "country"})

    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).rename(
        {"name": "country"}, axis=1
    )

    transmission_projects = snakemake.params.transmission_projects
    projects = [
        project
        for project, include in transmission_projects["include"].items()
        if include
    ]
    paths = snakemake.input.transmission_projects
    for project in projects:
        path = list(filter(lambda path: project in path, paths))[0]
        new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df = (
            add_projects(
                n,
                new_lines_df,
                new_links_df,
                adjust_lines_df,
                adjust_links_df,
                new_buses_df,
                region_shape,
                country_shapes,
                offshore_shapes,
                path=path,
                plan=project,
                status=transmission_projects["status"],
                skip=transmission_projects["skip"],
            )
        )
    if "underground" in adjust_lines_df.columns:
        adjust_lines_df["underground"] = (
            adjust_lines_df["underground"].astype("bool").fillna(False)
        )
       # Patch upgraded lines with s_nom still zero
    mask_upgraded_lines = new_lines_df.index.str.contains("_upgraded")

    if mask_upgraded_lines.any():
        # Identify upgraded lines where s_nom is missing or zero
        missing_capacity = new_lines_df.loc[mask_upgraded_lines, "s_nom"].fillna(0) == 0
        affected = new_lines_df.loc[mask_upgraded_lines][missing_capacity]

        for idx, row in affected.iterrows():
            original_idx = idx.replace("_upgraded", "")
            fallback_type = row.get("type", "N2XS(FL)2Y 1x240 RM/35 64/110 kV")
            fallback_v_nom = row.get("v_nom", 115)
            fallback_parallel = row.get("num_parallel", 1)

            # If type is missing or unknown, fallback
            if pd.isna(fallback_type) or fallback_type not in n.line_types.index:
                logger.warning(
                    f"Upgraded line {idx} has missing or unknown type. "
                    f"Falling back to default type N2XS(FL)2Y 1x240 RM/35 64/110 kV"
                )
                fallback_type = "N2XS(FL)2Y 1x240 RM/35 64/110 kV"

            # Use fallback current rating
            i_nom = n.line_types.at[fallback_type, "i_nom"]

            # Calculate fallback capacity
            fallback_s_nom = np.sqrt(3) * i_nom * fallback_v_nom * fallback_parallel
            fallback_s_nom = round(fallback_s_nom, 2)

            new_lines_df.at[idx, "s_nom"] = fallback_s_nom

            logger.warning(
                f"Fallback s_nom for upgraded line {idx}: "
                f"type={fallback_type}, v_nom={fallback_v_nom}, "
                f"num_parallel={fallback_parallel} → s_nom={fallback_s_nom}"
            )

    if not new_lines_df.empty:
        line_type = "Al/St 240/40 4-bundle 380.0"

        # Add new line type for new lines
        # Only fill missing line types for NEW lines
        mask_new_lines = ~new_lines_df.index.str.contains("_upgraded")
        new_lines_df.loc[mask_new_lines, "type"] = new_lines_df.loc[mask_new_lines, "type"].fillna(line_type)

        new_lines_df["num_parallel"] = new_lines_df["num_parallel"].fillna(2)

        if "underground" in new_lines_df.columns:
            new_lines_df["underground"] = (
                new_lines_df["underground"].astype("bool").fillna(False)
            )

        # Add carrier types of lines
        new_lines_df["carrier"] = "AC"

        # Fill empty length values with length calculated from geometry
        new_lines_df["length"] = new_lines_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # Only recalculate s_nom for NEW lines (not upgrades)
        mask_new_lines = ~new_lines_df.index.str.contains("_upgraded")

        if mask_new_lines.any():
            # Check if all new types exist in line_types
            types_new = new_lines_df.loc[mask_new_lines, "type"]
            missing_types = types_new[~types_new.isin(n.line_types.index)]
            if not missing_types.empty:
                logger.warning(
                    f"Line types missing in network line_types table: "
                    f"{missing_types.unique()}"
                )
                # fallback: set s_nom for these rows to 100
                new_lines_df.loc[missing_types.index, "s_nom"] = 100.0

            # Compute s_nom only for new lines where type is valid
            valid_mask = mask_new_lines & new_lines_df["type"].isin(n.line_types.index)

            if valid_mask.any():
                new_lines_df.loc[valid_mask, "s_nom"] = (
                    np.sqrt(3)
                    * n.line_types.loc[new_lines_df.loc[valid_mask, "type"], "i_nom"].values
                    * new_lines_df.loc[valid_mask, "v_nom"]
                    * new_lines_df.loc[valid_mask, "num_parallel"]
                ).round(2)
        # set s_max_pu
        new_lines_df["s_max_pu"] = s_max_pu
    if not new_links_df.empty:
        # Add carrier types of lines and links
        new_links_df["carrier"] = "DC"
        # Fill empty length values with length calculated from geometry
        new_links_df["length"] = new_links_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # Whether to keep existing link capacity or set to zero
        not_upgraded = ~new_links_df.index.str.contains("upgraded")
        if transmission_projects["new_link_capacity"] == "keep":
            new_links_df.loc[not_upgraded, "p_nom"] = new_links_df["p_nom"].fillna(0)
        elif transmission_projects["new_link_capacity"] == "zero":
            new_links_df.loc[not_upgraded, "p_nom"] = 0
    # export csv files for new buses, lines, links and adjusted lines and links
    new_lines_df.to_csv(snakemake.output.new_lines)
    new_links_df.to_csv(snakemake.output.new_links)
    adjust_lines_df.to_csv(snakemake.output.adjust_lines)
    adjust_links_df.to_csv(snakemake.output.adjust_links)
    new_buses_df.to_csv(snakemake.output.new_buses)
