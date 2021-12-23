import logging
import math
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from _helpers import _sets_path_to_root
from _helpers import _to_csv_nafix
from _helpers import _read_geojson
from _helpers import configure_logging
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.ops import linemerge
from shapely.ops import unary_union
from tqdm import tqdm
from download_osm_data import create_country_list

logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


def line_endings_to_bus_conversion(lines):
    # Assign to every line a start and end point

    lines["bounds"] = lines["geometry"].boundary  # create start and end point

    lines["bus_0_coors"] = lines["bounds"].map(lambda p: p.geoms[0])
    lines["bus_1_coors"] = lines["bounds"].map(lambda p: p.geoms[1])

    # splits into coordinates
    lines["bus0_lon"] = lines["bus_0_coors"].x
    lines["bus0_lat"] = lines["bus_0_coors"].y
    lines["bus1_lon"] = lines["bus_1_coors"].x
    lines["bus1_lat"] = lines["bus_1_coors"].y

    return lines


def create_bus_df_from_lines(substations, lines):
    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns)
    bus_e = gpd.GeoDataFrame(columns=substations.columns)

    # Read information from line.csv
    bus_s[["voltage", "lon", "lat", "geometry", "country"]] = lines[[
        "voltage", "bus0_lon", "bus0_lat", "bus_0_coors", "country"
    ]]  # line start points
    bus_e[["voltage", "lon", "lat", "geometry", "country"]] = lines[[
        "voltage", "bus1_lon", "bus1_lat", "bus_1_coors", "country"
    ]]  # line end points
    bus_all = bus_s.append(bus_e).reset_index(drop=True)

    # Assign index to bus_id
    bus_all.loc[:, "bus_id"] = bus_all.index
    buses = bus_all

    # Removing the NaN
    buses["dc"] = "False"
    buses["symbol"] = "False"
    buses["under_construction"] = "False"
    buses["tag_substation"] = "False"
    buses["tag_area"] = "False"
    buses["substation_lv"] = True

    return buses


def add_line_endings_tosubstations(substations, lines):
    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns)
    bus_e = gpd.GeoDataFrame(columns=substations.columns)

    # Read information from line.csv
    bus_s[["voltage", "country"]] = lines[["voltage",
                                           "country"]]  # line start points
    bus_s["geometry"] = lines.geometry.boundary.map(lambda p: p.geoms[0])
    bus_s["lon"] = bus_s["geometry"].x
    bus_s["lat"] = bus_s["geometry"].y
    bus_s["bus_id"] = lines["line_id"].astype(str) + "_s"

    bus_e[["voltage", "country"]] = lines[["voltage",
                                           "country"]]  # line start points
    bus_e["geometry"] = lines.geometry.boundary.map(lambda p: p.geoms[1])
    bus_e["lon"] = bus_s["geometry"].x
    bus_e["lat"] = bus_s["geometry"].y
    bus_e["bus_id"] = lines["line_id"].astype(str) + "_e"

    bus_all = bus_s.append(bus_e).reset_index(drop=True)
    # Assign index to bus_id
    bus_all.loc[:, "bus_id"] = bus_all.index
    buses = bus_all

    # Add NaN as default
    bus_all["station_id"] = np.nan
    bus_all["dc"] = False  # np.nan
    # Assuming substations completed for installed lines
    bus_all["under_construction"] = False
    bus_all["tag_area"] = 0.0  # np.nan
    bus_all["symbol"] = "substation"
    # TODO: this tag may be improved, maybe depending on voltage levels
    bus_all["tag_substation"] = "transmission"

    buses = substations.append(bus_all).reset_index(drop=True)

    return buses


# tol in m
def set_substations_ids(buses, tol=2000):
    """
    Function to set substations ids to buses, accounting for location tolerance

    The algorithm is as follows:

    1. initialize all substation ids to -1
    2. if the current substation has been already visited [substation_id < 0], then skip the calculation
    3. otherwise:
        1. identify the substations within the specified tolerance (tol)
        2. when all the substations in tolerance have substation_id < 0, then specify a new substation_id
        3. otherwise, if one of the substation in tolerance has a substation_id >= 0, then set that substation_id to all the others;
           in case of multiple substations with substation_ids >= 0, the first value is picked for all

    """

    buses["station_id"] = -1

    # create temporary series to execute distance calculations using m as reference distances
    temp_bus_geom = buses.geometry.to_crs(3736)
    
    # set tqdm options for substation ids
    tqdm_kwargs_substation_ids = dict(
        ascii=False,
        unit=" buses",
        total=buses.shape[0],
        desc="Set substation ids ",
    )

    station_id = 0
    for i, row in tqdm(buses.iterrows(), **tqdm_kwargs_substation_ids):
        if buses.loc[i, "station_id"] >= 0:
            continue

        # get substations within tolerance
        close_nodes = np.flatnonzero(
            temp_bus_geom.distance(temp_bus_geom.loc[i]) <= tol)

        # print("set_substations_ids: ", i, "/", buses.shape[0])

        if len(close_nodes) == 1:
            # if only one substation is in tolerance, then the substation is the current one iì
            # Note that the node cannot be with substation_id >= 0, given the preliminary check
            # at the beginning of the for loop
            buses.loc[buses.index[i], "station_id"] = station_id
            # update station id
            station_id += 1
        else:
            # several substations in tolerance

            # get their ids
            subset_substation_ids = buses.loc[buses.index[close_nodes],
                                              "station_id"]
            # check if all substation_ids are negative (<0)
            all_neg = subset_substation_ids.max() < 0
            # check if at least a substation_id is negative (<0)
            some_neg = subset_substation_ids.min() < 0

            if all_neg:
                # when all substation_ids are negative, then this is a new substation id
                # set the current station_id and increment the counter
                buses.loc[buses.index[close_nodes], "station_id"] = station_id
                station_id += 1
            elif some_neg:
                # otherwise, when at least a substation_id is non-negative, then pick the first value
                # and set it to all the other substations within tolerance
                sub_id = -1
                for substation_id in subset_substation_ids:
                    if substation_id >= 0:
                        sub_id = substation_id
                        break
                buses.loc[buses.index[close_nodes], "station_id"] = sub_id


def set_lines_ids(lines, buses):
    """
    Function to set line buses ids to the closest bus in the list

    """
    
    # set tqdm options for set lines ids
    tqdm_kwargs_line_ids = dict(
        ascii=False,
        unit=" lines",
        total=lines.shape[0],
        desc="Set line bus ids ",
    )

    # initialization
    lines["bus0"] = -1
    lines["bus1"] = -1

    for i, row in tqdm(lines.iterrows(), **tqdm_kwargs_line_ids):

        # print("set_lines_ids: ", i, "/", lines.shape[0])

        # select buses having the voltage level of the current line
        buses_sel = buses[buses["voltage"] == row["voltage"]]

        # find the closest node of the bus0 of the line
        bus0_id = buses_sel.geometry.distance(row["bus_0_coors"]).idxmin()
        lines.loc[i, "bus0"] = buses.loc[bus0_id, "bus_id"]

        # check if the line starts exactly in the node, otherwise modify the linestring
        distance_bus0 = buses.loc[bus0_id, "geometry"].distance(
            row["bus_0_coors"])
        if distance_bus0 > 0.0:
            # the line does not start in the node, thus modify the linestring
            lines.loc[i, "geometry"] = linemerge([
                LineString([
                    buses.loc[bus0_id, "geometry"],
                    row["bus_0_coors"],
                ]),
                lines.loc[i, "geometry"],
            ])

        # find the closest node of the bus1 of the line
        bus1_id = buses_sel.geometry.distance(row["bus_1_coors"]).idxmin()
        lines.loc[i, "bus1"] = buses.loc[bus1_id, "bus_id"]

        # check if the line ends exactly in the node, otherwise modify the linestring
        distance_bus1 = buses.loc[bus1_id, "geometry"].distance(
            row["bus_1_coors"])
        if distance_bus1 > 0.0:
            # the line does not end in the node, thus modify the linestring
            lines.loc[i, "geometry"] = linemerge([
                lines.loc[i, "geometry"],
                LineString([row["bus_1_coors"], buses.loc[bus1_id, "geometry"]]),
            ])

    return lines, buses


def merge_stations_same_station_id(buses, delta_lon=0.001, delta_lat=0.001):
    """
    Function to merge buses with same voltage and station_id
    This function iterattes over all substation ids and creates a bus_id for every substation and voltage level.
    Therefore, a substation with multiple voltage levels is represented with different buses, one per voltage level
    """
    # initialize empty dataset of cleaned buses
    buses_clean = gpd.GeoDataFrame(columns=buses.columns)

    # initalize the number of buses
    n_buses = 0

    for g_name, g_value in buses.groupby(by=["station_id"]):

        # average location of the buses having the same station_id
        station_loc = g_value[["lon", "lat"]].mean()

        # loop for every voltage level in the bus
        # The location of the buses is averaged; in the case of multiple voltage levels for the same station_id,
        # each bus corresponding to a voltage level is located at a distanceregulated by delta_lon/delta_lat
        v_it = 0
        for v_name, bus_row in g_value.groupby(by=["voltage"]):

            # add the bus
            buses_clean.loc[n_buses] = {
                "bus_id":
                n_buses,
                "station_id":
                g_name,
                "voltage":
                v_name,
                "dc":
                bus_row["dc"].all(),
                "symbol":
                "|".join(bus_row["symbol"].unique()),
                "under_construction":
                bus_row["under_construction"].any(),
                "tag_substation":
                "|".join(bus_row["tag_substation"].unique()),
                "tag_area":
                bus_row["tag_area"].sum(),
                "lon":
                station_loc.lon + v_it * delta_lon,
                "lat":
                station_loc.lat + v_it * delta_lat,
                "country":
                bus_row["country"].iloc[0],
                "geometry":
                Point(
                    station_loc.lon + v_it * delta_lon,
                    station_loc.lat + v_it * delta_lat,
                ),
            }

            # increase counters
            v_it += 1
            n_buses += 1

    return buses_clean


def get_transformers(buses, lines):
    """
    Function to create fake transformer lines that connect buses of the same station_id at different voltage
    """

    df_transformers = gpd.GeoDataFrame(columns=lines.columns)

    for g_name, g_value in buses.sort_values(
            "voltage", ascending=True).groupby(by=["station_id"]):
        
        # print("get_transformers: ", g_name)

        # note: by construction there cannot be more that two nodes with the same station_id and same voltage
        n_voltages = len(g_value)

        if n_voltages > 1:

            for id in range(0, n_voltages - 1):
                # when g_value has more than one node, it means that there are multiple voltages for the same bus
                geom_trans = LineString(
                    [g_value.geometry.iloc[id], g_value.geometry.iloc[id + 1]])
                transf_id = lines.shape[0] + df_transformers.shape[0] + 1

                df_transformers.loc[transf_id] = {
                    "line_id": f"transf_{g_name}_{id}",
                    "bus0": g_value["bus_id"].iloc[id],
                    "bus1": g_value["bus_id"].iloc[id + 1],
                    "voltage": g_value.voltage.iloc[[id, id + 1]].max(),
                    "circuits": 1,
                    "length": 0.0,
                    "underground": False,
                    "under_construction": False,
                    "tag_type": "transmission",
                    "tag_frequency": 50,
                    "country": g_value.country.iloc[id],
                    "geometry": geom_trans,
                    "bounds": geom_trans.bounds,
                    "bus_0_coors": g_value.geometry.iloc[id],
                    "bus_1_coors": g_value.geometry.iloc[id + 1],
                    "bus0_lon": g_value.geometry.iloc[id].x,
                    "bus0_lat": g_value.geometry.iloc[id].y,
                    "bus1_lon": g_value.geometry.iloc[id + 1].x,
                    "bus1_lat": g_value.geometry.iloc[id + 1].y,
                }

    # update line endings
    df_transformers = line_endings_to_bus_conversion(df_transformers)

    return df_transformers


def connect_stations_same_station_id(lines, buses):
    """
    Function to create fake links between substations with the same substation_id
    """

    station_id_list = buses.station_id.unique()

    add_lines = []
    from shapely.geometry import LineString

    for s_id in station_id_list:
        buses_station_id = buses[buses.station_id == s_id]

        if len(buses_station_id) > 1:
            for b_it in range(1, len(buses_station_id)):
                add_lines.append([
                    f"link{buses_station_id}_{b_it}",
                    buses_station_id.index[0],
                    buses_station_id.index[b_it],
                    400000,
                    1,
                    0.0,
                    False,
                    False,
                    "transmission",
                    50,
                    buses_station_id.country.iloc[0],
                    LineString([
                        buses_station_id.geometry.iloc[0],
                        buses_station_id.geometry.iloc[b_it],
                    ]),
                    LineString([
                        buses_station_id.geometry.iloc[0],
                        buses_station_id.geometry.iloc[b_it],
                    ]).bounds,
                    buses_station_id.geometry.iloc[0],
                    buses_station_id.geometry.iloc[b_it],
                    buses_station_id.lon.iloc[0],
                    buses_station_id.lat.iloc[0],
                    buses_station_id.lon.iloc[b_it],
                    buses_station_id.lat.iloc[b_it],
                ])
    return lines.append(gpd.GeoDataFrame(add_lines, columns=lines.keys()),
                        ignore_index=True)


def set_lv_substations(buses):
    """
    Function to set what nodes are lv, thereby setting substation_lv
    The current methodology is to set lv nodes to buses where multiple voltage level are found,
    hence when the station_id is duplicated
    """

    # initialize column substation_lv to true
    buses["substation_lv"] = True

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    bus_with_stations_duplicates = buses[buses.station_id.duplicated(
        keep=False)].sort_values(by=["station_id", "voltage"])
    lv_bus_at_station_duplicates = (buses[buses.station_id.duplicated(
        keep=False)].sort_values(by=["station_id", "voltage"]).drop_duplicates(
            subset=["station_id"]))
    # Set all buses with station duplicates "False"
    buses.loc[bus_with_stations_duplicates.index, "substation_lv"] = False
    # Set lv_buses with station duplicates "True"
    buses.loc[lv_bus_at_station_duplicates.index, "substation_lv"] = True

    return buses


# Note tolerance = 0.01 means around 700m
# TODO: the current tolerance is high to avoid an issue in the Nigeria case where line 565939360-1
#       seems to be interconnected to both ends, but at the eastern one, the node is actually not connected
#       another line seems to be exactly touching the node, but from the data point of view it only fly over it.
#       There may be the need to split a line in several segments in the case the line is within tolerance with
#       respect to a node


def merge_stations_lines_by_station_id_and_voltage(lines, buses, tol=2000):
    """
    Function to merge close stations and adapt the line datasets to adhere to the merged dataset

    """
    
    logger.info("Stage 3a/4: Set substation ids with tolerance of %.2f km" % (tol/1000))

    # set substation ids
    set_substations_ids(buses, tol=tol)
    
    logger.info("Stage 3b/4: Merge substations with the same id")

    # merge buses with same station id and voltage
    buses = merge_stations_same_station_id(buses)
    
    logger.info("Stage 3c/4: Specify the bus ids of the line endings")

    # set the bus ids to the line dataset
    lines, buses = set_lines_ids(lines, buses)

    # drop lines starting and ending in the same node
    lines = lines[lines["bus0"] != lines["bus1"]]

    # update line endings
    lines = line_endings_to_bus_conversion(lines)

    # set substation_lv
    set_lv_substations(buses)
    
    logger.info("Stage 3d/4: Add transformers")

    # get transformers: modelled as lines connecting buses with different voltage
    transformers = get_transformers(buses, lines)

    # append transformer lines
    lines = lines.append(transformers)

    # reset index
    lines.reset_index(drop=True, inplace=True)

    return lines, buses


def create_station_at_equal_bus_locations(lines, buses, tol=2000):
    # V1. Create station_id at same bus location
    # - We saw that buses are not connected exactly at one point, they are
    #   usually connected to a substation "area" (analysed on maps)
    # - Create station_id at exactly the same location might therefore be not
    #   always correct
    # - Though as you can see below, it might be still sometime the case.
    #   Examples are **station 4** (2 lines with the same voltage connect at the
    #   same point) and **station 23** (4 lines with two different voltages connect
    #   at the same point)
    # TODO: Filter out the generator lines - defined as going from generator to
    #       the next station which is connected to a load. Excluding generator
    #       lines make proably sense because they are not transmission expansion
    #       relevant. For now we simplify and include generator lines.

    # If same location/geometry make station
    bus_all = buses

    # set substation ids
    set_substations_ids(buses, tol=tol)

    # set the bus ids to the line dataset
    lines, buses = set_lines_ids(lines, buses)

    # update line endings
    lines = line_endings_to_bus_conversion(lines)

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    set_lv_substations(bus_all)

    # # Add station_id to line dataframe.
    # # Note: by construction, the first half of bus_all is "bus0" and the rest is "bus1"
    # n_row = int(bus_all.shape[0] / 2)  # row length
    # lines = lines.reset_index(drop=True)
    # lines["bus0"] = bus_all.loc[:(n_row - 1), ["bus_id"]]
    # lines["bus1"] = bus_all.loc[n_row:, ["bus_id"]].reset_index(drop=True)

    # TRY: Keep only buses that are not duplicated & lv_substation = True
    # TODO: Check if this is necessary. What effetc do duplicates have?
    bus_all = bus_all[bus_all["substation_lv"] == True]

    lines = connect_stations_same_station_id(lines, buses)

    return lines, buses


def built_network(inputs, outputs):
    # ----------- LOAD DATA -----------
    logger.info("Stage 1/4: Read input data")

    substations = gpd.read_file(inputs["substations"]).set_crs(epsg=4326,
                                                               inplace=True)
    lines = gpd.read_file(inputs["lines"]).set_crs(epsg=4326, inplace=True)
    generators = _read_geojson(inputs["generators"]).set_crs(epsg=4326,
                                                             inplace=True)

    
    logger.info("Stage 2/4: Add line endings to the substation datasets")

    # Use lines and create bus/line df
    lines = line_endings_to_bus_conversion(lines)
    buses = add_line_endings_tosubstations(substations, lines)
    # buses = create_bus_df_from_lines(substations, lines)

    # Methods to build network
    # METHOD 1 - Form station at exact same location
    # TODO: Select at clean_data(method = ...) the applied method
    # lines, buses = create_station_at_equal_bus_locations(lines, buses)

    # METHOD 2 - Use interference method
    # TODO: Add and test other method (ready in jupyter)
    
    logger.info("Stage 3/4: Aggregate close substations")

    # METHOD to merge buses with same voltage and within tolerance
    if snakemake.config["clean_osm_data_options"]["group_close_buses"]:
        lines, buses = merge_stations_lines_by_station_id_and_voltage(lines,
                                                                  buses,
                                                                  tol=4000)

    logger.info("Stage 4/4: Add augmented substation to country with no data")
    
    country_shapes_fn = snakemake.input.country_shapes
    country_shapes = (
        gpd.read_file(country_shapes_fn)
        .set_index("name")["geometry"]
        .set_crs(4326)
    )
    input = snakemake.config["countries"]
    country_list = create_country_list(input)
    bus_country_list = buses["country"].unique().tolist() 

    if len(bus_country_list) != len(country_list):
        no_data_countries = set(country_list).difference(set(bus_country_list))
        no_data_countries_shape = country_shapes[
            country_shapes.index.isin(no_data_countries)==True
            ].reset_index().set_crs(4326)
        length = len(no_data_countries)
        df = gpd.GeoDataFrame({
            'voltage': [220000]*length,
            'country': no_data_countries_shape["name"],
            'lon': no_data_countries_shape["geometry"].to_crs(epsg=4326).centroid.x,
            'lat': no_data_countries_shape["geometry"].to_crs(epsg=4326).centroid.y,
            'bus_id': np.arange(len(buses)+1, len(buses)+(length+1), 1),
            'station_id': [np.nan]*4,
            'dc': [False]*length,
            'under_construction': [False]*length,
            'tag_area': [0.0]*length,
            'symbol': ["substation"]*length,
            'tag_substation': ["transmission"]*length,
            'geometry': no_data_countries_shape["geometry"].to_crs(epsg=4326).centroid,
            "substation_lv": [True]*length,
        })
        buses = gpd.GeoDataFrame(pd.concat([buses, df], ignore_index=True), crs=buses.geometry.crs)

    logger.info("Save outputs")

    # create clean directory if not already exist
    if not os.path.exists(outputs["lines"]):
        os.makedirs(os.path.dirname(outputs["lines"]), exist_ok=True)

    _to_csv_nafix(lines, outputs["lines"])  # Generate CSV

    # Buses
    # Output file directory
    # outputfile_partial = os.path.join(
    #     os.getcwd(), "data", "base_network",
    #     "africa_all" + "_buses" + "_build_network")

    # create clean directory if not already exist
    if not os.path.exists(outputs["substations"]):
        os.makedirs(os.path.dirname(outputs["substations"]), exist_ok=True)
    # Generate CSV
    _to_csv_nafix(buses, outputs["substations"])

    # # save generators
    # if not os.path.exists(outputs["generators"]):
    #     os.makedirs(os.path.dirname(outputs["generators"]), exist_ok=True)
    # # Generate CSV
    # _to_csv_nafix(generators, outputs["generators"])

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake("build_osm_network")
    configure_logging(snakemake)

    _sets_path_to_root("pypsa-africa")
    # sys.path.append("./../../scripts")  ## alternative

    built_network(snakemake.input, snakemake.output)
