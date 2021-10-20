import logging
import os
import sys
import math

import geopandas as gpd
import numpy as np
from _helpers import _sets_path_to_root
from _helpers import _to_csv_nafix
from _helpers import configure_logging

logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


def line_endings_to_bus_conversion(lines):
    # Assign to every line a start and end point

    lines["bounds"] = lines["geometry"].boundary  # create start and end point
    lines["bus_0_coors"] = lines["bounds"].map(lambda p: p.geoms[0]
                                               if len(p.geoms) >= 2 else None)
    lines["bus_1_coors"] = lines["bounds"].map(lambda p: p.geoms[1]
                                               if len(p.geoms) >= 2 else None)

    # splits into coordinates
    lines["bus0_lon"] = lines["bus_0_coors"].map(lambda p: p.x
                                                 if p != None else None)
    lines["bus0_lat"] = lines["bus_0_coors"].map(lambda p: p.y
                                                 if p != None else None)
    lines["bus1_lon"] = lines["bus_1_coors"].map(lambda p: p.x
                                                 if p != None else None)
    lines["bus1_lat"] = lines["bus_1_coors"].map(lambda p: p.y
                                                 if p != None else None)

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


def set_substations_ids(buses, tol=0.01):  # tol=0.01, around 700m at latitude 44.
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

    station_id = 0
    for i, row in buses.iterrows():
        if buses.loc[i, "station_id"] >= 0:
            continue

        # get substations within tolerance
        close_nodes = np.where(buses.apply(
            lambda x: math.dist([row["lat"], row["lon"]], [x["lat"], x["lon"]]) <= tol,
            axis=1
        ))[0]

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
            subset_substation_ids = buses.loc[buses.index[close_nodes],"station_id"]
            all_neg = subset_substation_ids.max() < 0  # check if all substation_ids are negative (<0)
            some_neg = subset_substation_ids.min() < 0  # check if at least a substation_id is negative (<0)

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
                    f"link{buses_station_id}_{b_it}", buses_station_id.index[0], buses_station_id.index[b_it], 400000, 1, 0.0,
                    False, False, "transmission", 50, buses_station_id.country.iloc[0],
                    LineString([buses_station_id.geometry.iloc[0], buses_station_id.geometry.iloc[b_it]]), LineString([buses_station_id.geometry.iloc[0], buses_station_id.geometry.iloc[b_it]]).bounds,
                    buses_station_id.geometry.iloc[0], buses_station_id.geometry.iloc[b_it], buses_station_id.lon.iloc[0], buses_station_id.lat.iloc[0], buses_station_id.lon.iloc[b_it], buses_station_id.lat.iloc[b_it]
                ])
    return lines.append(gpd.GeoDataFrame(add_lines, columns=lines.keys()), ignore_index=True)
    


def create_station_at_equal_bus_locations(lines, buses):
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

    set_substations_ids(bus_all)
    # bus_all["station_id"] = bus_all.groupby(["lon", "lat"]).ngroup()

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    bus_with_stations_duplicates = bus_all[bus_all.station_id.duplicated(
        keep=False)].sort_values(by=["station_id", "voltage"])
    lv_bus_at_station_duplicates = (bus_all[bus_all.station_id.duplicated(
        keep=False)].sort_values(by=["station_id", "voltage"]).drop_duplicates(
            subset=["station_id"]))
    # Set all buses with station duplicates "False"
    bus_all.loc[bus_with_stations_duplicates.index, "substation_lv"] = False
    # Set lv_buses with station duplicates "True"
    bus_all.loc[lv_bus_at_station_duplicates.index, "substation_lv"] = True

    # Add station_id to line dataframe.
    # Note: by construction, the first half of bus_all is "bus0" and the rest is "bus1"
    n_row = int(bus_all.shape[0] / 2)  # row length
    lines = lines.reset_index(drop=True)
    lines["bus0"] = bus_all.loc[:(n_row - 1), ["bus_id"]]
    lines["bus1"] = bus_all.loc[n_row:, ["bus_id"]].reset_index(drop=True)

    # TRY: Keep only buses that are not duplicated & lv_substation = True
    # TODO: Check if this is necessary. What effetc do duplicates have?
    bus_all = bus_all[bus_all["substation_lv"] == True]

    lines = connect_stations_same_station_id(lines, buses)

    return lines, buses


def built_network():
    # ----------- LOAD DATA -----------
    paths = os.path.realpath("data/clean") + "/africa_all_substations.geojson"
    pathl = os.path.realpath("data/clean") + "/africa_all_lines.geojson"
    substations = gpd.read_file(paths).set_crs(epsg=4326, inplace=True)
    lines = gpd.read_file(pathl).set_crs(epsg=4326, inplace=True)

    # Filter only Nigeria
    # lines = lines[lines.loc[:,"country"] == "nigeria"].copy()
    # substations = substations[substations.loc[:,"country"] == "nigeria"].copy()

    # Use lines and create bus/line df
    lines = line_endings_to_bus_conversion(lines)
    buses = create_bus_df_from_lines(substations, lines)

    # Methods to build network
    # METHOD 1 - Form station at exact same location
    # TODO: Select at clean_data(method = ...) the applied method
    lines, buses = create_station_at_equal_bus_locations(lines, buses)
    # METHOD 2 - Use interference method
    # TODO: Add and test other method (ready in jupyter)

    # Export data
    # Lines
    # Output file directory
    outputfile_partial = os.path.join(
        os.getcwd(), "data", "base_network",
        "africa_all" + "_lines" + "_build_network")

    # create clean directory if not already exist
    if not os.path.exists(outputfile_partial):
        os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True)

    _to_csv_nafix(lines, outputfile_partial + ".csv")  # Generate CSV

    # Buses
    # Output file directory
    outputfile_partial = os.path.join(
        os.getcwd(), "data", "base_network",
        "africa_all" + "_buses" + "_build_network")
    # create clean directory if not already exist
    if not os.path.exists(outputfile_partial):
        os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True)
    # Generate CSV
    _to_csv_nafix(buses, outputfile_partial + ".csv")

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake("build_osm_network")
    configure_logging(snakemake)

    _sets_path_to_root("pypsa-africa")
    # sys.path.append("./../../scripts")  ## alternative

    built_network()
