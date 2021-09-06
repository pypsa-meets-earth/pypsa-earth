""" OSM build network."""

import os
import sys
import logging

import geopandas as gpd
import numpy as np

from _helpers import configure_logging, _sets_path_to_root
logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


def line_endings_to_bus_conversion(lines):
    # Assign to every line a start and end point

    lines["bounds"] = lines["geometry"].boundary  # create start and end point
    # splits into coordinates
    lines["bus0_lon"] = lines["bounds"].bounds.iloc[:, 0]
    lines["bus0_lat"] = lines["bounds"].bounds.iloc[:, 1]
    lines["bus1_lon"] = lines["bounds"].bounds.iloc[:, 2]
    lines["bus1_lat"] = lines["bounds"].bounds.iloc[:, 3]

    # TODO: FIX bug. Extract the coordinates is working for Nigeria bt not for
    # the whole African continent
    # lines_ng["bus_0_coors"]=lines_ng["bounds"].apply(lambda mp: mp[0])
    # lines_ng["bus_1_coors"]=lines_ng["bounds"].apply(lambda mp: mp[1])
    lines["bus_0_coors"] = np.nan
    lines["bus_1_coors"] = np.nan

    return lines


def create_bus_df_from_lines(substations, lines):
    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns)
    bus_e = gpd.GeoDataFrame(columns=substations.columns)

    # Read information from line.csv
    bus_s[["voltage", "lon", "lat", "geometry", "country"]] = lines[
        ["voltage", "bus0_lon", "bus0_lat", "bus_0_coors", "country"]
    ]  # line start points
    bus_e[["voltage", "lon", "lat", "geometry", "country"]] = lines[
        ["voltage", "bus1_lon", "bus1_lat", "bus_1_coors", "country"]
    ]  # line end points
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
    bus_all["station_id"] = bus_all.groupby(["lon", "lat"]).ngroup()

    # For each station number with multiple buses make lowest voltage `substation_lv = TRUE`
    bus_with_stations_duplicates = bus_all[bus_all.station_id.duplicated(keep=False)].sort_values(by=["station_id","voltage"])
    lv_bus_at_station_duplicates = bus_all[bus_all.station_id.duplicated(keep=False)].sort_values(by=["station_id","voltage"]).drop_duplicates(subset=["station_id"])
    # Set all buses with station duplicates "False"
    bus_all.loc[bus_with_stations_duplicates.index, "substation_lv"] = False
    # Set lv_buses with station duplicates "True"
    bus_all.loc[lv_bus_at_station_duplicates.index, "substation_lv"] = True

    # Add station_id to line dataframe
    n_row = int(bus_all.shape[0] / 2)  # row length
    lines = lines.reset_index(drop=True)
    lines["bus0"] = bus_all.loc[: (n_row - 1), ["bus_id"]]
    lines["bus1"] = bus_all.loc[n_row:, ["bus_id"]].reset_index(drop=True)

    # TRY: Keep only buses that are not duplicated & lv_substation = True
    # TODO: Check if this is necessary. What effetc do duplicates have?
    bus_all = bus_all[bus_all["substation_lv"] == True]

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
        os.getcwd(), "data", "base_network", "africa_all" + "_lines" + "_build_network"
    )

    # create clean directory if not already exist
    if not os.path.exists(outputfile_partial):
        os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True)

    lines.to_csv(outputfile_partial + ".csv")  # Generate CSV

    # Buses
    # Output file directory
    outputfile_partial = os.path.join(
        os.getcwd(), "data", "base_network", "africa_all" + "_buses" + "_build_network"
    )
    # create clean directory if not already exist
    if not os.path.exists(outputfile_partial):
        os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True)
    # Generate CSV
    buses.to_csv(outputfile_partial + ".csv")

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():    
        from _helpers import mock_snakemake
        snakemake = mock_snakemake("build_osm_network")
    configure_logging(snakemake)

    _sets_path_to_root("pypsa-africa")
    #sys.path.append("./../../scripts")  ## alternative

    built_network()
