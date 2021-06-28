# This script does the following
# 1. Downloads OSM files for specified countries from Geofabrik
# 2. Filters files for substations, lines and generators
# 3. Process and clean data
# 4. Exports to CSV
# 5. Exports to GeoJson

"""
OSM extraction scrpt
"""

import os
import sys

# IMPORTANT: RUN SCRIPT FROM THIS SCRIPTS DIRECTORY i.e data_exploration/ TODO: make more robust
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("../../scripts")

import logging
import shutil

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from esy.osmfilter import run_filter
from esy.osmfilter import Node, Relation, Way
from esy.osmfilter import osm_info as osm_info
from esy.osmfilter import osm_pickle as osm_pickle
from iso_country_codes import AFRICA_CC
from shapely.geometry import LineString, Point, Polygon

logger = logging.getLogger(__name__)

# https://gitlab.com/dlr-ve-esy/esy-osmfilter/-/tree/master/


# import logging
# logging.basicConfig()
# logger=logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
# logger.setLevel(logging.WARNING)

# Downloads PBF File for given Country Code


def download_pbf(country_code, update):
    """
    Downloads the pbf file from geofabrik for a given country code (see scripts/iso_country_codes.py).

    Parameters
    ----------
    country_code : str
    update : bool
        name of the network component
        update = true forces re-download of files
    """
    country_name = AFRICA_CC[country_code]
    # Filename for geofabrik
    geofabrik_filename = f"{country_name}-latest.osm.pbf"
    # https://download.geofabrik.de/africa/nigeria-latest.osm.pbf
    geofabrik_url = f"https://download.geofabrik.de/africa/{geofabrik_filename}"
    PBF_inputfile = os.path.join(
        os.getcwd(), "data", "osm", "pbf", geofabrik_filename
    )  # Input filepath

    if not os.path.exists(PBF_inputfile) or update is True:
        print(f"{geofabrik_filename} does not exist, downloading to {PBF_inputfile}")
        #  create data/osm directory
        os.makedirs(os.path.dirname(PBF_inputfile), exist_ok=True)
        with requests.get(geofabrik_url, stream=True) as r:
            with open(PBF_inputfile, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return PBF_inputfile


def download_and_filter(country_code, update=False):
    PBF_inputfile = download_pbf(country_code, update)

    filter_file_exists = False
    # json file for the Data dictionary
    JSON_outputfile = os.path.join(
        os.getcwd(), "data", "osm", country_code + "_power.json"
    )  # json file for the Elements dictionary is automatically written to "data/osm/Elements"+filename)

    if os.path.exists(JSON_outputfile):
        filter_file_exists = True

    # Load Previously Pre-Filtered Files
    if update is False and filter_file_exists is True:
        create_elements = False  # Do not create elements again
        new_prefilter_data = False  # Do not pre-filter data again
        # HACKY: esy.osmfilter code to re-create Data.pickle
        Data = osm_info.ReadJason(JSON_outputfile, verbose="no")
        DataDict = {"Data": Data}
        osm_pickle.picklesave(
            DataDict,
            os.path.realpath(
                os.path.join(os.getcwd(), os.path.dirname(JSON_outputfile))
            ),
        )
        print(f"Loading Pickle for {AFRICA_CC[country_code]}")  # TODO: Change to Logger
    else:
        create_elements = True
        new_prefilter_data = True
        print(
            f"Creating  New Elements for {AFRICA_CC[country_code]}"
        )  # TODO: Change to Logger

    prefilter = {
        Node: {"power": ["substation", "line", "generator"]},
        Way: {"power": ["substation", "line", "generator"]},
        Relation: {"power": ["substation", "line", "generator"]},
    }  # see https://dlr-ve-esy.gitlab.io/esy-osmfilter/filter.html for filter structures
    # HACKY: due to esy.osmfilter validation

    blackfilter = [
        ("", ""),
    ]

    for feature in ["substation", "line", "generator"]:
        whitefilter = [
            [
                ("power", feature),
            ],
        ]
        elementname = f"{country_code}_{feature}s"

        feature_data = run_filter(
            elementname,
            PBF_inputfile,
            JSON_outputfile,
            prefilter,
            whitefilter,
            blackfilter,
            NewPreFilterData=new_prefilter_data,
            CreateElements=create_elements,
            LoadElements=True,
            verbose=False,
            multiprocess=True,
        )
        # For better performance yield feature_data here
        if feature == "substation":
            substation_data = feature_data
        if feature == "line":
            line_data = feature_data
        if feature == "generator":
            generator_data = feature_data

    return (substation_data, line_data, generator_data)

# Convert Filtered Data, Elements to Pandas Dataframes


def convert_filtered_data_to_dfs(country_code, feature_data, feature):
    [Data, Elements] = feature_data
    elementname = f"{country_code}_{feature}s"
    df_way = pd.json_normalize(Elements[elementname]["Way"].values())
    df_node = pd.json_normalize(Elements[elementname]["Node"].values())
    return (df_node, df_way, Data)


# Lookup refs and convert to list of longlats


def lonlat_lookup(df_way, Data):
    lonlat_list = []

    col = "refs"
    if col not in df_way.columns:
        print ("refs column not found")
        df_way[col] = pd.Series().astype(float) # create empty "refs" if not in dataframe
      
    for ref in df_way["refs"]:
        lonlat_row = []
        for r in ref:
            lonlat = tuple(Data["Node"][str(r)]["lonlat"])
            lonlat_row.append(lonlat)
        lonlat_list.append(lonlat_row)
    return lonlat_list


# Convert Ways to Point Coordinates


def convert_ways_points(df_way, Data):
    lonlat_list = lonlat_lookup(df_way, Data)
    lonlat_column = []
    area_column = []
    for lonlat in lonlat_list:
        way_polygon = Polygon(lonlat)
        polygon_area = int(round(gpd.GeoSeries(way_polygon).set_crs("EPSG:4326").to_crs("EPSG:3857").area, -1)) # nearest tens m2
        # print('{:g}'.format(float('{:.3g}'.format(float(polygon_area))))) # For significant numbers
        area_column.append(polygon_area)
        center_point = way_polygon.centroid
        lonlat_column.append(list((center_point.x, center_point.y)))

    # df_way.drop("refs", axis=1, inplace=True, errors="ignore")
    df_way.insert(2, "Area", area_column)
    df_way.insert(0, "lonlat", lonlat_column)


# Convert Ways to Line Coordinates


def convert_ways_lines(df_way, Data):
    lonlat_column = lonlat_lookup(df_way, Data)
    df_way.insert(0, "lonlat", lonlat_column)


# Convert Points Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf(df_way):
    gdf = gpd.GeoDataFrame(df_way, geometry=[Point(x, y) for x, y in df_way.lonlat], crs="EPSG:4326")
    gdf.drop(columns=["lonlat"], inplace=True)
    return gdf


# Convert Lines Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf_lines(df_way, simplified=False):
    df_way["geometry"] = df_way["lonlat"].apply(lambda x: LineString(x))
    if simplified is True:
        df_way["geometry"] = df_way["geometry"].apply(
            lambda x: x.simplify(0.005, preserve_topology=False)
        )
    gdf = gpd.GeoDataFrame(df_way, geometry="geometry", crs="EPSG:4326")
    gdf.drop(columns=["lonlat"], inplace=True)

    return gdf


def process_bus_data(country_code, feature_data, feature):

    df_node, df_way, Data = convert_filtered_data_to_dfs(
        country_code, feature_data, feature
    )

    convert_ways_points(df_way, Data)
    # Add Type Column
    df_node["Type"] = "Node"
    df_way["Type"] = "Way"

    df_combined = pd.concat([df_node, df_way], axis=0)
    # Add Country Column
    df_combined["Country"] = AFRICA_CC[country_code]

    return df_combined


def process_line_data(country_code, line_data):
    df_node, df_way, Data = convert_filtered_data_to_dfs(
        country_code, line_data, "line"
    )
    convert_ways_lines(df_way, Data)
    # Add Type Column
    df_way["Type"] = "Way"

    # Add Country Column
    df_way["Country"] = AFRICA_CC[country_code]
    return df_way



def process_data():
    df_all_substations = pd.DataFrame()
    df_all_lines = pd.DataFrame()
    df_all_generators = pd.DataFrame()
    test_CC = {"NG": "nigeria"}
    for country_code in test_CC.keys():
        substation_data, line_data, generator_data = download_and_filter(country_code)
        for feature in ["substation", "line", "generator"]:
            if feature == "substation":
                df_substation = process_bus_data(country_code, substation_data, feature)
                df_all_substations = pd.concat([df_all_substations, df_substation])
            if feature == "line":
                df_line = process_line_data(country_code, line_data)
                df_all_lines = pd.concat([df_all_lines, df_line])
            if feature == "generator":
                df_generator = process_bus_data(country_code, generator_data, feature)
                df_all_generators = pd.concat([df_all_generators, df_generator])

    # ----------- SUBSTATIONS -----------

    # Columns of interest
    df_all_substations = df_all_substations[
        df_all_substations.columns &
        [
            "id",
            "lonlat",
            "Area",
            "tags.power",
            "tags.substation",
            "tags.voltage",
            "tags.frequency",
            "Type",
            # "refs",
            "Country",
        ]
    ]

    # Clean Data
    df_all_substations.drop(df_all_substations.loc[df_all_substations['tags.substation']=='industrial'].index, inplace=True) # Drop industrial substations
    df_all_substations.drop(df_all_substations.loc[df_all_substations['tags.substation']=='distribution'].index, inplace=True) # Drop distribution substations
    # df_all_substations.dropna(subset=['tags.voltage'], inplace = True) # Drop any substations with Voltage = N/A
    # df_all_substations.dropna(thresh=len(df_all_substations)*0.25, axis=1, how='all', inplace = True) #Drop Columns with 75% values as N/A
    # df_all_substations.reset_index(drop=True, inplace=True)

    # Generate Files
    outputfile_partial = os.path.join(
        os.getcwd(), "data", "africa_all" + "_substations."
    )
    df_all_substations.to_csv(outputfile_partial + "csv")  # Generate CSV
    gdf_substations = convert_pd_to_gdf(df_all_substations)
    gdf_substations.to_file(
        outputfile_partial + "geojson", driver="GeoJSON"
    )  # Generate GeoJson

    # ----------- LINES -----------

    # Columns of interest
    df_all_lines = df_all_lines[
        df_all_lines.columns &
        [
            "id",
            "lonlat",
            "tags.power",
            "tags.cables",
            "tags.voltage",
            "tags.circuits",
            "tags.frequency",
            "Type",
            # "refs",
            "Country"
        ]
    ]

    # Clean Data
    # TODO: Fix Cleaning by uncommenting and verifying results
    # df_all_lines = df_all_lines.dropna(subset=['tags.voltage']) # Drop any lines with Voltage = N/A
    # df_all_lines = df_all_lines.rename(columns = {'tags.voltage':"voltage_V"}) 
    # df_all_lines['voltage_V'] = df_all_lines['voltage_V'].str.split('*').str[0] #just keeps the 
    # df_all_lines['voltage_V'] = df_all_lines['voltage_V'].str.split(';').str[0]
    # df_all_lines['voltage_V'] = df_all_lines['voltage_V'].apply(lambda x: pd.to_numeric(x, errors='coerce')).dropna() ## if cell can't converted to float -> drop
    # df_all_lines = df_all_lines[df_all_lines.voltage_V > 10000]
    df_all_lines = df_all_lines.reset_index(drop=True)


    # Generate Files
    outputfile_partial = os.path.join(os.getcwd(), "data", "africa_all" + "_lines.")
    df_all_lines.to_csv(outputfile_partial + "csv")  # Generate CSV
    gdf_lines = convert_pd_to_gdf_lines(df_all_lines, simplified=False) # Set simplified = True to simplify lines
    gdf_lines.to_file(
        outputfile_partial + "geojson", driver="GeoJSON"
    )  # Generate GeoJson

    # ----------- Generator -----------

    # Columns of interest
    df_all_generators = df_all_generators[
        df_all_generators.columns &
        [
            "id",
            "lonlat",
            # "Area",
            "tags.power",
            "tags.generator:type",
            "tags.generator:method",
            "tags.generator:source",
            "tags.generator:output:electricity",
            "Type",
            # "refs",
            "Country",
        ]
    ]

    # Clean Data
    # TODO: Fix Cleaning by uncommenting and verifying results
    # TODO: Should Probably drop columns (as in the line below) or find a better solution for the above fix for selecting columns
    # df_all_generators.drop(columns = ["tags.fixme","tags.name:ar","tags.building","tags.barrier"], inplace = True, errors='ignore') #
    # df_all_generators = df_all_generators[df_all_generators['tags.generator:output:electricity'].astype(str).str.contains('MW')] #removes boolean 
    # df_all_generators['tags.generator:output:electricity'] = df_all_generators['tags.generator:output:electricity'].str.extract('(\d+)').astype(float)
    # df_all_generators.rename(columns = {'tags.generator:output:electricity':"power_output_MW"}, inplace = True)
    # df_all_generators.reset_index(drop=True, inplace=True)

    # Generate Files
    outputfile_partial = os.path.join(
        os.getcwd(), "data", "africa_all" + "_generators."
    )
    df_all_generators.to_csv(outputfile_partial + "csv")  # Generate CSV
    gdf_generators = convert_pd_to_gdf(df_all_generators)
    gdf_generators.to_file(
        outputfile_partial + "geojson", driver="GeoJSON"
    )  # Generate GeoJson


if __name__ == "__main__":
    process_data()
