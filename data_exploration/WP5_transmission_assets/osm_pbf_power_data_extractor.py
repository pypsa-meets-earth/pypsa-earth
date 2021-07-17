# This script does the following
# 1. Downloads OSM files for specified countries from Geofabrik
# 2. Filters files for substations, lines and generators
# 3. Process and clean data
# 4. Exports to CSV
# 5. Exports to GeoJson

""" OSM extraction script."""

import os
import sys

# IMPORTANT: RUN SCRIPT FROM THIS SCRIPTS DIRECTORY i.e data_exploration/ TODO: make more robust
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("./../../scripts")

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
#from ..scripts.iso_country_codes import AFRICA_CC

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
    Download pbf file from geofabrik for a given country code

    Parameters
    ----------
    country_code : str 
        Three letter country codes of the downloaded files 
    update : bool 
        Name of the network component 
        Update = true, forces re-download of files

    Returns
    -------
    Pbf file per country

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
    """
    Download OpenStreetMap raw file for selected tag.

    Apply pbf download and filter with esy.osmfilter selected OpenStreetMap
    tags or data. Examples of possible tags are listed at `OpenStreetMap wiki <https://wiki.openstreetmap.org/wiki/Key:power>`_.
    More information on esy.osmfilter `here <https://gitlab.com/dlr-ve-esy/esy-osmfilter>`_.

    Parameters
    ----------
    country_code : str 
        Three letter country codes of the downloaded files 
    update : bool 
        Name of the network component 
        Update = true, forces re-download of files
        Update = false, uses existing or previously downloaded files to safe time

    Returns
    -------
    substation_data : Data, Elements
    line_data : Data, Elements
    generator_data : Data, Elements
        Nested dictionary with all OpenStreetMap keys of specific component.
        Example of lines. See https://wiki.openstreetmap.org/wiki/Tag:power%3Dline
    """
    PBF_inputfile = download_pbf(country_code, update)

    filter_file_exists = False
    # json file for the Data dictionary
    JSON_outputfile = os.path.join(
        os.getcwd(), "data", "osm", country_code + "_power.json"
    ) # json file for the Elements dictionary is automatically written to "data/osm/Elements"+filename)

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
        Node: {"power": ["substation", "line", "generator", "cable"]},
        Way: {"power": ["substation", "line", "generator", "cable"]},
        Relation: {"power": ["substation", "line", "generator", "cable"]},
    }  # see https://dlr-ve-esy.gitlab.io/esy-osmfilter/filter.html for filter structures
    # HACKY: due to esy.osmfilter validation

    blackfilter = [
        ("", ""),
    ]

    for feature in ["substation", "line", "generator", "cable"]:
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
        if feature == "cable":
            cable_data = feature_data

    return (substation_data, line_data, generator_data, cable_data)

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
        print ("refs column not found") # TODO : Change to logger and do not create a hacky empty ref col
        df_way[col] = pd.Series([],dtype=pd.StringDtype()).astype(float) # create empty "refs" if not in dataframe
      
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
        if len(lonlat) >= 3: #Minimum for a triangle
            way_polygon = Polygon(lonlat)
            # TODO : Do set_crs and to_crs for whole lonlat_list and not indivudually, expected big performance boost.
            polygon_area = int(round(gpd.GeoSeries(way_polygon).set_crs("EPSG:4326").to_crs("EPSG:3857").area, -1)) # nearest tens m2
            # print('{:g}'.format(float('{:.3g}'.format(float(polygon_area))))) # For significant numbers
            area_column.append(polygon_area)
            center_point = way_polygon.centroid
            lonlat_column.append(list((center_point.x, center_point.y)))
        else:
            area_column.append(0)
            center_point = lonlat[0]
            lonlat_column.append(list(center_point))

    # df_way.drop("refs", axis=1, inplace=True, errors="ignore")
    df_way.insert(0, "Area", area_column)
    df_way.insert(0, "lonlat", lonlat_column)


# Convert Ways to Line Coordinates


def convert_ways_lines(df_way, Data):
    lonlat_list = lonlat_lookup(df_way, Data)
    lonlat_column = lonlat_list
    length_column = []
    for lonlat in lonlat_list:
        way_linestring = LineString(lonlat)
        line_length = gpd.GeoSeries(way_linestring).set_crs("EPSG:4326").to_crs("EPSG:3857").length
        length_column.append(float(line_length))

    df_way.insert(0, "Length", length_column)
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


def process_node_data(country_code, feature_data, feature):
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


def process_line_data(country_code, feature_data, feature):
    df_node, df_way, Data = convert_filtered_data_to_dfs(
        country_code, feature_data, feature
    )
    convert_ways_lines(df_way, Data)
    
    # Add Type Column
    df_way["Type"] = "Way"

    # Add Country Column
    df_way["Country"] = AFRICA_CC[country_code]
    return df_way



def process_data():
    columns_substation = [ 
            "id",
            "lonlat",
            "Area",
            "tags.power",
            "tags.substation",
            "tags.voltage",
            "Type",
            # "refs",
            "Country",
        ]

    columns_line = [
            "id",
            "lonlat",
            "Length",
            "tags.power",
            "tags.cables",
            "tags.voltage",
            "tags.circuits",
            "tags.frequency",
            "Type",
            # "refs",
            "Country"
        ]
    
    columns_generator= [
            "id",
            "lonlat",
            "Area",
            "tags.power",
            "tags.generator:type",
            "tags.generator:method",
            "tags.generator:source",
            "tags.generator:output:electricity",
            "Type",
            # "refs",
            "Country",
        ]

    columns_cable= [
            "id",
            "lonlat",
            "Length",
            "tags.power",
            "tags.cables",
            "tags.voltage",
            "tags.circuits",
            "tags.frequency",
            "tags.location",
            "Type",
            # "refs",
            "Country"
        ]

    df_all_substations = pd.DataFrame()
    df_all_lines = pd.DataFrame()
    df_all_generators = pd.DataFrame()
    df_all_cables = pd.DataFrame()
    # test_CC = {"NG": "nigeria"}
    # test_CC = {"ZA": "SOUTH AFRICA"} # or any other country
    # Africa_CC = {list of all African countries that are imported from script -> iso_countries_codes}
    for country_code in AFRICA_CC.keys(): # replace Africa_CC by test_CC to only download data for one country
        substation_data, line_data, generator_data, cable_data = download_and_filter(country_code)
        for feature in ["substation", "line", "generator", "cable"]:
            if feature == "substation":
                df_substation = process_node_data(country_code, substation_data, feature)
                df_all_substations = pd.concat([df_all_substations, df_substation])
            if feature == "line":
                df_line = process_line_data(country_code, line_data, feature)
                df_all_lines = pd.concat([df_all_lines, df_line])
            if feature == "generator":
                df_generator = process_node_data(country_code, generator_data, feature)
                df_all_generators = pd.concat([df_all_generators, df_generator])
            if feature == "cable":
                df_cable = process_line_data(country_code, cable_data, feature)
                df_all_cables = pd.concat([df_all_cables, df_cable])


    outputfile_partial = os.path.join(os.getcwd(), "data", "raw", "africa_all" + "_raw") # Output file directory

    if not os.path.exists(outputfile_partial):
        os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True) #  create raw directory


    # ----------- SUBSTATIONS -----------

    df_all_substations = df_all_substations[df_all_substations.columns.intersection(set(columns_substation))]
    df_all_substations.reset_index(drop=True, inplace=True)

    # Generate Files
    df_all_substations.to_csv(outputfile_partial + "_substations" + ".csv")  # Generate CSV
    gdf_substations = convert_pd_to_gdf(df_all_substations)
    gdf_substations.to_file(outputfile_partial + "_substations" + ".geojson", driver="GeoJSON")  # Generate GeoJson

    # ----------- LINES -----------

    df_all_lines = df_all_lines[df_all_lines.columns.intersection(set(columns_line))]
    df_all_lines.reset_index(drop=True, inplace=True)


    # Generate Files
    df_all_lines.to_csv(outputfile_partial + "_lines"+ ".csv")  # Generate CSV
    gdf_lines = convert_pd_to_gdf_lines(df_all_lines, simplified=False) # Set simplified = True to simplify lines
    gdf_lines.to_file(outputfile_partial + "_lines"+ ".geojson", driver="GeoJSON")  # Generate GeoJson

    # ----------- Generator -----------

    df_all_generators = df_all_generators[df_all_generators.columns.intersection(set(columns_generator))]
    df_all_generators.reset_index(drop=True, inplace=True)
    # df_all_generators.drop(columns = ["tags.fixme","tags.name:ar","tags.building","tags.barrier"], inplace = True, errors='ignore') # TODO: Should Probably drop columns inatead of selecting

    # Generate Files
    df_all_generators.to_csv(outputfile_partial + "_generators" + ".csv")  # Generate CSV
    gdf_generators = convert_pd_to_gdf(df_all_generators)
    gdf_generators.to_file(outputfile_partial +"_generators"+ ".geojson", driver="GeoJSON")  # Generate GeoJson
    
    # ----------- Cables -----------

    df_all_cables = df_all_cables[df_all_cables.columns.intersection(set(columns_cable))]
    df_all_cables.reset_index(drop=True, inplace=True)


    # Generate Files
    df_all_cables.to_csv(outputfile_partial + "_cables"+ ".csv")  # Generate CSV
    gdf_cables = convert_pd_to_gdf_lines(df_all_cables, simplified=False) # Set simplified = True to simplify lines
    gdf_cables.to_file(outputfile_partial + "_cables"+ ".geojson", driver="GeoJSON")  # Generate GeoJson


if __name__ == "__main__":
    process_data()
