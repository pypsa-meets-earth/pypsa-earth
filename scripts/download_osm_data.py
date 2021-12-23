# This script does the following
# 1. Downloads OSM files for specified countries from Geofabrik
# 2. Filters files for substations, lines and generators
# 3. Process and clean data
# 4. Exports to CSV
# 5. Exports to GeoJson
# Disables pylint problem in this scripts
# pylint: disable=E1120
""" OSM extraction script."""
import hashlib
import json
import logging
import os
import pickle
import shutil
import sys

import geopandas as gpd
import pandas as pd
import requests
import urllib3
from _helpers import _sets_path_to_root
from _helpers import _to_csv_nafix
from _helpers import configure_logging
from esy.osmfilter import Node
from esy.osmfilter import osm_info as osm_info
from esy.osmfilter import osm_pickle as osm_pickle
from esy.osmfilter import Relation
from esy.osmfilter import run_filter
from esy.osmfilter import Way
from osm_data_config import continent_regions
from osm_data_config import continents
from osm_data_config import feature_category
from osm_data_config import feature_columns
from osm_data_config import iso_to_geofk_dict
from osm_data_config import world
from shapely import geometry
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon

# https://gitlab.com/dlr-ve-esy/esy-osmfilter/-/tree/master/

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# logging.basicConfig()
_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)
# logger.setLevel(logging.WARNING)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Downloads PBF File for given Country Code


def getContinentCountry(code):
    for continent in world:
        country = world[continent].get(code, 0)
        if country:
            return continent, country
    return continent, country


def download_pbf(country_code, update, verify):
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
    continent, country_name = getContinentCountry(country_code)
    # Filename for geofabrik
    geofabrik_filename = f"{country_name}-latest.osm.pbf"
    # https://download.geofabrik.de/africa/nigeria-latest.osm.pbf
    geofabrik_url = f"https://download.geofabrik.de/{continent}/{geofabrik_filename}"
    PBF_inputfile = os.path.join(os.getcwd(), "data", "osm", continent, "pbf",
                                 geofabrik_filename)  # Input filepath

    if not os.path.exists(PBF_inputfile):
        _logger.info(f"{geofabrik_filename} downloading to {PBF_inputfile}")
        #  create data/osm directory
        os.makedirs(os.path.dirname(PBF_inputfile), exist_ok=True)
        with requests.get(geofabrik_url, stream=True, verify=False) as r:
            with open(PBF_inputfile, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    if verify is True:
        if verify_pbf(PBF_inputfile, geofabrik_url, update) is False:
            _logger.warning(f"md5 mismatch, deleting {geofabrik_filename}")
            if os.path.exists(PBF_inputfile):
                os.remove(PBF_inputfile)

            # Only try downloading once
            download_pbf(country_code, update=False, verify=False)

    return PBF_inputfile


verified_pbf = []


def verify_pbf(PBF_inputfile, geofabrik_url, update):
    if PBF_inputfile in verified_pbf:
        return True

    geofabrik_md5_url = geofabrik_url + ".md5"
    PBF_md5file = PBF_inputfile + ".md5"

    def calculate_md5(fname):
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    if update is True or not os.path.exists(PBF_md5file):
        with requests.get(geofabrik_md5_url, stream=True, verify=False) as r:
            with open(PBF_md5file, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    local_md5 = calculate_md5(PBF_inputfile)

    with open(PBF_md5file) as f:
        contents = f.read()
        remote_md5 = contents.split()[0]

    if local_md5 == remote_md5:
        verified_pbf.append(PBF_inputfile)
        return True
    else:
        # print(local_md5, remote_md5)
        return False


pre_filtered = []


def download_and_filter(feature, country_code, update=False, verify=False):
    """
    Download OpenStreetMap raw file for selected tag.

    Apply pbf download and filter with esy.osmfilter selected OpenStreetMap
    tags or data. Examples of possible tags are listed at `OpenStreetMap wiki
    <https://wiki.openstreetmap.org/wiki/Key:power>`_.
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
    PBF_inputfile = download_pbf(country_code, update, verify)

    continent, country_name = getContinentCountry(country_code)

    # folder path
    folder_path = os.path.join(os.getcwd(), "data", "osm", continent)
    
    # path of the Data.pickle used by run_filter
    file_pickle = os.path.join(folder_path, "Data.pickle")

    # path of the backup file
    file_stored_pickle = os.path.join(folder_path, f"Data_{country_code}.pickle")


    # json file for the Data dictionary
    JSON_outputfile = os.path.join(folder_path, country_code + "_power.json")
    # json file for the Elements dictionary is automatically written to "data/osm/Elements"+filename)

    # check if data have already been processed or not
    filter_file_exists = False
    if os.path.exists(file_stored_pickle):
        filter_file_exists = True
    else:
        _logger.warning("Element file not found so pre-filtering")
        filter_file_exists = False

    elementname = f"{country_code}_{feature}s"

    new_prefilter_data = True

    # Load Previously Pre-Filtered Files
    if update is False and verify is False and filter_file_exists is True:
        create_elements = True  # Do not create elements again
        new_prefilter_data = False

        # copy backup pickle into Data.pickle
        shutil.copyfile(file_stored_pickle, file_pickle)

        _logger.info(f"Loading {feature} Pickle for {country_name}, iso-code: {country_code}")
        # feature_data = Data, Elements
        # return feature_data

    else:
        create_elements = True
        if country_code not in pre_filtered:  # Ensures pre-filter is not run everytime
            new_prefilter_data = True
            _logger.info(f"Pre-filtering {country_name} ")
            pre_filtered.append(country_code)
        
        _logger.info(f"Creating new {feature} Elements for {country_name}, iso-code: {country_code}")

    prefilter = {
        Node: {
            "power": feature_list
        },
        Way: {
            "power": feature_list
        },
        Relation: {
            "power": feature_list
        },
    }  # see https://dlr-ve-esy.gitlab.io/esy-osmfilter/filter.html for filter structures

    blackfilter = [
        ("", ""),
    ]

    whitefilter = [
        [
            ("power", feature),
        ],
    ]

    Data, Elements = run_filter(
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

    # if new_prefilter_data, the prefiltering is performed and the Data.pickle for the country is created;
    # save a backup file
    if new_prefilter_data:
        if os.path.exists(file_stored_pickle):
            os.remove(file_pickle)
        else:
            # rename and store pickle country
            os.rename(file_pickle, file_stored_pickle)

    logging.disable(logging.NOTSET
                    )  # Re-enable logging as run_filter disables logging.INFO
    _logger.info(
        f"Pre: {new_prefilter_data}, Elem: {create_elements}, for {feature} in {country_code}"
    )

    feature_data = Data, Elements

    return feature_data


# Convert Filtered Data, Elements to Pandas Dataframes


def convert_filtered_data_to_dfs(country_code, feature_data, feature):
    Data, Elements = feature_data
    elementname = f"{country_code}_{feature}s"
    df_way = pd.json_normalize(Elements[elementname]["Way"].values())
    df_node = pd.json_normalize(Elements[elementname]["Node"].values())
    return (df_node, df_way, Data)


# Lookup refs and convert to list of longlats


def lonlat_lookup(df_way, Data):

    if "refs" not in df_way.columns:
        _logger.warning("refs column not found")
        print(df_way.columns)
        # df_way[col] = pd.Series([], dtype=pd.StringDtype()).astype(float)  # create empty "refs" if not in dataframe

    def look(ref):
        lonlat_row = list(
            map(lambda r: tuple(Data["Node"][str(r)]["lonlat"]), ref))
        return lonlat_row

    lonlat_list = df_way["refs"].apply(look)

    return lonlat_list


# Convert Ways to Point Coordinates


def convert_ways_points(df_way, Data):
    lonlat_list = lonlat_lookup(df_way, Data)
    way_polygon = list(
        map(
            lambda lonlat: Polygon(lonlat)
            if len(lonlat) >= 3 else Point(lonlat[0]),
            lonlat_list,
        ))
    area_column = list(
        map(
            int,
            round(
                gpd.GeoSeries(way_polygon).set_crs("EPSG:4326").to_crs(
                    "EPSG:3857").area,
                -1,
            ),
        ))  # TODO: Rounding should be down in cleaning scripts

    def find_center_point(p):
        if p.geom_type == "Polygon":
            center_point = p.centroid
        else:
            center_point = p
        return list((center_point.x, center_point.y))

    lonlat_column = list(map(find_center_point, way_polygon))

    # df_way.drop("refs", axis=1, inplace=True, errors="ignore")
    df_way.insert(0, "Area", area_column)
    df_way.insert(0, "lonlat", lonlat_column)


# Convert Ways to Line Coordinates


def convert_ways_lines(df_way, Data):
    lonlat_list = lonlat_lookup(df_way, Data)
    lonlat_column = lonlat_list
    df_way.insert(0, "lonlat", lonlat_column)

    way_linestring = map(lambda lonlats: LineString(lonlats), lonlat_list)
    length_column = (gpd.GeoSeries(way_linestring).set_crs("EPSG:4326").to_crs(
        "EPSG:3857").length)

    df_way.insert(0, "Length", length_column)


# Convert Points Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf_nodes(df_way):
    gdf = gpd.GeoDataFrame(df_way,
                           geometry=[Point(x, y) for x, y in df_way.lonlat],
                           crs="EPSG:4326")
    gdf.drop(columns=["lonlat"], inplace=True)
    return gdf


# Convert Lines Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf_lines(df_way, simplified=False):
    # df_way["geometry"] = df_way["lonlat"].apply(lambda x: LineString(x))
    if simplified is True:
        df_way["geometry"] = df_way["geometry"].apply(
            lambda x: x.simplify(0.005, preserve_topology=False))
    gdf = gpd.GeoDataFrame(df_way,
                           geometry=[LineString(x) for x in df_way.lonlat],
                           crs="EPSG:4326")
    gdf.drop(columns=["lonlat"], inplace=True)

    return gdf


def convert_iso_to_geofk(iso_code,
                         iso_coding=True,
                         convert_dict=iso_to_geofk_dict):
    """Function to convert the iso code name of a country into the corresponding geofabrik"""
    if iso_code in convert_dict:
        if not iso_coding:
            _logger.error(
                f"Unexpected iso code {iso_code}: expected only geofabrik codes"
            )
        return convert_dict[iso_code]
    else:
        return iso_code


def output_csv_geojson(output_files, country_code, df_all_feature, columns_feature, feature):
    "Function to save the feature as csv and geojson"

    continent, country_name = getContinentCountry(country_code)

    path_file_geojson = output_files[feature + "s"]  # get path from snakemake; expected geojson
    if not path_file_geojson.endswith(".geojson"):
        _logger.error(
            f"Output file feature {feature} is not a geojson file"
        )
    path_file_csv = path_file_geojson.replace(".geojson", ".csv")  # get csv file
    
    # outputfile_partial = os.path.join(os.getcwd(), "data", "raw",
    #                                   continent + "_all"
    #                                   "_raw")  # Output file directory

    if not os.path.exists(path_file_geojson):
        os.makedirs(os.path.dirname(path_file_geojson),
                    exist_ok=True)  # create raw directory

    df_all_feature = df_all_feature[df_all_feature.columns.intersection(
        set(columns_feature))]
    df_all_feature.reset_index(drop=True, inplace=True)

    # Generate Files

    # _to_csv_nafix(df_all_feature,
    #               outputfile_partial + f"_{feature}s" + ".csv")  # Generate CSV
    _to_csv_nafix(df_all_feature, path_file_csv)  # Generate CSV

    if df_all_feature.empty:
        _logger.warning(f"Store empty Dataframe for {feature}.")
        gdf_feature = []

        # create empty file to avoid issues with snakemake
        with open(path_file_geojson, "w") as fp:
            pass

        return None

    if feature_category[feature] == "way":
        gdf_feature = convert_pd_to_gdf_lines(df_all_feature)
    else:
        gdf_feature = convert_pd_to_gdf_nodes(df_all_feature)

    _logger.info("Writing GeoJSON file")
    gdf_feature.to_file(path_file_geojson,
                        driver="GeoJSON")  # Generate GeoJson


def process_data(feature_list, country_list, output_files,
                iso_coding=True, update=False, verify=False):
    """
    Download the features in feature_list for each country of the country_list
    """

    # loop the request for each feature

    for feature in feature_list:  # feature dataframe
        
        df_list = []

        # initialize dataframe
        # df_all_feature = pd.DataFrame()

        for country_code_isogeofk in country_list:

            country_code = convert_iso_to_geofk(country_code_isogeofk,
                                                iso_coding)

            feature_data = download_and_filter(feature, country_code, update,
                                               verify)

            df_node, df_way, Data = convert_filtered_data_to_dfs(
                country_code, feature_data, feature)

            if feature_category[feature] == "way":
                convert_ways_lines(
                    df_way, Data) if not df_way.empty else _logger.warning(
                        f"Empty Way Dataframe for {feature} in {country_code}")
                if not df_node.empty:
                    _logger.warning(
                        f"Node dataframe not empty for {feature} in {country_code}"
                    )

            if feature_category[feature] == "node":
                convert_ways_points(df_way, Data) if not df_way.empty else None

            # Add Type Column
            df_node["Type"] = "Node"
            df_way["Type"] = "Way"

            # Concatinate Nodes and Ways
            df_feature = pd.concat([df_node, df_way], axis=0)

            # Add Country Column with GeoFabrik coding
            df_feature["Country"] = country_code

            df_list.append(df_feature)

            # df_all_feature = pd.concat([df_all_feature, df_feature], ignore_index=True)
        
        df_all_feature = pd.concat(df_list, ignore_index=True)

        output_csv_geojson(output_files, country_code, df_all_feature,
                           feature_columns[feature], feature)


def create_country_list(input, iso_coding=True):
    """
    Create a country list for defined regions in osm_data_config.py

    Parameters
    ----------
    input : str
        Any two-letter country name, regional name, or continent given in osm_data_config.py
        Country name duplications won't distort the result.
        Examples are:
        ["NG","ZA"], downloading osm data for Nigeria and South Africa
        ["africa"], downloading data for Africa
        ["NAR"], downloading data for the North African Power Pool
        ["TEST"], downloading data for a customized test set.
        ["NG","ZA","NG"], won't distort result.

    Returns
    -------
    full_codes_list : list
        Example ["NG","ZA"]
    """
    def filter_codes(c_list, iso_coding=True):
        """
        Filter list according to the specified coding.
        When iso code are implemented (iso_coding=True), then remove the geofabrik-specific ones.
        When geofabrik codes are selected(iso_coding=False), ignore iso-specific names.
        """
        if iso_coding:
            return [c for c in c_list if len(c) == 2]
        else:
            return [c for c in c_list if c not in iso_to_geofk_dict]

    full_codes_list = []

    for value1 in input:

        codes_list = []

        # extract countries in world
        if value1 == "world":
            for continent in world.keys():
                codes_list.extend(list(world[continent]))

        # extract countries in continent
        elif value1 in world.keys():
            codes_list = list(world[value1])

        # extract countries in regions
        elif value1 in continent_regions.keys():
            codes_list = continent_regions[value1]

        # extract countries
        else:
            codes_list.extend([value1])

        # create a list with all countries
        full_codes_list.extend(codes_list)

    # Removing duplicates and filter outputs by coding
    full_codes_list = filter_codes(set(full_codes_list), iso_coding=iso_coding)

    return full_codes_list


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake("download_osm_data")
    configure_logging(snakemake)

    # Required to set path to pypsa-africa
    _sets_path_to_root("pypsa-africa")

    # ["substation", "generator", "line", "cable", "tower"]
    feature_list = ["substation", "generator", "line", "cable"]

    input = snakemake.config["countries"]
    output_files = snakemake.output

    country_list = create_country_list(input)

    # Set update # Verify = True checks local md5s and pre-filters data again
    process_data(feature_list, country_list, output_files,
                iso_coding=True, update=False, verify=False)
