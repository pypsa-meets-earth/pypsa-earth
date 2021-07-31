# This script does the following
# 1. Downloads OSM files for specified countries from Geofabrik
# 2. Filters files for substations, lines and generators
# 3. Process and clean data
# 4. Exports to CSV
# 5. Exports to GeoJson

""" OSM extraction script."""

import logging
import os
import shutil
import sys

import geopandas as gpd
from numpy import append
import pandas as pd
import requests
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
from esy.osmfilter import Node, Relation, Way # https://gitlab.com/dlr-ve-esy/esy-osmfilter/-/tree/master/
from esy.osmfilter import osm_info as osm_info
from esy.osmfilter import osm_pickle as osm_pickle
from esy.osmfilter import run_filter
from _helpers import _sets_path_to_root
from iso_country_codes import AFRICA_CC
from shapely.geometry import LineString, Point, Polygon
import hashlib

import logging
logging.basicConfig()
logger=logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# logger.setLevel(logging.WARNING)

_sets_path_to_root("pypsa-africa")
# os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) # move up to root directory

# Downloads PBF File for given Country Code

# feature_list = ["substation", "tower", "line", "generator", "cable"]
feature_list = ["substation", "line", "generator", "cable"]
# feature_list = ["substation", "line"]
# feature_list = ["line"]
# feature_list = ["substation"]

feature_category = {
    "substation": "node",
    "generator": "node",
    "line": "way",
    "tower": "node",
    "cable": "way"
}

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

    if not os.path.exists(PBF_inputfile):
        logger.info(f"{geofabrik_filename} downloading to {PBF_inputfile}")
        #  create data/osm directory
        os.makedirs(os.path.dirname(PBF_inputfile), exist_ok=True)
        with requests.get(geofabrik_url, stream=True, verify=False) as r:
            with open(PBF_inputfile, "wb") as f:
                shutil.copyfileobj(r.raw, f)
    
    if update is True:
        if verify_pbf(PBF_inputfile, geofabrik_url) is False:
            logger.warning(f"md5 mismatch, deleting {geofabrik_filename}")
            if os.path.exists(PBF_inputfile):
                os.remove(PBF_inputfile)

            download_pbf(country_code, update=False) # Only try downloading once


        
    return PBF_inputfile



verified_pbf =[]
def verify_pbf(PBF_inputfile, geofabrik_url):
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
def download_and_filter(feature, country_code, update=False):
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
    PBF_inputfile = download_pbf(country_code, update)

    filter_file_exists = False
    # json file for the Data dictionary
    JSON_outputfile = os.path.join(
        os.getcwd(), "data", "osm", country_code + "_power.json"
    )  # json file for the Elements dictionary is automatically written to "data/osm/Elements"+filename)

    if os.path.exists(JSON_outputfile):
        filter_file_exists = True

    if not os.path.exists(os.path.join(os.getcwd(), "data", "osm", "Elements", country_code + f"_{feature}s.json")):
        logger.warning("Element file not found so pre-filtering")
        filter_file_exists = False

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
        logger.info(f"Loading Pickle for {AFRICA_CC[country_code]}")
    else:
        create_elements = True
        if country_code not in pre_filtered: # Ensures pre-filter is not run everytime
            new_prefilter_data = True
            logger.info(f"pre-filtering {AFRICA_CC[country_code]} ")
            pre_filtered.append(country_code)
        else:
            new_prefilter_data = False
        logger.info(
            f"Creating  New {feature} Elements for {AFRICA_CC[country_code]}"
        ) 

    prefilter = {
        Node: {"power": feature_list},
        Way: {"power": feature_list},
        Relation: {"power": feature_list},
    }  # see https://dlr-ve-esy.gitlab.io/esy-osmfilter/filter.html for filter structures
    # HACKY: due to esy.osmfilter validation

    blackfilter = [
        ("", ""),
    ]

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

    print(new_prefilter_data, create_elements, feature, country_code)

    return feature_data

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
        logger.warning("refs column not found")  # TODO : do not create a hacky empty ref col
        print(df_way.columns)
        df_way[col] = pd.Series([], dtype=pd.StringDtype()).astype(
            float
        )  # create empty "refs" if not in dataframe

    # def look(ref):
    #     lonlat_row = []
    #     for r in ref:
    #         lonlat = tuple(Data["Node"][str(r)]["lonlat"])
    #         lonlat_row.append(lonlat)
    #     return lonlat_row

    # lonlat_list = df_way["refs"].apply(look)

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
    # lonlat_column = []
    # area_column = []
    way_polygon = list(map(lambda lonlat: Polygon(lonlat) if len(lonlat)>=3 else Point(lonlat[0]), lonlat_list))
    area_column = list(map(int,round(gpd.GeoSeries(way_polygon).set_crs("EPSG:4326").to_crs("EPSG:3857").area,-1)))

    def find_center_point(p):
        if p.geom_type == 'Polygon':
            center_point = p.centroid
        else:
            center_point = p
        return list((center_point.x, center_point.y))

    lonlat_column = list(map(find_center_point, way_polygon))

    # for lonlat in lonlat_list:
    #     if len(lonlat) >= 3:  # Minimum for a triangle
    #         way_polygon = Polygon(lonlat)
    #         # TODO : Do set_crs and to_crs for whole lonlat_list and not indivudually, expected big performance boost.
    #         polygon_area = int(
    #             round(
    #                 gpd.GeoSeries(way_polygon)
    #                 .set_crs("EPSG:4326")
    #                 .to_crs("EPSG:3857")
    #                 .area,
    #                 -1,
    #             )
    #         )  # nearest tens m2
    #         # print('{:g}'.format(float('{:.3g}'.format(float(polygon_area))))) # For significant numbers
    #         area_column.append(polygon_area)
    #         center_point = way_polygon.centroid
    #         lonlat_column.append(list((center_point.x, center_point.y)))
    #     else:
    #         area_column.append(0)
    #         center_point = lonlat[0]
    #         lonlat_column.append(list(center_point))

    # df_way.drop("refs", axis=1, inplace=True, errors="ignore")
    df_way.insert(0, "Area", area_column)
    df_way.insert(0, "lonlat", lonlat_column)


# Convert Ways to Line Coordinates


def convert_ways_lines(df_way, Data):
    lonlat_list = lonlat_lookup(df_way, Data)
    lonlat_column = lonlat_list
    df_way.insert(0, "lonlat", lonlat_column)

    way_linestring = map(lambda lonlats: LineString(lonlats), lonlat_list)
    length_column = gpd.GeoSeries(way_linestring).set_crs("EPSG:4326").to_crs("EPSG:3857").length
    # for lonlat in lonlat_list:
    #     way_linestring = LineString(lonlat)
    #     line_length = (
    #         gpd.GeoSeries(way_linestring)
    #         .set_crs("EPSG:4326")
    #         .to_crs("EPSG:3857")
    #         .length
    #     )
    #     length_column.append(float(line_length))

    df_way.insert(0, "Length", length_column)
    


# Convert Points Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf_nodes(df_way):
    gdf = gpd.GeoDataFrame(
        df_way, geometry=[Point(x, y) for x, y in df_way.lonlat], crs="EPSG:4326"
    )
    gdf.drop(columns=["lonlat"], inplace=True)
    return gdf


# Convert Lines Pandas Dataframe to GeoPandas Dataframe


def convert_pd_to_gdf_lines(df_way, simplified=False):
    # df_way["geometry"] = df_way["lonlat"].apply(lambda x: LineString(x))
    if simplified is True:
        df_way["geometry"] = df_way["geometry"].apply(
            lambda x: x.simplify(0.005, preserve_topology=False)
        )
    gdf = gpd.GeoDataFrame(df_way, geometry=[LineString(x) for x in df_way.lonlat], crs="EPSG:4326")
    gdf.drop(columns=["lonlat"], inplace=True)

    return gdf


def process_node_data(country_code, feature_data, feature):
    df_node, df_way, Data = convert_filtered_data_to_dfs(
        country_code, feature_data, feature
    )
    if not df_way.empty:
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

    if not df_way.empty:
        convert_ways_lines(df_way, Data)
    
    # Add Type Column
    df_way["Type"] = "Way"

    # Add Country Column
    df_way["Country"] = AFRICA_CC[country_code]
    return df_way


def process_data(update):
    columns_tower = [
        "id",
        "lonlat",
        "Area",
        "tags.power",
        "tags.tower",
        "tags.material",
        "tags.structure",
        "tags.operator",
        "tags.line_attachment",
        "tags.line_management",
        "tags.ref",
        "tags.height",
        "Type",
        "Country",
    ]

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
        "Country",
    ]

    columns_generator = [
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

    columns_cable = [
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
        "Country",
    ]

    def output_csv_geojson(df_all_feature, columns_feature, feature):
        outputfile_partial = os.path.join(os.getcwd(), "data", "raw", "africa_all" + "_raw")  # Output file directory

        if not os.path.exists(outputfile_partial):
            os.makedirs(os.path.dirname(outputfile_partial), exist_ok=True)  # create raw directory
        
        df_all_feature = df_all_feature[df_all_feature.columns.intersection(set(columns_feature))]
        df_all_feature.reset_index(drop=True, inplace=True)

        # Generate Files

        if df_all_feature.empty:
            logger.warning(f"All feature data frame empty for {feature}")
            return None

        df_all_feature.to_csv(outputfile_partial + f"_{feature}s" + ".csv")  # Generate CSV

        if feature_category[feature] == "way":
            gdf_feature = convert_pd_to_gdf_lines(df_all_feature)
        else:
            gdf_feature = convert_pd_to_gdf_nodes(df_all_feature)

        
        gdf_feature.to_file(
            outputfile_partial + f"_{feature}s" + ".geojson", driver="GeoJSON"
        )  # Generate GeoJson   

    # test_CC = {"NG": "nigeria"}
    # test_CC = {"DZ": "algeria"}
    # test_CC = {"BJ": "benin"}
    # AFRICA_CC = {list of all African countries that are imported from script -> iso_countries_codes}        
    for feature in feature_list:
        df_all_feature = pd.DataFrame()
        for country_code in AFRICA_CC.keys():
            # for country_code in AFRICA_CC.keys():
            # replace Africa_CC by test_CC to only download data for one country
            feature_data = download_and_filter(feature, country_code, update)
            
            if feature_category[feature] == "way":
                df_feature = process_line_data(country_code, feature_data, feature)
            else:
                df_feature = process_node_data(country_code, feature_data, feature)


            df_all_feature = pd.concat([df_all_feature, df_feature])

        if feature == "tower":
            output_csv_geojson(df_all_feature, columns_tower, "tower")
        if feature == "substation":
            output_csv_geojson(df_all_feature, columns_substation, "substation")
        if feature == "line":
            output_csv_geojson(df_all_feature, columns_line, "line")
        if feature == "generator":
            output_csv_geojson(df_all_feature, columns_generator, "generator")
        if feature == "cable":
            output_csv_geojson(df_all_feature, columns_cable, "cable")


if __name__ == "__main__":
    process_data(update=False) # Set update
