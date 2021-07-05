from shapely.geometry import LineString, Point, Polygon
from iso_country_codes import AFRICA_CC
import pandas as pd
import numpy as np
import geopandas as gpd
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("../../scripts")


def clean_data():

    outputfile_partial = os.path.join(os.getcwd(), "data", "africa_all") # Output file directory

    raw_outputfile_partial = os.path.join(os.getcwd(), "data", "raw", "africa_all" + "_raw") # Output file directory

    # ----------- SUBSTATIONS -----------

    df_all_substations = gpd.read_file(raw_outputfile_partial + "_substations" + ".geojson").set_crs(epsg=4326, inplace=True)

    # TODO : Cleaning goes here

    df_all_substations.to_file(outputfile_partial + "_substations"+ ".geojson", driver="GeoJSON")


    # ----------- LINES -----------

    # Load raw data
    df_all_lines = gpd.read_file(raw_outputfile_partial + "_lines" + ".geojson").set_crs(epsg=4326, inplace=True)

    # TODO : Cleaning Goes here

    # Modification - create final dataframe layout
    df_all_lines = df_all_lines.rename(
        columns={
            "id": "line_id",
            "tags.voltage": "voltage",
            "tags.circuits": "circuits",
            "tags.cables": "cables",
            "tags.frequency": "tag_frequency",
            "tags.power": "tag_type",
            "lonlat": "geometry",
        }
    )

    # Add NaN as default
    df_all_lines["bus0"] = np.nan
    df_all_lines["bus1"] = np.nan
    df_all_lines["length"] = np.nan
    df_all_lines["underground"] = np.nan
    df_all_lines["under_construction"] = np.nan

    #Rearrange columns
    clist = ["line_id","bus0","bus1","voltage","circuits","length","underground",
         "under_construction","tag_type","tag_frequency","geometry"]

    df_all_lines = df_all_lines[clist]

    df_all_lines.to_file(outputfile_partial + "_lines"+ ".geojson", driver="GeoJSON")


    # ----------- Generator -----------


    df_all_generators = gpd.read_file(raw_outputfile_partial + "_generators" + ".geojson").set_crs(epsg=4326, inplace=True)


    # TODO : Cleaning goes here

    df_all_generators.to_file(outputfile_partial + "_generators"+ ".geojson", driver="GeoJSON")



    return None


if __name__ == "__main__":
    clean_data()
