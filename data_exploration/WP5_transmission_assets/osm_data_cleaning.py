from shapely.geometry import LineString, Point, Polygon
from iso_country_codes import AFRICA_CC
import pandas as pd
import numpy as np
import geopandas as gpd
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("../../scripts")

def prepare_substation_df(df_all_substations):
    """
    Prepare raw substations dataframe to the structure compatible with PyPSA-Eur

    Parameters
    ----------
    df_all_substations : dataframe
        Raw substations dataframe as downloaded from OpenStreetMap
    """

    # Modify the naming of the DataFrame columns to adapt to the PyPSA-Eur-like format
    df_all_substations = df_all_substations.rename(
        columns = {
            "id": "bus_id",
            "tags.voltage": "voltage",
            # "dc", will be added below
            "tags.power": "symbol",
            # "under_construction", will be added below     
            "tags.substation": "tag_substation",
            "Country": "country",  # new/different to PyPSA-Eur
            "Area": "tag_area",
            "lonlat": "geometry",
        }
    )

    # Add longitute (lon) and latitude (lat) coordinates in the dataset
    df_all_substations["lon"] = df_all_substations["geometry"].x
    df_all_substations["lat"] = df_all_substations["geometry"].y

    # Add NaN as default
    df_all_substations["station_id"] = np.nan
    df_all_substations["dc"] = np.nan
    df_all_substations["under_construction"] = np.nan

    #Rearrange columns
    clist = ["bus_id","station_id","voltage","dc","symbol","under_construction","tag_substation",
            "tag_area","lon", "lat", "geometry","country"]
    df_all_substations = df_all_substations[clist]

    # add the under_construction and dc
    df_all_substations["under_construction"] = True
    df_all_substations["dc"] = False


def set_unique_bus_id(df_all_substations):
    """
    Create unique bus id's
    The steps below create unique bus id's without loosing the original OSM bus_id 

    Unique bus_id are created by simply adding -1,-2,-3 to the original bus_id
    Every unique id gets a -1 
    If a bus_id exist i.e. three times it it will the counted by cumcount -1,-2,-3 making the id unique
    """
    if df_all_substations["bus_id"].count() != df_all_substations["bus_id"].nunique(): # operate only if line_id is not already unique (nunique counts unique values)
        df_all_substations["cumcount"] = df_all_substations.groupby(["bus_id"]).cumcount() # create cumcount column. Cumcount counts 0,1,2,3 the number of duplicates
        df_all_substations["cumcount"] = df_all_substations["cumcount"] + 1 # avoid 0 value for better understanding
        df_all_substations["bus_id"] = df_all_substations["bus_id"].astype(str) + "-" + df_all_substations["cumcount"].values.astype(str) # add cumcount to line_id to make line_id unique
        df_all_substations.drop(columns = "cumcount", inplace=True) # remove cumcount column


def filter_substations(df_all_substations, tag_substation = "transmission", threshold_voltage = 110000):
    """
    Filter substations according to substation type and voltage.
    Substations corresponding to non-numeric voltages or nan are excluded
    
    """

    # filter only entries corresponding to transmission elements ("transmission")
    df_all_substations = df_all_substations[df_all_substations["tag_substation"] == tag_substation]

    # Covnert voltage to float, or set to nan the value
    df_all_substations['voltage'] = df_all_substations['voltage'].apply(lambda x: pd.to_numeric(x, errors='coerce')).astype(float)

    # Drop any row with N/A voltage
    df_all_substations = df_all_substations.dropna(subset=['voltage']) 

    # keep only lines with a voltage no lower than than threshold_voltage
    df_all_substations = df_all_substations[df_all_substations.voltage >= threshold_voltage]


def finalize_substation_types(df_all_substations):
    """
    Specify bus_id and voltage columns as integer
    """
    # make float to integer
    df_all_substations["bus_id"] = df_all_substations["bus_id"].astype(int)
    df_all_substations.loc[:,"voltage"]  = df_all_substations['voltage'].astype(int)
    

def clean_data(tag_substation = "transmission", threshold_voltage = 110000):

    outputfile_partial = os.path.join(os.getcwd(), "data", "africa_all") # Output file directory

    raw_outputfile_partial = os.path.join(os.getcwd(), "data", "raw", "africa_all" + "_raw") # Output file directory

    # ----------- SUBSTATIONS -----------

    df_all_substations = gpd.read_file(raw_outputfile_partial + "_substations" + ".geojson").set_crs(epsg=4326, inplace=True)

    # TODO : Cleaning goes here

    # prepare dataset for substations
    prepare_substation_df(df_all_substations)

    # filter substations
    filter_substations(df_all_substations, tag_substation, threshold_voltage)

    # set unique bus ids
    set_unique_bus_id(df_all_substations)

    # finalize dataframe types
    finalize_substation_types(df_all_substations)

    # save to csv file
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
