
import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("../../scripts")

# from shapely.geometry import LineString, Point, Polygon
# from iso_country_codes import AFRICA_CC
import pandas as pd
import numpy as np
import geopandas as gpd


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
            "lonlat": "geometry"
        }
    )

    # Add longitute (lon) and latitude (lat) coordinates in the dataset
    df_all_substations["lon"] = df_all_substations["geometry"].x
    df_all_substations["lat"] = df_all_substations["geometry"].y

    # Add NaN as default
    df_all_substations["station_id"] = np.nan
    df_all_substations["dc"] = np.nan
    df_all_substations["under_construction"] = np.nan

    # Rearrange columns
    clist = ["bus_id", "station_id", "voltage", "dc",
            "symbol", "under_construction", "tag_substation",
            "tag_area", "lon", "lat", "geometry", "country"]
    df_all_substations = df_all_substations[clist]

    # add the under_construction and dc
    df_all_substations["under_construction"] = True
    df_all_substations["dc"] = False

    return df_all_substations


def set_unique_id(df, col):
    """
    Create unique id's, where id is specified by the column "col"
    The steps below create unique bus id's without loosing the original OSM bus_id 

    Unique bus_id are created by simply adding -1,-2,-3 to the original bus_id
    Every unique id gets a -1 
    If a bus_id exist i.e. three times it it will the counted by cumcount -1,-2,-3 making the id unique
    
    Parameters
    ----------
    df : dataframe
        Dataframe considered for the analysis
    col : str
        Column name for the analyses; examples: "bus_id" for substations or "line_id" for lines
    """
    
    if df[col].count() != df[col].nunique():  # operate only if id is not already unique (nunique counts unique values)
        df["cumcount"] = df.groupby([col]).cumcount()  # create cumcount column. Cumcount counts 0,1,2,3 the number of duplicates
        df["cumcount"] = df["cumcount"] + 1  # avoid 0 value for better understanding
        df[col] = df[col].astype(str) + "-" + df["cumcount"].values.astype(str)  # add cumcount to id to make id unique
        df.drop(columns = "cumcount", inplace=True)  # remove cumcount column
    
    return df


def split_cells(df, lst_col = 'voltage'):
    """
    Split semicolon separated cells i.e. [66000;220000] and create new identical rows
    
        Parameters
    ----------
    df : dataframe
        Dataframe under analysis
    lst_col : str
        Target column over which to perform the analysis
    """
    x = df.assign(**{lst_col:df[lst_col].str.split(';')})
    x = pd.DataFrame({
        col:np.repeat(x[col].values, x[lst_col].str.len())
        for col in x.columns.difference([lst_col])
        }).assign(**{lst_col:np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return x


def filter_voltage(df, threshold_voltage = 110000):
    
   # Drop any row with N/A voltage
    df = df.dropna(subset=['voltage']) 

    #Split semicolon separated cells i.e. [66000;220000] and create new identical rows
    df = split_cells(df)

     # Convert voltage to float, if impossible, discard row
    df['voltage'] = df['voltage'].apply(lambda x: pd.to_numeric(x, errors='coerce')).astype(float)
    df = df.dropna(subset=['voltage'])  # Drop any row with Voltage = N/A

    # convert voltage to int
    df.loc[:,"voltage"]  = df['voltage'].astype(int)

    # keep only lines with a voltage no lower than than threshold_voltage
    df = df[df.voltage >= threshold_voltage]

    return df


def finalize_substation_types(df_all_substations):
    """
    Specify bus_id and voltage columns as integer
    """

    # make float to integer
    df_all_substations["bus_id"] = df_all_substations["bus_id"].astype(int)
    df_all_substations.loc[:, "voltage"]  = df_all_substations['voltage'].astype(int)

    return df_all_substations


def prepare_lines_df(df_lines):
    """
    This function prepares the dataframe for lines and cables
    
    Parameters
    ----------
    df_lines : dataframe
        Raw lines or cables dataframe as downloaded from OpenStreetMap
    """

    # Modification - create final dataframe layout
    df_lines = df_lines.rename(
        columns = {
            "id": "line_id",
            "tags.voltage": "voltage",
            "tags.circuits": "circuits",
            "tags.cables": "cables",
            "tags.frequency": "tag_frequency",
            "tags.power": "tag_type",
            "lonlat": "geometry",
            "Country": "country",  # new/different to PyPSA-Eur
            "Length": "length"
        }
    )

    # Add NaN as default
    df_lines["bus0"] = np.nan
    df_lines["bus1"] = np.nan
    df_lines["underground"] = np.nan
    df_lines["under_construction"] = np.nan

    # Rearrange columns
    clist = ["line_id", "bus0", "bus1", "voltage", "circuits", "length", "underground",
            "under_construction", "tag_type", "tag_frequency", "cables", "geometry", "country"]
    df_lines = df_lines[clist]

    return df_lines


def finalize_lines_type(df_lines):
    """
    This function is aimed at finalizing the type of the columns of the dataframe
    """
    df_lines["line_id"] = df_lines["line_id"].astype(int)

    return df_lines


def integrate_lines_df(df_all_lines):
    """
    Function to add underground, under_construction, frequency and circuits
    """

    # Add under construction info
    df_all_lines["under_construction"] = False # default. Not more information atm available

    # Add underground flag to check whether the line (cable) is underground
    df_all_lines["underground"] = (df_all_lines["tag_type"] == "cable") # Simplified. If tag_type cable then underground is True. 

    # More information extractable for "underground" by looking at "tag_location".
    if 'tag_location' in df_all_lines: # drop column if exist
        df_all_lines.drop(columns = "tag_location", inplace=True)

    # Add frequency column
    df_all_lines["tag_frequency"] = 50

    # Add circuits information
    if df_all_lines["cables"].dtype != int: # if not int make int
        df_all_lines.loc[(df_all_lines["cables"] < "3") | df_all_lines["cables"].isna(), "cables"] = "0" #HERE. "0" if cables "None", "nan" or "1"
        df_all_lines["cables"] = df_all_lines["cables"].astype("int")
    if 4 or 5 in df_all_lines["cables"].values: # downgrade 4 and 5 cables to 3... 
        # Reason: 4 cables have 1 lighting protection cables, 5 cables has 2 LP cables - not transferring energy; 
        # see https://hackaday.com/2019/06/11/a-field-guide-to-transmission-lines/
        df_all_lines.loc[(df_all_lines["cables"] == 4) | (df_all_lines["cables"] == 5), "cables"] = 3 # where circuits are "0" make "1"
    df_all_lines.loc[df_all_lines["circuits"].isna(), "circuits"] = df_all_lines.loc[df_all_lines['circuits'].isna(), "cables"] / 3 # one circuit contains 3 cables
    df_all_lines["circuits"] = df_all_lines["circuits"].astype(int)
    df_all_lines.loc[(df_all_lines["circuits"] == "0") | (df_all_lines["circuits"] == 0), "circuits"] = 1 # where circuits are "0" make "1"

    if 'cables' in df_all_lines: # drop column if exist
        df_all_lines.drop(columns = "cables", inplace=True)

    return df_all_lines


def prepare_generators_df(df_all_generators):
    """
    Prepare the dataframe for generators
    """

    # reset index
    df_all_generators = df_all_generators.reset_index(drop=True)

    # rename columns
    df_all_generators = df_all_generators.rename(columns = {'tags.generator:output:electricity':"power_output_MW"})

    # convert electricity column from string to float value
    df_all_generators = df_all_generators[df_all_generators['power_output_MW'].astype(str).str.contains('MW')]
    df_all_generators['power_output_MW'] = df_all_generators['power_output_MW'].str.extract('(\d+)').astype(float)

    return df_all_generators


def clean_data(tag_substation = "transmission", threshold_voltage = 110000):

    outputfile_partial = os.path.join(os.getcwd(), "data", "clean", "africa_all")  # Output file directory

    raw_outputfile_partial = os.path.join(os.getcwd(), "data", "raw", "africa_all" + "_raw")  # Output file directory

    # ----------- SUBSTATIONS -----------

    df_all_substations = gpd.read_file(raw_outputfile_partial + "_substations" + ".geojson").set_crs(epsg=4326, inplace=True)
    
    # prepare dataset for substations
    df_all_substations = prepare_substation_df(df_all_substations)

    # filter substations by tag
    df_all_substations = df_all_substations[df_all_substations["tag_substation"] == tag_substation]

    # filter substation by voltage
    df_all_substations = filter_voltage(df_all_substations, threshold_voltage)

    # finalize dataframe types
    df_all_substations = finalize_substation_types(df_all_substations)

    # set unique bus ids
    df_all_substations = set_unique_id(df_all_substations, "bus_id")

    # save to csv file
    df_all_substations = gpd.GeoDataFrame(df_all_substations, geometry="geometry",crs="EPSG:4326")
    df_all_substations.to_file(outputfile_partial + "_substations"+ ".geojson", driver="GeoJSON") 

    # ----------- LINES AND CABLES -----------

    # Load raw data lines
    df_lines = gpd.read_file(raw_outputfile_partial + "_lines" + ".geojson").set_crs(epsg=4326, inplace=True)
    
    # prepare lines dataframe and data types
    df_lines = prepare_lines_df(df_lines)
    df_lines = finalize_lines_type(df_lines)
    
    # Load raw data lines
    df_cables = gpd.read_file(raw_outputfile_partial + "_cables" + ".geojson").set_crs(epsg=4326, inplace=True)

    # prepare cables dataframe and data types
    df_cables = prepare_lines_df(df_cables)
    df_cables = finalize_lines_type(df_cables)

    # concatenate lines and cables in a single dataframe
    df_all_lines = pd.concat([df_lines, df_cables])
    
    # Add underground, under_construction, frequency and circuits columns to the dataframe
    # and drop corresponding unused columns 
    df_all_lines = integrate_lines_df(df_all_lines)

    # filter lines by voltage
    df_all_lines = filter_voltage(df_all_lines, threshold_voltage)
    
    # set unique line ids
    df_all_lines = set_unique_id(df_all_lines, "line_id")

    df_all_lines = gpd.GeoDataFrame(df_all_lines, geometry="geometry",crs="EPSG:4326")
    df_all_lines.to_file(outputfile_partial + "_lines" + ".geojson", driver="GeoJSON")

    # ----------- Generator -----------


    df_all_generators = gpd.read_file(raw_outputfile_partial + "_generators" + ".geojson").set_crs(epsg=4326, inplace=True)

    # prepare the generator dataset
    df_all_generators = prepare_generators_df(df_all_generators)

    df_all_generators.to_file(outputfile_partial + "_generators"+ ".geojson", driver="GeoJSON")

    return None


if __name__ == "__main__":
    clean_data()
