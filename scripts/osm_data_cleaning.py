import logging
import math
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.ops import unary_union
from _helpers import _sets_path_to_root
from _helpers import _to_csv_nafix
from _helpers import configure_logging
from _helpers import _save_to_geojson

# from shapely.geometry import LineString, Point, Polygon
# from osm_data_config import AFRICA_CC

logger = logging.getLogger(__name__)

# Requirement to set path to filepath for execution
# os.chdir(os.path.dirname(os.path.abspath(__file__)))


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
        columns={
            "id": "bus_id",
            "tags.voltage": "voltage",
            # "dc", will be added below
            "tags.power": "symbol",
            # "under_construction", will be added below
            "tags.substation": "tag_substation",
            "Country": "country",  # new/different to PyPSA-Eur
            "Area": "tag_area",
            "lonlat": "geometry",
        })

    # Add longitute (lon) and latitude (lat) coordinates in the dataset
    df_all_substations["lon"] = df_all_substations["geometry"].x
    df_all_substations["lat"] = df_all_substations["geometry"].y

    # Add NaN as default
    df_all_substations["station_id"] = np.nan
    df_all_substations["dc"] = np.nan
    df_all_substations["under_construction"] = np.nan

    # Rearrange columns
    clist = [
        "bus_id",
        "station_id",
        "voltage",
        "dc",
        "symbol",
        "under_construction",
        "tag_substation",
        "tag_area",
        "lon",
        "lat",
        "geometry",
        "country",
    ]
    df_all_substations = df_all_substations[clist]

    # add the under_construction and dc
    df_all_substations["under_construction"] = False
    df_all_substations["dc"] = False

    return df_all_substations


def add_line_endings_tosubstations(substations, lines):
    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns)
    bus_e = gpd.GeoDataFrame(columns=substations.columns)

    # Read information from line.csv
    bus_s[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_s["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[0] if len(p.geoms) >= 2 else None)
    bus_s["lon"] = bus_s["geometry"].map(lambda p: p.x if p != None else None)
    bus_s["lat"] = bus_s["geometry"].map(lambda p: p.y if p != None else None)
    bus_s["bus_id"] = (
        (substations["bus_id"].max() if "bus_id" in substations else 0) + 1 +
        bus_s.index)

    bus_e[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_e["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[1] if len(p.geoms) >= 2 else None)
    bus_e["lon"] = bus_e["geometry"].map(lambda p: p.x if p != None else None)
    bus_e["lat"] = bus_e["geometry"].map(lambda p: p.y if p != None else None)
    bus_e["bus_id"] = bus_s["bus_id"].max() + 1 + bus_e.index

    bus_all = bus_s.append(bus_e).reset_index(drop=True)

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

    # Assign index to bus_id
    buses.loc[:, "bus_id"] = buses.index

    return buses


# # tol=0.01, around 700m at latitude 44.
# def set_substations_ids(buses, tol=0.02):
#     """
#     Function to set substations ids to buses, accounting for location tolerance

#     The algorithm is as follows:

#     1. initialize all substation ids to -1
#     2. if the current substation has been already visited [substation_id < 0], then skip the calculation
#     3. otherwise:
#         1. identify the substations within the specified tolerance (tol)
#         2. when all the substations in tolerance have substation_id < 0, then specify a new substation_id
#         3. otherwise, if one of the substation in tolerance has a substation_id >= 0, then set that substation_id to all the others;
#            in case of multiple substations with substation_ids >= 0, the first value is picked for all

#     """

#     buses["station_id"] = -1

#     station_id = 0
#     for i, row in buses.iterrows():
#         if buses.loc[i, "station_id"] >= 0:
#             continue

#         # get substations within tolerance
#         close_nodes = np.where(
#             buses.apply(
#                 lambda x: math.dist([row["lat"], row["lon"]],
#                                     [x["lat"], x["lon"]]) <= tol,
#                 axis=1,
#             ))[0]

#         if len(close_nodes) == 1:
#             # if only one substation is in tolerance, then the substation is the current one iì
#             # Note that the node cannot be with substation_id >= 0, given the preliminary check
#             # at the beginning of the for loop
#             buses.loc[buses.index[i], "station_id"] = station_id
#             # update station id
#             station_id += 1
#         else:
#             # several substations in tolerance

#             # get their ids
#             subset_substation_ids = buses.loc[buses.index[close_nodes],
#                                               "station_id"]
#             # check if all substation_ids are negative (<0)
#             all_neg = subset_substation_ids.max() < 0
#             # check if at least a substation_id is negative (<0)
#             some_neg = subset_substation_ids.min() < 0

#             if all_neg:
#                 # when all substation_ids are negative, then this is a new substation id
#                 # set the current station_id and increment the counter
#                 buses.loc[buses.index[close_nodes], "station_id"] = station_id
#                 station_id += 1
#             elif some_neg:
#                 # otherwise, when at least a substation_id is non-negative, then pick the first value
#                 # and set it to all the other substations within tolerance
#                 sub_id = -1
#                 for substation_id in subset_substation_ids:
#                     if substation_id >= 0:
#                         sub_id = substation_id
#                         break
#                 buses.loc[buses.index[close_nodes], "station_id"] = sub_id


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
    # operate only if id is not already unique (nunique counts unique values)
    if df[col].count() != df[col].nunique():
        # create cumcount column. Cumcount counts 0,1,2,3 the number of duplicates
        df["cumcount"] = df.groupby([col]).cumcount()
        # avoid 0 value for better understanding
        df["cumcount"] = df["cumcount"] + 1
        # add cumcount to id to make id unique
        df[col] = df[col].astype(str) + "-" + df["cumcount"].values.astype(str)
        # remove cumcount column
        df.drop(columns="cumcount", inplace=True)

    return df


def split_cells(df, lst_col="voltage"):
    """
    Split semicolon separated cells i.e. [66000;220000] and create new identical rows

    Parameters
    ----------
    df : dataframe
        Dataframe under analysis
    lst_col : str
        Target column over which to perform the analysis
    """
    x = df.assign(**{lst_col: df[lst_col].str.split(";")})
    x = pd.DataFrame({
        col: np.repeat(x[col].values, x[lst_col].str.len())
        for col in x.columns.difference([lst_col])
    }).assign(
        **{lst_col: np.concatenate(x[lst_col].values)})[x.columns.tolist()]
    return x


def filter_voltage(df, threshold_voltage=35000):

    # Drop any row with N/A voltage
    df = df.dropna(subset=["voltage"])

    # Split semicolon separated cells i.e. [66000;220000] and create new identical rows
    df = split_cells(df)

    # Convert voltage to float, if impossible, discard row
    df["voltage"] = (df["voltage"].apply(
        lambda x: pd.to_numeric(x, errors="coerce")).astype(float))
    df = df.dropna(subset=["voltage"])  # Drop any row with Voltage = N/A

    # convert voltage to int
    df.loc[:, "voltage"] = df["voltage"].astype(int)

    # keep only lines with a voltage no lower than than threshold_voltage
    df = df[df.voltage >= threshold_voltage]

    return df


def finalize_substation_types(df_all_substations):
    """
    Specify bus_id and voltage columns as integer
    """

    # make float to integer
    df_all_substations["bus_id"] = df_all_substations["bus_id"].astype(int)
    df_all_substations.loc[:,
                           "voltage"] = df_all_substations["voltage"].astype(
                               int)

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
        columns={
            "id": "line_id",
            "tags.voltage": "voltage",
            "tags.circuits": "circuits",
            "tags.cables": "cables",
            "tags.frequency": "tag_frequency",
            "tags.power": "tag_type",
            "lonlat": "geometry",
            "Country": "country",  # new/different to PyPSA-Eur
            "Length": "length",
        })

    # Add NaN as default
    df_lines["bus0"] = np.nan
    df_lines["bus1"] = np.nan
    df_lines["underground"] = np.nan
    df_lines["under_construction"] = np.nan

    # Rearrange columns
    clist = [
        "line_id",
        "bus0",
        "bus1",
        "voltage",
        "circuits",
        "length",
        "underground",
        "under_construction",
        "tag_type",
        "tag_frequency",
        "cables",
        "geometry",
        "country",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in df_lines:
            df_lines[c] = np.nan

    df_lines = df_lines[clist]

    return df_lines


def finalize_lines_type(df_lines):
    """
    This function is aimed at finalizing the type of the columns of the dataframe
    """
    df_lines["line_id"] = df_lines["line_id"].astype(int)

    return df_lines


# split function for cables and voltage
def split_cells_multiple(df, list_col=["cables", "voltage"]):
    # TODO: split multiple cell need fix!
    # One line with 'cable':{3,6} and 'voltage':{220000,110000} or
    # only 'voltage':{110000,220000, 3330000} needs to split in several lines
    # with new line_ID
    n_rows = df.shape[0]
    for i in range(n_rows):
        sub = df[list_col].iloc[i]  # for each cables and voltage
        if sub.notnull().all() == True:  # check not both empty
            # check both contain ";"
            if [";" in s for s in sub].count(True) == len(list_col):
                d = [s.split(";") for s in sub]  # split them
                r = df.loc[i, :].copy()
                df.loc[i, list_col[0]] = d[0][0]  # first split  [0]
                df.loc[i, list_col[1]] = d[1][0]
                r[list_col[0]] = d[0][1]  # second split [1]
                r[list_col[1]] = d[1][1]
                df = df.append(r)

    # if some columns still contain ";" then sum the values
    for cl_name in list_col:
        sel_ids = (df[cl_name].str.contains(";") == True)
        # TODO: split multiple cell need fix!
        df = df.drop(df.loc[sel_ids, cl_name].index, )
        # Original fix breaks with one input like [30;t], a string in the element
        # df.loc[sel_ids, cl_name] = (
        #   df.loc[sel_ids, cl_name].apply(
        #         lambda x: str(sum(map(int, x.split(";")))) if ";" in str(x) else x)
        # )
    return df  # return new frame


def integrate_lines_df(df_all_lines):
    """
    Function to add underground, under_construction, frequency and circuits
    """

    # Add under construction info
    # Default = False. No more information available atm
    df_all_lines["under_construction"] = False

    # Add underground flag to check whether the line (cable) is underground
    # Simplified. If tag_type cable then underground is True
    df_all_lines["underground"] = df_all_lines["tag_type"] == "cable"

    # More information extractable for "underground" by looking at "tag_location".
    if "tag_location" in df_all_lines:  # drop column if exist
        df_all_lines.drop(columns="tag_location", inplace=True)

    # Add frequency column
    df_all_lines["tag_frequency"] = 50

    df_all_lines = split_cells_multiple(df_all_lines)
    # Add circuits information
    # if not int make int
    if df_all_lines["cables"].dtype != int:
        # HERE. "0" if cables "None", "nan" or "1"
        df_all_lines.loc[(df_all_lines["cables"] < "3")
                         | df_all_lines["cables"].isna(), "cables"] = "0"
        df_all_lines["cables"] = df_all_lines["cables"].astype("int")

    # downgrade 4 and 5 cables to 3...
    if 4 or 5 in df_all_lines["cables"].values:
        # Reason: 4 cables have 1 lighting protection cables, 5 cables has 2 LP cables - not transferring energy;
        # see https://hackaday.com/2019/06/11/a-field-guide-to-transmission-lines/
        # where circuits are "0" make "1"
        df_all_lines.loc[(df_all_lines["cables"] == 4) |
                         (df_all_lines["cables"] == 5), "cables"] = 3

    # one circuit contains 3 cable
    df_all_lines.loc[df_all_lines["circuits"].isna(), "circuits"] = (
        df_all_lines.loc[df_all_lines["circuits"].isna(), "cables"] / 3)
    df_all_lines["circuits"] = df_all_lines["circuits"].astype(int)

    # where circuits are "0" make "1"
    df_all_lines.loc[(df_all_lines["circuits"] == "0") |
                     (df_all_lines["circuits"] == 0), "circuits"] = 1

    # drop column if exist
    if "cables" in df_all_lines:
        df_all_lines.drop(columns="cables", inplace=True)

    return df_all_lines


def filter_lines_by_geometry(df_all_lines):

    # drop None geometries
    df_all_lines.dropna(subset=["geometry"], axis=0, inplace=True)

    # remove lines without endings (Temporary fix for a Tanzanian line TODO: reformulation?)
    df_all_lines = df_all_lines[df_all_lines["geometry"].map(
        lambda g: len(g.boundary.geoms) >= 2)]

    return df_all_lines


def prepare_generators_df(df_all_generators):
    """
    Prepare the dataframe for generators
    """

    # reset index
    df_all_generators = df_all_generators.reset_index(drop=True)

    # rename columns
    df_all_generators = df_all_generators.rename(
        columns={
            "tags.generator:output:electricity": "power_output_MW",
            "tags.name": "name"})

    # convert electricity column from string to float value
    df_all_generators = df_all_generators[
        df_all_generators["power_output_MW"].astype(str).str.contains("MW")]
    df_all_generators["power_output_MW"] = (
        df_all_generators["power_output_MW"].str.extract("(\\d+)").astype(
            float))

    return df_all_generators


def find_first_overlap(geom, country_geoms, default_name):
    """Return the first country whose shape intersects the geometry"""
    for c_name, c_geom in country_geoms.items():
        if not geom.disjoint(c_geom):
            return c_name
    return default_name


def set_countryname_by_shape(df,
                             ext_country_shapes,
                             names_by_shapes=True,
                             exclude_external=True,
                             col_country="country"):
    "Set the country name by the name shape"
    if names_by_shapes:
        df[col_country] = [
            find_first_overlap(
                row["geometry"],
                ext_country_shapes,
                None if exclude_external else row[col_country],
            ) for id, row in df.iterrows()
        ]
        df.dropna(subset=[col_country], inplace=True)
    return df


def create_extended_country_shapes(country_shapes, offshore_shapes):
    """Obtain the extended country shape by merging on- and off-shore shapes"""

    merged_shapes = (gpd.GeoDataFrame({
        "name":
        list(country_shapes.index),
        "geometry": [
            c_geom.unary_union(offshore_shapes[c_code])
            if c_code in offshore_shapes else c_geom
            for c_code, c_geom in country_shapes.items()
        ],
    }).set_index("name")["geometry"].set_crs(4326))

    return merged_shapes


def clean_data(
    input_files,
    output_files,
    africa_shape,
    ext_country_shapes=None,
    names_by_shapes=True,
    tag_substation="transmission",
    threshold_voltage=35000,
    add_line_endings=True,
):

    # ----------- LINES AND CABLES -----------

    # Load raw data lines
    df_lines = gpd.read_file(input_files["lines"]).set_crs(epsg=4326,
                                                           inplace=True)

    # prepare lines dataframe and data types
    df_lines = prepare_lines_df(df_lines)
    df_lines = finalize_lines_type(df_lines)

    # initialize name of the final dataframe
    df_all_lines = df_lines

    # load cables only if data are stored
    if os.path.getsize(input_files["cables"]) > 0:
        # Load raw data lines
        df_cables = gpd.read_file(input_files["cables"]).set_crs(epsg=4326,
                                                                 inplace=True)

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

    # filter lines to make sure the geometry is appropriate
    df_all_lines = filter_lines_by_geometry(df_all_lines)

    # drop lines crossing regions with and without the region under interest
    df_all_lines = df_all_lines[
        df_all_lines.apply(lambda x: africa_shape.contains(
            x.geometry.boundary), axis=1)
    ]

    # set unique line ids
    df_all_lines = set_unique_id(df_all_lines, "line_id")

    df_all_lines = gpd.GeoDataFrame(df_all_lines,
                                    geometry="geometry",
                                    crs="EPSG:4326")

    # set the country name by the shape
    df_all_lines = set_countryname_by_shape(df_all_lines,
                                            ext_country_shapes,
                                            names_by_shapes=names_by_shapes)

    _save_to_geojson(df_all_lines, output_files["lines"])

    # ----------- SUBSTATIONS -----------

    df_all_substations = gpd.read_file(input_files["substations"]).set_crs(
        epsg=4326, inplace=True)

    # prepare dataset for substations
    df_all_substations = prepare_substation_df(df_all_substations)

    # add line endings if option is enabled
    if add_line_endings:
        df_all_substations = add_line_endings_tosubstations(
            gpd.GeoDataFrame(), df_all_lines)

    # filter substations by tag
    if tag_substation:  # if the string is not empty check it
        df_all_substations = df_all_substations[
            df_all_substations["tag_substation"] == tag_substation]

    # filter substation by voltage
    df_all_substations = filter_voltage(df_all_substations, threshold_voltage)

    # # set substation_id
    # set_substations_ids(df_all_substations, tol=0.02)

    # finalize dataframe types
    df_all_substations = finalize_substation_types(df_all_substations)

    # set unique bus ids
    df_all_substations = set_unique_id(df_all_substations, "bus_id")

    # save to geojson file
    df_all_substations = gpd.GeoDataFrame(df_all_substations,
                                          geometry="geometry",
                                          crs="EPSG:4326")

    # set the country name by the shape
    df_all_substations = set_countryname_by_shape(
        df_all_substations,
        ext_country_shapes,
        names_by_shapes=names_by_shapes,
        col_country="Country")

    _save_to_geojson(df_all_substations, output_files["substations"])

    # ----------- GENERATORS -----------

    df_all_generators = gpd.read_file(input_files["generators"]).set_crs(
        epsg=4326, inplace=True)

    # prepare the generator dataset
    df_all_generators = prepare_generators_df(df_all_generators)

    # set the country name by the shape
    df_all_generators = set_countryname_by_shape(
        df_all_generators,
        ext_country_shapes,
        names_by_shapes=names_by_shapes)

    # save to csv
    _to_csv_nafix(df_all_generators, output_files["generators_csv"])

    # save to geojson
    _save_to_geojson(df_all_generators, output_files["generators"])

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        # needed to run mock_snakemake
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake("clean_osm_data")
    configure_logging(snakemake)

    # Required to set path to pypsa-africa
    # _sets_path_to_root("pypsa-africa")

    tag_substation = snakemake.config["osm_data_cleaning_options"][
        "tag_substation"]
    threshold_voltage = snakemake.config["osm_data_cleaning_options"][
        "threshold_voltage"]
    names_by_shapes = snakemake.config["osm_data_cleaning_options"][
        "names_by_shapes"]
    add_line_endings = snakemake.config["osm_data_cleaning_options"][
        "add_line_endings"]

    input_files = snakemake.input
    output_files = snakemake.output

    africa_shape = (gpd.read_file(
        snakemake.input.africa_shape).set_crs(4326)["geometry"].iloc[0])

    # only when country names are defined by shapes, load the info
    if names_by_shapes:
        country_shapes = (gpd.read_file(
            snakemake.input.country_shapes).set_index("name")
            ["geometry"].set_crs(4326))
        offshore_shapes = (gpd.read_file(
            snakemake.input.offshore_shapes).set_index("name")
            ["geometry"].set_crs(4326))
        ext_country_shapes = create_extended_country_shapes(
            country_shapes, offshore_shapes)
    else:
        ext_country_shapes = None

    clean_data(
        input_files,
        output_files,
        africa_shape,
        ext_country_shapes=ext_country_shapes,
        names_by_shapes=names_by_shapes,
        tag_substation=tag_substation,
        threshold_voltage=threshold_voltage,
        add_line_endings=add_line_endings,
    )
