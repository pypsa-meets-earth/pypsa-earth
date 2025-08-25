# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import os

import geopandas as gpd
import numpy as np
import pandas as pd
import reverse_geocode as rg
from _helpers import (
    REGION_COLS,
    configure_logging,
    create_logger,
    save_to_geojson,
    to_csv_nafix,
)

logger = create_logger(__name__)


def prepare_substation_df(df_all_substations):
    """
    Prepare raw substations dataframe to the structure compatible with PyPSA-
    Eur.

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
            "Country": "country",
            "Area": "tag_area",
            "lonlat": "geometry",
        }
    )

    # Convert polygons to points
    df_all_substations["geometry"] = df_all_substations["geometry"].centroid

    # Add longitude (lon) and latitude (lat) coordinates in the dataset
    df_all_substations["lon"] = df_all_substations["geometry"].x
    df_all_substations["lat"] = df_all_substations["geometry"].y

    # Initialize columns to default value
    df_all_substations["dc"] = False
    df_all_substations["under_construction"] = False

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

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in df_all_substations:
            df_all_substations[c] = np.nan

    df_all_substations.drop(
        df_all_substations.columns[~df_all_substations.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    return df_all_substations


def add_line_endings_tosubstations(substations, lines):
    if lines.empty:
        return substations

    # extract columns from substation df
    bus_s = gpd.GeoDataFrame(columns=substations.columns, crs=substations.crs)
    bus_e = gpd.GeoDataFrame(columns=substations.columns, crs=substations.crs)

    # Read information from line.csv
    bus_s[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_s["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[0] if len(p.geoms) >= 2 else None
    )
    bus_s["lon"] = bus_s["geometry"].map(lambda p: p.x if p != None else None)
    bus_s["lat"] = bus_s["geometry"].map(lambda p: p.y if p != None else None)
    bus_s["bus_id"] = (
        (substations["bus_id"].max() if "bus_id" in substations else 0)
        + 1
        + bus_s.index
    )
    bus_s["dc"] = lines["dc"]

    bus_e[["voltage", "country"]] = lines[["voltage", "country"]].astype(str)
    bus_e["geometry"] = lines.geometry.boundary.map(
        lambda p: p.geoms[1] if len(p.geoms) >= 2 else None
    )
    bus_e["lon"] = bus_e["geometry"].map(lambda p: p.x if p != None else None)
    bus_e["lat"] = bus_e["geometry"].map(lambda p: p.y if p != None else None)
    bus_e["bus_id"] = bus_s["bus_id"].max() + 1 + bus_e.index
    bus_e["dc"] = lines["dc"]

    bus_all = pd.concat([bus_s, bus_e], ignore_index=True)

    # Initialize default values
    bus_all["station_id"] = np.nan
    # Assuming substations completed for installed lines
    bus_all["under_construction"] = False
    bus_all["tag_area"] = 0.0
    bus_all["symbol"] = "substation"
    # TODO: this tag may be improved, maybe depending on voltage levels
    bus_all["tag_substation"] = "transmission"

    buses = pd.concat([substations, bus_all], ignore_index=True)

    # Assign index to bus_id
    buses["bus_id"] = buses.index

    return buses


def set_unique_id(df, col):
    """
    Create unique id's, where id is specified by the column "col" The steps
    below create unique bus id's without losing the original OSM bus_id.

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


def split_cells(df, cols=["voltage"]):
    """
    Split semicolon separated cells i.e. [66000;220000] and create new
    identical rows.

    Parameters
    ----------
    df : dataframe
        Dataframe under analysis
    cols : list
        List of target columns over which to perform the analysis

    Example
    -------
    Original data:
    row 1: '66000;220000', '50'

    After applying split_cells():
    row 1, '66000', '50'
    row 2, '220000', '50'
    """
    if df.empty:
        return df

    x = df.assign(**{col: df[col].str.split(";") for col in cols})

    return x.explode(cols, ignore_index=True)


def filter_voltage(df, threshold_voltage=35000):
    """
    Filters df to contain only lines with voltage above threshold_voltage.
    """
    # Convert to numeric and drop any row with N/A voltage
    df["voltage"] = pd.to_numeric(df["voltage"], errors="coerce").astype(float)
    df.dropna(subset=["voltage"], inplace=True)

    # convert voltage to int
    df["voltage"] = df["voltage"].astype(int)

    # drop lines with a voltage lower than than threshold_voltage
    df.drop(
        df[df.voltage < threshold_voltage].index,
        axis=0,
        inplace=True,
        errors="ignore",
    )

    return df


def filter_frequency(df, accepted_values=[50, 60, 0], threshold=0.1):
    """
    Filters df to contain only lines with frequency with accepted_values.
    """
    df["tag_frequency"] = pd.to_numeric(df["tag_frequency"], errors="coerce").astype(
        float
    )
    df.dropna(subset=["tag_frequency"], inplace=True)

    accepted_rows = pd.concat(
        [(df["tag_frequency"] - f_val).abs() <= threshold for f_val in accepted_values],
        axis=1,
    ).any(axis="columns")

    df.drop(df[~accepted_rows].index, inplace=True)

    df["dc"] = df["tag_frequency"].abs() <= threshold

    return df


def filter_circuits(df, min_value_circuit=0.1):
    """
    Filters df to contain only lines with circuit value above
    min_value_circuit.
    """
    df["circuits"] = pd.to_numeric(df["circuits"], errors="coerce").astype(float)
    df.dropna(subset=["circuits"], inplace=True)

    accepted_rows = df["circuits"] >= min_value_circuit

    df.drop(df[~accepted_rows].index, inplace=True)

    return df


def finalize_substation_types(df_all_substations):
    """
    Specify bus_id and voltage columns as integer.
    """
    df_all_substations["bus_id"] = df_all_substations["bus_id"].astype(int)
    df_all_substations["voltage"] = df_all_substations["voltage"].astype(int)

    return df_all_substations


def prepare_lines_df(df_lines):
    """
    This function prepares the dataframe for lines and cables.

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
            "Country": "country",
            "Length": "length",
        }
    )

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
        "dc",
        "cables",
        "geometry",
        "country",
    ]

    # Check. If column is not in df create an empty one.
    for c in clist:
        if c not in df_lines:
            df_lines[c] = np.nan

    df_lines.drop(
        df_lines.columns[~df_lines.columns.isin(clist)],
        axis=1,
        inplace=True,
        errors="ignore",
    )

    return df_lines


def finalize_lines_type(df_lines):
    """
    This function is aimed at finalizing the type of the columns of the
    dataframe.
    """
    df_lines["line_id"] = df_lines["line_id"].astype(int)

    return df_lines


def clean_frequency(df, default_frequency="50"):
    """
    Function to clean raw frequency column: manual fixing and fill nan values
    """
    # replace raw values
    repl_freq = {
        "16.67": "16.7",
        "50;50;16.716.7": "50;50;16.7;16.7",
        "50;16.7?": "50;16.7",
        "50.0": "50",
        # "24 kHz": "24000",
    }

    # TODO: default frequency may be by country
    df["tag_frequency"] = (
        df["tag_frequency"].fillna(default_frequency).astype(str).replace(repl_freq)
    )

    return df


def clean_voltage(df):
    """
    Function to clean the raw voltage column: manual fixing and drop nan values
    """
    # replace raw values
    repl_voltage = {
        "medium": "33000",
        "19.1 kV": "19100",
        "high": "220000",
        "240 VAC": "240",
        "2*220000": "220000;220000",
        "KV30": "30kV",
    }

    df.dropna(subset=["voltage"], inplace=True)

    df["voltage"] = (
        df["voltage"]
        .astype(str)
        .replace(repl_voltage)
        .str.lower()
        .str.replace(" ", "")
        .str.replace("_", "")
        .str.replace("kv", "000")
        .str.replace("v", "")
        # .str.replace("/", ";")  # few OSM entries are separated by / instead of ;
        # this line can be a fix for that if relevant
    )

    return df


def clean_circuits(df):
    """
    Function to clean the raw circuits column: manual fixing and clean nan values
    """
    # replace raw values
    repl_circuits = {
        "1/3": "1",
        "2/3": "2",
        # assumption of two lines and one grounding wire
        "2-1": "2",
        "single": "1",
        "partial": "1",
        "1;1 disused": "1;0",
        # assuming that in case of a typo there is at least one line
        "`": "1",
        "^1": "1",
        "e": "1",
        "d": "1",
        "1.": "1",
    }

    # note: no string conversion for all entries in clean_circuits! it is performed later on
    df["circuits"] = (
        df["circuits"]
        .replace(repl_circuits)
        .map(lambda x: x.replace(" ", "") if isinstance(x, str) else x)
    )

    # Convert numbers in different dtypes to string while preserving NaN or other strings.
    is_numeric = ~pd.to_numeric(df["circuits"], errors="coerce").isna()
    df["circuits"] = df["circuits"].mask(is_numeric, df["circuits"].astype(str))

    # Report non-numeric and non-NaN values, which should be added to repl_circuits.
    if df.loc[~is_numeric, "circuits"].notna().any():
        logger.warning(
            "Non-numeric and non-NaN values found in circuits column, consider replacement: "
            + str(df.loc[~is_numeric, "circuits"].unique())
        )

    return df


def clean_cables(df):
    """
    Function to clean the raw cables column: manual fixing and drop undesired values
    """
    # replace raw values
    repl_cables = {
        "1 disused": "0",
        "ground": "0",
        "single": "1",
        "triple": "3",
        "3;3 disused": "3;0",
        "1 (Looped - Haul & Return) + 1 power wire": "1",
        "2-1": "3",
        "3+3": "6",
        "6+1": "6",
        "2x3": "6",
        "3x2": "6",
        "2x2": "4",
        # assuming that in case of a typo there is at least one line
        "partial": "1",
        "`": "1",
        "^1": "1",
        "e": "1",
        "d": "1",
        "line": "1",
    }

    df["cables"] = (
        df["cables"]
        .replace(repl_cables)
        .map(lambda x: x.replace(" ", "") if isinstance(x, str) else x)
    )

    # Convert numbers in different dtypes to string while preserving NaN or other strings.
    is_numeric = ~pd.to_numeric(df["cables"], errors="coerce").isna()
    df["cables"] = df["cables"].mask(is_numeric, df["cables"].astype(str))

    # Report non-numeric and non-NaN values, which should be added to repl_cables.
    if df.loc[~is_numeric, "cables"].notna().any():
        logger.warning(
            "Non-numeric and non-NaN values found in cables column, consider replacement: "
            + str(df.loc[~is_numeric, "cables"].unique())
        )

    return df


def split_and_match_voltage_frequency_size(df):
    """
    Function to match the length of the columns in subset by duplicating the
    last value in the column.

    The function does as follows:

    1. First, it splits voltage and frequency columns by semicolon
       For example, the following lines
       row 1: '50', '220000
       row 2: '50;50;50', '220000;380000'

       become:
       row 1: ['50'], ['220000']
       row 2: ['50','50','50'], ['220000','380000']

    2. Then, it harmonize each row to match the length of the lists
       by filling the missing values with the last elements of each list.
       In agreement to the example of before, after the cleaning:

       row 1: ['50'], ['220000']
       row 2: ['50','50','50'], ['220000','380000','380000']
    """
    for col in ["tag_frequency", "voltage"]:
        df[col] = df[col].str.split(";")

    len_freq = df["tag_frequency"].map(len)
    len_voltage = df["voltage"].map(len)

    def _fill_by_last(row, col_to_fill, size_col):
        """
        This functions takes a series and checks two elements in locations
        col_to_fill and size_col that are lists.

        The list of col_to_fill has less elements than of size_col. This
        function extends the col_to_fill element to match the size of
        size_col by replicating the last element as necessary.
        """
        size_to_fill = len(row[size_col])
        if not row[col_to_fill]:
            return None
        n_missing = size_to_fill - len(row[col_to_fill])
        fill_val = row[col_to_fill][-1]

        return row[col_to_fill] + [fill_val] * n_missing

    df_fhv = df[len_freq > len_voltage].apply(
        lambda row: _fill_by_last(row, "voltage", "tag_frequency"),
        axis=1,
    )
    df.loc[df_fhv.index, "voltage"] = df_fhv

    df_vhf = df[len_freq < len_voltage].apply(
        lambda row: _fill_by_last(row, "tag_frequency", "voltage"),
        axis=1,
    )
    df.loc[df_vhf.index, "tag_frequency"] = df_vhf

    return df


def fill_circuits(df):
    """
    This function fills the rows circuits column so that the size of each list
    element matches the size of the list in the frequency column.

    Multiple procedure are adopted:

    1. In the rows of circuits where the number of elements matches
       the number of the frequency column, nothing is done
    2. Where the number of elements in the cables column match the ones
       in the frequency column, then the values of cables are used.
    3. Where the number of elements in cables exceed those in frequency,
       the cables elements are downscaled and the last values of cables
       are summed.
       Let's assume that cables is [3,3,3] but frequency is [50,50].
       With this procedure, cables is treated as [3,6] and used for
       calculating the circuits
    4. Where the number in cables has an unique number, e.g. ['6'],
       but frequency does not, e.g. ['50', '50'],
       then distribute the cables proportionally across the values.
       Note: the distribution accounts for the frequency type;
       when the frequency is 50 or 60, then a circuit requires 3 cables,
       when DC (0 frequency) is used, a circuit requires 2 cables.
    5. Where no information of cables or circuits is available,
       a circuit is assumed for every frequency entry.
    """

    def _get_circuits_status(df):
        len_f = df["tag_frequency"].map(len)
        len_c = df["circuits"].map(
            lambda x: x.count(";") + 1 if isinstance(x, str) else np.nan
        )
        isna_c = df["circuits"].isna()
        len_cab = df["cables"].map(
            lambda x: x.count(";") + 1 if isinstance(x, str) else np.nan
        )
        isna_cab = df["cables"].isna()
        return len_f, len_c, isna_c, len_cab, isna_cab

    def _parse_float(x, ret_def=0.0):
        try:
            return float(x)
        except:
            return ret_def

    # cables requirement for circuits calculation
    cables_req = {"50": 3, "60": 3, "16.7": 2, "0": 2}

    def _basic_cables(f_val, cables_req=cables_req, def_circ=2):
        return cables_req[f_val] if f_val in cables_req.keys() else def_circ

    len_f, len_c, isna_c, len_cab, isna_cab = _get_circuits_status(df)

    is_numeric_cables = ~pd.to_numeric(df["cables"], errors="coerce").isna()

    to_fill = isna_c | (len_f != len_c)
    to_fill_direct = to_fill & (len_cab == len_f)
    to_fill_merge = to_fill & (len_cab > len_f)
    to_fill_indirect = to_fill & ~to_fill_direct & is_numeric_cables
    to_fill_default = to_fill & ~to_fill_merge & ~to_fill_direct & ~to_fill_indirect

    # length of cables match the frequency one
    # matching uses directly only the cables series
    df_match_by_cables = df[to_fill_direct][["tag_frequency", "cables"]].copy()
    df_match_by_cables["cables"] = (
        df_match_by_cables["cables"].astype(str).str.split(";")
    )

    def _filter_cables(row):
        return ";".join(
            [
                str(_parse_float(vc) / _basic_cables(vf))
                for (vc, vf) in zip(row["cables"], row["tag_frequency"])
            ]
        )

    df.loc[df_match_by_cables.index, "circuits"] = df_match_by_cables.apply(
        _filter_cables, axis=1
    )

    # length of cables elements is larger than frequency; the last cable data are merged to match
    df_merge_by_cables = df[to_fill_merge][["tag_frequency", "cables"]].copy()
    df_merge_by_cables["cables"] = (
        df_merge_by_cables["cables"].astype(str).str.split(";")
    )

    def _parse_cables_to_len(row):
        lf = len(row["tag_frequency"])
        float_cable = [_parse_float(vc) for vc in row["cables"]]
        parsed_cable_list = float_cable[0 : lf - 1] + [sum(float_cable[lf - 1 :])]

        return ";".join(
            [
                str(vc / _basic_cables(vf))
                for (vc, vf) in zip(parsed_cable_list, row["tag_frequency"])
            ]
        )

    df.loc[df_merge_by_cables.index, "circuits"] = df_merge_by_cables.apply(
        _parse_cables_to_len, axis=1
    )

    # indirect matching exploiting the total numeric value in cables
    # the minimum requirement of cables by frequency value is calculated
    # then using the total numeric cables number, the values are scaled proportionally
    df_indirect = df[to_fill_indirect]

    basic_cables = (
        df_indirect["tag_frequency"]
        .map(lambda x: [_basic_cables(v) for v in x])
        .rename("basic_cables")
    )
    min_cables = basic_cables.map(sum).rename("min_cables")
    multiplier = (
        pd.to_numeric(df_indirect["cables"], errors="coerce") / min_cables
    ).rename("multiplier")
    filled_values = pd.concat([basic_cables, multiplier], axis=1).apply(
        lambda x: ";".join([str(x["multiplier"] * v) for v in x["basic_cables"]]),
        axis=1,
    )
    df["circuits"] = df["circuits"].astype(str)
    df.loc[filled_values.index, "circuits"] = filled_values

    # otherwise assume a circuit per element
    df_fill_default = df[to_fill_default]
    df.loc[df_fill_default.index, "circuits"] = len_f.loc[df_fill_default.index].map(
        lambda x: ";".join(["1"] * x)
    )

    # explode column
    df["circuits"] = df["circuits"].astype(str).str.split(";")

    return df


def explode_rows(df, cols):
    """
    Function that explodes the rows as specified in cols, including warning
    alerts for unexpected values.

    Example
    --------
    row 1: [50,50], [33000, 110000]

    after explode_rows applied on the two columns becomes
    row 1: 50, 33000
    row 2: 50, 110000
    """
    # check if all row elements are list
    is_all_list = df[cols].map(lambda x: isinstance(x, list)).all(axis=1)
    if not is_all_list.all():
        df_nonlist = df[~is_all_list]
        logger.warning(
            f"Unexpected non-list values in dataframe; dropping rows. Needed fix for entries:\n{df_nonlist}"
        )
        df.drop(df_nonlist.index, inplace=True)

    # check if errors in the columns
    nunique_values = df[cols].map(len).nunique(axis=1)
    df_nunique = df[nunique_values != 1]
    if not df_nunique.empty:
        logger.warning(
            f"Improper explosion of dataframe entries; dropping rows. Needed fix for entries:\n{df_nunique}"
        )
        df.drop(df_nunique.index, inplace=True)

    df = df.explode(cols, ignore_index=True)

    return df


def integrate_lines_df(df_all_lines, distance_crs):
    """
    Function to add underground, under_construction, frequency and circuits.
    """
    # explode frequency and columns
    df = pd.DataFrame(df_all_lines)

    # preliminary raw parsing of raw columns
    clean_frequency(df)
    clean_voltage(df)
    clean_circuits(df)
    clean_cables(df)

    # analyse each row of voltage and frequency and match their content
    split_and_match_voltage_frequency_size(df)

    # fill the circuits column for explode
    fill_circuits(df)

    # Add under construction info
    # Default = False. No more information available atm
    df["under_construction"] = False

    # Add underground flag to check whether the line (cable) is underground
    # Simplified. If tag_type cable then underground is True
    df["underground"] = df["tag_type"] == "cable"

    # drop columns
    df.drop(columns=["tag_location", "cables"], errors="ignore", inplace=True)

    # explode rows
    df = explode_rows(df, ["tag_frequency", "voltage", "circuits"])

    return gpd.GeoDataFrame(df, crs=df_all_lines.crs)


def filter_lines_by_geometry(df_all_lines):
    if df_all_lines.empty:
        return df_all_lines
    # drop None geometries
    df_all_lines.dropna(subset=["geometry"], axis=0, inplace=True)

    # remove lines represented as Polygons
    df_all_lines = df_all_lines[df_all_lines.geometry.geom_type == "LineString"]

    return df_all_lines


def prepare_generators_df(df_all_generators):
    """
    Prepare the dataframe for generators.
    """
    # reset index
    df_all_generators = df_all_generators.reset_index(drop=True)

    check_fields_for_generators = ["tags.generator:output:electricity"]
    for field_to_add in check_fields_for_generators:
        if field_to_add not in df_all_generators.columns.tolist():
            df_all_generators[field_to_add] = ""

    df_all_generators = df_all_generators.rename(
        columns={
            "tags.generator:output:electricity": "power_output_MW",
            "tags.name": "name",
        }
    )

    # convert electricity column from string to float value
    # TODO: this filtering can be improved
    df_all_generators = df_all_generators[
        df_all_generators["power_output_MW"].astype(str).str.contains("MW")
    ]
    df_all_generators["power_output_MW"] = (
        df_all_generators["power_output_MW"]
        .astype(str)
        .str.extract("(\\d+)")
        .astype(float)
    )

    return df_all_generators


def find_first_overlap(geom, country_geoms, default_name):
    """
    Return the first index whose shape intersects the geometry.
    """
    for c_name, c_geom in country_geoms.items():
        if not geom.disjoint(c_geom):
            return c_name
    return default_name


def set_countryname_by_shape(
    df,
    ext_country_shapes,
    exclude_external=True,
    col_country="country",
):
    "Set the country name by the name shape"
    df[col_country] = [
        find_first_overlap(
            row["geometry"],
            ext_country_shapes,
            None if exclude_external else row[col_country],
        )
        for id, row in df.iterrows()
    ]
    df.dropna(subset=[col_country], inplace=True)
    return df


def create_extended_country_shapes(country_shapes, offshore_shapes, tolerance=0.01):
    """
    Obtain the extended country shape by merging on- and off-shore shapes.
    """

    merged_shapes = (
        gpd.GeoDataFrame(
            {
                "name": list(country_shapes.index),
                "geometry": [
                    (
                        c_geom.unary_union(offshore_shapes[c_code])
                        if c_code in offshore_shapes
                        else c_geom
                    )
                    for c_code, c_geom in country_shapes.items()
                ],
            },
            crs=country_shapes.crs,
        )
        .set_index("name")["geometry"]
        .buffer(tolerance)
    )

    return merged_shapes


def set_name_by_closestcity(df_all_generators, colname="name"):
    """
    Function to set the name column equal to the name of the closest city.
    """

    # get cities name
    list_cities = rg.search([g.coords[0] for g in df_all_generators.geometry])

    # replace name
    df_all_generators[colname] = [
        l["city"] + "_" + str(id) + " - " + c_code
        for (l, c_code, id) in zip(
            list_cities, df_all_generators.country, df_all_generators.index
        )
    ]

    return df_all_generators


def load_network_data(network_asset, data_options):
    """
    Function to check if OSM or custom data should be considered.

    The network_asset should be a string named "lines", "cables" or
    "substations".
    """

    # checks the options for loading data to be used based on the network_asset defined (lines/cables/substations)
    try:
        cleaning_data_options = data_options[f"use_custom_{network_asset}"]
        custom_path = data_options[f"path_custom_{network_asset}"]

    except:
        logger.error(
            f"Missing use_custom_{network_asset} or path_custom_{network_asset} options in the config file"
        )

    # creates a dataframe for the network_asset defined
    if cleaning_data_options == "custom_only":
        loaded_df = gpd.read_file(custom_path)

    elif cleaning_data_options == "add_custom":
        loaded_df1 = gpd.read_file(input_files[network_asset])
        loaded_df2 = gpd.read_file(custom_path)
        loaded_df = pd.concat([loaded_df1, loaded_df2], ignore_index=True)

    else:
        if cleaning_data_options != "OSM_only":
            logger.warning(
                f"Unrecognized option {data_options} for handling custom data of {network_asset}."
                + "Default OSM_only option used. Options available in clean_OSM_data_options configtable"
            )

        loaded_df = gpd.read_file(input_files[network_asset])

    # returns dataframe to be read in each section of the code depending on the component type (lines, substations or cables)
    return loaded_df


def clean_data(
    input_files,
    output_files,
    africa_shape,
    geo_crs,
    distance_crs,
    data_options,
    ext_country_shapes=None,
    names_by_shapes=True,
    tag_substation="transmission",
    threshold_voltage=35000,
    add_line_endings=True,
    generator_name_method="OSM",
):
    logger.info("Process OSM lines")

    if os.path.getsize(input_files["lines"]) > 0:
        # Load raw data lines
        df_lines = load_network_data("lines", data_options)

        # prepare lines dataframe and data types
        df_lines = prepare_lines_df(df_lines)
        df_lines = finalize_lines_type(df_lines)
    else:
        logger.info("No OSM lines")
        df_lines = gpd.GeoDataFrame()

    # initialize name of the final dataframe
    df_all_lines = df_lines

    # load cables only if data are stored
    if os.path.getsize(input_files["cables"]) > 0:
        logger.info("Add OSM cables to data")
        # Load raw data lines
        df_cables = load_network_data("cables", data_options)

        # prepare cables dataframe and data types
        df_cables = prepare_lines_df(df_cables)
        df_cables = finalize_lines_type(df_cables)

        # concatenate lines and cables in a single dataframe
        df_all_lines = pd.concat([df_lines, df_cables], ignore_index=True)
    else:
        logger.info("No OSM cables to add: skipping")

    if not df_all_lines.empty:
        # Add underground, under_construction, frequency and circuits columns to the dataframe
        # and drop corresponding unused columns
        df_all_lines = integrate_lines_df(df_all_lines, distance_crs)

        logger.info("Filter lines by voltage, frequency, circuits and geometry")

        # filter lines
        df_all_lines = filter_voltage(df_all_lines, threshold_voltage)
        df_all_lines = filter_frequency(df_all_lines)
        df_all_lines = filter_circuits(df_all_lines)
        df_all_lines = filter_lines_by_geometry(df_all_lines)

        logger.info("Select lines and cables in the region of interest")

        # drop lines crossing regions with and without the region under interest
        df_all_lines = df_all_lines[df_all_lines.geometry.boundary.within(africa_shape)]

        df_all_lines = gpd.GeoDataFrame(df_all_lines, geometry="geometry")

        # set the country name by the shape
        if names_by_shapes:
            logger.info("Setting lines country name using the GADM shapes")
            df_all_lines = set_countryname_by_shape(df_all_lines, ext_country_shapes)

        # set unique line ids
        df_all_lines = set_unique_id(df_all_lines, "line_id")

    # save lines output
    logger.info("Saving lines output")
    save_to_geojson(df_all_lines, output_files["lines"])

    # ----------- SUBSTATIONS -----------

    logger.info("Process OSM substations")

    if os.path.getsize(input_files["substations"]) > 0:
        df_all_substations = load_network_data("substations", data_options)

        # prepare dataset for substations
        df_all_substations = prepare_substation_df(df_all_substations)

        # filter substations by tag
        if tag_substation:  # if the string is not empty check it
            df_all_substations = df_all_substations[
                df_all_substations["tag_substation"] == tag_substation
            ]

        # clean voltage and make sure it is string
        df_all_substations = clean_voltage(df_all_substations)

        df_all_substations = gpd.GeoDataFrame(
            split_cells(pd.DataFrame(df_all_substations)),
            crs=df_all_substations.crs,
        )

        # add line endings if option is enabled
        if add_line_endings:
            df_all_substations = add_line_endings_tosubstations(
                df_all_substations, df_all_lines
            )

        # drop substations with nan geometry
        df_all_substations.dropna(subset=["geometry"], axis=0, inplace=True)

        # filter substation by voltage
        df_all_substations = filter_voltage(df_all_substations, threshold_voltage)

        # finalize dataframe types
        df_all_substations = finalize_substation_types(df_all_substations)

        # save to geojson file
        df_all_substations = gpd.GeoDataFrame(df_all_substations, geometry="geometry")

        if names_by_shapes:
            # set the country name by the shape
            logger.info("Setting substations country name using the GADM shapes")
            df_all_substations = set_countryname_by_shape(
                df_all_substations,
                ext_country_shapes,
            )

        # set unique bus ids
        df_all_substations = set_unique_id(df_all_substations, "bus_id")
    else:
        logger.info("No OSM substations")
        df_all_substations = gpd.GeoDataFrame()

    # save substations output
    logger.info("Saving substations output")
    save_to_geojson(df_all_substations, output_files["substations"])

    # ----------- GENERATORS -----------

    logger.info("Process OSM generators")

    if os.path.getsize(input_files["generators"]) > 0:
        df_all_generators = gpd.read_file(input_files["generators"])

        # prepare the generator dataset
        df_all_generators = prepare_generators_df(df_all_generators)

        if names_by_shapes:
            # set the country name by the shape
            logger.info("Setting generators country name using the GADM shapes")
            df_all_generators = set_countryname_by_shape(
                df_all_generators,
                ext_country_shapes,
                col_country="Country",
            )

        # set name tag by closest city when the value is nan
        if generator_name_method == "closest_city":
            logger.info("Setting unknown generators name using the closest city")
            df_all_generators = set_name_by_closestcity(df_all_generators)
    else:
        logger.info("No OSM generators")
        df_all_generators = gpd.GeoDataFrame()

    # save to csv
    to_csv_nafix(df_all_generators, output_files["generators_csv"])

    # save to geojson
    logger.info("Saving generators output")
    save_to_geojson(df_all_generators, output_files["generators"])

    return None


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("clean_osm_data")
    configure_logging(snakemake)

    tag_substation = snakemake.params.clean_osm_data_options["tag_substation"]
    threshold_voltage = snakemake.params.clean_osm_data_options["threshold_voltage"]
    names_by_shapes = snakemake.params.clean_osm_data_options["names_by_shapes"]
    add_line_endings = snakemake.params.clean_osm_data_options["add_line_endings"]
    generator_name_method = snakemake.params.clean_osm_data_options.get(
        "generator_name_method", "OSM"
    )
    offshore_shape_path = snakemake.input.offshore_shapes
    onshore_shape_path = snakemake.input.country_shapes
    geo_crs = snakemake.params.crs["geo_crs"]
    distance_crs = snakemake.params.crs["distance_crs"]
    data_options = snakemake.params["clean_osm_data_options"]

    input_files = snakemake.input
    output_files = snakemake.output

    africa_shape = gpd.read_file(snakemake.input.africa_shape)["geometry"].iloc[0]

    # only when country names are defined by shapes, load the info
    if names_by_shapes:
        country_shapes = gpd.read_file(onshore_shape_path).set_index("name")["geometry"]

        offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes)

        offshore_shapes = offshore_shapes.reindex(columns=REGION_COLS).set_index(
            "name"
        )["geometry"]

        ext_country_shapes = create_extended_country_shapes(
            country_shapes, offshore_shapes
        )

    else:
        ext_country_shapes = None

    clean_data(
        input_files,
        output_files,
        africa_shape,
        geo_crs,
        distance_crs,
        data_options,
        ext_country_shapes=ext_country_shapes,
        names_by_shapes=names_by_shapes,
        tag_substation=tag_substation,
        threshold_voltage=threshold_voltage,
        add_line_endings=add_line_endings,
        generator_name_method=generator_name_method,
    )
