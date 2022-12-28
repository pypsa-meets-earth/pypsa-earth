from earth_osm import eo
import os
from pathlib import Path
from _helpers import configure_logging
import glob
import pandas as pd
import geopandas as gpd

from config_osm_data import (
    iso_to_geofk_dict,
)


def country_list_to_geofk(country_list):
    """
    Convert the requested country list into geofk norm

    Parameters
    ----------
    input : str
        Any two-letter country name or aggregation of countries given in config_osm_data.py
        Country name duplications won't distort the result.
        Examples are:
        ["NG","ZA"], downloading osm data for Nigeria and South Africa
        ["SNGM"], downloading data for Senegal&Gambia shape
        ["NG","ZA","NG"], won't distort result.

    Returns
    -------
    full_codes_list : list
        Example ["NG","ZA"]
    """

    full_codes_list = [convert_iso_to_geofk(c_code) for c_code in set(country_list)]

    return full_codes_list


def convert_iso_to_geofk(iso_code, iso_coding=True, convert_dict=iso_to_geofk_dict):
    """
    Function to convert the iso code name of a country into the corresponding geofabrik
    In Geofabrik, some countries are aggregated, thus if a single country is requested,
    then all the agglomeration shall be downloaded
    For example, Senegal (SN) and Gambia (GM) cannot be found alone in geofabrik,
    but they can be downloaded as a whole SNGM

    The conversion directory, initialized to iso_to_geofk_dict is used to perform such conversion
    When a two-letter code country is found in convert_dict, and iso_coding is enabled,
    then that two-letter code is converted into the corresponding value of the dictionary

    Parameters
    ----------
    iso_code : str
        Two-code country code to be converted
    iso_coding : bool
        When true, the iso to geofk is performed
    convert_dict : dict
        Dictionary used to apply the conversion iso to geofk
        The keys correspond to the countries iso codes that need a different region to be downloaded
    """
    if iso_coding and iso_code in convert_dict:
        return convert_dict[iso_code]
    else:
        return iso_code


def concat_outputs(store_path_data, store_path_resources):
  names = ["generators", "cables", "lines", "substations"]
  # create output folder if it does not exist
  os.makedirs(str(store_path_resources), exist_ok=True)
  for n in names:  
    filter = Path.joinpath(store_path_data, "out", f"*{n}.csv")
    filenames = [i for i in glob.glob(str(filter))]
    path = Path.joinpath(store_path_resources, f"africa_all_raw_{n}.csv")
    if len(filenames) > 0:
      concat_csv(filenames).to_csv(path, index=False)
    else:
      open(path, "a").close()

    filter = Path.joinpath(store_path_data, "out", f"*{n}.geojson")
    filenames = [i for i in glob.glob(str(filter))]
    path = Path.joinpath(store_path_resources, f"africa_all_raw_{n}.geojson")
    if len(filenames) > 0:
      concat_geojson(filenames).to_file(path, driver="GeoJSON")
    else:
      open(path, "a").close()


def concat_csv(filenames):
    # Read the files into a list of DataFrames
    data = list(map(pd.read_csv, filenames))
    # Concatenate the DataFrames into a single DataFrame
    result = pd.concat(data)
    return result


def concat_geojson(filenames):
    # Read the files into a list of GeoDataFrames\
    assert len(filenames) > 0, "No files to concatenate"
    data = list(map(gpd.read_file, filenames))
    # Concatenate the GeoDataFrames into a single GeoDataFrame
    result = gpd.GeoDataFrame(pd.concat(data, ignore_index=True))
    return result


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake, sets_path_to_root
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("download_osm_data") 
        sets_path_to_root("pypsa-earth")   
    configure_logging(snakemake)

    store_path_data = Path.joinpath(Path().cwd(), "data", "osm")
    store_path_resources = Path.joinpath(Path().cwd(), "resources", "osm", "raw")

    # get list of countries into geofabrik convention; expected iso norm in input
    country_list = country_list_to_geofk(snakemake.config["countries"])

    eo.get_osm_data(
      primary_name = 'power',
      region_list = country_list,
      feature_list = ['substation', 'line', 'cable', 'generator'],
      update = False,
      mp = True,
      data_dir = store_path_data,
    )

    concat_outputs(store_path_data, store_path_resources)
