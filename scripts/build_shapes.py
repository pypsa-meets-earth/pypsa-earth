# SPDX-FileCopyrightText: : 2021 PyPSA-Africa authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
import logging
import multiprocessing as mp
import os
import shutil
import zipfile
from itertools import takewhile
from operator import attrgetter

import fiona
import geopandas as gpd
import numpy as np
import rasterio
import requests
import rioxarray as rx
import xarray as xr
from _helpers import sets_path_to_root
from _helpers import three_2_two_digits_country
from _helpers import two_2_three_digits_country
from _helpers import two_digits_2_name_country
from _helpers import configure_logging
from rasterio.mask import mask
from shapely.geometry import LineString
from shapely.geometry import MultiPolygon
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from shapely.validation import make_valid
from tqdm import tqdm

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)

sets_path_to_root("pypsa-africa")


def download_GADM(country_code, update=False, out_logging=False):
    """
    Download gpkg file from GADM for a given country code

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files
    update : bool
        Update = true, forces re-download of files

    Returns
    -------
    gpkg file per country

    """

    GADM_filename = f"gadm36_{two_2_three_digits_country(country_code)}"
    GADM_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/{GADM_filename}_gpkg.zip"

    GADM_inputfile_zip = os.path.join(
        os.getcwd(),
        "data",
        "raw",
        "gadm",
        GADM_filename,
        GADM_filename + ".zip",
    )  # Input filepath zip

    GADM_inputfile_gpkg = os.path.join(
        os.getcwd(),
        "data",
        "raw",
        "gadm",
        GADM_filename,
        GADM_filename + ".gpkg",
    )  # Input filepath gpkg

    if not os.path.exists(GADM_inputfile_gpkg) or update is True:
        if out_logging:
            _logger.warning(
                f"Stage 4/4: {GADM_filename} of country {two_digits_2_name_country(country_code)} does not exist, downloading to {GADM_inputfile_zip}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_zip), exist_ok=True)

        with requests.get(GADM_url, stream=True) as r:
            with open(GADM_inputfile_zip, "wb") as f:
                shutil.copyfileobj(r.raw, f)

        with zipfile.ZipFile(GADM_inputfile_zip, "r") as zip_ref:
            zip_ref.extractall(os.path.dirname(GADM_inputfile_zip))

    return GADM_inputfile_gpkg, GADM_filename


def get_GADM_layer(country_list, layer_id, update=False, outlogging=False):
    """
    Function to retrive a specific layer id of a geopackage for a selection of countries

    Parameters
    ----------
    country_list : str
        List of the countries
    layer_id : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested

    """
    # initialization of the geoDataFrame
    geodf_GADM = gpd.GeoDataFrame()

    for country_code in country_list:
        # download file gpkg
        file_gpkg, name_file = download_GADM(country_code, update, outlogging)

        # get layers of a geopackage
        list_layers = fiona.listlayers(file_gpkg)

        # get layer name
        if layer_id < 0 | layer_id >= len(list_layers):
            # when layer id is negative or larger than the number of layers, select the last layer
            layer_id = len(list_layers) - 1
        code_layer = np.mod(layer_id, len(list_layers))
        layer_name = (
            f"gadm36_{two_2_three_digits_country(country_code).upper()}_{code_layer}"
        )

        # read gpkg file
        geodf_temp = gpd.read_file(file_gpkg, layer=layer_name)

        # convert country name representation of the main country (GID_0 column)
        geodf_temp["GID_0"] = [
            three_2_two_digits_country(twoD_c)
            for twoD_c in geodf_temp["GID_0"]
        ]

        # create a subindex column that is useful
        # in the GADM processing of sub-national zones
        geodf_temp["GADM_ID"] = geodf_temp[f"GID_{code_layer}"]

        # append geodataframes
        geodf_GADM = geodf_GADM.append(geodf_temp)

    geodf_GADM.reset_index(drop=True, inplace=True)

    return geodf_GADM


def _simplify_polys(polys,
                    minarea=0.0001,
                    tolerance=0.008,
                    filterremote=False):
    "Function to simplify the shape polygons"
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon([
                p for p in takewhile(lambda p: p.area > minarea, polys)
                if not filterremote or (mainpoly.distance(p) < mainlength)
            ])
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)


def countries(countries, update=False, out_logging=False):
    "Create country shapes"

    if out_logging:
        _logger.info("Stage 1 of 4: Create country shapes")

    # download data if needed and get the layer id 0, corresponding to the countries
    df_countries = get_GADM_layer(countries, 0, update, out_logging)

    # select and rename columns
    df_countries = df_countries[["GID_0", "geometry"]].copy()
    df_countries.rename(columns={"GID_0": "name"}, inplace=True)

    # set index and simplify polygons
    ret_df = df_countries.set_index("name")["geometry"].map(_simplify_polys)

    return ret_df


def country_cover(country_shapes,
                  eez_shapes=None,
                  out_logging=False,
                  distance=0.1):

    if out_logging:
        _logger.info(
            "Stage 3 of 4: Merge country shapes to create continent shape")

    shapes = country_shapes.apply(lambda x: x.buffer(distance))
    shapes_list = list(shapes)
    if eez_shapes is not None:
        shapes_list += list(eez_shapes)

    africa_shape = unary_union(shapes_list)

    return africa_shape


def save_to_geojson(df, fn):
    if os.path.exists(fn):
        os.unlink(fn)  # remove file if it exists
    if not isinstance(df, gpd.GeoDataFrame):
        df = gpd.GeoDataFrame(dict(geometry=df))

    # save file if the GeoDataFrame is non-empty
    if df.shape[0] > 0:
        df = df.reset_index()
        schema = {**gpd.io.file.infer_schema(df), "geometry": "Unknown"}
        df.to_file(fn, driver="GeoJSON", schema=schema)
    else:
        # create empty file to avoid issues with snakemake
        with open(fn, "w") as fp:
            pass


def load_EEZ(countries_codes, EEZ_gpkg="./data/raw/eez/eez_v11.gpkg"):
    """
    Function to load the database of the Exclusive Economic Zones.
    The dataset shall be downloaded independently by the user (see guide) or
    together with pypsa-africa package.
    """
    if not os.path.exists(EEZ_gpkg):
        raise Exception(
            f"File EEZ {EEZ_gpkg} not found, please download it from https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip and copy it in {os.path.dirname(EEZ_gpkg)}"
        )

    geodf_EEZ = gpd.read_file(EEZ_gpkg)
    geodf_EEZ.dropna(axis=0, how="any", subset=["ISO_TER1"], inplace=True)
    # [["ISO_TER1", "TERRITORY1", "ISO_SOV1", "ISO_SOV2", "ISO_SOV3", "geometry"]]
    geodf_EEZ = geodf_EEZ[["ISO_TER1", "geometry"]]
    selected_countries_codes_3D = [
        two_2_three_digits_country(x) for x in countries_codes
    ]
    geodf_EEZ = geodf_EEZ[[
        any([x in selected_countries_codes_3D]) for x in geodf_EEZ["ISO_TER1"]
    ]]
    geodf_EEZ["ISO_TER1"] = geodf_EEZ["ISO_TER1"].map(
        lambda x: three_2_two_digits_country(x))
    geodf_EEZ.reset_index(drop=True, inplace=True)

    geodf_EEZ.rename(columns={"ISO_TER1": "name"}, inplace=True)

    return geodf_EEZ


def eez(countries, country_shapes, EEZ_gpkg, out_logging=False, distance=0.01):
    """
    Creates offshore shapes by
    - buffer smooth countryshape (=offset country shape)
    - and differ that with the offshore shape
    Leads to for instance a 100m non-build coastline

    """

    if out_logging:
        _logger.info("Stage 2 of 4: Create offshore shapes")

    # load data
    df_eez = load_EEZ(countries, EEZ_gpkg)

    ret_df = df_eez[["name", "geometry"]]
    # create unique shape if country is described by multiple shapes
    for c_code in countries:
        selection = ret_df.name == c_code
        n_offshore_shapes = selection.sum()

        if n_offshore_shapes > 1:
            # when multiple shapes per country, then merge polygons
            geom = ret_df[selection].geometry.unary_union
            ret_df.drop(ret_df[selection].index, inplace=True)
            ret_df = ret_df.append({
                "name": c_code,
                "geometry": geom
            },
                                   ignore_index=True)

    ret_df = ret_df.set_index("name")["geometry"].map(
        lambda x: _simplify_polys(x, minarea=0.001, tolerance=0.0001))

    ret_df = ret_df.apply(lambda x: make_valid(x))
    country_shapes = country_shapes.apply(lambda x: make_valid(x))

    country_shapes_with_buffer = country_shapes.buffer(distance)
    ret_df_new = ret_df.difference(country_shapes_with_buffer)

    # repeat to simplify after the buffer correction
    ret_df_new = ret_df_new.map(lambda x: x if x is None else _simplify_polys(
        x, minarea=0.001, tolerance=0.0001))
    ret_df_new = ret_df_new.apply(lambda x: x if x is None else make_valid(x))

    # Drops empty geometry
    ret_df = ret_df_new.dropna()

    return ret_df


def download_WorldPop(country_code,
                      year=2020,
                      update=False,
                      out_logging=False,
                      size_min=300):
    """
    Download tiff file for each country code

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download

    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    WorldPop_filename : str
        Name of the file

    """
    if out_logging:
        _logger.info("Stage 3/4: Download WorldPop datasets")

    WorldPop_filename = f"{two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
    # Urls used to possibly download the file
    WorldPop_urls = [
        f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
    ]
    WorldPop_inputfile = os.path.join(os.getcwd(), "data", "raw", "WorldPop",
                                      WorldPop_filename)  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            _logger.warning(
                f"Stage 4/4: {WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break
        if not loaded:
            _logger.error(
                f"Stage 4/4: Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def convert_GDP(name_file_nc, year=2015, out_logging=False):
    """
    Function to convert the nc database of the GDP to tif, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or toghether with pypsa-africa package.
    """

    if out_logging:
        _logger.info("Stage 4/4: Access to GDP raster data")

    # tif namefile
    name_file_tif = name_file_nc[:-2] + "tif"

    # path of the nc file
    GDP_nc = os.path.join(os.getcwd(), "data", "raw", "GDP",
                          name_file_nc)  # Input filepath nc

    # path of the tif file
    GDP_tif = os.path.join(os.getcwd(), "data", "raw", "GDP",
                           name_file_tif)  # Input filepath nc

    # Check if file exists, otherwise throw exception
    if not os.path.exists(GDP_nc):
        raise Exception(
            f"File {name_file_nc} not found, please download it from https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0 and copy it in {os.path.dirname(GDP_nc)}"
        )

    # open nc dataset
    GDP_dataset = xr.open_dataset(GDP_nc)

    # get the requested year of data or its closest one
    list_years = GDP_dataset["time"]
    if year not in list_years:
        if out_logging:
            _logger.warning(
                f"Stage 3/4 GDP data of year {year} not found, selected the most recent data ({int(list_years[-1])})"
            )
        year = float(list_years[-1])

    # subset of the database and conversion to dataframe
    GDP_dataset = GDP_dataset.sel(time=year).drop("time")
    GDP_dataset.rio.to_raster(GDP_tif)

    return GDP_tif, name_file_tif


def load_GDP(
    countries_codes,
    year=2015,
    update=False,
    out_logging=False,
    name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
):
    """
    Function to load the database of the GDP, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or toghether with pypsa-africa package.
    """

    if out_logging:
        _logger.info("Stage 4/4: Access to GDP raster data")

    # path of the nc file
    name_file_tif = name_file_nc[:-2] + "tif"
    GDP_tif = os.path.join(os.getcwd(), "data", "raw", "GDP",
                           name_file_tif)  # Input filepath tif

    if update | (not os.path.exists(GDP_tif)):
        if out_logging:
            _logger.warning(
                f"Stage 4/4: File {name_file_tif} not found, the file will be produced by processing {name_file_nc}"
            )
        convert_GDP(name_file_nc, year, out_logging)

    return GDP_tif, name_file_tif


def generalized_mask(src, geom, **kwargs):
    "Generalize mask function to account for Polygon and MultiPolygon"
    if geom.geom_type == "Polygon":
        return mask(src, [geom], **kwargs)
    else:
        return mask(src, geom, **kwargs)


def _sum_raster_over_mask(shape, img):
    """
    Function to sum the raster value within a shape
    """
    # select the desired area of the raster corresponding to each polygon
    # Approximation: the population is measured including the pixels
    #   where the border of the shape lays. This leads to slightly overestimate
    #   the output, but the error is limited and it enables halving the
    #   computational time
    out_image, out_transform = generalized_mask(img,
                                                shape,
                                                all_touched=True,
                                                invert=False,
                                                nodata=0.0)
    # calculate total output in the selected geometry
    out_sum = out_image.sum()
    # out_sum = out_image.sum()/2 + out_image_int.sum()/2

    return out_sum


def add_gdp_data(
    df_gadm,
    year=2020,
    update=False,
    out_logging=False,
    name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
    nprocesses=2,
    disable_progressbar=False,
):
    """
    Function to add gdp data to arbitrary number of shapes in a country

    Inputs:
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Outputs:
    --------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["gdp"]
    """
    if out_logging:
        _logger.info("Stage 4/4: Add gdp data to GADM GeoDataFrame")

    # initialize new gdp column
    df_gadm["gdp"] = 0.0

    GDP_tif, name_tif = load_GDP(year, update, out_logging, name_file_nc)

    with rasterio.open(GDP_tif) as src:
        # resample data to target shape
        tqdm_kwargs = dict(
            ascii=False,
            unit=" geometries",
            total=df_gadm.shape[0],
            desc="Compute GDP ",
        )
        for i in tqdm(df_gadm.index, **tqdm_kwargs):
            df_gadm.loc[i,
                        "gdp"] = _sum_raster_over_mask(df_gadm.geometry.loc[i],
                                                       src)
    return df_gadm

    # for index, row in tqdm(df_gadm.iterrows(), index=df_gadm.shape[0]):
    #     # select the desired area of the raster corresponding to each polygon
    #     # Approximation: the gdp is measured excluding the pixels
    #     #   where the border of the shape lays. This may affect the computation
    #     #   but it is conservative and avoids considering multiple times the same
    #     #   pixels
    #     out_image, out_transform = generalized_mask(src,
    #                                                 row["geometry"],
    #                                                 all_touched=True,
    #                                                 invert=False,
    #                                                 nodata=0.0)
    #     # out_image_int, out_transform = mask(src,
    #     #                                row["geometry"],
    #     #                                all_touched=False,
    #     #                                invert=False,
    #     #                                nodata=0.0)

    #     # calculate total gdp in the selected geometry
    #     gdp_by_geom = np.nansum(out_image)
    #     # gdp_by_geom = out_image.sum()/2 + out_image_int.sum()/2

    #     if out_logging == True:
    #         _logger.info("Stage 4/4 GDP: shape: " + str(index) +
    #                      " out of " + str(df_gadm.shape[0]))

    #     # update the gdp data in the dataset
    #     df_gadm.loc[index, "gdp"] = gdp_by_geom


def _init_process_pop(df_gadm_, year_):
    global df_gadm, year
    df_gadm, year = df_gadm_, year_


def _process_func_pop(c_code):

    # get subset by country code
    country_rows = df_gadm.loc[df_gadm["country"] == c_code].copy()

    # get worldpop image
    WorldPop_inputfile, WorldPop_filename = download_WorldPop(
        c_code, year, False, False)

    with rasterio.open(WorldPop_inputfile) as src:

        for i, row in country_rows.iterrows():
            country_rows.loc[i,
                             "pop"] = _sum_raster_over_mask(row.geometry, src)

    return country_rows


def add_population_data(
    df_gadm,
    country_codes,
    year=2020,
    update=False,
    out_logging=False,
    nprocesses=2,
    disable_progressbar=False,
):
    """
    Function to add population data to arbitrary number of shapes in a country

    Inputs:
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Outputs:
    --------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["pop"]
    """

    if out_logging:
        _logger.info("Stage 4/4 POP: Add population data to GADM GeoDataFrame")

    # initialize new population column
    df_gadm["pop"] = 0.0

    tqdm_kwargs = dict(
        ascii=False,
        unit=" countries",
        # total=len(country_codes),
        desc="Compute population ",
    )
    if (nprocesses is None) or (nprocesses == 1):

        with tqdm(total=df_gadm.shape[0], **tqdm_kwargs) as pbar:

            for c_code in country_codes:
                # get subset by country code
                country_rows = df_gadm.loc[df_gadm["country"] == c_code]

                # get worldpop image
                WorldPop_inputfile, WorldPop_filename = download_WorldPop(
                    c_code, year, update, out_logging)

                with rasterio.open(WorldPop_inputfile) as src:

                    for i, row in country_rows.iterrows():
                        df_gadm.loc[i, "pop"] = _sum_raster_over_mask(
                            row.geometry, src)
                        pbar.update(1)

    else:

        kwargs = {
            "initializer": _init_process_pop,
            "initargs": (df_gadm, year),
            "processes": nprocesses,
        }
        with mp.get_context("spawn").Pool(**kwargs) as pool:
            if disable_progressbar:
                _ = list(pool.map(_process_func_pop, country_codes))
                for elem in _:
                    df_gadm.loc[elem.index, "pop"] = elem["pop"]
            else:
                _ = list(
                    tqdm(
                        pool.imap(_process_func_pop, country_codes),
                        total=len(country_codes),
                        **tqdm_kwargs,
                    ))
                for elem in _:
                    df_gadm.loc[elem.index, "pop"] = elem["pop"]


def gadm(countries,
         layer_id=2,
         update=False,
         out_logging=False,
         year=2020,
         nprocesses=None):

    if out_logging:
        _logger.info("Stage 4/4: Creation GADM GeoDataFrame")

    # download data if needed and get the desired layer_id
    df_gadm = get_GADM_layer(countries, layer_id, update)

    # select and rename columns
    df_gadm.rename(columns={"GID_0": "country"}, inplace=True)

    # drop useless columns
    df_gadm = df_gadm[["country", "GADM_ID", "geometry"]]

    # add the population data to the dataset
    add_population_data(df_gadm,
                        countries,
                        year,
                        update,
                        out_logging,
                        nprocesses=nprocesses)

    # add the gdp data to the dataset
    add_gdp_data(
        df_gadm,
        year,
        update,
        out_logging,
        name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
    )

    # set index and simplify polygons
    df_gadm.set_index("GADM_ID", inplace=True)
    df_gadm["geometry"] = df_gadm["geometry"].map(_simplify_polys)

    return df_gadm


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_shapes")
        sets_path_to_root("pypsa-africa")
    configure_logging(snakemake)

    out = snakemake.output

    countries_list = snakemake.config["countries"]
    layer_id = snakemake.config["build_shape_options"]["gadm_layer_id"]
    update = snakemake.config["build_shape_options"]["update_file"]
    out_logging = snakemake.config["build_shape_options"]["out_logging"]
    year = snakemake.config["build_shape_options"]["year"]
    nprocesses = snakemake.config["build_shape_options"]["nprocesses"]
    EEZ_gpkg = snakemake.input["eez"]

    country_shapes = countries(countries_list, update, out_logging)
    save_to_geojson(country_shapes, out.country_shapes)

    offshore_shapes = eez(countries_list, country_shapes, EEZ_gpkg,
                          out_logging)
    save_to_geojson(offshore_shapes, out.offshore_shapes)

    africa_shape = country_cover(country_shapes, offshore_shapes, out_logging)
    save_to_geojson(gpd.GeoSeries(africa_shape), out.africa_shape)

    gadm_shapes = gadm(countries_list,
                       layer_id,
                       update,
                       out_logging,
                       year,
                       nprocesses=nprocesses)
    save_to_geojson(gadm_shapes, out.gadm_shapes)
