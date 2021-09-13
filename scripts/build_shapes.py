# import
import logging
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
#import rioxarray
import rioxarray as rx
from shapely.geometry.base import BaseGeometry
import xarray as xr
from _helpers import _sets_path_to_root
from _helpers import _three_2_two_digits_country
from _helpers import _two_2_three_digits_country
from _helpers import _two_digits_2_name_country
from _helpers import configure_logging
from osm_pbf_power_data_extractor import create_country_list
from rasterio.mask import mask
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon, Point, LineString
from shapely.ops import cascaded_union
#from osm_data_config import AFRICA_CC


_logger = logging.getLogger("build_shapes")
_logger.setLevel(logging.INFO)

# IMPORTANT: RUN SCRIPT FROM THIS SCRIPTS DIRECTORY i.e data_exploration/ TODO: make more robust
#os.chdir(os.path.dirname(os.path.abspath(__file__)))

# import sys

# from ..scripts.iso_country_codes import AFRICA_CC

_sets_path_to_root("pypsa-africa")


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

    GADM_filename = f"gadm36_{_two_2_three_digits_country(country_code)}"
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
                f"Stage 4/4: {GADM_filename} of country {_two_digits_2_name_country(country_code)} does not exist, downloading to {GADM_inputfile_zip}"
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
            f"gadm36_{_two_2_three_digits_country(country_code).upper()}_{code_layer}"
        )

        # read gpkg file
        geodf_temp = gpd.read_file(file_gpkg, layer=layer_name)

        # convert country name representation of the main country (GID_0 column)
        geodf_temp["GID_0"] = [
            _three_2_two_digits_country(twoD_c)
            for twoD_c in geodf_temp["GID_0"]
        ]

        # create a subindex column that is useful
        # in the GADM processing of sub-national zones
        geodf_temp["GADM_ID"] = geodf_temp[f"GID_{code_layer}"]

        # append geodataframes
        geodf_GADM = geodf_GADM.append(geodf_temp)

    geodf_GADM.reset_index(drop=True, inplace=True)

    return geodf_GADM


def _simplify_polys(polys, minarea=0.0001, tolerance=0.008, filterremote=False):
    "Function to simplify the shape polygons"
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys, key=attrgetter("area"), reverse=True) # here deprecation warning: Iteration over multi-part geometries is deprecated and will be removed in Shapely 2.0. Use the `geoms` property to access the constituent parts of a multi-part geometry.
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


def country_cover(country_shapes, eez_shapes=None, out_logging=False):

    if out_logging:
        _logger.info("Stage 3 of 4: Merge country shapes to create continent shape")

    shapes = list(country_shapes)
    if eez_shapes is not None:
        shapes += list(eez_shapes)

    africa_shape = cascaded_union(shapes)
    if isinstance(africa_shape, MultiPolygon):
        africa_shape = max(africa_shape, key=attrgetter("area"))
    return Polygon(shell=africa_shape.exterior)


def save_to_geojson(df, fn): # error occurs here: ERROR:shapely.geos:IllegalArgumentException: Geometry must be a Point or LineString
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
        with os.open(fn, "w") as fp:
            pass


def load_EEZ(countries_codes, name_file="eez_v11.gpkg"):
    """
    Function to load the database of the Exclusive Economic Zones.
    The dataset shall be downloaded independently by the user (see guide) or toghether with pypsa-africa package.
    """

    EEZ_gpkg = os.path.join(os.getcwd(), "data", "raw", "eez",
                            name_file)  # Input filepath gpkg

    if not os.path.exists(EEZ_gpkg):
        raise Exception(
            f"File EEZ {name_file} not found, please download it from https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip and copy it in {os.path.dirname(EEZ_gpkg)}"
        )

    geodf_EEZ = gpd.read_file(EEZ_gpkg)
    geodf_EEZ.dropna(axis=0, how="any", subset=["ISO_TER1"], inplace=True)
    # [["ISO_TER1", "TERRITORY1", "ISO_SOV1", "ISO_SOV2", "ISO_SOV3", "geometry"]]
    geodf_EEZ = geodf_EEZ[["ISO_TER1", "geometry"]]
    selected_countries_codes_3D = [
        _two_2_three_digits_country(x) for x in countries_codes
    ]
    geodf_EEZ = geodf_EEZ[[
        any([x in selected_countries_codes_3D]) for x in geodf_EEZ["ISO_TER1"]
    ]]
    geodf_EEZ["ISO_TER1"] = geodf_EEZ["ISO_TER1"].map(
        lambda x: _three_2_two_digits_country(x))
    geodf_EEZ.reset_index(drop=True, inplace=True)

    geodf_EEZ.rename(columns={"ISO_TER1": "name"}, inplace=True)

    return geodf_EEZ


def eez(countries, country_shapes, update=False, out_logging=False, tol=1e-3):

    if out_logging:
        _logger.info("Stage 2bis of 4: Create offshore shapes - old method")

    # load data
    df_eez = load_EEZ(countries)

    ret_df = df_eez[["name", "geometry"]]
    # create unique shape if country is described by multiple shapes
    for c_code in countries:
        selection = (ret_df.name == c_code)
        n_offshore_shapes = selection.sum()

        if n_offshore_shapes > 1:
            # when multiple shapes per country, then merge polygons
            
            geom = ret_df[selection].geometry.unary_union
            ret_df.drop(ret_df[selection].index, inplace=True)
            ret_df.loc[ret_df.shape[0]] = [c_code, geom]

    ret_df = ret_df.set_index("name")["geometry"].map(
        lambda x: _simplify_polys(x, minarea=0.001, tolerance=0.0001))

    ret_df = gpd.GeoSeries({
        k: v
        for k, v in ret_df.iteritems() if v.distance(country_shapes[k]) < tol
    })
    ret_df.index.name = "name"

    return ret_df


# def eez2(countries, country_shapes, update=False, out_logging=False, tol=1e-3):

#     if out_logging:
#         print("Create offshore shapes")

#     # load data
#     df_eez = load_EEZ(countries)

#     # set index and simplify polygons
#     # ret_df = df_eez.set_index("name").geometry.map(
#     #     lambda x: _simplify_polys(x, filterremote=False))

#     # ret_df = gpd.GeoSeries({
#     #     k: v
#     #     for k, v in ret_df.iteritems() if v.distance(country_shapes[k]) < tol
#     # })
#     # ret_df.index.name = "name"
#     ret_df = df_eez.set_index("name")["geometry"].map(_simplify_polys)

#     return ret_df


def eez_new(countries, country_shapes, out_logging=False, distance=0.01):
    """
    Creates offshore shapes by 
    - buffer smooth countryshape (=offset country shape)
    - and differ that with the offshore shape
    Leads to for instance a 100m non-build coastline

    """
    from shapely.validation import make_valid

    
    if out_logging:
        _logger.info("Stage 2 of 4: Create offshore shapes")

    # load data
    df_eez = load_EEZ(countries)

    # simplified offshore_shape
    # ret_df = df_eez.set_index("name")["geometry"].map(
    #     lambda x: _simplify_polys(x, minarea=0.001, tolerance=0.0001))
    
    ret_df = df_eez[["name", "geometry"]]
    # create unique shape if country is described by multiple shapes
    for c_code in countries:
        selection = (ret_df.name == c_code)
        n_offshore_shapes = selection.sum()

        if n_offshore_shapes > 1:
            # when multiple shapes per country, then merge polygons
            
            geom = ret_df[selection].geometry.unary_union
            ret_df.drop(ret_df[selection].index, inplace=True)
            ret_df.loc[ret_df.shape[0]] = [c_code, geom]

    ret_df = ret_df.set_index("name")["geometry"].map(
        lambda x: _simplify_polys(x, minarea=0.001, tolerance=0.0001))


    ret_df = ret_df.apply(lambda x: make_valid(x))  # hole lies outside occurs here
    country_shapes = country_shapes.apply(lambda x: make_valid(x))

    country_shapes_with_buffer = country_shapes.buffer(distance)
    ret_df_new = ret_df.difference(country_shapes_with_buffer)

    # repeat to simplify after the buffer correction
    ret_df_new = ret_df_new.map(
        lambda x: x if x is None else _simplify_polys(x, minarea=0.001, tolerance=0.0001))
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

    # UN not adjusted
    # WorldPop_filename = f"{_two_2_three_digits_country(country_code).lower()}_ppp_{year}_constrained.tif"
    # WorldPop_url = f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{_two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}"

    WorldPop_filename = f"{_two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
    # Urls used to possibly download the file
    WorldPop_urls = [
        f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{_two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{_two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
    ]

    WorldPop_inputfile = os.path.join(os.getcwd(), "data",
                                      "raw", "WorldPop",
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
            _logger.error(f"Stage 4/4: Impossible to download {WorldPop_filename}")

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
            f"File EEZ {name_file_nc} not found, please download it from https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0 and copy it in {os.path.dirname(GDP_nc)}"
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


def add_gdp_data(
    df_gadm,
    year=2020,
    update=False,
    out_logging=False,
    name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
):
    """Function to add the population info for each country shape in the gadm dataset"""

    if out_logging:
        _logger.info("Stage 4/4: Add population data to GADM GeoDataFrame")

    # initialize new population column
    df_gadm["gdp"] = 0.0

    GDP_tif, name_tif = load_GDP(year, update, out_logging, name_file_nc)

    with rasterio.open(GDP_tif) as src:
        # data_GDP = src.read(1)
        # resample data to target shape

        for index, row in df_gadm.iterrows():
            # select the desired area of the raster corresponding to each polygon
            # Approximation: the gdp is measured excluding the pixels
            #   where the border of the shape lays. This may affect the computation
            #   but it is conservative and avoids considering multiple times the same
            #   pixels
            if row["geometry"].geom_type == 'Polygon':
                out_image, out_transform = mask(src,
                                            [row["geometry"]],
                                            all_touched=True,
                                            invert=False,
                                            nodata=0.0)
            else:
                out_image, out_transform = mask(src,
                                            row["geometry"],
                                            all_touched=True,
                                            invert=False,
                                            nodata=0.0)
            # out_image_int, out_transform = mask(src,
            #                                row["geometry"],
            #                                all_touched=False,
            #                                invert=False,
            #                                nodata=0.0)

            # calculate total gdp in the selected geometry
            gdp_by_geom = np.nansum(out_image)
            # gdp_by_geom = out_image.sum()/2 + out_image_int.sum()/2

            if out_logging == True:
                _logger.info("Stage 4/4 GDP: shape: " + str(index) + " out of " + str(df_gadm.shape[0]))

            # update the gdp data in the dataset
            df_gadm.loc[index, "gdp"] = gdp_by_geom
    return df_gadm



def add_population_data(df_gadm,
                        country_codes,
                        year=2020,
                        update=False,
                        out_logging=False):
    """Function to add the population info for each country shape in the gadm dataset"""

    if out_logging:
        _logger.info("Stage 4/4 POP: Add population data to GADM GeoDataFrame")

    # initialize new population column
    df_gadm["pop"] = 0.0

    # count elements
    count = 0

    for c_code in country_codes:
        WorldPop_inputfile, WorldPop_filename = download_WorldPop(
            c_code, year, update, out_logging)

        with rasterio.open(WorldPop_inputfile) as src:
            country_rows = df_gadm.loc[df_gadm["country"] == c_code]

            for index, row in country_rows.iterrows():
                # select the desired area of the raster corresponding to each polygon
                # Approximation: the population is measured including the pixels
                #   where the border of the shape lays. This leads to slightly overestimate
                #   the population, but the error is limited and it enables halving the
                #   computational time
                if row["geometry"].geom_type == 'Polygon':
                    out_image, out_transform = mask(src,
                                                [row["geometry"]],
                                                all_touched=True,
                                                invert=False,
                                                nodata=0.0)
                else:
                    out_image, out_transform = mask(src,
                                                row["geometry"],
                                                all_touched=True,
                                                invert=False,
                                                nodata=0.0)
                # out_image_int, out_transform = mask(src,
                #                                row["geometry"],
                #                                all_touched=False,
                #                                invert=False,
                #                                nodata=0.0)

                # calculate total population in the selected geometry
                pop_by_geom = out_image.sum()
                # pop_by_geom = out_image.sum()/2 + out_image_int.sum()/2

                # update the population data in the dataset
                df_gadm.loc[index, "pop"] = pop_by_geom

                count += 1

                if out_logging == True:
                    _logger.info("Stage 4/4 POP: " + str(count) + " out of " + str(df_gadm.shape[0]) + " [" + c_code + "]")
                    # print(c_code, ": ", index, " out of ",
                    #      country_rows.shape[0])


def gadm(countries, layer_id=2, update=False, out_logging=False, year=2020):

    if out_logging:
        _logger.info("Stage 4/4: Creation GADM GeoDataFrame")

    # download data if needed and get the desired layer_id
    df_gadm = get_GADM_layer(countries, layer_id, update)

    # select and rename columns
    df_gadm.rename(columns={"GID_0": "country"}, inplace=True)

    # drop useless columns
    df_gadm = df_gadm[["country", "GADM_ID", "geometry"]]

    # # add the population data to the dataset
    # add_population_data(df_gadm, countries, year, update, out_logging)

    # # add the gdp data to the dataset
    # add_gdp_data(
    #     df_gadm,
    #     year,
    #     update,
    #     out_logging,
    #     name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
    # )

    # set index and simplify polygons
    df_gadm.set_index("GADM_ID", inplace=True)
    df_gadm["geometry"] = df_gadm["geometry"].map(_simplify_polys)

    return df_gadm


if __name__ == "__main__":
    # print(os.path.dirname(os.path.abspath(__file__)))
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake("build_shapes")
    configure_logging(snakemake)

    out = snakemake.output

    # Parameters to be later initialized through snakemake
    # countries_list = list(AFRICA_CC) # TODO: implement some coding to automatically process all africa by config.yaml
    countries_list = create_country_list(snakemake.config["countries"])
    layer_id = snakemake.config['build_shape_options']['gadm_layer_id']
    update = snakemake.config['build_shape_options']['update_file']
    out_logging = snakemake.config['build_shape_options']['out_logging']
    year = snakemake.config['build_shape_options']['year']

    # print(snakemake.config)

    country_shapes = countries(countries_list, update, out_logging)
    save_to_geojson(country_shapes, out.country_shapes)

    offshore_shapes = eez_new(countries_list, country_shapes, out_logging)
    save_to_geojson(offshore_shapes, out.offshore_shapes)

    offshore_shapes_old = eez(countries_list, country_shapes, update, out_logging)
    save_to_geojson(offshore_shapes_old, out.offshore_shapes_old)

    africa_shape = country_cover(country_shapes, offshore_shapes, out_logging)
    save_to_geojson(gpd.GeoSeries(africa_shape), out.africa_shape)

    gadm_shapes = gadm(countries_list, layer_id, update, out_logging, year)
    save_to_geojson(gadm_shapes, out.gadm_shapes)
