# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import multiprocessing as mp
import os
import shutil
from itertools import takewhile
from operator import attrgetter

import fiona
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import requests
import shapely
import xarray as xr
from _helpers import (
    BASE_DIR,
    configure_logging,
    create_logger,
    save_to_geojson,
    three_2_two_digits_country,
    two_2_three_digits_country,
    two_digits_2_name_country,
)
from numba import njit
from numba.core import types
from numba.typed import Dict
from pyproj import Transformer
from rasterio.mask import mask
from rasterio.windows import Window
from shapely.geometry import MultiPolygon
from shapely.ops import transform, unary_union
from shapely.validation import make_valid
from tqdm import tqdm

logger = create_logger(__name__)


def get_GADM_filename(country_code: str) -> str:
    """
    Function to get the GADM filename given the country code.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files

    Returns
    -------
    str
        GADM filename corresponding to the country code
    """
    special_codes_GADM = {
        "XK": "XKO",  # kosovo
        "CP": "XCL",  # clipperton island
        "SX": "MAF",  # sint maartin
        "TF": "ATF",  # french southern territories
        "AX": "ALA",  # aland
        "IO": "IOT",  # british indian ocean territory
        "CC": "CCK",  # cocos island
        "NF": "NFK",  # norfolk
        "PN": "PCN",  # pitcairn islands
        "JE": "JEY",  # jersey
        "XS": "XSP",  # spratly
        "GG": "GGY",  # guernsey
        "UM": "UMI",  # united states minor outlying islands
        "SJ": "SJM",  # svalbard
        "CX": "CXR",  # Christmas island
    }

    if country_code in special_codes_GADM:
        return f"gadm41_{special_codes_GADM[country_code]}"
    else:
        return f"gadm41_{two_2_three_digits_country(country_code)}"


def download_GADM(
    country_code: str, update: bool = False, out_logging: bool = False
) -> tuple[str, str]:
    """
    Download gpkg file from GADM for a given country code.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files
    update : bool
        Update = true, forces re-download of files

    Returns
    -------
    GADM_inputfile_gpkg : str
        Path of the downloaded gpkg file
    GADM_filename : str
        Name of the gpkg file per country
    """
    GADM_filename = get_GADM_filename(country_code)
    GADM_url = f"https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/{GADM_filename}.gpkg"

    GADM_inputfile_gpkg = os.path.join(
        BASE_DIR,
        "data",
        "gadm",
        GADM_filename,
        GADM_filename + ".gpkg",
    )  # Input filepath gpkg

    if not os.path.exists(GADM_inputfile_gpkg) or update is True:
        if out_logging:
            logger.warning(
                f"Stage 5 of 5: {GADM_filename} of country {two_digits_2_name_country(country_code)} does not exist, downloading to {GADM_inputfile_gpkg}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_gpkg), exist_ok=True)

        try:
            r = requests.get(GADM_url, stream=True, timeout=300)
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout):
            raise Exception(
                f"GADM server is down at {GADM_url}. Data needed for building shapes can't be extracted.\n\r"
            )
        except Exception as exception:
            raise Exception(
                f"An error happened when trying to load GADM data by {GADM_url}.\n\r"
                + str(exception)
                + "\n\r"
            )
        else:
            with open(GADM_inputfile_gpkg, "wb") as f:
                shutil.copyfileobj(r.raw, f)

    return GADM_inputfile_gpkg, GADM_filename


def filter_gadm(
    geodf: gpd.GeoDataFrame,
    layer: int,
    cc: str,
    contended_flag: str,
    output_nonstd_to_csv: bool = False,
) -> gpd.GeoDataFrame:
    """
    Function to filter the GADM geodataframe according to the contented_flag option.

    Parameters
    ----------
    geodf : gpd.GeoDataFrame
        Geodataframe to filter
    layer : int
        Layer of the geodataframe to filter
    cc : str
        Country code to filter
    contended_flag : str
        Option to treat contended areas, i.e. areas that are not assigned to the country in the GADM layer
        but are part of the country according to the country code (GID_0) of the geodataframe. The options are:
        - "drop": drop contended areas from the geodataframe
        - "set_by_country": set GID_0 of contended areas to the country code of the country (default)
    output_nonstd_to_csv : bool
        If True, outputs the non-standard rows to a csv file for debugging purposes (default False)

    Returns
    -------
    gpd.GeoDataFrame
        Filtered GADM geodataframe
    """
    # identify non standard geodf rows
    geodf_non_std = geodf[geodf["GID_0"] != two_2_three_digits_country(cc)].copy()

    if not geodf_non_std.empty:
        logger.info(
            f"Contended areas have been found for gadm layer {layer}. They will be treated according to {contended_flag} option"
        )

        # NOTE: in these options GID_0 is not changed because it is modified below
        if contended_flag == "drop":
            geodf.drop(geodf_non_std.index, inplace=True)
        elif contended_flag != "set_by_country":
            # "set_by_country" option is the default; if this elif applies, the desired option falls back to the default
            logger.warning(
                f"Value '{contended_flag}' for option contented_flag is not recognized.\n"
                + "Fallback to 'set_by_country'"
            )

    # force GID_0 to be the country code for the relevant countries
    geodf["GID_0"] = cc

    # country shape should have a single geometry
    if (layer == 0) and (geodf.shape[0] > 1):
        logger.warning(
            f"Country shape is composed by multiple shapes that are being merged in agreement to contented_flag option '{contended_flag}'"
        )
        # take the first row only to re-define geometry keeping other columns
        geodf = geodf.iloc[[0]].set_geometry([geodf.union_all()])

    # debug output to file
    if output_nonstd_to_csv and not geodf_non_std.empty:
        geodf_non_std.to_csv(
            f"resources/non_standard_gadm{layer}_{cc}_raw.csv", index=False
        )

    return geodf


def get_GADM_layer(
    country_list: list[str],
    layer_id: int,
    geo_crs: str = "EPSG:4326",
    contended_flag: str = "set_by_country",
    update: bool = False,
    outlogging: bool = False,
) -> gpd.GeoDataFrame:
    """
    Function to retrieve a specific layer id of a geopackage for a selection of
    countries.

    Parameters
    ----------
    country_list : list[str]
        List of the countries
    layer_id : int
        Layer to consider in the format GID_{layer_id}.
        When the requested layer_id is greater than the last available layer, then the last layer is selected.
        When a negative value is requested, then, the last layer is requested
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g. "EPSG:4326")
    contended_flag : str
        Option to treat contended areas, i.e. areas that are not assigned to the country in the GADM layer
        but are part of the country according to the country code (GID_0) of the geodataframe. The options are:
        - "drop": drop contended areas from the geodataframe
        - "set_by_country": set GID_0 of contended areas to the country code of the country (default)
    update : bool
        Update = true, forces re-download of files
    outlogging : bool
        If True, emits progress information via the module logger.

    Returns
    -------
    gpd.GeoDataFrame
        Geodataframe with the requested GADM layer for the selected countries
    """
    # initialization of the geoDataFrame
    geodf_list = []

    for country_code in country_list:
        # Set the current layer id (cur_layer_id) to global layer_id
        cur_layer_id = layer_id

        # download file gpkg
        file_gpkg, name_file = download_GADM(country_code, update, outlogging)

        # get layers of a geopackage
        list_layers = fiona.listlayers(file_gpkg)

        # get layer name
        if (cur_layer_id < 0) or (cur_layer_id >= len(list_layers)):
            # when layer id is negative or larger than the number of layers, select the last layer
            cur_layer_id = len(list_layers) - 1

        # read gpkg file
        geodf_temp = gpd.read_file(
            file_gpkg, layer="ADM_ADM_" + str(cur_layer_id), engine="pyogrio"
        ).to_crs(geo_crs)

        geodf_temp = filter_gadm(
            geodf=geodf_temp,
            layer=cur_layer_id,
            cc=country_code,
            contended_flag=contended_flag,
            output_nonstd_to_csv=False,
        )

        # create a subindex column that is useful
        # in the GADM processing of sub-national zones
        geodf_temp["GADM_ID"] = geodf_temp[f"GID_{cur_layer_id}"]

        # from pypsa-earth-sec
        # if layer_id == 0:
        #     geodf_temp["GADM_ID"] = geodf_temp[f"GID_{cur_layer_id}"].apply(
        #         lambda x: two_2_three_digits_country(x[:2])
        #     ) + pd.Series(range(1, geodf_temp.shape[0] + 1)).astype(str)
        # else:
        #     # create a subindex column that is useful
        #     # in the GADM processing of sub-national zones
        #     # Fix issues with missing "." in selected cases
        #     geodf_temp["GADM_ID"] = geodf_temp[f"GID_{cur_layer_id}"].apply(
        #         lambda x: x if x[3] == "." else x[:3] + "." + x[3:]
        #     )

        # append geodataframes
        geodf_list.append(geodf_temp)

    geodf_GADM = gpd.GeoDataFrame(pd.concat(geodf_list, ignore_index=True))
    geodf_GADM.set_crs(geo_crs)

    return geodf_GADM


def _simplify_polys(
    polys: gpd.GeoDataFrame,
    minarea: float = 0.01,
    tolerance: float = 0.01,
    filterremote: bool = False,
) -> gpd.GeoDataFrame:
    """
    Function to simplify the shape polygons

    Parameters
    ----------
    polys : gpd.GeoDataFrame
        Geodataframe with the polygons to simplify
    minarea : float
        Minimum area of the polygons to keep (default 0.01)
    tolerance : float
        Tolerance for the simplification (default 0.01)
    filterremote : bool
        If True, filters out polygons that are remote from the main polygon (default False).
        Remote is defined as having a distance from the main polygon larger than the main polygon length.

    Returns
    -------
    gpd.GeoDataFrame
        Geodataframe with the simplified polygons
    """
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys.geoms, key=attrgetter("area"), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area / (2.0 * np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon(
                [
                    p
                    for p in takewhile(lambda p: p.area > minarea, polys)
                    if not filterremote or (mainpoly.distance(p) < mainlength)
                ]
            )
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)


def countries(
    countries: list[str],
    geo_crs: str,
    contended_flag: str,
    update: bool = False,
    out_logging: bool = False,
    tolerance: float = 0.01,
    minarea: float = 0.01,
) -> gpd.GeoSeries:
    """
    Create country shapes

    Parameters
    ----------
    countries : list[str]
        List of the countries
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g. "EPSG:4326")
    contended_flag : str
        Flag indicating whether to include contended areas
    update : bool
        Update = true, forces re-download of files
    out_logging : bool
        If True, emits progress information via the module logger.
    tolerance : float
        Tolerance for the simplification (default 0.01)
    minarea : float
        Minimum area of the polygons to keep (default 0.01)

    Returns
    -------
    gpd.GeoSeries
        Geoseries with the country shapes
    """

    if out_logging:
        logger.info("Stage 1 of 5: Create country shapes")

    # download data if needed and get the layer id 0, corresponding to the countries
    df_countries = get_GADM_layer(
        countries,
        0,
        geo_crs,
        contended_flag,
        update,
        out_logging,
    )

    # select and rename columns
    df_countries = df_countries[["GID_0", "geometry"]].copy()
    df_countries.rename(columns={"GID_0": "name"}, inplace=True)

    # set index and simplify polygons
    ret_df = df_countries.set_index("name")["geometry"].map(
        lambda x: _simplify_polys(x, tolerance=tolerance, minarea=minarea)
    )

    # there may be "holes" in the countries geometry which cause troubles along the workflow
    # e.g. that is the case for enclaves like Dahagram–Angarpota for IN/BD
    ret_df = ret_df.make_valid()

    return ret_df


def country_cover(
    country_shapes: gpd.GeoSeries,
    eez_shapes: gpd.GeoSeries = None,
    out_logging: bool = False,
    distance: float = 0.02,
) -> gpd.GeoSeries:
    """
    Create a continent shape by merging the country shapes and the EEZ shapes (if provided) with a buffer.

    Parameters
    ----------
    country_shapes : gpd.GeoSeries
        Geoseries with the country shapes
    eez_shapes : gpd.GeoSeries, optional
        Geoseries with the EEZ shapes
    out_logging : bool
        If True, emits progress information via the module logger.
    distance : float
        Distance for the buffer (default 0.02)

    Returns
    -------
    gpd.GeoSeries
        Geoseries with the continent shape
    """
    if out_logging:
        logger.info("Stage 3 of 5: Merge country shapes to create continent shape")

    shapes = country_shapes.apply(lambda x: x.buffer(distance))
    shapes_list = list(shapes)
    if eez_shapes is not None:
        shapes_list += list(eez_shapes)

    africa_shape = make_valid(unary_union(shapes_list))

    return africa_shape


def load_EEZ(
    countries_codes: list[str], geo_crs: str, EEZ_gpkg: str = "./data/eez/eez_v11.gpkg"
) -> gpd.GeoDataFrame:
    """
    Function to load the database of the Exclusive Economic Zones.

    The dataset shall be downloaded independently by the user (see
    guide) or together with pypsa-earth package.

    Parameters
    ----------
    countries_codes : list[str]
        List of two-letter ISO country codes.
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g. "EPSG:4326").
    EEZ_gpkg : str, default "./data/eez/eez_v11.gpkg"
        Path to the Marine Regions *World EEZ v11* geopackage.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with the EEZ geometries.
    """
    if not os.path.exists(EEZ_gpkg):
        raise Exception(
            f"File EEZ {EEZ_gpkg} not found, please download it from https://www.marineregions.org/download_file.php?name=World_EEZ_v11_20191118_gpkg.zip and copy it in {os.path.dirname(EEZ_gpkg)}"
        )

    geodf_EEZ = gpd.read_file(EEZ_gpkg, engine="pyogrio").to_crs(geo_crs)
    geodf_EEZ.dropna(axis=0, how="any", subset=["ISO_TER1"], inplace=True)
    # [["ISO_TER1", "TERRITORY1", "ISO_SOV1", "ISO_SOV2", "ISO_SOV3", "geometry"]]
    geodf_EEZ = geodf_EEZ[["ISO_TER1", "geometry"]]
    selected_countries_codes_3D = [
        two_2_three_digits_country(x) for x in countries_codes
    ]
    geodf_EEZ = geodf_EEZ[
        [any([x in selected_countries_codes_3D]) for x in geodf_EEZ["ISO_TER1"]]
    ]
    geodf_EEZ["ISO_TER1"] = geodf_EEZ["ISO_TER1"].map(
        lambda x: three_2_two_digits_country(x)
    )
    geodf_EEZ.reset_index(drop=True, inplace=True)

    geodf_EEZ.rename(columns={"ISO_TER1": "name"}, inplace=True)

    return geodf_EEZ


def eez(
    countries: list[str],
    geo_crs: str,
    country_shapes: gpd.GeoSeries | gpd.GeoDataFrame,
    EEZ_gpkg: str,
    out_logging: bool = False,
    distance: float = 0.0,
    minarea: float = 0.01,
    tolerance: float = 0.01,
    simplify_gadm: bool = True,
) -> gpd.GeoDataFrame:
    """
    Build offshore (EEZ) shapes for the requested countries.

    The function loads EEZ geometries, unions them per country, and optionally
    simplifies them. If ``distance`` is non-zero, the onshore country
    shapes are buffered and subtracted from EEZ geometries to create an
    "offshore-only" region (e.g. to enforce a non-build coastal strip).

    Parameters
    ----------
    countries : list[str]
        Two-letter ISO country codes to process (e.g. ``["DE", "FR"]``).
    geo_crs : str
        CRS used for geographic projection, passed to GeoPandas (e.g."EPSG:4326").
        Note: buffering distances depend on the CRS units.
    country_shapes : geopandas.GeoSeries or geopandas.GeoDataFrame
        Country geometries indexed by the same two-letter ISO codes. Must be in ``geo_crs``.
    EEZ_gpkg : str
        Path to the Marine Regions *World EEZ v11* geopackage.
    out_logging : bool, default False
        If True, emits progress information via the module logger.
    distance : float, default 0.0
        Buffer distance applied to ``country_shapes`` before subtraction from
        the EEZ geometries. Units are CRS-dependent.
    minarea : float, default 0.01
        Minimum polygon area threshold used by the simplification routine.
    tolerance : float, default 0.01
        Tolerance passed to Shapely ``simplify`` via ``_simplify_polys``.
    simplify_gadm : bool, default True
        If True, simplify and validate geometries before and after the coastal
        subtraction.

    Returns
    -------
    geopandas.GeoDataFrame
        Offshore EEZ geometries indexed by country code (column ``name`` as
        index, geometry column ``geometry``). Empty/invalid geometries are
        removed.
    """
    if out_logging:
        logger.info("Stage 2 of 5: Create offshore shapes")

    # load data
    df_eez = load_EEZ(countries, geo_crs, EEZ_gpkg)

    eez_countries = [cc for cc in countries if df_eez.name.str.contains(cc).any()]
    ret_df = gpd.GeoDataFrame(
        {
            "name": eez_countries,
            "geometry": [
                df_eez.geometry.loc[df_eez.name == cc].geometry.union_all()
                for cc in eez_countries
            ],
        }
    ).set_index("name")

    if simplify_gadm:
        ret_df = ret_df.geometry.map(
            lambda x: _simplify_polys(x, minarea=minarea, tolerance=tolerance)
        )

        ret_df = ret_df.apply(lambda x: make_valid(x))

    country_shapes_with_buffer = country_shapes.buffer(distance)
    ret_df_new = ret_df.difference(country_shapes_with_buffer)

    if simplify_gadm:
        # repeat to simplify after the buffer correction
        ret_df_new = ret_df_new.map(
            lambda x: (
                x
                if x is None
                else _simplify_polys(x, minarea=minarea, tolerance=tolerance)
            )
        )
        ret_df_new = ret_df_new.apply(lambda x: x if x is None else make_valid(x))

    # Drops empty geometry
    ret_df = ret_df_new.dropna()
    ret_df = ret_df[ret_df.geometry.is_valid & ~ret_df.geometry.is_empty]

    return ret_df


def download_WorldPop(
    country_code: str,
    worldpop_method: str,
    year: int = 2020,
    update: bool = False,
    out_logging: bool = False,
    size_min: int = 300,
) -> tuple[str, str]:
    """
    Download Worldpop using either the standard method or the API method.

    Parameters
    ----------
    worldpop_method: str
         worldpop_method = "api" will use the API method to access the WorldPop 100mx100m dataset.  worldpop_method = "standard" will use the standard method to access the WorldPop 1KMx1KM dataset.
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
    if worldpop_method == "api":
        return download_WorldPop_API(country_code, year, update, out_logging, size_min)

    elif worldpop_method == "standard":
        return download_WorldPop_standard(
            country_code, year, update, out_logging, size_min
        )


def download_WorldPop_standard(
    country_code: str,
    year: int = 2020,
    update: bool = False,
    out_logging: bool = False,
    size_min: int = 300,
) -> tuple[str, str]:
    """
    Download tiff file for each country code using the standard method from
    worldpop datastore with 1kmx1km resolution.

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

    if country_code == "XK":
        WorldPop_filename = f"kos_ppp_{year}_constrained.tif"
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/KOS/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/KOS/{WorldPop_filename}",
        ]
    else:
        WorldPop_filename = f"{two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
        # Urls used to possibly download the file
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        ]

    WorldPop_inputfile = os.path.join(
        BASE_DIR, "data", "WorldPop", WorldPop_filename
    )  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            logger.warning(
                f"Stage 3 of 5: {WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
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
            logger.error(f"Stage 3 of 5: Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def download_WorldPop_API(
    country_code: str,
    year: int = 2020,
    update: bool = False,
    out_logging: bool = False,
    size_min: int = 300,
) -> tuple[str, str]:
    """
    Download tiff file for each country code using the api method from worldpop
    API with 100mx100m resolution.

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

    WorldPop_filename = f"{two_2_three_digits_country(country_code).lower()}_ppp_{year}_UNadj_constrained.tif"
    # Request to get the file
    WorldPop_inputfile = os.path.join(
        BASE_DIR, "data", "WorldPop", WorldPop_filename
    )  # Input filepath tif
    os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)
    year_api = int(str(year)[2:])
    loaded = False
    WorldPop_api_urls = [
        f"https://www.worldpop.org/rest/data/pop/wpgp?iso3={two_2_three_digits_country(country_code)}",
    ]
    for WorldPop_api_url in WorldPop_api_urls:
        with requests.get(WorldPop_api_url, stream=True) as r:
            WorldPop_tif_url = r.json()["data"][year_api]["files"][0]

        with requests.get(WorldPop_tif_url, stream=True) as r:
            with open(WorldPop_inputfile, "wb") as f:
                if float(r.headers["Content-length"]) > size_min:
                    shutil.copyfileobj(r.raw, f)
                    loaded = True
                    break
    if not loaded:
        logger.error(f"Stage 3 of 5: Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def convert_GDP(
    name_file_nc: str, year: int = 2015, out_logging: bool = False
) -> tuple[str, str]:
    """
    Function to convert the nc database of the GDP to tif, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or together with pypsa-earth package.

    Parameters
    ----------
    name_file_nc : str
        Name of the nc file containing the GDP data (e.g. "GDP_PPP_1990_2015_5arcmin_v2.nc")
    year : int
        Year of the data to convert
    out_logging : bool
        If True, emits progress information via the module logger.

    Returns
    -------
    GDP_tif : str
        Path of the converted tif file
    name_file_tif : str
        Name of the converted tif file
    """

    if out_logging:
        logger.info("Stage 5 of 5: Access to GDP raster data")

    # tif namefile
    name_file_tif = name_file_nc[:-2] + "tif"

    # path of the nc file
    GDP_nc = os.path.join(BASE_DIR, "data", "GDP", name_file_nc)  # Input filepath nc

    # path of the tif file
    GDP_tif = os.path.join(BASE_DIR, "data", "GDP", name_file_tif)  # Input filepath nc

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
            logger.warning(
                f"Stage 5 of 5 GDP data of year {year} not found, selected the most recent data ({int(list_years[-1])})"
            )
        year = float(list_years[-1])

    # subset of the database and conversion to dataframe
    GDP_dataset = GDP_dataset.sel(time=year).drop("time")
    GDP_dataset.rio.to_raster(GDP_tif)

    return GDP_tif, name_file_tif


def load_GDP(
    year: int = 2015,
    update: bool = False,
    out_logging: bool = False,
    name_file_nc: str = "GDP_PPP_1990_2015_5arcmin_v2.nc",
) -> tuple[str, str]:
    """
    Function to load the database of the GDP, based on the work at https://doi.org/10.1038/sdata.2018.4.
    The dataset shall be downloaded independently by the user (see guide) or together with pypsa-earth package.

    Parameters
    ----------
    year : int
        Year of the data to load
    update : bool
        Update = true, forces re-download of files
    out_logging : bool
        If True, emits progress information via the module logger.
    name_file_nc : str
        Name of the nc file containing the GDP data (e.g. "GDP_PPP_1990_2015_5arcmin_v2.nc")

    Returns
    -------
    GDP_tif : str
        Path of the converted tif file
    name_file_tif : str
        Name of the converted tif file
    """

    if out_logging:
        logger.info("Stage 5 of 5: Access to GDP raster data")

    # path of the nc file
    name_file_tif = name_file_nc[:-2] + "tif"
    GDP_tif = os.path.join(BASE_DIR, "data", "GDP", name_file_tif)  # Input filepath tif

    if update | (not os.path.exists(GDP_tif)):
        if out_logging:
            logger.warning(
                f"Stage 5 of 5: File {name_file_tif} not found, the file will be produced by processing {name_file_nc}"
            )
        convert_GDP(name_file_nc, year, out_logging)

    return GDP_tif, name_file_tif


def generalized_mask(
    src: rasterio.io.DatasetReader, geom: shapely.geometry.base.BaseGeometry, **kwargs
) -> rasterio.io.DatasetReader:
    """
    Generalize mask function to account for Polygon and MultiPolygon

    Parameters
    ----------
    src : rasterio.io.DatasetReader
        Rasterio dataset reader object
    geom : shapely.geometry.base.BaseGeometry
        Geometry to mask the raster with. Can be a Polygon or a MultiPolygon.
    **kwargs : dict
        Additional keyword arguments to pass to rasterio.mask.mask function.

    Returns
    -------
    rasterio.io.DatasetReader
        Masked rasterio dataset reader object
    """
    if geom.geom_type == "Polygon":
        return mask(src, [geom], **kwargs)
    elif geom.geom_type == "MultiPolygon":
        return mask(src, geom.geoms, **kwargs)
    else:
        return mask(src, geom, **kwargs)


def _sum_raster_over_mask(
    shape: shapely.geometry.base.BaseGeometry, img: rasterio.io.DatasetReader
) -> float:
    """
    Function to sum the raster value within a shape.

    The raster is masked using the provided geometry, and all pixels touched
    by the geometry boundary are included (`all_touched=True`). This approach
    slightly overestimates the result, as boundary pixels are fully counted,
    but it significantly reduces computational cost.

    Parameters
    ----------
    shape : shapely.geometry.base.BaseGeometry
        Geometry to mask the raster with. Can be a Polygon or a MultiPolygon.
    img : rasterio.io.DatasetReader
        Rasterio dataset reader object

    Returns
    -------
    float
        Sum of the raster values within the shape
    """
    out_image, out_transform = generalized_mask(
        img, shape, all_touched=True, invert=False, nodata=0.0
    )
    # calculate total output in the selected geometry
    out_image[np.isnan(out_image)] = 0
    out_sum = out_image.sum()
    # out_sum = out_image.sum()/2 + out_image_int.sum()/2

    return out_sum


def add_gdp_data(
    df_gadm: gpd.GeoDataFrame,
    year: int = 2020,
    update: bool = False,
    out_logging: bool = False,
    name_file_nc: str = "GDP_PPP_1990_2015_5arcmin_v2.nc",
    nprocesses: int = 2,
    disable_progressbar: bool = False,
) -> gpd.GeoDataFrame:
    """
    Function to add gdp data to arbitrary number of shapes in a country.

    Parameters
    ----------
    df_gadm: Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]

    Returns
    -------
    df_gadm: Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["gdp"]
    """
    if out_logging:
        logger.info("Stage 5 of 5: Add gdp data to GADM GeoDataFrame")

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
            df_gadm.loc[i, "gdp"] = _sum_raster_over_mask(df_gadm.geometry.loc[i], src)
    return df_gadm


def _init_process_pop(
    df_gadm_: gpd.GeoDataFrame,
    df_tasks_: pd.DataFrame,
    dict_worldpop_file_locations_: dict,
) -> None:
    """
    Initialize the global variables for the population processing function.
    This function is called by the multiprocessing Pool to set up the global variables
    for the population processing function.

    Parameters
    ----------
    df_gadm_: gpd.GeoDataFrame
        GeoDataFrame containing the administrative boundaries
    df_tasks_: pd.DataFrame
        DataFrame containing the tasks for population processing
    dict_worldpop_file_locations_: dict
        Dictionary mapping country codes to worldpop file locations
    """
    global df_gadm, df_tasks
    df_gadm, df_tasks = df_gadm_, df_tasks_

    global dict_worldpop_file_locations
    dict_worldpop_file_locations = dict_worldpop_file_locations_


def process_function_population(row_id: int) -> pd.DataFrame:
    """
    Function that reads the task from df_tasks and executes all the methods.

    to obtain population values for the specified region

    Parameters
    ----------
    row_id: int
        integer which indicates a specific row of df_tasks

    Returns
    -------
    windowed_pop_count: pd.DataFrame
        Dataframe containing "GADM_ID" and "pop" columns
        It represents the amount of population per region (GADM_ID),
        for the settings given by the row in df_tasks
    """
    # Get the current task values
    current_row = df_tasks.iloc[row_id]

    c_code = current_row["c_code"]
    window_dimensions = current_row["window_dimensions"]
    transform = current_row["affine_transform"]
    latlong_topleft = current_row["latlong_coordinate_topleft"]
    latlong_botright = current_row["latlong_coordinate_botright"]

    # Obtain the inputfile location from dict
    WorldPop_inputfile = dict_worldpop_file_locations[c_code]
    # get subset by country code
    country_rows = df_gadm.loc[df_gadm["country"] == c_code]

    np_pop_val, np_pop_xy = get_worldpop_val_xy(WorldPop_inputfile, window_dimensions)

    # If no values are present in the current window skip the remaining steps
    if len(np_pop_val) == 0:
        return []

    # get the geomask with id mappings
    region_geomask, id_mapping = compute_geomask_region(
        country_rows, transform, window_dimensions, latlong_topleft, latlong_botright
    )

    # If no values are present in the id_mapping skip the remaining steps
    if len(id_mapping) == 0:
        return []

    # Calculate the population for each region
    windowed_pop_count = sum_values_using_geomask(
        np_pop_val, np_pop_xy, region_geomask, id_mapping
    )
    return windowed_pop_count


def get_worldpop_val_xy(
    WorldPop_inputfile: str, window_dimensions: tuple
) -> tuple[np.ndarray, np.ndarray]:
    """
    Function to extract data from .tif input file.

    Parameters
    ----------
    WorldPop_inputfile: str
        file location of worldpop file
    window_dimensions: tuple
        dimensions of window used when reading file

    Returns
    -------
    np_pop_valid: np.ndarray
        array filled with values for each nonzero pixel in the worldpop file
    np_pop_xy: np.ndarray
        array with [x,y] coordinates of the corresponding nonzero values in np_pop_valid
    """
    col_offset, row_offset, width, height = window_dimensions

    current_window = Window(col_offset, row_offset, width, height)

    # Open the file using rasterio
    with rasterio.open(WorldPop_inputfile) as src:
        # --- Process the pixels in the image for population data ---

        # Read the gray layer (1) to get an np.array of this band
        # Rasterio doesn't support lower than float32 readout
        # Hence np_pop_raster will have nbytes = 4 * width * height
        np_pop_raster = src.read(1, window=current_window)

        # Set 'nodata' values to 0
        np_pop_raster[np_pop_raster == src.nodata] = 0

        # Set np_pop_xy to pixel locations of non zero values
        np_pop_xy = np_pop_raster.nonzero()

        # Transform to get [ x, y ] array
        # np_pop_xy as 'I' (uintc), see
        # https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.uintc
        np_pop_xy = np.array([np_pop_xy[0], np_pop_xy[1]]).T.astype("I")

        # Extract the values from the locations of non zero pixels
        np_pop_valid = np_pop_raster[np_pop_xy.T[0], np_pop_xy.T[1]]

    return np_pop_valid, np_pop_xy


def compute_geomask_region(
    country_rows: gpd.GeoDataFrame,
    affine_transform: rasterio.transform.Affine,
    window_dimensions: tuple,
    latlong_topleft: list,
    latlong_botright: list,
) -> tuple[np.ndarray, pd.DataFrame]:
    """
    Function to mask geometries into np_map_ID using an incrementing counter.

    Parameters
    ----------
    country_rows: gpd.GeoDataFrame
        geoDataFrame filled with geometries and their GADM_ID
    affine_transform: rasterio.transform.Affine
        affine transform of current window
    window_dimensions: tuple
        dimensions of window used when reading file
    latlong_topleft: list
        [latitude, longitude] of top left corner of the window
    latlong_botright: list
        [latitude, longitude] of bottom right corner of the window

    Returns
    -------
    np_map_ID.astype("H"): np.ndarray
        np_map_ID contains an ID for each location (undefined is 0)
        dimensions are taken from window_dimensions, .astype("H") for memory savings
    id_result: pd.DataFrame
        DataFrame of the mapping from id (from counter) to GADM_ID
    """
    col_offset, row_offset, x_axis_len, y_axis_len = window_dimensions

    # Set an empty numpy array with the dimensions of the country .tif file
    # np_map_ID will contain an ID for each location (undefined is 0)
    # ID corresponds to a specific geometry in country_rows
    np_map_ID = np.zeros((y_axis_len, x_axis_len))

    # List to contain the mappings of id to GADM_ID
    id_to_GADM_ID = []

    # Loop the country_rows geoDataFrame
    for i in range(len(country_rows)):
        # Set the current geometry
        cur_geometry = country_rows.iloc[i]["geometry"]

        latitude_min = cur_geometry.bounds[1]
        latitude_max = cur_geometry.bounds[3]

        # Check if bounds of geometry overlap the window
        # In the following cases we can skip the rest of the loop
        # If the geometry is above the window
        if latitude_min > latlong_topleft[0]:
            continue
        # If the geometry is below the window
        if latitude_max < latlong_botright[0]:
            continue

        # Generate a mask for the specific geometry
        temp_mask = rasterio.features.geometry_mask(
            [cur_geometry],
            (y_axis_len, x_axis_len),
            transform=affine_transform,
            all_touched=True,
            invert=True,
        )

        # Map the values of counter value to np_map_ID
        np_map_ID[temp_mask] = i + 1

        # Store the id -> GADM_ID mapping
        id_to_GADM_ID.append([i + 1, country_rows.iloc[i]["GADM_ID"]])

    if len(id_to_GADM_ID) > 0:
        id_result = pd.DataFrame(id_to_GADM_ID).set_index(0)
    else:
        id_result = pd.DataFrame()

    # Return np_map_ID as type 'H' np.ushort
    # 'H' -> https://numpy.org/doc/stable/reference/arrays.scalars.html#numpy.ushort
    # This lowers memory usage, note: ID has to be within the range [0,65535]
    return np_map_ID.astype("H"), id_result


def sum_values_using_geomask(
    np_pop_val: np.ndarray,
    np_pop_xy: np.ndarray,
    region_geomask: np.ndarray,
    id_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """
    Function that sums all the population values in np_pop_val into the correct
    GADM_ID It uses np_pop_xy to access the key stored in region_geomask[x][y]

    The relation of this key to GADM_ID is stored in id_mapping

    Parameters
    ----------
    np_pop_val: np.ndarray
        array filled with values for each nonzero pixel in the worldpop file
    np_pop_xy: np.ndarray
        array with [x,y] coordinates of the corresponding nonzero values in np_pop_valid
    region_geomask: np.ndarray
        array with dimensions of window, values are keys that map to GADM_ID using id_mapping
    id_mapping: pd.DataFrame
        Dataframe that contains mappings of region_geomask values to GADM_IDs

    Returns
    -------
    df_pop_count: pd.DataFrame
        Dataframe with columns
            - "GADM_ID"
            - "pop" containing population of GADM_ID region
    """
    # Initialize a dictionary
    dict_id = Dict.empty(
        key_type=types.int64,
        value_type=types.int64,
    )
    dict_id[0] = 0
    counter = 1
    # Loop over ip mapping and add indices to the dictionary
    for ID_index in np.array(id_mapping.index):
        dict_id[ID_index] = counter
        counter += 1

    # Declare an array to contain population counts
    np_pop_count = np.zeros(len(id_mapping) + 1)

    # Calculate population count of region using a numba njit compiled function
    np_pop_count = loop_and_extact_val_x_y(
        np_pop_count, np_pop_val, np_pop_xy, region_geomask, dict_id
    )

    df_pop_count = pd.DataFrame(np_pop_count, columns=["pop"])
    df_pop_count["GADM_ID"] = np.append(np.array("NaN"), id_mapping.values)
    df_pop_count = df_pop_count[["GADM_ID", "pop"]]

    return df_pop_count


@njit
def loop_and_extact_val_x_y(
    np_pop_count: np.ndarray,
    np_pop_val: np.ndarray,
    np_pop_xy: np.ndarray,
    region_geomask: np.ndarray,
    dict_id: dict,
) -> np.ndarray:
    """
    Function that will be compiled using @njit (numba) It takes all the
    population values from np_pop_val and stores them in np_pop_count.

    where each location in np_pop_count is mapped to a GADM_ID through dict_id (id_mapping by extension)

    Parameters
    ----------
    np_pop_count: np.ndarray
        np.zeros array, which will store population counts
    np_pop_val: np.ndarray
        array filled with values for each nonzero pixel in the worldpop file
    np_pop_xy: np.ndarray
        array with [x,y] coordinates of the corresponding nonzero values in np_pop_valid
    region_geomask: np.ndarray
        array with dimensions of window, values are keys that map to GADM_ID using id_mapping
    dict_id: dict
        numba typed.dict containing id_mapping.index -> location in np_pop_count

    Returns
    -------
    np_pop_count: np.ndarray
        np.array containing population counts
    """
    # Loop the population data
    for i in range(len(np_pop_val)):
        cur_value = np_pop_val[i]
        cur_x, cur_y = np_pop_xy[i]

        # Set the current id to the id at the same coordinate of the geomask
        cur_id = region_geomask[int(cur_x)][int(cur_y)]

        # Add the current value to the population
        np_pop_count[dict_id[cur_id]] += cur_value

    return np_pop_count


def calculate_transform_and_coords_for_window(
    current_transform: rasterio.transform.Affine,
    window_dimensions: tuple,
    original_window=False,
) -> list:
    """
    Function which calculates the [lat,long] corners of the window given
    window_dimensions, if not(original_window) it also changes the affine
    transform to match the window.

    Parameters
    ----------
    current_transform: rasterio.transform.Affine
        affine transform of source image
    window_dimensions: tuple
        dimensions of window used when reading file
    original_window: bool
        boolean to track if window covers entire country

    Returns
    -------
    list
        A list of: [
            adjusted_transform: affine transform adjusted to window
            coordinate_topleft: [latitude, longitude] of top left corner of the window
            coordinate_botright: [latitude, longitude] of bottom right corner of the window ]
    """
    col_offset, row_offset, x_axis_len, y_axis_len = window_dimensions

    # Declare a affine transformer with given current_transform
    transformer = rasterio.transform.AffineTransformer(current_transform)

    # Obtain the coordinates of the upper left corner of window
    window_topleft_longitude, window_topleft_latitude = transformer.xy(
        row_offset, col_offset
    )

    # Obtain the coordinates of the bottom right corner of window
    window_botright_longitude, window_botright_latitude = transformer.xy(
        row_offset + y_axis_len, col_offset + x_axis_len
    )

    if original_window:
        # If the window covers the entire country then use original affine transform
        adjusted_transform = current_transform
    else:
        # Set the adjusted transform to the correct lat and long
        adjusted_transform = rasterio.Affine(
            current_transform[0],
            current_transform[1],
            window_topleft_longitude,
            current_transform[3],
            current_transform[4],
            window_topleft_latitude,
        )

    coordinate_topleft = [window_topleft_latitude, window_topleft_longitude]
    coordinate_botright = [window_botright_latitude, window_botright_longitude]

    return [adjusted_transform, coordinate_topleft, coordinate_botright]


def generate_df_tasks(
    c_code: str, mem_read_limit_per_process: int, WorldPop_inputfile: str
) -> pd.DataFrame:
    """
    Function to generate a list of tasks based on the memory constraints.

    One task represents a single window of the image

    Parameters
    ----------
    c_code: str
        country code
    mem_read_limit_per_process: int
        memory limit for src.read() operation
    WorldPop_inputfile: str
        file location of worldpop file

    Returns
    -------
    pd.DataFrame
        Dataframe of task_list
    """
    task_list = []

    # Read out dimensions and transform of image
    with rasterio.open(WorldPop_inputfile) as src:
        worldpop_y_dim, worldpop_x_dim = src.shape
        transform = src.meta["transform"]
        block_y_dim, block_x_dim = [int(i) for i in src.block_shapes[0]]

    # Rasterio doesn't support lower than float32 readout
    # Hence reading the file will take up: nbytes = 4 * y_dim * x_dim
    expected_bytes_input_read = 4 * worldpop_y_dim * worldpop_x_dim

    # Introduce a limit that avoids overfilling RAM during readout
    # Ensure worldpop_byte_limit >= 883 * 10**6 (minimum memory for 'US')
    worldpop_byte_limit = max(883, mem_read_limit_per_process) * 10**6

    if expected_bytes_input_read < worldpop_byte_limit:
        # If the rasterio read will be within byte limit
        bool_original_window = True
        # Set the window to entire image dimensions
        current_window = [0, 0, worldpop_x_dim, worldpop_y_dim]
        # Get latlong coordinates of window
        transform_and_coords = calculate_transform_and_coords_for_window(
            transform, current_window, bool_original_window
        )
        task_list.append(
            [c_code, current_window, bool_original_window] + transform_and_coords
        )
    else:
        # Reading operation has to be split up into multiple windows
        bool_original_window = False
        # Set the windows x dimension to the input x dimension (width)
        window_x_dim = worldpop_x_dim
        # From testing we can assume max x dimension will always fit in memory:
        #   Largest memory requirement is 'US' at ~882.1 MB = 4 * 512 * window_x_dim
        #   Hence worldpop_byte_limit has to be greater than 883 MB

        # As the window spans the x dimension, set column offset to 0
        window_col_offset = 0

        # Calculate the bytes for reading the window using window_x_dim (readout is float32)
        read_block_size = 4 * block_y_dim * window_x_dim

        # Calculate the amount of blocks that fit into the memory budget
        # Using the calculated x dimension
        window_block_count = int(worldpop_byte_limit // read_block_size)

        # Multiply the y_dimension by the amount of blocks in the window
        # window_y_dim will be height of the window
        window_y_dim = window_block_count * block_y_dim

        # Calculate row offsets for all windows
        window_row_offset = np.arange(0, worldpop_y_dim, window_y_dim)

        # Loop the windows and add task for each one to task_list
        for row_offset in window_row_offset:
            current_window = [window_col_offset, row_offset, window_x_dim, window_y_dim]

            transform_and_coords = calculate_transform_and_coords_for_window(
                transform, current_window, bool_original_window
            )

            task_list.append(
                [c_code, current_window, bool_original_window] + transform_and_coords
            )

    return pd.DataFrame(task_list)


def add_population_data(
    df_gadm: gpd.GeoDataFrame,
    country_codes: list,
    worldpop_method: str,
    year: int = 2020,
    update: bool = False,
    out_logging: bool = False,
    mem_read_limit_per_process: int = 1024,
    nprocesses: int = 2,
    disable_progressbar: bool = False,
) -> gpd.GeoDataFrame:
    """
    Function to add population data to arbitrary number of shapes in a country.
    It loads data from WorldPop raster files where each pixel represents the
    population in that square region. Each square polygon (or pixel) is then
    mapped into the corresponding GADM shape. Then the population in a GADM
    shape is identified by summing over all pixels mapped to that region.

    This is performed with an iterative approach:

    1. All necessary WorldPop data tiff file are downloaded
    2. The so-called windows are created to handle RAM limitations related to large WorldPop files.
       Large WorldPop files require significant RAM to handle, which may not be available,
       hence, the entire activity is decomposed into multiple windows (or tasks).
       Each window represents a subset of a raster file on which the following algorithm is applied.
       Note: when enough RAM is available only a window is created for efficiency purposes.
    3. Execute all tasks by summing the values of the pixels mapped into each GADM shape.
       Parallelization applies in this task.

    Parameters
    ----------
    df_gadm: gpd.GeoDataFrame
        Geodataframe with one Multipolygon per row
        - Essential column ["country", "geometry"]
        - Non-essential column ["GADM_ID"]
    country_codes: list
        List of country codes to download and process
    worldpop_method: str
        Method to download worldpop data, either "api" or "ftp"
    year: int
        Year of the data to download
    update: bool
        Update = true, forces re-download of files
    out_logging: bool
        If True, emits progress information via the module logger.
    mem_read_limit_per_process: int
        Memory limit for src.read() operation in MB, used to determine window size
    nprocesses: int
        Number of processes to use for parallelization, default is 2
    disable_progressbar: bool
        If True, disables the progress bar, default is False

    Returns
    -------
    df_gadm: gpd.GeoDataFrame
        Geodataframe with one Multipolygon per row
        - Same columns as input
        - Includes a new column ["pop"]
    """

    # Initialize new population column
    df_gadm["pop"] = 0.0

    # Initialize new dict to hold worldpop_inputfile strings
    dict_worldpop_file_locations = {}

    # Initialize new dataframe to hold tasks for processing
    df_tasks = pd.DataFrame()

    if out_logging:
        logger.info(
            "Stage 3 of 5: Download WorldPop datasets, method: " + str(worldpop_method)
        )

    tqdm_kwargs_download = dict(
        ascii=False,
        desc="Downloading worldpop file per country and assigning tasks",
    )
    with tqdm(total=len(country_codes), **tqdm_kwargs_download) as pbar:
        for c_code in country_codes:
            # Download worldpop image (if required) and get file location
            WorldPop_inputfile, WorldPop_filename = download_WorldPop(
                c_code, worldpop_method, year, update, out_logging
            )
            dict_worldpop_file_locations[c_code] = WorldPop_inputfile

            # Given the memory constraint, generate a task for all windows to cover the country (c_code)
            df_new_tasks = generate_df_tasks(
                c_code, mem_read_limit_per_process, WorldPop_inputfile
            )

            df_tasks = pd.concat([df_tasks, df_new_tasks], ignore_index=True)

            pbar.update(1)

    df_tasks.columns = [
        "c_code",
        "window_dimensions",
        "is_original_window",
        "affine_transform",
        "latlong_coordinate_topleft",
        "latlong_coordinate_botright",
    ]

    if out_logging:
        logger.info("Stage 4 of 5: Add population data to GADM GeoDataFrame")

    if out_logging:
        logger.info("Stage 4 of 5: Starting multiprocessing")

    tqdm_kwargs_compute = dict(
        ascii=False, desc="Compute population per window", unit=" window"
    )
    lock = mp.Lock()
    kwargs = {
        "initializer": _init_process_pop,
        "initargs": (
            df_gadm,
            df_tasks,
            dict_worldpop_file_locations,
        ),
        "processes": nprocesses,
    }
    # Spawn processes with the parameters from kwargs
    with mp.get_context("spawn").Pool(**kwargs) as pool:
        # Create a progress bar
        with tqdm(total=len(df_tasks), **tqdm_kwargs_compute) as pbar:
            # Give the pool a workload
            for df_pop_count in pool.imap_unordered(
                process_function_population, range(len(df_tasks))
            ):
                # Acquire the lock before accessing df_gadm and pbar
                with lock:
                    # Loop the regions and write population to df_gadm
                    for i in range(len(df_pop_count)):
                        gadm_id, pop_count = df_pop_count.iloc[i]
                        # Select the row with the same "GADM_ID" and set the population count
                        df_gadm.loc[df_gadm["GADM_ID"] == gadm_id, "pop"] += pop_count

                    # update bar
                    pbar.update(1)


def gadm(
    worldpop_method: str,
    gdp_method: str,
    countries: list,
    geo_crs: str,
    contended_flag: bool,
    mem_mb: int,
    layer_id: int = 2,
    update: bool = False,
    out_logging: bool = False,
    year: int = 2020,
    nprocesses: int = None,
    simplify_gadm: bool = True,
    tolerance: float = 0.01,
    minarea: float = 0.01,
) -> gpd.GeoDataFrame:
    """
    Function to create a GeoDataFrame with GADM shapes and population and gdp data.

    Parameters
    ----------
    worldpop_method: str
        Method to download worldpop data, either "api" or "ftp"
    gdp_method: str
        Method to download gdp data, either "api" or "ftp"
    countries: list
        List of country codes to download and process
    geo_crs: str
        CRS used for geographic projection, passed to GeoPandas
    contended_flag: bool
        If True, includes contended territories in the GADM layer, which may have overlapping geometries.
        If False, contended territories are excluded.
    mem_mb: int
        Memory limit in megabytes for processing WorldPop data, used to determine window size for population data processing
    layer_id: int
        GADM layer ID to download, default is 2 (admin level 2)
    update: bool
        Update = true, forces re-download of files
    out_logging: bool
        If True, emits progress information via the module logger.
    year: int
        Year of the data to download
    nprocesses: int
        Number of processes to use for parallelization, default is None (uses os.cpu_count())
    simplify_gadm: bool
        If True, simplifies the geometries in the GADM GeoDataFrame to reduce file size and speed up processing, default is True
    tolerance: float
        Tolerance parameter for geometry simplification, default is 0.01 (units depend on the CRS of the geometries)
    minarea: float
        Minimum area threshold for geometry simplification, geometries smaller than this area will be removed after simplification,
        default is 0.01 (units depend on the CRS of the geometries)

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with one Multipolygon per row, indexed by "GADM_ID"
        - Essential columns: ["country", "geometry"]
        - Optional columns: ["gdp", "pop"] (if gdp_method and worldpop_method are not False, respectively)
    """
    if out_logging:
        logger.info("Stage 3 of 5: Creation GADM GeoDataFrame")

    # download data if needed and get the desired layer_id
    df_gadm = get_GADM_layer(countries, layer_id, geo_crs, contended_flag, update)

    # select and rename columns
    df_gadm.rename(columns={"GID_0": "country"}, inplace=True)

    # drop useless columns
    df_gadm.drop(
        df_gadm.columns.difference(["country", "GADM_ID", "geometry"]),
        axis=1,
        inplace=True,
        errors="ignore",
    )

    if worldpop_method != False:
        mem_read_limit_per_process = mem_mb / nprocesses
        # add the population data to the dataset
        add_population_data(
            df_gadm,
            countries,
            worldpop_method,
            year,
            update,
            out_logging,
            mem_read_limit_per_process,
            nprocesses=nprocesses,
        )

    if gdp_method != False:
        # add the gdp data to the dataset
        add_gdp_data(
            df_gadm,
            year,
            update,
            out_logging,
            name_file_nc="GDP_PPP_1990_2015_5arcmin_v2.nc",
        )

    # renaming 3 letter to 2 letter ISO code before saving GADM file
    # In the case of a contested territory in the form 'Z00.00_0', save 'AA.00_0'
    # Include bugfix for the case of 'XXX00_0' where the "." is missing, such as for Ghana
    df_gadm["GADM_ID"] = df_gadm["country"] + df_gadm["GADM_ID"].str[3:].apply(
        lambda x: x if x.find(".") == 0 else "." + x
    )
    df_gadm.set_index("GADM_ID", inplace=True)

    if simplify_gadm:
        df_gadm["geometry"] = df_gadm["geometry"].map(
            lambda x: _simplify_polys(x, tolerance=tolerance, minarea=minarea)
        )
    df_gadm.geometry = df_gadm.geometry.apply(
        lambda r: make_valid(r) if not r.is_valid else r
    )
    df_gadm = df_gadm[df_gadm.geometry.is_valid & ~df_gadm.geometry.is_empty]

    return df_gadm


def crop_country(
    gadm_shapes: gpd.GeoDataFrame, subregion_config: dict[str, list[str]]
) -> gpd.GeoDataFrame:
    """
    Merge GADM administrative units into custom subregions and aggregate remaining units at the country level,
    returning combined geometries.

    Parameters
    ----------
    gadm_shapes : GeoDataFrame
        GeoDataFrame of GADM units indexed by administrative ID, containing a ``country`` column and a valid CRS.
    subregion_config : dict[str, list[str]]
        Mapping of subregion names to lists of GADM index values to be merged into each subregion.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame indexed by name with merged subregion and country geometries,
        preserving the input CRS.
    """

    country_shapes_new = gpd.GeoDataFrame(columns=["name", "geometry"]).set_index(
        "name"
    )
    remain_gadm_shapes = gadm_shapes.copy(deep=True)

    for sub_region in subregion_config:

        region_GADM = subregion_config[sub_region]
        sub_gadm = remain_gadm_shapes.loc[remain_gadm_shapes.index.isin(region_GADM)]
        sub_geometry = sub_gadm.union_all()

        if not sub_geometry:
            logger.warning(f"No subregion shape generated for {sub_region}")
            continue

        logger.info(
            f"Created the {sub_region} subregion based on {list(sub_gadm.index)}"
        )

        sub_country_shapes = gpd.GeoDataFrame(
            {
                "name": sub_region,
                "geometry": [sub_geometry],
            }
        ).set_index("name")

        country_shapes_new = pd.concat([country_shapes_new, sub_country_shapes])

        remain_gadm_shapes = remain_gadm_shapes.loc[
            ~remain_gadm_shapes.index.isin(region_GADM)
        ]

    for country in remain_gadm_shapes.country.unique():
        country_geometry_new = remain_gadm_shapes.query(
            "country == @country"
        ).union_all()
        country_shapes_country = gpd.GeoDataFrame(
            {
                "name": country,
                "geometry": [country_geometry_new],
            }
        ).set_index("name")

        country_shapes_new = pd.concat([country_shapes_new, country_shapes_country])

    return gpd.GeoDataFrame(
        country_shapes_new,
        crs=gadm_shapes.crs,
        geometry=country_shapes_new.geometry,
    )


def generate_points_every_km(
    gdf: gpd.GeoDataFrame,
    distance_crs: str = "EPSG:3857",
    interval_km: float = 100,
) -> gpd.GeoDataFrame:
    """
    Generate perimeter points for each Polygon/MultiPolygon in a GeoDataFrame

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame containing Polygon or MultiPolygon geometries.
    distance_crs : str
        Projected CRS in meters for distance calculation.
    interval_km : float
        Target spacing between consecutive points in kilometers.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame of sampled points with a column 'name'
        preserving the original index of the input geometry.
    """

    if gdf.crs is None:
        raise ValueError("Input GeoDataFrame must have a CRS defined.")

    interval_m = interval_km * 1000
    geo_crs = gdf.crs

    # Transformer objects
    forward_transformer = Transformer.from_crs(geo_crs, distance_crs, always_xy=True)
    backward_transformer = Transformer.from_crs(distance_crs, geo_crs, always_xy=True)

    output_rows = []

    for idx, row in gdf.iterrows():
        geometry = row.geometry

        # Project geometry to metric CRS
        projected_geom = transform(forward_transformer.transform, geometry)

        if projected_geom.geom_type == "Polygon":
            polygons = [projected_geom]
        elif projected_geom.geom_type == "MultiPolygon":
            polygons = list(projected_geom.geoms)
        else:
            continue  # skip non-polygon geometries

        for poly in polygons:
            ring = poly.exterior
            total_length = ring.length

            current_distance_along = 0
            current_point = ring.interpolate(0)

            output_rows.append(
                {
                    "name": idx,
                    "geometry": transform(
                        backward_transformer.transform, current_point
                    ),
                }
            )

            step_resolution = 1000  # meters

            while current_distance_along < total_length:
                search_distance = current_distance_along + step_resolution
                found_next = False

                while search_distance <= total_length:
                    candidate = ring.interpolate(search_distance)

                    if current_point.distance(candidate) >= interval_m:
                        current_point = candidate
                        current_distance_along = search_distance

                        output_rows.append(
                            {
                                "name": idx,
                                "geometry": transform(
                                    backward_transformer.transform, candidate
                                ),
                            }
                        )

                        found_next = True
                        break

                    search_distance += step_resolution

                if not found_next:
                    break

    # Create output GeoDataFrame
    points_gdf = gpd.GeoDataFrame(
        output_rows,
        geometry="geometry",
        crs=geo_crs,
    )

    return points_gdf


def determine_subregion_country(
    subregion_shapes: gpd.GeoDataFrame,
    country_shapes: gpd.GeoDataFrame,
    distance_crs: str,
) -> dict[str, str]:
    """
    Assign each subregion to the country with which it has the largest spatial overlap.

    Parameters
    ----------
    subregion_shapes : GeoDataFrame
        GeoDataFrame indexed by subregion name containing onshore geometries with a defined CRS.
    country_shapes : GeoDataFrame
        GeoDataFrame of country geometries used to determine which country each subregion belongs to, with a defined CRS.
    distance_crs : str
        Projected CRS used for accurate area calculations, passed to GeoPandas for area computation.

    Returns
    -------
    dict[str, str]
        Mapping of subregion name to country name, where each subregion is assigned to the country with which it has the largest spatial overlap.
    """
    intersections = gpd.overlay(
        subregion_shapes.reset_index(),
        country_shapes.reset_index(),
        how="intersection",
        keep_geom_type=False,
    )

    intersections["overlap_area"] = intersections.to_crs(distance_crs).geometry.area
    idx = intersections.groupby("name_1")["overlap_area"].idxmax()
    largest = intersections.loc[idx]

    return largest.set_index("name_1")["name_2"].to_dict()


def crop_offshore(
    subregion_shapes: gpd.GeoDataFrame,
    country_shapes: gpd.GeoDataFrame,
    offshore_shapes: gpd.GeoSeries,
    distance_crs: str = "EPSG:3857",
) -> gpd.GeoDataFrame:
    """
    Split offshore (EEZ) geometries among subregions of the same country,
    assigning each subregion a portion of the country's offshore area.

    Parameters
    ----------
    subregion_shapes : GeoDataFrame
        GeoDataFrame indexed by subregion name containing onshore geometries with a defined CRS.
    country_shapes : GeoDataFrame
        GeoDataFrame of country geometries used to determine which country each subregion belongs to.
    offshore_shapes : GeoSeries
        GeoSeries indexed by name containing offshore (e.g., EEZ) geometries to be partitioned.
    distance_crs : str, default "EPSG:3857"
        Projected CRS used for distance-based operations and Voronoi partitioning.

    Returns
    -------
    GeoDataFrame
        GeoDataFrame indexed by subregion name containing offshore geometries,
        either equal to the full country offshore area (single subregion) or Voronoi-partitioned among multiple subregions.
    """

    from build_bus_regions import custom_voronoi_partition_pts

    # Determine country for each subregion
    subregion_dict = determine_subregion_country(
        subregion_shapes,
        country_shapes,
        distance_crs=distance_crs,
    )

    subregion_shapes = subregion_shapes.copy()
    subregion_shapes["country"] = subregion_shapes.index.map(subregion_dict)

    results = []

    # Group once instead of value_counts + filtering
    for country, outline in offshore_shapes.items():
        sub_shape = subregion_shapes[subregion_shapes["country"] == country].copy()

        # If only one subregion → just use outline
        if len(sub_shape) == 1:
            new_gdf = gpd.GeoDataFrame(
                geometry=[outline],
                index=sub_shape.index,
                crs=subregion_shapes.crs,
            )

        # Multiple subregions → Voronoi partition
        else:
            points_gdf = generate_points_every_km(
                sub_shape,
                distance_crs=distance_crs,
                interval_km=100,
            )

            coords = points_gdf.geometry.get_coordinates()[["x", "y"]]

            if len(coords) > 1000:
                logger.warning(
                    f"Voronoi partitioning for {country} with {len(coords)} points may be slow. "
                    "Consider increasing the `interval_km` parameter to reduce the number of points, "
                    "or set `build_shape_options: simplify_gadm:` as True."
                )

            voronoi_geoms = custom_voronoi_partition_pts(
                coords,
                outline,
                add_bounds_shape=True,
                multiplier=5,
            )

            points_gdf = points_gdf.copy()
            points_gdf["geometry"] = voronoi_geoms

            new_gdf = points_gdf.dissolve(by="name")
            new_gdf = new_gdf[~new_gdf.geometry.is_empty]

        results.append(new_gdf)

    if not results:
        return gpd.GeoDataFrame(
            columns=["name", "geometry"], geometry="geometry"
        ).set_index("name")

    return pd.concat(results)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_shapes")
    configure_logging(snakemake)

    out = snakemake.output

    EEZ_gpkg = snakemake.input["eez"]
    mem_mb = snakemake.resources["mem_mb"]

    countries_list = snakemake.params.countries
    geo_crs = snakemake.params.crs["geo_crs"]
    distance_crs = snakemake.params.crs["distance_crs"]

    layer_id = snakemake.params.build_shape_options["gadm_layer_id"]
    update = snakemake.params.build_shape_options["update_file"]
    out_logging = snakemake.params.build_shape_options["out_logging"]
    year = snakemake.params.build_shape_options["year"]
    nprocesses = snakemake.params.build_shape_options["nprocesses"]
    contended_flag = snakemake.params.build_shape_options["contended_flag"]
    worldpop_method = snakemake.params.build_shape_options["worldpop_method"]
    gdp_method = snakemake.params.build_shape_options["gdp_method"]
    tolerance = snakemake.params.build_shape_options["simplify_tolerance"]
    simplify_gadm = snakemake.params.build_shape_options["simplify_gadm"]
    minarea = snakemake.params.build_shape_options["minarea"]

    country_shapes = countries(
        countries_list,
        geo_crs,
        contended_flag,
        update,
        out_logging,
        tolerance=tolerance,
    )
    country_shapes.to_file(out.country_shapes)

    offshore_shapes = eez(
        countries_list,
        geo_crs,
        country_shapes,
        EEZ_gpkg,
        out_logging=out_logging,
        tolerance=tolerance,
        minarea=minarea,
        simplify_gadm=simplify_gadm,
    )

    offshore_shapes.reset_index().to_file(out.offshore_shapes)

    africa_shape = gpd.GeoDataFrame(
        geometry=[country_cover(country_shapes, offshore_shapes.geometry)]
    )
    africa_shape.reset_index().to_file(out.africa_shape)

    gadm_shapes = gadm(
        worldpop_method,
        gdp_method,
        countries_list,
        geo_crs,
        contended_flag,
        mem_mb,
        layer_id,
        update,
        out_logging,
        year,
        nprocesses=nprocesses,
        simplify_gadm=simplify_gadm,
        tolerance=tolerance,
        minarea=minarea,
    )

    save_to_geojson(gadm_shapes, out.gadm_shapes)

    subregion_config = snakemake.params.subregion
    subregion_method = subregion_config.get("method")

    # Prepare subregion shapes based on GADM, otherwise place empty GeoDataFrame.
    if subregion_method == "gadm":
        define_by_gadm = subregion_config["define_by_gadm"]
        subregion_shapes = crop_country(gadm_shapes, define_by_gadm)

    else:
        subregion_shapes = gpd.GeoDataFrame(
            columns=["name", "geometry"],
            geometry="geometry",
            crs=gadm_shapes.crs,
        ).set_index("name")

    save_to_geojson(subregion_shapes, out.subregion_shapes)

    # If the custom method is selected and no custom offshore path is provided,
    # generate offshore shapes from the custom subregion shapes.
    has_custom_offshore = subregion_config["path_custom_offshore"]
    if subregion_method == "custom" and not has_custom_offshore:
        custom_path = subregion_config["path_custom_shapes"]
        subregion_shapes = gpd.read_file(custom_path).set_index("name")

    if not subregion_shapes.empty:
        subregion_offshore = crop_offshore(
            subregion_shapes,
            country_shapes,
            offshore_shapes,
            distance_crs=distance_crs,
        )

    else:
        subregion_offshore = subregion_shapes  # empty

    subregion_offshore.reset_index().to_file(out.subregion_offshore)
