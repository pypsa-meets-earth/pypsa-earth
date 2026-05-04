# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-

import os
import pathlib
import sys

import fiona
import geopandas as gpd
import numpy as np

sys.path.append("./scripts")

from build_shapes import (
    _simplify_polys,
    add_population_data,
    convert_GDP,
    countries,
    country_cover,
    download_GADM,
    download_WorldPop_API,
    download_WorldPop_standard,
    gadm,
    get_GADM_filename,
    get_GADM_layer,
    save_to_geojson,
)

path_cwd = str(pathlib.Path.cwd())


def test_simplify_polys(get_config_dict):
    """
    Verify what is returned by _simplify_polys.
    """

    config_dict = get_config_dict

    countries_list = ["NG"]
    geo_crs = config_dict["crs"]["geo_crs"]

    update = config_dict["build_shape_options"]["update_file"]
    out_logging = config_dict["build_shape_options"]["out_logging"]
    contended_flag = config_dict["build_shape_options"]["contended_flag"]

    country_shapes_df = countries(
        countries_list,
        geo_crs,
        contended_flag,
        update,
        out_logging,
    )

    simplified_poly = _simplify_polys(country_shapes_df)

    simplified_poly_df = gpd.GeoDataFrame(
        geometry=[
            country_cover(
                simplified_poly, eez_shapes=None, out_logging=False, distance=0.02
            )
        ]
    )
    simplified_poly_df["area"] = simplified_poly_df.area
    simplified_poly_df["centroid"] = simplified_poly_df.centroid
    simplified_poly_df["perimeter"] = simplified_poly_df.length
    print(simplified_poly_df["perimeter"][0])
    assert np.round(simplified_poly_df.area[0], 2) == 75.75
    assert np.round(simplified_poly_df["perimeter"][0], 2) == 42.51


def test_get_countries_shapes(get_config_dict):
    """
    Verify what is returned by get_countries_shapes.
    """

    config_dict = get_config_dict

    countries_list = ["XK"]
    geo_crs = config_dict["crs"]["geo_crs"]

    update = config_dict["build_shape_options"]["update_file"]
    out_logging = config_dict["build_shape_options"]["out_logging"]
    contended_flag = config_dict["build_shape_options"]["contended_flag"]

    country_shapes_df = countries(
        countries_list,
        geo_crs,
        contended_flag,
        update,
        out_logging,
    )

    assert country_shapes_df.shape == (1,)
    assert country_shapes_df.index.unique().tolist() == ["XK"]


def test_country_cover(get_config_dict):
    """
    Verify what is returned by country_cover.
    """

    config_dict = get_config_dict

    countries_list = ["NG"]
    geo_crs = config_dict["crs"]["geo_crs"]

    update = config_dict["build_shape_options"]["update_file"]
    out_logging = config_dict["build_shape_options"]["out_logging"]
    contended_flag = config_dict["build_shape_options"]["contended_flag"]

    country_shapes_df = countries(
        countries_list,
        geo_crs,
        contended_flag,
        update,
        out_logging,
    )

    africa_shapes_df = gpd.GeoDataFrame(
        geometry=[
            country_cover(
                country_shapes_df, eez_shapes=None, out_logging=False, distance=0.02
            )
        ]
    )
    africa_shapes_df["area"] = africa_shapes_df.area
    africa_shapes_df["centroid"] = africa_shapes_df.centroid
    africa_shapes_df["perimeter"] = africa_shapes_df.length
    print(africa_shapes_df["perimeter"])
    assert np.round(africa_shapes_df.area[0], 2) == 75.75
    assert np.round(africa_shapes_df["perimeter"][0], 2) == 42.51


def test_download_world_pop_standard(get_config_dict):
    """
    Verify what is returned by download_WorldPop_standard.
    """

    config_dict = get_config_dict
    update_val = config_dict["build_shape_options"]["update_file"]
    out_logging_val = config_dict["build_shape_options"]["out_logging"]

    world_pop_input_file, world_pop_file_name = download_WorldPop_standard(
        "NG",
        year=2020,
        update=update_val,
        out_logging=out_logging_val,
        size_min=300,
    )
    assert world_pop_file_name == "nga_ppp_2020_UNadj_constrained.tif"


def test_download_world_pop_api():
    """
    Verify what is returned by download_WorldPop_API.
    """
    world_pop_input_file, world_pop_file_name = download_WorldPop_API(
        "NG", year=2020, size_min=300
    )
    assert world_pop_file_name == "nga_ppp_2020_UNadj_constrained.tif"


def test_get_gadm_shapes(get_config_dict):
    """
    Verify what is returned by get_gadm_shapes.
    """
    config_dict = get_config_dict

    mem_mb = 3096

    countries_list = ["XK"]
    geo_crs = config_dict["crs"]["geo_crs"]

    layer_id = config_dict["build_shape_options"]["gadm_layer_id"]
    update = config_dict["build_shape_options"]["update_file"]
    out_logging = config_dict["build_shape_options"]["out_logging"]
    year = config_dict["build_shape_options"]["year"]
    nprocesses = config_dict["build_shape_options"]["nprocesses"]
    contended_flag = config_dict["build_shape_options"]["contended_flag"]
    worldpop_method = config_dict["build_shape_options"]["worldpop_method"]
    gdp_method = False  # config_dict["build_shape_options"]["gdp_method"]  # TODO: enable with a tutorial dataset

    gadm_shapes_df = gadm(
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
    )

    assert gadm_shapes_df.shape == (7, 3)
    assert gadm_shapes_df.index.unique().tolist() == [f"XK.{x}_1" for x in range(1, 8)]
    assert round(gadm_shapes_df.loc["XK.1_1"]["pop"], -2) == 207500


def test_add_population_data(get_config_dict):
    """
    Verify what is returned by add_population_data.
    """
    config_dict = get_config_dict

    mem_mb = 3096

    countries_list = ["XK"]
    geo_crs = config_dict["crs"]["geo_crs"]

    layer_id = config_dict["build_shape_options"]["gadm_layer_id"]
    update = config_dict["build_shape_options"]["update_file"]
    out_logging = config_dict["build_shape_options"]["out_logging"]
    year = config_dict["build_shape_options"]["year"]
    nprocesses = config_dict["build_shape_options"]["nprocesses"]
    contended_flag = config_dict["build_shape_options"]["contended_flag"]
    worldpop_method = config_dict["build_shape_options"]["worldpop_method"]

    mem_read_limit_per_process = mem_mb / nprocesses

    df_gadm = get_GADM_layer(
        countries_list,
        layer_id,
        geo_crs,
        contended_flag,
        update,
        out_logging,
    )

    # select and rename columns
    df_gadm.rename(columns={"GID_0": "country"}, inplace=True)

    # drop useless columns
    df_gadm.drop(
        df_gadm.columns.difference(["country", "GADM_ID", "geometry"]),
        axis=1,
        inplace=True,
        errors="ignore",
    )

    add_population_data(
        df_gadm,
        countries_list,
        worldpop_method,
        year,
        update,
        out_logging,
        mem_read_limit_per_process,
        nprocesses=nprocesses,
    )

    assert np.round(df_gadm["pop"].values[0], 0) == 207474.0
    assert np.round(df_gadm["pop"].values[1], 0) == 208332.0
    assert np.round(df_gadm["pop"].values[2], 0) == 257191.0
    assert np.round(df_gadm["pop"].values[3], 0) == 215703.0
    assert np.round(df_gadm["pop"].values[4], 0) == 610695.0
    assert np.round(df_gadm["pop"].values[5], 0) == 420344.0
    assert np.round(df_gadm["pop"].values[6], 0) == 215316.0


# TODO: re-enable with a tutorial dataset or better data handling
# def test_load_gdp(get_config_dict):
#     """
#     Verify what is returned by load_gdp.
#     """
#     config_dict = get_config_dict

#     out_logging = config_dict["build_shape_options"]["out_logging"]
#     year = config_dict["build_shape_options"]["year"]
#     name_file_nc = "GDP_PPP_1990_2015_5arcmin_v2.nc"
#     GDP_tif, name_tif = convert_GDP(name_file_nc, year, out_logging)
#     assert name_tif == "GDP_PPP_1990_2015_5arcmin_v2.tif"


def test_get_gadm_filename():
    """
    Verify what is returned by get_gadm_filename.
    """
    # Kosovo
    assert get_GADM_filename("XK") == "gadm41_XKO"
    # Clipperton island
    assert get_GADM_filename("CP") == "gadm41_XCL"
    # Saint-Martin
    assert get_GADM_filename("SX") == "gadm41_MAF"
    # French Southern Territories
    assert get_GADM_filename("TF") == "gadm41_ATF"
    # Micronesia (Federated States of) with different file_prefix
    assert get_GADM_filename("US") == "gadm41_USA"


def test_download_GADM():
    """
    Verify what is returned by download_GADM.
    """
    gadm_input_file_gpkg_41, gadm_filename_41 = download_GADM(
        "XK",
        update=True,
    )
    assert gadm_input_file_gpkg_41 == os.path.join(
        path_cwd, "data/gadm/gadm41_XKO/gadm41_XKO.gpkg"
    )
    assert gadm_filename_41 == "gadm41_XKO"
    list_layers_41 = fiona.listlayers(gadm_input_file_gpkg_41)
    assert list_layers_41[0] == "ADM_ADM_0"
    assert list_layers_41[1] == "ADM_ADM_1"
    assert list_layers_41[2] == "ADM_ADM_2"
