
# import
import os
# import sys

# IMPORTANT: RUN SCRIPT FROM THIS SCRIPTS DIRECTORY i.e data_exploration/ TODO: make more robust
# os.chdir(os.path.dirname(os.path.abspath(__file__)))

from _helpers import configure_logging, _sets_path_to_root, _two_2_three_digits_country, _three_2_two_digits_country

import logging
import shutil

import geopandas as gpd
import fiona

import numpy as np
from itertools import takewhile

import requests
import zipfile
from operator import attrgetter
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import cascaded_union

from iso_country_codes import AFRICA_CC
#from ..scripts.iso_country_codes import AFRICA_CC

from shapely.geometry import LineString, Point, Polygon

logger = logging.getLogger(__name__)
_sets_path_to_root("pypsa-africa")


def _get_country(target, **keys):
    """
    Function to convert country codes using pycountry
    
    Parameters
    ----------
    target: str
        Desired type of country code.
        Examples:
            'alpha_3' for 3-digit
            'alpha_2' for 2-digit
            'name' for full country name
    keys: 
        Specification of the country name and reference system.
        Examples:
            alpha_3="ZAF" for 3-digit
            alpha_2="ZA" for 2-digit
            name="South Africa" for full country name
    
    Returns
    -------
    country code as requested in keys or np.nan, when country code is not recognized

    Example of usage
    -------
    _get_country('alpha_3', alpha_2="ZA")
    _get_country('alpha_2', alpha_3="ZAF")
    _get_country('name', alpha_2="ZA")
    
    """
    assert len(keys) == 1
    try:
        return getattr(pyc.countries.get(**keys), target)
    except (KeyError, AttributeError):
        return np.nan


def download_GADM(country_code, update=False):
    """
    Download gpkg file from GADM for a given country code

    Parameters
    ----------
    country_code : str 
        Three letter country codes of the downloaded files 
    update : bool 
        Name of the network component 
        Update = true, forces re-download of files

    Returns
    -------
    gpkg file per country

    """

    GADM_filename = f"gadm36_{_two_2_three_digits_country(country_code)}"
    GADM_url = f"https://biogeo.ucdavis.edu/data/gadm3.6/gpkg/{GADM_filename}_gpkg.zip"

    GADM_inputfile_zip = os.path.join(
        os.getcwd(), "data", "gadm", GADM_filename, GADM_filename + ".zip"
    )  # Input filepath zip

    GADM_inputfile_gpkg = os.path.join(
        os.getcwd(), "data", "gadm", GADM_filename, GADM_filename + ".gpkg"
    )  # Input filepath gpkg

    if not os.path.exists(GADM_inputfile_zip) or update is True:
        print(f"{GADM_filename} does not exist, downloading to {GADM_inputfile_zip}")
        #  create data/osm directory
        os.makedirs(os.path.dirname(GADM_inputfile_zip), exist_ok=True)

        with requests.get(GADM_url, stream=True) as r:
            with open(GADM_inputfile_zip, "wb") as f:
                shutil.copyfileobj(r.raw, f)

        with zipfile.ZipFile(GADM_inputfile_zip,"r") as zip_ref:
            zip_ref.extractall(os.path.dirname(GADM_inputfile_zip))
    
    return GADM_inputfile_gpkg, GADM_filename


def get_GADM_layer(country_list, layer_id, update=False):
    """
    Function to retrive a specific layer id of a geopackage for a selection of countries
    """

    # initialization of the geoDataFrame
    geodf_GADM = gpd.GeoDataFrame()

    for country_code in country_list:
        # download file gpkg
        file_gpkg, name_file = download_GADM(country_code, False)

        # get layers of a geopackage
        list_layers = fiona.listlayers(file_gpkg)

        # read gpkg file
        geodf_temp = gpd.read_file(file_gpkg, layer = list_layers[layer_id])

        # append geodataframes
        geodf_GADM = geodf_GADM.append(geodf_temp)
    
    geodf_GADM.reset_index(drop=True, inplace=True)

    return geodf_GADM


def _simplify_polys(polys, minarea=0.1, tolerance=0.01, filterremote=True):
    if isinstance(polys, MultiPolygon):
        polys = sorted(polys, key=attrgetter('area'), reverse=True)
        mainpoly = polys[0]
        mainlength = np.sqrt(mainpoly.area/(2.*np.pi))
        if mainpoly.area > minarea:
            polys = MultiPolygon([p
                                  for p in takewhile(lambda p: p.area > minarea, polys)
                                  if not filterremote or (mainpoly.distance(p) < mainlength)])
        else:
            polys = mainpoly
    return polys.simplify(tolerance=tolerance)



def countries(update = False):
    countries = snakemake.config['countries']

    # download data if needed and get the layer id 0, corresponding to the countries
    df_countries = get_GADM_layer(countries, 0, update)

    # select and rename columns
    df_countries = df_countries[["GID_0", "geometry"]].copy()
    df_countries.rename(columns = {"GID_0": "name"}, inplace=True)

    # set index and simplify polygons
    ret_df = df_countries.set_index('name')['geometry'].map(_simplify_polys)

    return ret_df


def country_cover(country_shapes, eez_shapes=None):
    shapes = list(country_shapes)
    if eez_shapes is not None:
        shapes += list(eez_shapes)

    europe_shape = cascaded_union(shapes)
    if isinstance(europe_shape, MultiPolygon):
        europe_shape = max(europe_shape, key=attrgetter('area'))
    return Polygon(shell=europe_shape.exterior)


def save_to_geojson(df, fn):
    if os.path.exists(fn):
        os.unlink(fn) # remove file if it exists
    if not isinstance(df, gpd.GeoDataFrame):
        df = gpd.GeoDataFrame(dict(geometry=df))
    df = df.reset_index()
    schema = {**gpd.io.file.infer_schema(df), 'geometry': 'Unknown'}
    df.to_file(fn, driver='GeoJSON', schema=schema)


def eez(update = False, tol=1e-3):
    countries = snakemake.config['countries']

    # download data if needed and get the last layer id [-1], corresponding to the highest resolution
    df_eez = get_GADM_layer(countries, -1, update)

    # select and rename columns
    # df_eez = df_eez[["GID_0", "geometry"]].copy()
    df_eez.rename(columns = {"GID_0": "name"}, inplace=True)

    # set index and simplify polygons
    ret_df = df_eez.set_index('name')['geometry'].map(lambda s: _simplify_polys(s, filterremote=False))

    ret_df = gpd.GeoSeries({k:v for k,v in ret_df.iteritems() if v.distance(country_shapes[k]) < tol})

    df = gpd.read_file(snakemake.input.eez)
    df = df.loc[df['ISO_3digit'].isin([_get_country('alpha_3', alpha_2=c) for c in snakemake.config['countries']])]
    df['name'] = df['ISO_3digit'].map(lambda c: _get_country('alpha_2', alpha_3=c))
    s = df.set_index('name').geometry.map(lambda s: _simplify_polys(s, filterremote=False))
    s = gpd.GeoSeries({k:v for k,v in s.iteritems() if v.distance(country_shapes[k]) < 1e-3})
    s.index.name = "name"
    return s


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('build_shapes')
    configure_logging(snakemake)

    out = snakemake.output

    country_shapes = countries()
    save_to_geojson(country_shapes, out.country_shapes)

    offshore_shapes = eez()
    save_to_geojson(offshore_shapes, out.offshore_shapes)

    africa_shape = country_cover(country_shapes, offshore_shapes)
    save_to_geojson(gpd.GeoSeries(africa_shape), out.africa_shape)

    nuts3_shapes = nuts3(country_shapes)
    save_to_geojson(nuts3_shapes, out.nuts3_shapes)