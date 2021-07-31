# import
import logging
import os
import shutil
import sys
import zipfile

import fiona
import geopandas as gpd
import geoplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from esy.osmfilter import Node
from esy.osmfilter import osm_info as osm_info
from esy.osmfilter import osm_pickle as osm_pickle
from esy.osmfilter import Relation
from esy.osmfilter import run_filter
from esy.osmfilter import Way
from iso_country_codes import AFRICA_CC
from shapely.geometry import LineString
from shapely.geometry import Point
from shapely.geometry import Polygon

# IMPORTANT: RUN SCRIPT FROM THIS SCRIPTS DIRECTORY i.e data_exploration/ TODO: make more robust
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("./../../scripts")

# from ..scripts.iso_country_codes import AFRICA_CC

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_shapes")
    configure_logging(snakemake)

    out = snakemake.output

    country_shapes = countries()
    save_to_geojson(country_shapes, out.country_shapes)

    offshore_shapes = eez(country_shapes)
    save_to_geojson(offshore_shapes, out.offshore_shapes)

    europe_shape = country_cover(country_shapes, offshore_shapes)
    save_to_geojson(gpd.GeoSeries(europe_shape), out.europe_shape)

    nuts3_shapes = nuts3(country_shapes)
    save_to_geojson(nuts3_shapes, out.nuts3_shapes)
