import os
import sys

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.append("../../scripts")

import geopandas as gpd
import numpy as np
import pandas as pd
from iso_country_codes import AFRICA_CC
from shapely.geometry import LineString, Point, Polygon


def clean_data():
# add cleaning here
    return None



if __name__ == "__main__":
    clean_data()