import atlite
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import rasterio as rio
import shapely
import xarray as xr

from rasterio.plot import show
from rasterio.warp import calculate_default_transform, reproject, Resampling, transform_bounds

from shapely.geometry import box, Polygon
from shapely.ops import transform

natura_tiff_path = "./resources/natura.tiff"
# please change according to the region of interest
cutout_path = "./cutouts/africa-2013-era5.nc"

# work with nc data ------------------------------------------------------------

nc_data = atlite.Cutout(cutout_path)
nc_geom = shapely.wkt.loads(box(*nc_data.bounds).wkt)

# work with natura.tiff --------------------------------------------------------

data = rio.open(natura_tiff_path)
 
print("*--------* data.crs *----------*")
print(data.crs)

# transform the bounds range coordinates
geo_bounds = transform_bounds(data.crs, "EPSG:4326", 
        data.bounds[0], data.bounds[1], data.bounds[2], data.bounds[3]
    )
# build a polygon first to reproject after
tiff_geom_orig = shapely.wkt.loads(box(*data.bounds).wkt)

# transform coordinates --------------------------------------------------------

# NB strictly speaking, EPSG:4326 is an unprojected CRS
era5_crs = pyproj.CRS("EPSG:4326")
tiff_crs = pyproj.CRS(data.crs)

project_geo_tiff = pyproj.Transformer.from_crs(tiff_crs, era5_crs, always_xy=True).transform
tiff_geo = transform(project_geo_tiff, tiff_geom_orig)

# draw pictures ----------------------------------------------------------------

fig, ax = plt.subplots()
ax2 = plt.axes(projection = ccrs.PlateCarree())
ax2.coastlines()
plt.plot(*nc_geom.exterior.xy)
plt.plot(*tiff_geo.exterior.xy)
ax2.add_feature(cfeature.RIVERS)
plt.savefig("natura_cutout_overlay.pdf")



