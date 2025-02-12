
import numpy as np
from pyproj import Proj, Transformer


def coords_to_relative_utm(coords):
    """
    Transforms a list of longitude and latitude coordinates to UTM coordinates relative to the centroid.

    Parameters:
    - coords: list of tuples
        List containing (latitude, longitude) tuples.

    Returns:
    - relative_coords_km: numpy.ndarray
        Array of transformed coordinates in kilometers relative to the centroid.
    """
    coords_array = np.array(coords)
    lons = coords_array[:, 1]
    lats = coords_array[:, 0]

    centroid_lon = np.mean(lons)
    centroid_lat = np.mean(lats)

    utm_zone = int(np.floor((centroid_lon + 180) / 6) % 60) + 1
    hemisphere = 'north' if centroid_lat >= 0 else 'south'

    utm_proj = Proj(proj='utm', zone=utm_zone, hemisphere=hemisphere)

    transformer = Transformer.from_proj(
        proj_from='epsg:4326',  # WGS84 Latitude and Longitude
        proj_to=utm_proj,
        always_xy=True
    )

    eastings, northings = transformer.transform(lons, lats)
    centroid_easting, centroid_northing = transformer.transform(centroid_lon, centroid_lat)

    relative_eastings = eastings - centroid_easting
    relative_northings = northings - centroid_northing

    relative_coords = np.column_stack((relative_eastings, relative_northings))

    return relative_coords / 1000  # Convert to kilometers
