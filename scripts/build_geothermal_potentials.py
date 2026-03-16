
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from pathlib import Path


# --------------------------------------------------
# Pfade relativ zum Skript
# --------------------------------------------------



INPUT_FILE = snakemake.input.geothermal
REGIONS_FILE = snakemake.input.regions
OUTPUT_FILE = snakemake.output[0]



# --------------------------------------------------
# Geothermie CSV laden
# --------------------------------------------------

def load_geothermal_data(path):

    df = pd.read_csv(path, sep=";")
    print(df.columns)

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.lon, df.lat),
        crs="EPSG:4326",
    )

    return gdf





# --------------------------------------------------
# 500 m × 500 m Rasterzellen erzeugen
# --------------------------------------------------

def build_cells(points):

    points = points.to_crs(3857)

    half = 250  # Meter

    cells = points.copy()

    cells["geometry"] = points.geometry.apply(
        lambda p: box(
            p.x - half,
            p.y - half,
            p.x + half,
            p.y + half,
        )
    )

    return cells


# --------------------------------------------------
# Potenzial pro Region berechnen
# --------------------------------------------------

def geothermal_potential_by_region(cells, regions):

    cells = cells.to_crs(regions.crs)

    joined = gpd.sjoin(
        cells,
        regions[["name", "geometry"]],
        how="inner",
        predicate="intersects",
    )

    potentials = (
        joined.groupby("name")["potential_kw"]
        .sum()
        .div(1000) #umrechnung in mw
        .to_frame(name="potential [MWel]") # in mw
    )
    
    # potentials = potentials.reindex(regions["name"]).fillna(0)
    # potentials = potentials.to_frame(name="potential [MWel]")
    potentials.index.name = "name"
    

    return potentials


# --------------------------------------------------
# Hauptworkflow
# --------------------------------------------------

def build_geothermal_potentials():

    print("Loading geothermal data...")
    points = load_geothermal_data(INPUT_FILE)

    print("Building raster cells...")
    cells = build_cells(points)
    #cells = points.to_crs(3857)

    print("Loading regions...")
    regions = gpd.read_file(REGIONS_FILE)
    regions = regions.to_crs(3857)

   

    print("Calculating geothermal potentials per region...")
    potentials = geothermal_potential_by_region(cells, regions)

    print("Saving results...")
    potentials.to_csv(OUTPUT_FILE)

    print("Finished.")


# --------------------------------------------------
# Skript ausführen
# --------------------------------------------------

if __name__ == "__main__":
    build_geothermal_potentials()