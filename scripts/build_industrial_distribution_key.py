"""Build industrial distribution keys from hotmaps database."""

import uuid
import pandas as pd
import geopandas as gpd

from itertools import product
from distutils.version import StrictVersion
from helpers import locate_bus
from helpers import three_2_two_digits_country

gpd_version = StrictVersion(gpd.__version__)


def locate_missing_industrial_sites(df):
    """
    Locate industrial sites without valid locations based on
    city and countries. Should only be used if the model's
    spatial resolution is coarser than individual cities.
    """

    try:
        from geopy.geocoders import Nominatim
        from geopy.extra.rate_limiter import RateLimiter
    except:
        raise ModuleNotFoundError("Optional dependency 'geopy' not found."
                                  "Install via 'conda install -c conda-forge geopy'"
                                  "or set 'industry: hotmaps_locate_missing: false'.")

    locator = Nominatim(user_agent=str(uuid.uuid4()))
    geocode = RateLimiter(locator.geocode, min_delay_seconds=2)

    def locate_missing(s):

        if pd.isna(s.City) or s.City == "CONFIDENTIAL":
            return None

        loc = geocode([s.City, s.Country], geometry='wkt')
        if loc is not None:
            print(f"Found:\t{loc}\nFor:\t{s['City']}, {s['Country']}\n")
            return f"POINT({loc.longitude} {loc.latitude})"
        else:
            return None

    missing = df.index[df.geom.isna()]
    df.loc[missing, 'coordinates'] = df.loc[missing].apply(locate_missing, axis=1)

    # report stats
    num_still_missing = df.coordinates.isna().sum()
    num_found = len(missing) - num_still_missing
    share_missing = len(missing) / len(df) * 100
    share_still_missing = num_still_missing / len(df) * 100
    print(f"Found {num_found} missing locations.",
          f"Share of missing locations reduced from {share_missing:.2f}% to {share_still_missing:.2f}%.")

    return df


def prepare_GID_database(regions):
    """
    Load hotmaps database of industrial sites and map onto bus regions.
    Build industrial demand... Change name and add other functions. 
    Function similar to aviation/shipping. Use functions to disaggregate. 
    Only cement not steel - proof of concept. 
    Change hotmaps to more descriptive name, etc.
    """

#<<<<<<< HEAD
    df = pd.read_csv(snakemake.input.GID_industrial_database, sep=",", header=0, encoding=('latin1'))
#=======
    #df = pd.read_csv(snakemake.input.hotmaps_industrial_database, sep=",", header=0)
#>>>>>>> ffded8f (industrial_distribution_key)

    df["gadm_{}".format(gadm_level)] = df[["x", "y", "country"]].apply(
        lambda site: locate_bus(site[["x", "y"]].astype("float"), site["country"], gadm_level),
        axis=1,
    )
    
    df["gadm_{}".format(gadm_level)] = df["gadm_{}".format(gadm_level)].apply(
     lambda cocode: three_2_two_digits_country(cocode[:3]) + cocode[3:])

    return df 


def build_nodal_distribution_key(GID_data, regions): #returns percentage of co2 emissions
    """Build nodal distribution keys for each sector."""

    sectors = GID_data.Sector.unique() #TODO add more than just cement to data

    countries = regions.index.str[:2].unique() 

    keys = pd.DataFrame(index=regions.index, columns=sectors, dtype=float)

    pop = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    
    pop['country'] = pop.index.str[:2]
    keys['population'] = pop['total'].values

    for sector, country in product(sectors, countries):

        regions_ct = regions.index[regions.index.str.contains(country)]

        facilities = GID_data.query("country == @country and Sector == @sector")

        if not facilities.empty:
            emissions = facilities["Total CO2 emission"]
            if emissions.sum() == 0:
                key = pd.Series(1 / len(facilities), facilities.index)
            else:
                #TODO BEWARE: this is a strong assumption
                emissions = emissions.fillna(0)
                key = emissions / emissions.sum()
            key = key.groupby(facilities["gadm_{}".format(gadm_level)]).sum().reindex(regions_ct, fill_value=0.)
        else:
            key = keys.loc[regions_ct, 'population']

        keys.loc[regions_ct, sector] = key

    return keys


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helpers import mock_snakemake
        snakemake = mock_snakemake(
            'build_industrial_distribution_key',
            simpl='',
            clusters=12,
        )

    options = snakemake.config["sector"]
    gadm_level = options["gadm_level"]
    
    regions = gpd.read_file(snakemake.input.regions_onshore)
    
    
    #regions["name"] = regions["name"].apply(
    #lambda name: three_2_two_digits_country(name[:3]) + name[3:])
    
    regions = regions.set_index('name')
    
    GID_data = prepare_GID_database(regions)

    keys = build_nodal_distribution_key(GID_data, regions)

    keys.to_csv(snakemake.output.industrial_distribution_key)