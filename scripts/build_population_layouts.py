# -*- coding: utf-8 -*-
"""Build mapping between grid cells and population (total, urban, rural)"""
import multiprocessing as mp
import os

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from vresutils import shapes as vshapes

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake, sets_path_to_root

        snakemake = mock_snakemake("build_population_layouts")
        sets_path_to_root("pypsa-earth-sec")

    cutout_path = (
        snakemake.input.cutout
    )  # os.path.abspath(snakemake.config["atlite"]["cutout"])
    cutout = atlite.Cutout(cutout_path)

    grid_cells = cutout.grid_cells()

    # nuts3 has columns country, gdp, pop, geometry
    nuts3 = gpd.read_file(snakemake.input.nuts3_shapes).set_index("GADM_ID")

    # Set value of population to same dimension as in PyPSA-Eur-Sec, where the value is given in 1e3
    nuts3["pop"] = nuts3["pop"] / 1000

    # Indicator matrix NUTS3 -> grid cells
    I = atlite.cutout.compute_indicatormatrix(nuts3.geometry, grid_cells)

    # Indicator matrix grid_cells -> NUTS3; inprinciple Iinv*I is identity
    # but imprecisions mean not perfect
    Iinv = cutout.indicatormatrix(nuts3.geometry)

    countries = np.sort(nuts3.country.unique())
    # countries = np.array(["MA"])
    urban_fraction = (
        pd.read_csv(
            snakemake.input.urban_percent,
            header=None,
            index_col=0,
            names=["fraction"],
            squeeze=True,
        )
        / 100.0
    )

    # fill missing Balkans values
    # missing = ["AL", "ME", "MK"]
    # reference = ["RS", "BA"]
    # average = urban_fraction[reference].mean()
    # fill_values = pd.Series({ct: average for ct in missing})
    # urban_fraction = urban_fraction.append(fill_values)

    # population in each grid cell
    pop_cells = pd.Series(I.dot(nuts3["pop"]))

    # in km^2
    with mp.Pool(processes=snakemake.threads) as pool:
        cell_areas = pd.Series(pool.map(vshapes.area, grid_cells)) / 1e6

    # pop per km^2
    density_cells = pop_cells / cell_areas

    # rural or urban population in grid cell
    pop_rural = pd.Series(0.0, density_cells.index)
    pop_urban = pd.Series(0.0, density_cells.index)

    for ct in countries:

        indicator_nuts3_ct = nuts3.country.apply(lambda x: 1.0 if x == ct else 0.0)

        indicator_cells_ct = pd.Series(Iinv.T.dot(indicator_nuts3_ct))

        density_cells_ct = indicator_cells_ct * density_cells

        pop_cells_ct = indicator_cells_ct * pop_cells

        # correct for imprecision of Iinv*I
        pop_ct = nuts3.loc[nuts3.country == ct, "pop"].sum()
        pop_cells_ct *= pop_ct / pop_cells_ct.sum()

        # The first low density grid cells to reach rural fraction are rural
        asc_density_i = density_cells_ct.sort_values().index
        asc_density_cumsum = pop_cells_ct[asc_density_i].cumsum() / pop_cells_ct.sum()
        rural_fraction_ct = 1 - urban_fraction[ct]
        pop_ct_rural_b = asc_density_cumsum < rural_fraction_ct
        pop_ct_urban_b = ~pop_ct_rural_b

        pop_ct_rural_b[indicator_cells_ct == 0.0] = False
        pop_ct_urban_b[indicator_cells_ct == 0.0] = False

        pop_rural += pop_cells_ct.where(pop_ct_rural_b, 0.0)
        pop_urban += pop_cells_ct.where(pop_ct_urban_b, 0.0)

    pop_cells = {"total": pop_cells}
    pop_cells["rural"] = pop_rural
    pop_cells["urban"] = pop_urban

    for key, pop in pop_cells.items():

        ycoords = ("y", cutout.coords["y"].data)
        xcoords = ("x", cutout.coords["x"].data)
        values = pop.values.reshape(cutout.shape)
        layout = xr.DataArray(values, [ycoords, xcoords])

        layout.to_netcdf(snakemake.output[f"pop_layout_{key}"])
