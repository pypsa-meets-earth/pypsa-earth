# -*- coding: utf-8 -*-

import os
from itertools import dropwhile
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from helpers import mock_snakemake, override_component_attrs, sets_path_to_root


def override_values(tech, year, dr):
    buses = list(n.buses[n.buses.carrier == "AC"].index)

    custom_res_t = pd.read_csv(
        snakemake.input["custom_res_pot_{0}_{1}_{2}".format(tech, year, dr)],
        index_col=0,
        parse_dates=True,
    ).filter(buses, axis=1)

    custom_res = (
        pd.read_csv(
            snakemake.input["custom_res_ins_{0}_{1}_{2}".format(tech, year, dr)],
            index_col=0,
        )
        .filter(buses, axis=0)
        .reset_index()
    )

    custom_res["Generator"] = custom_res["Generator"].apply(lambda x: x + " " + tech)
    custom_res = custom_res.set_index("Generator")

    if tech.replace("-", " ") in n.generators.carrier.unique():

        to_drop = n.generators[n.generators.carrier == tech].index
        n.mremove("Generator", to_drop)

    print(n.generators)
    print(n.generators)
    print(n.generators)

    n.madd(
        "Generator",
        buses,
        " " + tech,
        bus=buses,
        carrier=tech,
        p_nom_extendable=True,
        p_nom_max=custom_res["p_nom_max"].values,
        # weight=ds["weight"].to_pandas(),
        # marginal_cost=custom_res["fixedomEuroPKW"].values * 1000,
        capital_cost=custom_res["annualcostEuroPMW"].values,
        efficiency=1.0,
        p_max_pu=custom_res_t,
        lifetime=custom_res["lifetime"][0],
        p_nom_min=custom_res["installedcapacity"].values,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "override_respot",
            simpl="",
            clusters="931",
            ll="c1.0",
            opts="Co2L",
            planning_horizons="2030",
            sopts="73H",
            discountrate=0.069,
        )
        sets_path_to_root("pypsa-earth-sec")

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    m = n.copy()
    if snakemake.config["custom_data"]["renewables"]:
        techs = snakemake.config["custom_data"]["renewables"]
        year = snakemake.wildcards["planning_horizons"]
        dr = snakemake.wildcards["discountrate"]

        m = n.copy()

        for tech in techs:
            override_values(tech, year, dr)

        n.export_to_netcdf(snakemake.output[0])
    else:
        print("No RES potential techs to override...")
        n.export_to_netcdf(snakemake.output[0])
