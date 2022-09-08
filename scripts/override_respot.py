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

# def annuity_factor(v):
#     return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

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

        n.generators_t.p_max_pu.update(custom_res_t)
        n.generators.update(pd.Series(custom_res["p_nom_max"]))

    else:
        n.madd(
            "Generator",
            buses,
            " " + tech,
            bus=buses,
            carrier=tech,
            p_nom_extendable=True,
            p_nom_max=custom_res["p_nom_max"].values,
            # weight=ds["weight"].to_pandas(),
            marginal_cost=custom_res["fixedomEuroPKW"].values * 1000,
            capital_cost=custom_res["annualcostEuroPMW"].values * 1000,
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
            clusters="1004",
            ll="c1.0",
            opts="Co2L",
            planning_horizons="2030",
            sopts="3H",
            discountrate=0.071,
        )
        sets_path_to_root("pypsa-earth-sec")

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

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

# Nyears = n.snapshot_weightings.generators.sum() / 8760

# costs["fixed"] = [
#     annuity_factor(v) * v["investment"] * Nyears for i, v in costs.iterrows()
# ]
