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

    if snakemake.wildcards["planning_horizons"] == 2050:
        directory = "results/" + snakemake.config.run.replace("2050", "2030")
        n_name = snakemake.input.network.split("/")[-1].replace(
            n.config["scenario"]["clusters"], ""
        )
        df = pd.read_csv(directory + "/res_caps_" + n_name, index_col=0)
        # df = pd.read_csv(snakemake.config["custom_data"]["existing_renewables"], index_col=0)
        existing_res = df.loc[tech]
        existing_res.index = existing_res.index.str.apply(lambda x: x + tech)
    else:
        existing_res = custom_res["installedcapacity"].values

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
        p_nom_min=existing_res,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        snakemake = mock_snakemake(
            "override_respot",
            simpl="",
            clusters="171",
            ll="c1.0",
            opts="Co2L",
            planning_horizons="2030",
            sopts="3H",
            demand="AP",
            discountrate=0.071,
        )
        sets_path_to_root("pypsa-earth-sec")

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)
    m = n.copy()
    buses = list(n.buses[n.buses.carrier == "AC"].index)
    energy_totals = pd.read_csv(snakemake.input.energy_totals, index_col=0)
    countries = snakemake.config["countries"]
    if snakemake.config["custom_data"]["renewables"]:
        techs = snakemake.config["custom_data"]["renewables"]
        year = snakemake.wildcards["planning_horizons"]
        dr = snakemake.wildcards["discountrate"]

        m = n.copy()

        for tech in techs:
            override_values(tech, year, dr)

    else:
        print("No RES potential techs to override...")

    if (
        snakemake.config["custom_data"]["add_existing"]
        and snakemake.wildcards["planning_horizons"] == "2030"
    ):
        n.generators.loc["MAR010003 onwind", "p_nom_min"] = 50
        n.generators.loc["MAR010005 offwind", "p_nom_min"] = 120
        n.generators.loc["MAR010005 offwind", "p_nom_max"] = 120
        n.generators.loc["MAR010001 onwind", "p_nom_min"] = 154
        n.generators.loc["MAR004005 onwind", "p_nom_min"] = 150
        n.generators.loc["MAR003003 onwind", "p_nom_min"] = 180
        n.generators.loc["MAR006003 onwind", "p_nom_min"] = 60
        n.generators.loc["MAR011001 onwind", "p_nom_min"] = 208

        n.generators.loc["MAR003001 csp", "p_nom_min"] = 500

        n.generators.loc["MAR003001 solar", "p_nom_min"] = 72

    if snakemake.config["custom_data"]["elec_demand"]:
        for country in countries:
            n.loads_t.p_set.filter(like=country)[buses] = (
                (
                    n.loads_t.p_set.filter(like=country)[buses]
                    / n.loads_t.p_set.filter(like=country)[buses].sum().sum()
                )
                * energy_totals.loc[country, "electricity residential"]
                * 1e6
            )

    n.export_to_netcdf(snakemake.output[0])
