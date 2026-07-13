# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Prepares brownfield data from previous planning horizon.

Relevant Settings
-----------------

```yaml

    sector:
        hydrogen:
            network:
            H2_retrofit_capacity_per_CH4:
            network_limit:
            network_routes:
            gas_network_repurposing:
            underground_storage:
            hydrogen_colors:
            set_color_shares:
            blue_share:
            pink_share:
            production_technologies:

    existing_capacities
        grouping_years_power:
        grouping_years_heat:
        threshold_capacity:
        default_heating_lifetime:
        conventional_carriers:

    snapshots:
        start:
        end:
        inclusive:

    electricity:
        renewable_carriers:
```
Inputs
------
- ``resources/{RDIR}/bus_regions/busmap_elec_s{simpl}.csv``: Busmap after simplifying the network
- ``resources/{RDIR}/bus_regions/busmap_elec_s{simpl}_{clusters}.csv``: Busmap after clustering the network
- ``{RESDIR}/prenetworks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc``: prenetwork file obtained prior to solving
- ``solved_previous_horizon``: Network solved at previous time step
- ``resources/{RDIR}/costs_{planning_horizons}_sec.csv``: Technology costs data
- ``resources/{SECDIR}/cops/cop_soil_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Ground/soil source heat pump COP time series aligned to the network snapshots
- ``resources/{SECDIR}/cops/cop_air_total_elec_s{simpl}_{clusters}_{planning_horizons}.nc``: Air source heat pump COP time series aligned to the network snapshots

Output
------
- ``{RESDIR}/prenetworks-brownfield/elec_s{simpl}_{clusters}_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}_{demand}_{h2export}export.nc``: Brownfield prenetwork file

Description
-----------
To prepare network for brownfield expansion

"""

import logging

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import read_csv_nafix, sanitize_carriers, sanitize_locations
from add_existing_baseyear import add_build_year_to_new_assets

# from pypsa.clustering.spatial import normed_or_uniform

logger = logging.getLogger(__name__)
idx = pd.IndexSlice


def add_brownfield(
    n: pypsa.Network,
    n_p: pypsa.Network,
    year: int,
    threshold_capacity: float,
    H2_retrofit=False,
    H2_retrofit_capacity_per_CH4=None,
) -> None:
    """
    Adds brownfield assets from the previous planning horizon to the network.

    Parameters
    ----------
    n : pypsa.Network
        The new PyPSA network to which brownfield assets will be added.
    n_p : pypsa.Network
        The previous PyPSA network from which brownfield assets will be sourced.
    year : int
        The planning horizon year for which brownfield assets are being prepared.
    threshold_capacity : float
        Minimum optimised capacity (MW) below which assets are dropped.
    H2_retrofit : bool or dict, optional
        H2 retrofit configuration (sector-coupled only). Falsy disables gas-network logic.
    H2_retrofit_capacity_per_CH4 : float, optional
        Capacity ratio for H2/CH4 retrofit (sector-coupled only).

    Returns
    -------
    None
    """
    logger.info(f"Preparing brownfield for the year {year}")

    # electric transmission grid set optimised capacities of previous as minimum
    n.lines.s_nom_min = n_p.lines.s_nom_opt
    dc_i = n.links[n.links.carrier == "DC"].index
    n.links.loc[dc_i, "p_nom_min"] = n_p.links.loc[dc_i, "p_nom_opt"]

    # Reset p_nom_min on extendable generators that was set from IRENA stats
    # in add_electricity to prevent double-counting.
    extendable_gens = n.generators.index[n.generators.p_nom_extendable]
    if not extendable_gens.empty:
        n.generators.loc[extendable_gens, "p_nom_min"] = 0.0
        n.generators.loc[extendable_gens, "p_nom"] = 0.0

    for c in n_p.iterate_components(["Link", "Generator", "Store"]):
        attr = "e" if c.name == "Store" else "p"

        # Remove generators, links and stores that track global values since they exist in n
        n_p.mremove(c.name, c.df.index[c.df.lifetime == np.inf])

        # Remove assets whose build_year + lifetime < year
        n_p.mremove(c.name, c.df.index[c.df.build_year + c.df.lifetime < year])

        # Remove assets that are not extendable since they already exist in n
        n_p.mremove(c.name, c.df.index[~c.df[f"{attr}_nom_extendable"]])

        # Remove assets if their optimized nominal capacity is lower than a threshold
        n_p.mremove(
            c.name,
            c.df.index[
                (
                    c.df[f"{attr}_nom_extendable"]
                    & (c.df[f"{attr}_nom_opt"] < threshold_capacity)
                )
            ],
        )

        # Copy optimized assets from previous horizon to current and fix their capacity
        c.df[f"{attr}_nom"] = c.df[f"{attr}_nom_opt"]
        c.df[f"{attr}_nom_extendable"] = False

        # Remove assets if name already exist in the new network
        n_p.mremove(c.name, c.df.index.intersection(getattr(n, c.list_name).index))

        n.import_components_from_dataframe(c.df, c.name)

        # Copy time-dependent parameters of the optimized assets from previous horizon to current
        selection = n.component_attrs[c.name].type.str.contains(
            "series"
        ) & n.component_attrs[c.name].status.str.contains("Input")
        for tattr in n.component_attrs[c.name].index[selection]:
            n.import_series_from_dataframe(c.pnl[tattr], c.name, tattr)

        # deal with gas network
        pipe_carrier = ["gas pipeline"]
        if H2_retrofit:
            # drop capacities of previous year to avoid duplicating
            to_drop = n.links.carrier.isin(pipe_carrier) & (n.links.build_year != year)
            n.mremove("Link", n.links.loc[to_drop].index)

            # subtract the already retrofitted from today's gas grid capacity
            h2_retrofitted_fixed_i = n.links[
                (n.links.carrier == "H2 pipeline retrofitted")
                & (n.links.build_year != year)
            ].index
            gas_pipes_i = n.links[n.links.carrier.isin(pipe_carrier)].index
            CH4_per_H2 = 1 / H2_retrofit_capacity_per_CH4
            fr = "H2 pipeline retrofitted"
            to = "gas pipeline"
            # today's pipe capacity
            pipe_capacity = n.links.loc[gas_pipes_i, "p_nom"]
            # already retrofitted capacity from gas -> H2
            already_retrofitted = (
                n.links.loc[h2_retrofitted_fixed_i, "p_nom"]
                .rename(lambda x: x.split("-2")[0].replace(fr, to))
                .groupby(level=0)
                .sum()
            )
            remaining_capacity = (
                pipe_capacity
                - CH4_per_H2
                * already_retrofitted.reindex(index=pipe_capacity.index).fillna(0)
            )
            n.links.loc[gas_pipes_i, "p_nom"] = remaining_capacity
        else:
            new_pipes = n.links.carrier.isin(pipe_carrier) & (
                n.links.build_year == year
            )
            n.links.loc[new_pipes, "p_nom"] = 0.0
            n.links.loc[new_pipes, "p_nom_min"] = 0.0


def disable_grid_expansion_if_limit_hit(n: pypsa.Network) -> None:
    """
    Check if transmission expansion limit is already reached; then turn off.

    In particular, this function checks if the total transmission
    capital cost or volume implied by s_nom_min and p_nom_min are
    numerically close to the respective global limit set in
    n.global_constraints. If so, the nominal capacities are set to the
    minimum and extendable is turned off; the corresponding global
    constraint is then dropped.

    Parameters
    ----------
    n : pypsa.Network
        The PyPSA network to check and adjust for transmission expansion limits.

    Returns
    -------
    None
    """
    cols = {"cost": "capital_cost", "volume": "length"}
    for limit_type in ["cost", "volume"]:
        glcs = n.global_constraints.query(
            f"type == 'transmission_expansion_{limit_type}_limit'"
        )

        for name, glc in glcs.iterrows():
            total_expansion = (
                (
                    n.lines.query("s_nom_extendable")
                    .eval(f"s_nom_min * {cols[limit_type]}")
                    .sum()
                )
                + (
                    n.links.query("carrier == 'DC' and p_nom_extendable")
                    .eval(f"p_nom_min * {cols[limit_type]}")
                    .sum()
                )
            ).sum()

            # Allow small numerical differences
            if np.abs(glc.constant - total_expansion) / glc.constant < 1e-6:
                logger.info(
                    f"Transmission expansion {limit_type} is already reached, disabling expansion and limit"
                )
                extendable_acs = n.lines.query("s_nom_extendable").index
                n.lines.loc[extendable_acs, "s_nom_extendable"] = False
                n.lines.loc[extendable_acs, "s_nom"] = n.lines.loc[
                    extendable_acs, "s_nom_min"
                ]

                extendable_dcs = n.links.query(
                    "carrier == 'DC' and p_nom_extendable"
                ).index
                n.links.loc[extendable_dcs, "p_nom_extendable"] = False
                n.links.loc[extendable_dcs, "p_nom"] = n.links.loc[
                    extendable_dcs, "p_nom_min"
                ]

                n.global_constraints.drop(name, inplace=True)


# def adjust_renewable_profiles(n, input_profiles, params, year):
#     """
#     Adjusts renewable profiles according to the renewable technology specified,
#     using the latest year below or equal to the selected year.
#     """

#     # spatial clustering
#     cluster_busmap = read_csv_nafix(snakemake.input.cluster_busmap, index_col=0).squeeze()
#     simplify_busmap = read_csv_nafix(
#         snakemake.input.simplify_busmap, index_col=0
#     ).squeeze()
#     clustermaps = simplify_busmap.map(cluster_busmap)
#     clustermaps.index = clustermaps.index.astype(str)

#     # temporal clustering
#     dr = pd.date_range(**params["snapshots"], freq="h")
#     snapshotmaps = (
#         pd.Series(dr, index=dr).where(lambda x: x.isin(n.snapshots), pd.NA).ffill()
#     )

#     for carrier in params["carriers"]:
#         if carrier == "hydro":
#             continue
#         with xr.open_dataset(getattr(input_profiles, "profile_" + carrier)) as ds:
#             if ds.indexes["bus"].empty or "year" not in ds.indexes:
#                 continue

#             closest_year = max(
#                 (y for y in ds.year.values if y <= year), default=min(ds.year.values)
#             )

#             p_max_pu = (
#                 ds["profile"]
#                 .sel(year=closest_year)
#                 .transpose("time", "bus")
#                 .to_pandas()
#             )

#             # spatial clustering
#             weight = ds["weight"].sel(year=closest_year).to_pandas()
#             weight = weight.groupby(clustermaps).transform(normed_or_uniform)
#             p_max_pu = (p_max_pu * weight).T.groupby(clustermaps).sum().T
#             p_max_pu.columns = p_max_pu.columns + f" {carrier}"

#             # temporal_clustering
#             p_max_pu = p_max_pu.groupby(snapshotmaps).mean()

#             # replace renewable time series
#             n.generators_t.p_max_pu.loc[:, p_max_pu.columns] = p_max_pu


if __name__ == "__main__":
    if "snakemake" not in globals():

        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "add_brownfield",
            simpl="",
            clusters="4",
            ll="c1",
            opts="Co2L-4H",
            planning_horizons="2030",
            sopts="144H",
            discountrate=0.071,
            demand="AB",
            h2export="120",
        )

    logger.info(f"Preparing brownfield from the file {snakemake.input.network_p}")

    year = int(snakemake.wildcards.planning_horizons)

    n = pypsa.Network(snakemake.input.network)

    # TODO
    # adjust_renewable_profiles(n, snakemake.input, snakemake.params, year)

    add_build_year_to_new_assets(n, year)

    n_p = pypsa.Network(snakemake.input.network_p)

    is_sector_coupled = "sopts" in snakemake.wildcards.keys()
    threshold_capacity = snakemake.params.threshold_capacity

    if is_sector_coupled:
        H2_retrofit = snakemake.params.H2_retrofit
        H2_retrofit_capacity_per_CH4 = snakemake.params.H2_retrofit_capacity_per_CH4
    else:
        H2_retrofit = False
        H2_retrofit_capacity_per_CH4 = None

    add_brownfield(
        n, n_p, year, threshold_capacity, H2_retrofit, H2_retrofit_capacity_per_CH4
    )

    disable_grid_expansion_if_limit_hit(n)

    sanitize_carriers(n, snakemake.config)
    sanitize_locations(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
