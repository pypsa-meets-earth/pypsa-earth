# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Cost helpers and the ``assign_costs`` Snakemake rule.

This module is the single source of truth for electricity cost assignment.
It owns:

- :func:`attach_dc_costs` — capital cost formula for DC lines/links.
- :func:`update_transmission_costs` — applies AC/DC capital costs to a network.
- :func:`calculate_renewable_capital_cost` — per-bus capital cost for a
  renewable carrier, including distance-dependent offshore connection costs.
- :func:`update_electricity_costs` — re-applies all cost-vintage-dependent
  attributes (capital_cost, marginal_cost, efficiency, lifetime) to an
  already-built electricity network without touching its structure.

The electricity build pipeline is split into two phases:

1. **Structure** (built once per scenario, shared across all horizons):
   ``add_electricity`` → ``simplify_network`` → ``cluster_network`` →
   ``add_extra_components`` → ``elec_s{simpl}_{clusters}_ec.nc``

2. **Costs** (assigned once per planning horizon, cheap):
   ``assign_costs`` → ``elec_s{simpl}_{clusters}_ec_{planning_horizons}.nc``

The ``__main__`` block at the bottom of this file implements phase 2.
"""

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


# ---------------------------------------------------------------------------
# Transmission helpers
# ---------------------------------------------------------------------------


def attach_dc_costs(lines_or_links, costs, length_factor=1.0, simple_hvdc_costs=False):
    """Apply capital costs to DC lines or links in-place."""
    if lines_or_links.empty:
        return

    if lines_or_links.loc[lines_or_links.carrier == "DC"].empty:
        return

    dc_b = lines_or_links.carrier == "DC"
    if simple_hvdc_costs:
        cap = (
            lines_or_links.loc[dc_b, "length"]
            * length_factor
            * costs.at["HVDC overhead", "capital_cost"]
        )
    else:
        cap = (
            lines_or_links.loc[dc_b, "length"]
            * length_factor
            * (
                (1.0 - lines_or_links.loc[dc_b, "underwater_fraction"])
                * costs.at["HVDC overhead", "capital_cost"]
                + lines_or_links.loc[dc_b, "underwater_fraction"]
                * costs.at["HVDC submarine", "capital_cost"]
            )
            + costs.at["HVDC inverter pair", "capital_cost"]
        )
    lines_or_links.loc[dc_b, "capital_cost"] = cap


def update_transmission_costs(n, costs, length_factor=1.0, simple_hvdc_costs=False):
    """Refresh AC line and DC link capital costs from the cost table."""
    n.lines["capital_cost"] = (
        n.lines["length"] * length_factor * costs.at["HVAC overhead", "capital_cost"]
    )

    attach_dc_costs(
        lines_or_links=n.links,
        costs=costs,
        length_factor=length_factor,
        simple_hvdc_costs=simple_hvdc_costs,
    )
    attach_dc_costs(
        lines_or_links=n.lines,
        costs=costs,
        length_factor=length_factor,
        simple_hvdc_costs=simple_hvdc_costs,
    )


# ---------------------------------------------------------------------------
# Renewable capital cost
# ---------------------------------------------------------------------------


def calculate_renewable_capital_cost(
    carrier: str,
    ds: xr.Dataset,
    costs: pd.DataFrame,
    line_length_factor: float,
    output_currency: str = "EUR",
):
    """
    Compute the per-bus capital cost of a renewable carrier from a cost table.

    For offshore wind the capital cost includes a distance-dependent grid
    connection cost (submarine/underground), which is why the renewable profile
    dataset ``ds`` (providing ``average_distance`` and ``underwater_fraction``)
    is required. For all other carriers the capital cost is the carrier-level
    value from ``costs``.

    This helper is the single source of truth for renewable capital costs and is
    used both at network build time (``attach_wind_and_solar``) and for
    per-horizon re-costing, so the two paths can never diverge.

    Parameters
    ----------
    carrier : str
        Renewable carrier name (e.g. ``"onwind"``, ``"offwind-ac"``).
    ds : xr.Dataset
        Renewable profile dataset for the carrier.
    costs : pd.DataFrame
        Technology cost table indexed by technology.
    line_length_factor : float
        Factor scaling connection distances for offshore wind.
    output_currency : str
        Currency label used only for logging.

    Returns
    -------
    capital_cost : float or pandas.Series
        Total carrier capital cost (technology + connection). A per-bus
        ``Series`` for offshore wind, a scalar otherwise.
    capital_cost_tech : float
        The technology-only component of the capital cost (a scalar per
        carrier), i.e. excluding the distance-dependent grid connection cost.
        For non-offshore carriers this equals ``capital_cost``.
    """
    supcarrier = carrier.split("-", 2)[0]

    if supcarrier == "offwind":
        capital_cost_tech = (
            costs.at["offwind", "capital_cost"]
            + costs.at[carrier + "-station", "capital_cost"]
        )
        underwater_fraction = ds["underwater_fraction"].to_pandas()
        connection_cost = (
            line_length_factor
            * ds["average_distance"].to_pandas()
            * (
                underwater_fraction
                * costs.at[carrier + "-connection-submarine", "capital_cost"]
                + (1.0 - underwater_fraction)
                * costs.at[carrier + "-connection-underground", "capital_cost"]
            )
        )
        capital_cost = capital_cost_tech + connection_cost

        logger.info(
            "Added connection cost of {:0.0f}-{:0.0f} {}/MW/a to {}".format(
                connection_cost.min(),
                connection_cost.max(),
                output_currency,
                carrier,
            )
        )
    else:
        capital_cost = costs.at[carrier, "capital_cost"]
        capital_cost_tech = capital_cost

    return capital_cost, capital_cost_tech


# ---------------------------------------------------------------------------
# Per-horizon re-costing
# ---------------------------------------------------------------------------


def update_electricity_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    renewable_carriers,
    length_factor: float = 1.0,
    hydro_capital_cost: bool = False,
) -> None:
    """
    Re-apply all cost-vintage-dependent attributes on an already-built
    electricity network, without rebuilding it.

    This is the single entry point for per-horizon re-costing. The structural
    network (topology, installed capacities, renewable profiles, build years,
    and the plant-specific efficiencies/lifetimes of *existing* units) is built
    once by ``add_electricity`` + ``add_extra_components``; this function then
    overwrites only the economic attributes that depend on the cost year, using
    a horizon-specific ``costs`` table.

    The formulas below mirror exactly those used at build time in
    ``add_electricity`` and ``add_extra_components``.  Attributes that do NOT
    depend on the cost vintage are deliberately left untouched:
      * ``p_nom``/``p_nom_min``/``p_nom_max`` and ``p_max_pu`` (capacities, profiles),
      * ``build_year`` and the per-plant ``efficiency``/``lifetime`` of existing
        (non-extendable) conventional generators,
      * transmission ``efficiency`` (set from ``sector.transmission_efficiency``).

    Conventional generators store a ``cost_key`` attribute at build time (set by
    ``attach_conventional_generators``).  The update path reads that tag instead
    of re-deriving the cost-table key from the carrier string, so the mapping
    between the network component and its cost row is explicit and
    change-tolerant.

    Parameters
    ----------
    n : pypsa.Network
        The already-built (clustered) network to re-cost in place.
    costs : pd.DataFrame
        Horizon-specific technology cost table indexed by technology.
    renewable_carriers : iterable of str
        Renewable generator carriers (used to route generator re-costing).
    length_factor : float
        Transmission/offshore connection length factor.
    hydro_capital_cost : bool
        Whether hydro reservoirs carry a capital cost.
    """
    renewable_carriers = set(renewable_carriers)

    def _assign(df, mask, **attrs):
        idx = df.index[mask] if not isinstance(mask, pd.Index) else mask
        if idx.empty:
            return
        for attr, value in attrs.items():
            df.loc[idx, attr] = value

    # --- Transmission lines and DC links ---
    update_transmission_costs(n, costs, length_factor=length_factor)

    gens = n.generators
    storage = n.storage_units
    stores = n.stores
    links = n.links

    # --- Generators ---
    for carrier in gens.carrier.dropna().unique():
        supcarrier = carrier.split("-", 2)[0]
        is_carrier = gens.carrier == carrier
        is_ext = is_carrier & gens.p_nom_extendable

        if carrier == "ror":
            _assign(
                n.generators,
                is_carrier,
                capital_cost=costs.at["ror", "capital_cost"],
                efficiency=costs.at["ror", "efficiency"],
            )

        elif carrier in renewable_carriers or supcarrier in renewable_carriers:
            # Wind/solar/csp: refresh all cost attributes from the horizon table.
            # For offshore wind, geometry (average_distance, underwater_fraction)
            # stored at build time lets us recompute the exact connection cost.
            if supcarrier == "offwind":
                lifetime = costs.at["offwind", "lifetime"]
            else:
                lifetime = costs.at[carrier, "lifetime"]

            _assign(
                n.generators,
                is_carrier,
                marginal_cost=costs.at[supcarrier, "marginal_cost"],
                efficiency=costs.at[supcarrier, "efficiency"],
                lifetime=lifetime,
            )

            sub = gens.index[is_carrier]
            has_geometry = (
                "average_distance" in n.generators.columns
                and "underwater_fraction" in n.generators.columns
                and n.generators.loc[sub, "average_distance"].notna().all()
            )

            if supcarrier == "offwind" and has_geometry:
                avg_dist = n.generators.loc[sub, "average_distance"]
                uw_frac = n.generators.loc[sub, "underwater_fraction"]
                connection_cost = (
                    length_factor
                    * avg_dist
                    * (
                        uw_frac
                        * costs.at[carrier + "-connection-submarine", "capital_cost"]
                        + (1.0 - uw_frac)
                        * costs.at[carrier + "-connection-underground", "capital_cost"]
                    )
                )
                n.generators.loc[sub, "capital_cost"] = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[carrier + "-station", "capital_cost"]
                    + connection_cost
                )
            elif supcarrier == "offwind":
                logger.warning(
                    f"No stored geometry for offshore carrier '{carrier}'; "
                    "capital_cost left unchanged during re-costing."
                )
            else:
                n.generators.loc[sub, "capital_cost"] = costs.at[
                    carrier, "capital_cost"
                ]

        elif carrier in costs.index:
            # Conventional generator. carrier is the cost-table row key.
            # marginal_cost is refreshed for all units; capital_cost, efficiency
            # and lifetime are only refreshed for extendable units (existing
            # plants keep their plant-specific values).
            _assign(
                n.generators,
                is_carrier,
                marginal_cost=costs.at[carrier, "marginal_cost"],
            )
            _assign(
                n.generators,
                is_ext,
                capital_cost=costs.at[carrier, "capital_cost"],
                efficiency=costs.at[carrier, "efficiency"],
                lifetime=costs.at[carrier, "lifetime"],
            )

    # --- Storage units ---
    if not storage.empty:
        _assign(
            n.storage_units,
            storage.carrier == "PHS",
            capital_cost=costs.at["PHS", "capital_cost"],
            efficiency_store=np.sqrt(costs.at["PHS", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["PHS", "efficiency"]),
        )

        _assign(
            n.storage_units,
            storage.carrier == "hydro",
            capital_cost=(
                costs.at["hydro", "capital_cost"] if hydro_capital_cost else 0.0
            ),
            marginal_cost=costs.at["hydro", "marginal_cost"],
            efficiency_dispatch=costs.at["hydro", "efficiency"],
        )

        is_batt = storage.carrier == "battery"
        is_batt_existing = is_batt & (~storage.p_nom_extendable)
        is_batt_ext = is_batt & storage.p_nom_extendable
        _assign(
            n.storage_units,
            is_batt_existing,
            capital_cost=costs.at["battery", "capital_cost"],
            marginal_cost=costs.at["battery", "marginal_cost"],
            efficiency_store=np.sqrt(costs.at["battery", "efficiency"]),
            efficiency_dispatch=np.sqrt(costs.at["battery", "efficiency"]),
        )
        _assign(
            n.storage_units,
            is_batt_ext,
            capital_cost=costs.at["battery", "capital_cost"],
            marginal_cost=costs.at["battery", "marginal_cost"],
            efficiency_store=costs.at["battery inverter", "efficiency"],
            efficiency_dispatch=costs.at["battery inverter", "efficiency"],
        )

        _assign(
            n.storage_units,
            storage.carrier == "H2",
            capital_cost=costs.at["H2", "capital_cost"],
            marginal_cost=costs.at["H2", "marginal_cost"],
            efficiency_store=costs.at["electrolysis", "efficiency"],
            efficiency_dispatch=costs.at["fuel cell", "efficiency"],
        )

    # --- Stores (Store-Link-Bus storage representation) ---
    if not stores.empty:
        _assign(
            n.stores,
            stores.carrier == "H2",
            capital_cost=costs.at["hydrogen storage tank", "capital_cost"],
        )
        _assign(
            n.stores,
            stores.carrier == "battery",
            capital_cost=costs.at["battery storage", "capital_cost"],
            marginal_cost=costs.at["battery", "marginal_cost"],
        )
        _assign(
            n.stores,
            stores.carrier == "csp",
            capital_cost=costs.at["csp-tower TES", "capital_cost"],
            marginal_cost=costs.at["csp-tower TES", "marginal_cost"],
        )

    # --- Extendable links (charging/discharging, pipelines, csp) ---
    if not links.empty:
        _assign(
            n.links,
            links.carrier == "H2 electrolysis",
            efficiency=costs.at["electrolysis", "efficiency"],
            capital_cost=costs.at["electrolysis", "capital_cost"],
            marginal_cost=costs.at["electrolysis", "marginal_cost"],
        )
        _assign(
            n.links,
            links.carrier == "H2 fuel cell",
            efficiency=costs.at["fuel cell", "efficiency"],
            capital_cost=costs.at["fuel cell", "capital_cost"]
            * costs.at["fuel cell", "efficiency"],
            marginal_cost=costs.at["fuel cell", "marginal_cost"],
        )
        _assign(
            n.links,
            links.carrier == "battery charger",
            efficiency=costs.at["battery inverter", "efficiency"],
            capital_cost=costs.at["battery inverter", "capital_cost"],
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )
        _assign(
            n.links,
            links.carrier == "battery discharger",
            efficiency=costs.at["battery inverter", "efficiency"],
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        )
        _assign(
            n.links,
            links.carrier == "csp",
            efficiency=costs.at["csp-tower", "efficiency"],
            capital_cost=costs.at["csp-tower", "capital_cost"],
            marginal_cost=costs.at["csp-tower", "marginal_cost"],
        )
        h2_pipe = links.carrier.str.startswith("H2 pipeline", na=False)
        if h2_pipe.any():
            n.links.loc[h2_pipe, "capital_cost"] = (
                costs.at["H2 pipeline", "capital_cost"] * n.links.loc[h2_pipe, "length"]
            )


# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "assign_costs",
            simpl="",
            clusters="4",
            planning_horizons=2030,
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)

    update_electricity_costs(
        n,
        costs,
        renewable_carriers=snakemake.params.electricity["renewable_carriers"],
        length_factor=snakemake.params.length_factor,
        hydro_capital_cost=snakemake.params.hydro_capital_cost,
    )

    logger.info(
        "Assigned costs from planning horizon %s to network %s.",
        snakemake.wildcards.planning_horizons,
        snakemake.input.network,
    )

    n.export_to_netcdf(snakemake.output.network)
