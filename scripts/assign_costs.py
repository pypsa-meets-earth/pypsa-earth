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

from typing import Any

import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import configure_logging, create_logger

logger = create_logger(__name__)


def attach_dc_costs(
    lines_or_links: pd.DataFrame,
    costs: pd.DataFrame,
    length_factor: float = 1.0,
    simple_hvdc_costs: bool = False,
) -> None:
    """
    Apply capital costs to DC lines or links in-place.

    Parameters
    ----------
    lines_or_links : pd.DataFrame
        DataFrame containing the lines or links to attach costs to.
    costs : pd.DataFrame
        DataFrame containing the costs to attach.
    length_factor : float
        Factor scaling connection distances for offshore wind.
    simple_hvdc_costs : bool
        Whether to use simple HVDC costs.

    Returns
    -------
    None
    """
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


def update_transmission_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    length_factor: float = 1.0,
    simple_hvdc_costs: bool = False,
) -> None:
    """
    Refresh AC line and DC link capital costs from the cost table.

    Parameters
    ----------
    n : pypsa.Network
        The network to update the transmission costs for.
    costs : pd.DataFrame
        The costs to update the transmission costs for.
    length_factor : float
        Factor scaling connection distances for offshore wind.
    simple_hvdc_costs : bool
        Whether to use simple HVDC costs.

    Returns
    -------
    None
    """
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


def _offwind_connection_cost(
    carrier: str,
    avg_dist: "pd.Series",
    uw_frac: "pd.Series",
    costs: pd.DataFrame,
    length_factor: float,
) -> pd.Series:
    """
    Compute the per-bus grid connection capital cost for an offshore wind carrier.

    This is the distance- and substrate-dependent part of the offshore wind
    capital cost (submarine + underground cable), separate from the turbine and
    substation costs. It is the single formula used both at network build time
    (from the renewable profile dataset) and during per-horizon re-costing
    (from attributes stored on generators).

    Parameters
    ----------
    carrier : str
        Offshore carrier name, e.g. ``"offwind-ac"`` or ``"offwind-dc"``.
    avg_dist : pd.Series
        Average connection distance per bus [km].
    uw_frac : pd.Series
        Fraction of the connection that is submarine (0–1) per bus.
    costs : pd.DataFrame
        Technology cost table indexed by technology.
    length_factor : float
        Scalar multiplier applied to distances.

    Returns
    -------
    pd.Series
        Per-bus connection capital cost [currency/MW/a].
    """
    return (
        length_factor
        * avg_dist
        * (
            uw_frac * costs.at[carrier + "-connection-submarine", "capital_cost"]
            + (1.0 - uw_frac)
            * costs.at[carrier + "-connection-underground", "capital_cost"]
        )
    )


def calculate_renewable_capital_cost(
    carrier: str,
    ds: xr.Dataset,
    costs: pd.DataFrame,
    line_length_factor: float,
    output_currency: str = "EUR",
) -> tuple:
    """
    Compute the per-bus capital cost of a renewable carrier from a cost table.

    For offshore wind the capital cost includes a distance-dependent grid
    connection cost (submarine/underground cable), computed via
    :func:`_offwind_connection_cost`. For all other carriers the capital cost
    is a scalar taken directly from the cost table.

    Used at network build time by ``attach_wind_and_solar``. The geometry
    (``average_distance``, ``underwater_fraction``) is stored on generators
    so that :func:`update_generator_costs` can recompute the same value per
    horizon without the profile dataset.

    Parameters
    ----------
    carrier : str
        Renewable carrier name, e.g. ``"onwind"``, ``"offwind-ac"``.
    ds : xr.Dataset
        Renewable profile dataset providing ``average_distance`` and
        ``underwater_fraction`` variables (buses as coordinate).
    costs : pd.DataFrame
        Technology cost table indexed by technology.
    line_length_factor : float
        Factor scaling connection distances for offshore wind.
    output_currency : str
        Currency label used only for logging.

    Returns
    -------
    capital_cost : float or pd.Series
        Total capital cost per bus (technology + connection for offshore).
    capital_cost_tech : float
        Technology-only component (turbine + substation), without connection.
        Equals ``capital_cost`` for non-offshore carriers.
    """
    supcarrier = carrier.split("-", 2)[0]

    if supcarrier == "offwind":
        capital_cost_tech = (
            costs.at["offwind", "capital_cost"]
            + costs.at[carrier + "-station", "capital_cost"]
        )
        conn_cost = _offwind_connection_cost(
            carrier,
            avg_dist=ds["average_distance"].to_pandas(),
            uw_frac=ds["underwater_fraction"].to_pandas(),
            costs=costs,
            length_factor=line_length_factor,
        )
        capital_cost = capital_cost_tech + conn_cost

        logger.info(
            "Added connection cost of {:0.0f}-{:0.0f} {}/MW/a to {}".format(
                conn_cost.min(),
                conn_cost.max(),
                output_currency,
                carrier,
            )
        )
    else:
        capital_cost = costs.at[carrier, "capital_cost"]
        capital_cost_tech = capital_cost

    return capital_cost, capital_cost_tech


def update_generator_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    renewable_carriers: set,
    length_factor: float,
) -> None:
    """
    Update cost-vintage-dependent attributes on all generators in-place.

    Two rules govern what gets updated:

    - **marginal_cost** is refreshed for **all** generators whose carrier is
      in the cost table.  Operating costs (fuel, VOM) change per horizon for
      both existing and new-build units.
    - **capital_cost, efficiency, lifetime** are refreshed **only for
      extendable** (new-build) generators.  For non-extendable (existing)
      plants these attributes are either irrelevant to the optimisation
      (capital_cost) or were set plant-specifically at build time and must
      not be overwritten (efficiency, lifetime).

    Offshore wind capital_cost is recomputed from ``average_distance`` and
    ``underwater_fraction`` stored on generators at build time, via
    :func:`_offwind_connection_cost`, keeping the formula consistent with
    :func:`calculate_renewable_capital_cost`.

    Carriers absent from the cost table are silently skipped.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    renewable_carriers : set of str
    length_factor : float
    """
    for carrier in n.generators.carrier.unique():
        supcarrier = carrier.split("-", 2)[0]
        cost_key = "offwind" if supcarrier == "offwind" else carrier

        if cost_key not in costs.index:
            continue

        is_carrier = n.generators.carrier == carrier
        is_ext = is_carrier & n.generators.p_nom_extendable

        # marginal_cost: refresh for all units — operating costs change per horizon
        _assign_component_attrs(
            n.generators,
            is_carrier,
            marginal_cost=costs.at[cost_key, "marginal_cost"],
        )

        # capital_cost / efficiency / lifetime: extendable units only
        if not is_ext.any():
            continue

        if carrier in renewable_carriers or supcarrier in renewable_carriers:
            _assign_component_attrs(
                n.generators,
                is_ext,
                efficiency=costs.at[cost_key, "efficiency"],
                lifetime=costs.at[cost_key, "lifetime"],
            )

            sub = n.generators.index[is_ext]
            has_geometry = (
                "average_distance" in n.generators.columns
                and "underwater_fraction" in n.generators.columns
                and n.generators.loc[sub, "average_distance"].notna().all()
            )

            if supcarrier == "offwind" and has_geometry:
                conn_cost = _offwind_connection_cost(
                    carrier,
                    avg_dist=n.generators.loc[sub, "average_distance"],
                    uw_frac=n.generators.loc[sub, "underwater_fraction"],
                    costs=costs,
                    length_factor=length_factor,
                )
                n.generators.loc[sub, "capital_cost"] = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[carrier + "-station", "capital_cost"]
                    + conn_cost
                )
            elif supcarrier == "offwind":
                logger.warning(
                    f"No stored geometry for offshore carrier '{carrier}'; "
                    "capital_cost left unchanged during re-costing."
                )
            else:
                _assign_component_attrs(
                    n.generators,
                    is_ext,
                    capital_cost=costs.at[carrier, "capital_cost"],
                )
        else:
            _assign_component_attrs(
                n.generators,
                is_ext,
                capital_cost=costs.at[carrier, "capital_cost"],
                efficiency=costs.at[carrier, "efficiency"],
                lifetime=costs.at[carrier, "lifetime"],
            )


def update_storage_unit_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    hydro_capital_cost: bool = False,
) -> None:
    """
    Update cost-vintage-dependent attributes on storage units in-place.

    Same rules as :func:`update_generator_costs`:

    - **marginal_cost** is refreshed for all units of a carrier.
    - **capital_cost** and **efficiency** are refreshed only for extendable
      (new-build) units. Existing PHS, hydro, and battery assets keep their
      build-time values.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    hydro_capital_cost : bool
        Whether hydro reservoirs carry a capital cost (mirrors
        ``renewable: hydro: hydro_capital_cost``).
    """
    # PHS — typically non-extendable; capital_cost irrelevant for existing units
    is_phs = n.storage_units.carrier == "PHS"
    is_phs_ext = is_phs & n.storage_units.p_nom_extendable
    _assign_component_attrs(
        n.storage_units,
        is_phs_ext,
        capital_cost=costs.at["PHS", "capital_cost"],
        efficiency_store=np.sqrt(costs.at["PHS", "efficiency"]),
        efficiency_dispatch=np.sqrt(costs.at["PHS", "efficiency"]),
    )

    # Hydro reservoir — marginal_cost for all; capital_cost/efficiency extendable only
    is_hydro = n.storage_units.carrier == "hydro"
    is_hydro_ext = is_hydro & n.storage_units.p_nom_extendable
    _assign_component_attrs(
        n.storage_units,
        is_hydro,
        marginal_cost=costs.at["hydro", "marginal_cost"],
    )
    _assign_component_attrs(
        n.storage_units,
        is_hydro_ext,
        capital_cost=(costs.at["hydro", "capital_cost"] if hydro_capital_cost else 0.0),
        efficiency_dispatch=costs.at["hydro", "efficiency"],
    )

    # Battery — marginal_cost for all; full update for extendable only
    is_batt = n.storage_units.carrier == "battery"
    is_batt_ext = is_batt & n.storage_units.p_nom_extendable
    _assign_component_attrs(
        n.storage_units,
        is_batt,
        marginal_cost=costs.at["battery", "marginal_cost"],
    )
    _assign_component_attrs(
        n.storage_units,
        is_batt_ext,
        capital_cost=costs.at["battery", "capital_cost"],
        efficiency_store=costs.at["battery inverter", "efficiency"],
        efficiency_dispatch=costs.at["battery inverter", "efficiency"],
    )

    # H2 StorageUnit — only present when config sets
    # extendable_carriers.StorageUnit: [H2]. When H2 is modelled as
    # Store+Link instead, this carrier appears in n.stores and is handled by
    # update_store_costs. The block below is skipped if H2 is absent.
    is_h2 = n.storage_units.carrier == "H2"
    is_h2_ext = is_h2 & n.storage_units.p_nom_extendable
    _assign_component_attrs(
        n.storage_units,
        is_h2,
        marginal_cost=costs.at["H2", "marginal_cost"],
    )
    _assign_component_attrs(
        n.storage_units,
        is_h2_ext,
        capital_cost=costs.at["H2", "capital_cost"],
        efficiency_store=costs.at["electrolysis", "efficiency"],
        efficiency_dispatch=costs.at["fuel cell", "efficiency"],
    )


def update_store_costs(n: pypsa.Network, costs: pd.DataFrame) -> None:
    """
    Update cost-vintage-dependent attributes on stores in-place.

    Stores represent energy capacity (MWh) in the Store-Link-Bus model (H2,
    battery, CSP thermal storage). They use ``e_nom_extendable`` rather than
    ``p_nom_extendable``.

    - **marginal_cost** is refreshed for all stores of a carrier.
    - **capital_cost** is refreshed only for extendable stores (new-build
      energy capacity). In the electricity workflow stores are added by
      ``add_extra_components`` and are typically all extendable.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    """
    is_h2 = n.stores.carrier == "H2"
    is_h2_ext = is_h2 & n.stores.e_nom_extendable
    _assign_component_attrs(
        n.stores,
        is_h2_ext,
        capital_cost=costs.at["hydrogen storage tank", "capital_cost"],
    )

    is_batt = n.stores.carrier == "battery"
    is_batt_ext = is_batt & n.stores.e_nom_extendable
    _assign_component_attrs(
        n.stores,
        is_batt,
        marginal_cost=costs.at["battery", "marginal_cost"],
    )
    _assign_component_attrs(
        n.stores,
        is_batt_ext,
        capital_cost=costs.at["battery storage", "capital_cost"],
    )

    is_csp = n.stores.carrier == "csp"
    is_csp_ext = is_csp & n.stores.e_nom_extendable
    _assign_component_attrs(
        n.stores,
        is_csp,
        marginal_cost=costs.at["csp-tower TES", "marginal_cost"],
    )
    _assign_component_attrs(
        n.stores,
        is_csp_ext,
        capital_cost=costs.at["csp-tower TES", "capital_cost"],
    )


def update_link_costs(n: pypsa.Network, costs: pd.DataFrame) -> None:
    """
    Update cost-vintage-dependent attributes on links in-place.

    Links in the electricity network are charging/discharging converters (H2,
    battery), CSP discharge links, and optionally H2 pipelines. They use
    ``p_nom_extendable``.

    - **marginal_cost** is refreshed for all links of a carrier (where defined).
    - **capital_cost** and **efficiency** are refreshed only for extendable
      links. Transmission efficiency set from config (not the cost table) is
      left untouched.

    H2 pipeline ``capital_cost`` scales with link length and is only updated
    for extendable pipeline links.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    """
    carrier_updates = {
        "H2 electrolysis": dict(
            efficiency=costs.at["electrolysis", "efficiency"],
            capital_cost=costs.at["electrolysis", "capital_cost"],
            marginal_cost=costs.at["electrolysis", "marginal_cost"],
        ),
        "H2 fuel cell": dict(
            efficiency=costs.at["fuel cell", "efficiency"],
            capital_cost=costs.at["fuel cell", "capital_cost"]
            * costs.at["fuel cell", "efficiency"],
            marginal_cost=costs.at["fuel cell", "marginal_cost"],
        ),
        "battery charger": dict(
            efficiency=costs.at["battery inverter", "efficiency"],
            capital_cost=costs.at["battery inverter", "capital_cost"],
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        ),
        "battery discharger": dict(
            efficiency=costs.at["battery inverter", "efficiency"],
            marginal_cost=costs.at["battery inverter", "marginal_cost"],
        ),
        "csp": dict(
            efficiency=costs.at["csp-tower", "efficiency"],
            capital_cost=costs.at["csp-tower", "capital_cost"],
            marginal_cost=costs.at["csp-tower", "marginal_cost"],
        ),
    }

    for carrier, attrs in carrier_updates.items():
        is_carrier = n.links.carrier == carrier
        is_ext = is_carrier & n.links.p_nom_extendable

        if "marginal_cost" in attrs:
            _assign_component_attrs(
                n.links, is_carrier, marginal_cost=attrs["marginal_cost"]
            )

        ext_attrs = {k: v for k, v in attrs.items() if k != "marginal_cost"}
        if ext_attrs:
            _assign_component_attrs(n.links, is_ext, **ext_attrs)

    h2_pipe = n.links.carrier.str.startswith("H2 pipeline", na=False)
    h2_pipe_ext = h2_pipe & n.links.p_nom_extendable
    if h2_pipe_ext.any():
        n.links.loc[h2_pipe_ext, "capital_cost"] = (
            costs.at["H2 pipeline", "capital_cost"] * n.links.loc[h2_pipe_ext, "length"]
        )


def _assign_component_attrs(
    df: pd.DataFrame,
    mask: pd.Series | pd.Index,
    **attrs: Any,
) -> None:
    """
    Assign one or more column values to rows of a PyPSA component table.

    The update is skipped when the selected rows are empty, so callers can
    safely pass carrier masks that may match no components.

    Parameters
    ----------
    df : pandas.DataFrame
        A PyPSA component table such as ``n.generators`` or ``n.links``.
    mask : pandas.Series or pandas.Index
        Boolean mask selecting rows, or an index of row labels.
    **attrs
        Column names and values to assign, e.g.
        ``capital_cost=1000, marginal_cost=5``.
    """
    idx = df.index[mask] if not isinstance(mask, pd.Index) else mask
    if idx.empty:
        return
    for attr, value in attrs.items():
        df.loc[idx, attr] = value


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
    a horizon-specific ``costs`` table. Attributes that do NOT
    depend on the cost vintage are deliberately left untouched.

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

    Returns
    -------
    None
    """
    renewable_carriers = set(renewable_carriers)

    # Updating transmission costs
    update_transmission_costs(n, costs, length_factor=length_factor)

    # Updating generator costs
    update_generator_costs(
        n,
        costs,
        renewable_carriers,
        length_factor,
    )

    # Updating storage unit costs
    update_storage_unit_costs(n, costs, hydro_capital_cost=hydro_capital_cost)

    # Updating store costs (Store-Link-Bus energy capacity)
    update_store_costs(n, costs)

    # Updating link costs (chargers, dischargers, pipelines)
    update_link_costs(n, costs)


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
