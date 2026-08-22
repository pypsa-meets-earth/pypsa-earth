# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Assign horizon-specific technology costs to an already-built electricity network.

The ``assign_costs`` Snakemake rule loads a clustered network that already
contains topology, capacities, and time series, and updates only the economic
attributes that depend on the planning horizon (``capital_cost``,
``marginal_cost``, ``efficiency``, ``lifetime``). Network structure and
existing-plant build years are left unchanged.

This module is also the shared implementation for cost assignment helpers used
elsewhere in the electricity workflow.

Relevant Settings
-----------------

.. code:: yaml

    electricity:
      renewable_carriers:
    lines:
      length_factor:
    renewable:
      hydro:
        hydro_capital_cost:
    storage_techs:

Inputs
------

- ``networks/elec_s{simpl}_{clusters}_ec.nc``: Clustered electricity network
  with components and profiles, but without horizon-specific costs.
- ``resources/costs_{planning_horizons}_elec.csv``: Technology cost table for
  the planning horizon.

Outputs
-------

- ``networks/elec_s{simpl}_{clusters}_ec_{planning_horizons}.nc``: Same network
  with costs assigned for the given planning horizon.
- ``resources/bus_regions/connection_costs_s{simpl}_{clusters}_{planning_horizons}.csv``:
  Horizon-specific displacement connection CAPEX for offshore generators moved
  during ``simplify_network`` (from stored path geometry).

Description
-----------

When multiple planning horizons are modelled, building the full network for
each year is expensive. The workflow therefore separates **structure**
(``add_electricity`` → ``simplify_network`` → ``cluster_network`` →
``add_extra_components``) from **cost assignment** (this rule). This
script reads the shared structure network and the horizon-specific cost
table, calls :func:`update_electricity_costs`, and writes the
horizon-tagged output network.

Main functions in this module:

- :func:`update_electricity_costs` — entry point; re-costs generators,
  storage units, stores, links, and transmission in place.
- :func:`update_generator_costs` — horizon costs for generators, including
  offshore farm and displacement connection CAPEX.
- :func:`update_storage_unit_costs` — horizon costs for StorageUnits
  (``storage_techs``, plus PHS/hydro special cases).
- :func:`update_store_costs` — horizon costs for Store energy capacity.
- :func:`update_link_costs` — horizon costs for charger/discharger links,
  CSP, and H2 pipelines.
- :func:`update_transmission_costs` — refresh AC/DC ``capital_cost`` from the
  cost table.
- :func:`calculate_renewable_capital_cost` — per-bus renewable capital cost
  at build time, including offshore connection costs.
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
        Factor scaling DC line/link lengths when computing capital costs.
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
        capital_costs = (
            lines_or_links.loc[dc_b, "length"]
            * length_factor
            * costs.at["HVDC overhead", "capital_cost"]
        )
    else:
        capital_costs = (
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
    lines_or_links.loc[dc_b, "capital_cost"] = capital_costs


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
        Factor scaling AC/DC line and link lengths when computing capital
        costs.
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
    output_currency: str = "EUR",
) -> pd.DataFrame:
    """
    Update horizon-specific cost attributes on generators in-place.

    - **marginal_cost** is refreshed for **all** generators whose carrier is
      in the cost table (operating costs change per horizon for existing and
      new-build units).
    - **capital_cost**, **efficiency**, and **lifetime** are refreshed **only
      for extendable** (new-build) generators. Non-extendable plants keep
      build-time efficiency and lifetime; their capital_cost is irrelevant to
      the optimisation.

    For offshore wind, extendable ``capital_cost`` is rebuilt as::

        offwind turbine + station + farm_connection + displacement_connection

    where farm and displacement connection costs come from
    :func:`_offwind_connection_cost` using ``average_distance`` /
    ``underwater_fraction`` and optional ``displacement_length`` from
    ``simplify_network``. If offshore geometry is missing, capital_cost is
    left unchanged and a warning is logged.

    Carriers absent from the cost table are skipped.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    renewable_carriers : set of str
        Renewable generator carriers (routes offshore vs standard re-costing).
    length_factor : float
        Scalar applied to offshore farm connection distances.
    output_currency : str
        Currency label used only for logging.

    Returns
    -------
    pd.DataFrame
        Displacement connection costs [currency/MW/a] indexed by generator,
        with one column per offshore carrier that had a positive surcharge.
        Empty if no displacement geometry is present.
    """
    displacement_connection_costs = {}

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
                farm_connection_cost = _offwind_connection_cost(
                    carrier,
                    avg_dist=n.generators.loc[sub, "average_distance"],
                    uw_frac=n.generators.loc[sub, "underwater_fraction"],
                    costs=costs,
                    length_factor=length_factor,
                )

                displacement_connection_cost = 0.0
                # displacement_length already includes length_factor from simplify_network
                if "displacement_length" in n.generators.columns:
                    displacement_length = n.generators.loc[
                        sub, "displacement_length"
                    ].fillna(0.0)
                    if displacement_length.gt(0).any():
                        displacement_connection_cost = _offwind_connection_cost(
                            carrier,
                            avg_dist=displacement_length,
                            uw_frac=n.generators.loc[sub, "underwater_fraction"],
                            costs=costs,
                            length_factor=1.0,
                        )
                        positive_displacement_costs = displacement_connection_cost[
                            displacement_connection_cost > 0
                        ]
                        if not positive_displacement_costs.empty:
                            displacement_connection_costs[carrier] = (
                                positive_displacement_costs
                            )
                            logger.info(
                                "Added displacement connection cost of {:.0f}-{:.0f} "
                                "{}/MW/a to {} generators of carrier {}".format(
                                    positive_displacement_costs.min(),
                                    positive_displacement_costs.max(),
                                    output_currency,
                                    len(positive_displacement_costs),
                                    carrier,
                                )
                            )
                total_connection_cost = (
                    farm_connection_cost + displacement_connection_cost
                )
                n.generators.loc[sub, "capital_cost"] = (
                    costs.at["offwind", "capital_cost"]
                    + costs.at[carrier + "-station", "capital_cost"]
                    + total_connection_cost
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

    return pd.DataFrame(displacement_connection_costs)


def update_storage_unit_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    hydro_capital_cost: bool = False,
    storage_techs: dict | None = None,
) -> None:
    """
    Update horizon-specific cost attributes on storage units in-place.

    Same rules as :func:`update_generator_costs`:

    - **marginal_cost** is refreshed for all units of a carrier.
    - **capital_cost**, charge/dispatch **efficiency**, and **lifetime** are
      refreshed only for extendable (new-build) units.

    Special cases outside ``storage_techs``:

    - **PHS** — efficiencies from the ``PHS`` cost row (sqrt split); capital
      cost only if extendable.
    - **hydro** — marginal cost for all; capital cost only if extendable and
      ``hydro_capital_cost`` is enabled.

    For carriers added via ``attach_storageunits`` (e.g. ``battery``, ``H2``),
    charger / discharger / bicharger names come from ``storage_techs``.
    Combined capital and marginal costs are keyed by carrier name in the
    processed cost table (see ``calculate_cost_for_storage_units``). Battery /
    li-ion efficiencies use the same round-trip correction as
    ``attach_storageunits``.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    hydro_capital_cost : bool
        Whether hydro reservoirs carry a capital cost (mirrors
        ``renewable: hydro: hydro_capital_cost``).
    storage_techs : dict, optional
        Mapping of storage carriers to ``store`` / ``charger`` /
        ``discharger`` / ``bicharger`` technology names. Same structure as
        ``config["storage_techs"]``.
    """
    storage_techs = storage_techs or {}

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

    # Configurable StorageUnit carriers from storage_techs (battery, H2, …)
    for carrier in sorted(set(n.storage_units.carrier) & set(storage_techs)):
        if carrier not in costs.index:
            logger.warning(
                f"No cost row for StorageUnit carrier '{carrier}'; skipping re-costing."
            )
            continue

        lookup = storage_techs[carrier]
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
        else:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]

        # Same round-trip split as attach_storageunits
        roundtrip_correction = 0.5 if carrier in ["battery", "li-ion"] else 1

        is_carrier = n.storage_units.carrier == carrier
        is_ext = is_carrier & n.storage_units.p_nom_extendable

        _assign_component_attrs(
            n.storage_units,
            is_carrier,
            marginal_cost=costs.at[carrier, "marginal_cost"],
        )
        if is_ext.any():
            _assign_component_attrs(
                n.storage_units,
                is_ext,
                capital_cost=costs.at[carrier, "capital_cost"],
                efficiency_store=costs.at[lookup_charge, "efficiency"]
                ** roundtrip_correction,
                efficiency_dispatch=costs.at[lookup_discharge, "efficiency"]
                ** roundtrip_correction,
                lifetime=costs.at[carrier, "lifetime"],
            )


def update_store_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    storage_techs: dict | None = None,
) -> None:
    """
    Update horizon-specific cost attributes on stores in-place.

    Stores hold energy capacity (MWh) in the Store-Link-Bus model and use
    ``e_nom_extendable``.

    Same rules as :func:`update_generator_costs`:

    - **marginal_cost** is refreshed for all stores of a carrier.
    - **capital_cost** and **lifetime** are refreshed only for extendable
      stores.

    Store technology names are resolved in this order:

    1. Special map (``csp`` → ``csp-tower TES`` from ``attach_advance_csp``)
    2. ``storage_techs[carrier]["store"]`` (same as ``attach_stores``)
    3. Carrier name itself if it exists in the cost table

    Carriers without a resolvable cost row are skipped with a warning.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    storage_techs : dict, optional
        Mapping of storage carriers to ``store`` / ``charger`` /
        ``discharger`` / ``bicharger`` technology names. Same structure as
        ``config["storage_techs"]``.
    """
    storage_techs = storage_techs or {}
    # Carriers not modelled via storage_techs (e.g. CSP TES from attach_advance_csp)
    special_store_techs = {
        "csp": "csp-tower TES",
    }

    for carrier in sorted(n.stores.carrier.unique()):
        if carrier in special_store_techs:
            store_tech = special_store_techs[carrier]
        elif carrier in storage_techs and "store" in storage_techs[carrier]:
            store_tech = storage_techs[carrier]["store"]
        elif carrier in costs.index:
            store_tech = carrier
        else:
            logger.warning(
                f"No store technology mapping for carrier '{carrier}'; "
                "skipping store re-costing."
            )
            continue

        if store_tech not in costs.index:
            logger.warning(
                f"No cost row '{store_tech}' for store carrier '{carrier}'; "
                "skipping store re-costing."
            )
            continue

        is_carrier = n.stores.carrier == carrier
        is_ext = is_carrier & n.stores.e_nom_extendable

        _assign_component_attrs(
            n.stores,
            is_carrier,
            marginal_cost=costs.at[store_tech, "marginal_cost"],
        )
        if is_ext.any():
            _assign_component_attrs(
                n.stores,
                is_ext,
                capital_cost=costs.at[store_tech, "capital_cost"],
                lifetime=costs.at[store_tech, "lifetime"],
            )


def update_link_costs(
    n: pypsa.Network,
    costs: pd.DataFrame,
    storage_techs: dict | None = None,
) -> None:
    """
    Update horizon-specific cost attributes on links in-place.

    Covers Store-Link converter links (chargers / dischargers), CSP discharge
    links, and optional H2 pipelines.

    Same rules as :func:`update_generator_costs`:

    - **marginal_cost** is refreshed for all links of a carrier.
    - **capital_cost**, **efficiency**, and **lifetime** are refreshed only
      for extendable (new-build) links.

    For converters added via ``attach_stores``, technology names and link
    carrier labels (e.g. ``H2 electrolysis``, ``battery charger``) come from
    ``storage_techs``, matching build-time naming. Battery / li-ion efficiencies
    use the same round-trip correction as ``attach_stores``. For bicharger
    techs, CAPEX is applied only on the charger link. Fuel-cell CAPEX is
    scaled by efficiency (cost is per MWel).

    Special cases outside ``storage_techs``:

    - **CSP** discharge links use the ``csp-tower`` cost row.
    - **H2 pipeline** CAPEX is ``costs["H2 pipeline"] * length`` for
      extendable pipeline links only. Length-based transmission efficiency
      from config is left unchanged.

    Parameters
    ----------
    n : pypsa.Network
    costs : pd.DataFrame
        Horizon-specific cost table indexed by technology.
    storage_techs : dict, optional
        Mapping of storage carriers to ``charger`` / ``discharger`` /
        ``bicharger`` technology names. Same structure as
        ``config["storage_techs"]``.
    """
    storage_techs = storage_techs or {}

    for carrier, lookup in sorted(storage_techs.items()):
        if "bicharger" in lookup:
            lookup_charge = lookup_discharge = lookup["bicharger"]
            is_bicharger = True
        elif "charger" in lookup and "discharger" in lookup:
            lookup_charge = lookup["charger"]
            lookup_discharge = lookup["discharger"]
            is_bicharger = False
        else:
            continue

        # Same link-carrier naming as attach_stores
        charge_name = "Electrolysis" if lookup_charge == "electrolysis" else "charger"
        discharge_name = (
            "Fuel Cell" if lookup_discharge == "fuel cell" else "discharger"
        )
        charge_carrier = f"{carrier} {charge_name.lower()}"
        discharge_carrier = f"{carrier} {discharge_name.lower()}"

        # storage_techs lists many carriers; skip those with no links on this network
        if not (
            (n.links.carrier == charge_carrier).any()
            or (n.links.carrier == discharge_carrier).any()
        ):
            continue

        roundtrip_correction = 0.5 if carrier in ["battery", "li-ion"] else 1

        for link_carrier, tech, role in (
            (charge_carrier, lookup_charge, "charger"),
            (discharge_carrier, lookup_discharge, "discharger"),
        ):
            if tech not in costs.index:
                logger.warning(
                    f"No cost row '{tech}' for link carrier '{link_carrier}'; "
                    "skipping link re-costing."
                )
                continue

            is_carrier = n.links.carrier == link_carrier
            if not is_carrier.any():
                continue
            is_ext = is_carrier & n.links.p_nom_extendable

            _assign_component_attrs(
                n.links,
                is_carrier,
                marginal_cost=costs.at[tech, "marginal_cost"],
            )

            if is_ext.any():
                efficiency = costs.at[tech, "efficiency"] ** roundtrip_correction
                _assign_component_attrs(
                    n.links,
                    is_ext,
                    efficiency=efficiency,
                    lifetime=costs.at[tech, "lifetime"],
                )
                # Bicharger CAPEX sits on the charger link only (discharger CAPEX = 0)
                if not (role == "discharger" and is_bicharger):
                    capital_cost = costs.at[tech, "capital_cost"]
                    # Fuel cell investment cost is per MWel (same as attach_stores)
                    if tech == "fuel cell":
                        capital_cost = capital_cost * costs.at[tech, "efficiency"]
                    _assign_component_attrs(
                        n.links,
                        is_ext,
                        capital_cost=capital_cost,
                    )

    # CSP TES discharge link (from attach_advance_csp; not in storage_techs)
    is_csp = n.links.carrier == "csp"
    if is_csp.any() and "csp-tower" in costs.index:
        is_csp_ext = is_csp & n.links.p_nom_extendable
        _assign_component_attrs(
            n.links,
            is_csp,
            marginal_cost=costs.at["csp-tower", "marginal_cost"],
        )
        _assign_component_attrs(
            n.links,
            is_csp_ext,
            efficiency=costs.at["csp-tower", "efficiency"],
            capital_cost=costs.at["csp-tower", "capital_cost"],
            lifetime=costs.at["csp-tower", "lifetime"],
        )

    h2_pipe = n.links.carrier.str.startswith("H2 pipeline", na=False)
    h2_pipe_ext = h2_pipe & n.links.p_nom_extendable
    if h2_pipe_ext.any():
        n.links.loc[h2_pipe_ext, "capital_cost"] = (
            costs.at["H2 pipeline", "capital_cost"] * n.links.loc[h2_pipe_ext, "length"]
        )
        n.links.loc[h2_pipe_ext, "lifetime"] = costs.at["H2 pipeline", "lifetime"]


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
    output_currency: str = "EUR",
    storage_techs: dict | None = None,
) -> pd.DataFrame:
    """
    Re-apply all horizon-specific cost attributes on an already-built
    electricity network, without rebuilding it.

    This is the single entry point for per-horizon re-costing. The structural
    network (topology, installed capacities, renewable profiles, build years,
    and the plant-specific efficiencies/lifetimes of *existing* units) is built
    once by ``add_electricity`` + ``add_extra_components``; this function then
    overwrites only the economic attributes that depend on the planning
    horizon, using a horizon-specific ``costs`` table. Attributes that do NOT
    depend on the planning horizon are deliberately left untouched.

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
    output_currency : str
        Currency label used only for logging.
    storage_techs : dict, optional
        Mapping of storage carriers to technology names
        (``config["storage_techs"]``).

    Returns
    -------
    pd.DataFrame
        Displacement connection costs written by :func:`update_generator_costs`.
    """
    renewable_carriers = set(renewable_carriers)

    # Updating transmission costs
    update_transmission_costs(n, costs, length_factor=length_factor)

    # Updating generator costs
    displacement_connection_costs = update_generator_costs(
        n,
        costs,
        renewable_carriers,
        length_factor,
        output_currency=output_currency,
    )

    # Updating storage unit costs
    update_storage_unit_costs(
        n,
        costs,
        hydro_capital_cost=hydro_capital_cost,
        storage_techs=storage_techs,
    )

    # Updating store costs (Store-Link-Bus energy capacity)
    update_store_costs(n, costs, storage_techs=storage_techs)

    # Updating link costs (chargers, dischargers, pipelines)
    update_link_costs(n, costs, storage_techs=storage_techs)

    return displacement_connection_costs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "assign_costs",
            simpl="",
            clusters="4",
            planning_horizons="2030",
        )

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input.network)

    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)

    displacement_connection_costs = update_electricity_costs(
        n,
        costs,
        renewable_carriers=snakemake.params.electricity["renewable_carriers"],
        length_factor=snakemake.params.length_factor,
        hydro_capital_cost=snakemake.params.hydro_capital_cost,
        output_currency=snakemake.params.output_currency,
        storage_techs=snakemake.params.storage_techs,
    )

    displacement_connection_costs.to_csv(snakemake.output.connection_costs)

    logger.info(
        "Assigned costs from planning horizon %s to network %s.",
        snakemake.wildcards.planning_horizons,
        snakemake.input.network,
    )

    n.export_to_netcdf(snakemake.output.network)
