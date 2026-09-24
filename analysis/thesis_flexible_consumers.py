"""
Common thesis-specific extraction functions for flexible consumers.

Purpose
-------
Provide one authoritative implementation for extracting BTC-mining
and PEM-electrolysis quantities from solved PyPSA-Earth networks.

Final BTC representation
------------------------
BTC mining is represented by:

    electricity bus
        -> THESIS BTC mining link
        -> THESIS BTC service bus
        -> THESIS BTC service accumulator

Electricity consumption is positive Link p0.

Legacy Generator representation
--------------------------------
The old Generator implementation is retained ONLY so that historical
validation networks remain readable. It must not be used for final
S2/S4/S6 thesis scenarios.

PEM representation
------------------
PEM electrolysis is represented by:

    electricity bus
        -> THESIS PEM electrolyser
        -> hydrogen product bus
        -> THESIS H2 annual product accumulator

All quantities are returned in explicit physical units.
"""

from __future__ import annotations

import pandas as pd

# =============================================================================
# Component names
# =============================================================================

BTC_LINK_NAME = "THESIS BTC mining link"
BTC_STORE_NAME = "THESIS BTC service accumulator"

LEGACY_BTC_GENERATOR_NAME = "THESIS BTC flexible mining"

PEM_LINK_NAME = "THESIS PEM electrolyser"
H2_PRODUCT_STORE_NAME = "THESIS H2 annual product accumulator"


# =============================================================================
# Snapshot weighting
# =============================================================================


def snapshot_weights(n, column: str = "generators") -> pd.Series:
    """
    Return snapshot weights aligned to the network snapshots.
    """

    if column not in n.snapshot_weightings.columns:
        raise KeyError(
            f"Snapshot-weighting column '{column}' "
            "not found. Available columns: "
            f"{list(n.snapshot_weightings.columns)}"
        )

    return n.snapshot_weightings[column].reindex(n.snapshots).astype(float)


def weighted_sum(
    series: pd.Series,
    n,
    column: str = "generators",
) -> float:
    """
    Weighted annual sum of an hourly PyPSA time series.
    """

    w = snapshot_weights(n, column)

    s = series.reindex(n.snapshots).astype(float)

    return float(s.mul(w).sum())


# =============================================================================
# BTC representation
# =============================================================================


def btc_representation(n) -> str:
    """
    Identify the BTC representation in a network.

    Returns
    -------
    "link"
        Final thesis representation.

    "legacy_generator"
        Historical validation representation.

    "none"
        No BTC component present.
    """

    has_link = BTC_LINK_NAME in n.links.index

    has_legacy_generator = LEGACY_BTC_GENERATOR_NAME in n.generators.index

    if has_link and has_legacy_generator:
        raise RuntimeError(
            "Network contains BOTH the final BTC Link "
            "and the legacy BTC Generator. "
            "This is not a valid thesis scenario."
        )

    if has_link:
        return "link"

    if has_legacy_generator:
        return "legacy_generator"

    return "none"


def btc_consumption_series_mw(n) -> pd.Series:
    """
    Positive BTC facility electricity consumption [MW].
    """

    representation = btc_representation(n)

    if representation == "link":

        return n.links_t.p0[BTC_LINK_NAME].astype(float).clip(lower=0.0)

    if representation == "legacy_generator":

        return (
            (-n.generators_t.p[LEGACY_BTC_GENERATOR_NAME]).astype(float).clip(lower=0.0)
        )

    return pd.Series(
        0.0,
        index=n.snapshots,
        dtype=float,
    )


def btc_consumption_mwh(n) -> float:
    """
    Annual BTC electricity consumption [MWh/a].
    """

    return weighted_sum(
        btc_consumption_series_mw(n),
        n,
        column="generators",
    )


def btc_consumption_twh(n) -> float:
    """
    Annual BTC electricity consumption [TWh/a].
    """

    return btc_consumption_mwh(n) / 1e6


def btc_capacity_mw(n) -> float:
    """
    Installed BTC electrical capacity [MW].
    """

    representation = btc_representation(n)

    if representation == "link":

        return float(
            n.links.at[
                BTC_LINK_NAME,
                "p_nom",
            ]
        )

    if representation == "legacy_generator":

        return float(
            n.generators.at[
                LEGACY_BTC_GENERATOR_NAME,
                "p_nom",
            ]
        )

    return 0.0


def btc_capacity_factor_pct(n) -> float:
    """
    BTC annual electrical capacity factor [%].
    """

    capacity_mw = btc_capacity_mw(n)

    if capacity_mw <= 0:
        return 0.0

    hours = float(
        snapshot_weights(
            n,
            "generators",
        ).sum()
    )

    return 100.0 * btc_consumption_mwh(n) / (capacity_mw * hours)


# =============================================================================
# BTC objective accounting
# =============================================================================


def btc_objective_contribution_eur(n) -> float:
    """
    BTC contribution to the solved PyPSA objective [EUR/a].

    Important
    ---------
    This uses the ACTUAL marginal-cost coefficient stored in the
    solved network after PyPSA-Earth's noisy_costs preparation.

    For BTC this quantity is negative because mining contributes
    economic utility/revenue to the optimization objective.
    """

    representation = btc_representation(n)

    if representation == "none":
        return 0.0

    w = snapshot_weights(
        n,
        "objective",
    )

    if representation == "link":

        dispatch = n.links_t.p0[BTC_LINK_NAME].reindex(n.snapshots).astype(float)

        marginal_cost = float(
            n.links.at[
                BTC_LINK_NAME,
                "marginal_cost",
            ]
        )

    else:

        dispatch = (
            n.generators_t.p[LEGACY_BTC_GENERATOR_NAME]
            .reindex(n.snapshots)
            .astype(float)
        )

        marginal_cost = float(
            n.generators.at[
                LEGACY_BTC_GENERATOR_NAME,
                "marginal_cost",
            ]
        )

    return float(dispatch.mul(w).sum() * marginal_cost)


def upstream_power_system_cost_eur(n) -> float:
    """
    Corrected upstream electricity-system cost proxy [EUR/a].

    Raw PyPSA objective in BTC scenarios contains the negative
    BTC dispatch-value term.

    Therefore:

        C_system = C_objective - C_BTC_objective

    where C_BTC_objective is negative.

    This removes BTC economic utility from the reported
    electricity-system expenditure proxy.
    """

    return float(n.objective) - btc_objective_contribution_eur(n)


# =============================================================================
# PEM / hydrogen
# =============================================================================


def has_pem(n) -> bool:
    """
    Return True when the thesis PEM Link is present.
    """

    return PEM_LINK_NAME in n.links.index


def pem_electricity_series_mw(n) -> pd.Series:
    """
    Positive PEM electricity consumption [MW].
    """

    if not has_pem(n):

        return pd.Series(
            0.0,
            index=n.snapshots,
            dtype=float,
        )

    return n.links_t.p0[PEM_LINK_NAME].astype(float).clip(lower=0.0)


def pem_electricity_mwh(n) -> float:
    """
    Annual PEM electricity consumption [MWh_el/a].
    """

    return weighted_sum(
        pem_electricity_series_mw(n),
        n,
        column="generators",
    )


def pem_electricity_twh(n) -> float:
    """
    Annual PEM electricity consumption [TWh_el/a].
    """

    return pem_electricity_mwh(n) / 1e6


def h2_product_series_mw_lhv(n) -> pd.Series:
    """
    Hydrogen product output [MW_H2,LHV].

    PyPSA Link p1 is negative when energy is injected into bus1,
    therefore hydrogen output is -p1.
    """

    if not has_pem(n):

        return pd.Series(
            0.0,
            index=n.snapshots,
            dtype=float,
        )

    return (-n.links_t.p1[PEM_LINK_NAME]).astype(float).clip(lower=0.0)


def h2_product_mwh_lhv(n) -> float:
    """
    Annual hydrogen product [MWh_H2,LHV/a].
    """

    return weighted_sum(
        h2_product_series_mw_lhv(n),
        n,
        column="generators",
    )


def h2_product_twh_lhv(n) -> float:
    """
    Annual hydrogen product [TWh_H2,LHV/a].
    """

    return h2_product_mwh_lhv(n) / 1e6


def h2_product_kt(n) -> float:
    """
    Annual hydrogen production [kt_H2/a].

    Thesis conversion:
        1 kt H2 = 33,330 MWh_H2,LHV
    """

    return h2_product_mwh_lhv(n) / 33330.0


# =============================================================================
# Combined flexible-consumer reporting
# =============================================================================


def flexible_electricity_twh(n) -> float:
    """
    Total electricity consumed by BTC + PEM [TWh/a].
    """

    return btc_consumption_twh(n) + pem_electricity_twh(n)


def summary(n) -> dict:
    """
    Compact thesis flexible-consumer KPI summary.
    """

    return {
        "btc_representation": btc_representation(n),
        "btc_capacity_mw": btc_capacity_mw(n),
        "btc_electricity_twh": btc_consumption_twh(n),
        "btc_capacity_factor_pct": btc_capacity_factor_pct(n),
        "btc_objective_contribution_eur": btc_objective_contribution_eur(n),
        "pem_present": has_pem(n),
        "pem_electricity_twh": pem_electricity_twh(n),
        "h2_product_twh_lhv": h2_product_twh_lhv(n),
        "h2_product_kt": h2_product_kt(n),
        "flexible_electricity_twh": flexible_electricity_twh(n),
        "raw_objective_eur": float(n.objective),
        "upstream_power_system_cost_eur": upstream_power_system_cost_eur(n),
    }
