"""
Add thesis-specific flexible consumers to a prepared PyPSA-Earth network.

Final thesis representations
----------------------------
Bitcoin mining
    Implemented exclusively through scripts/thesis_btc.py as

        electricity bus
            -> BTC mining Link
            -> BTC service bus
            -> BTC service accumulator Store

    The BTC Link has positive electricity withdrawal p0 and a
    negative marginal cost equal to the derived marginal value of
    mining electricity.

PEM hydrogen production
    Implemented as

        electricity bus
            -> PEM electrolyser Link
            -> H2 product bus
            -> annual H2 product accumulator Store

    The H2 Store is a bookkeeping accumulator for the annual product
    target. It is not physical hydrogen storage.

This script creates an UNSOLVED network.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import sys

import pandas as pd
import pypsa


# =============================================================================
# Make repository root importable when script is called as:
#
#     python analysis/add_flexible_consumer.py ...
# =============================================================================

REPO_ROOT = Path(__file__).resolve().parents[1]

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(
        0,
        str(REPO_ROOT),
    )


from scripts.thesis_btc import (  # noqa: E402
    BTC_LINK_NAME,
    BTC_BUS_NAME,
    BTC_STORE_NAME,
    DEFAULT_HASHPRICE_USD2026_PER_PH_DAY,
    DEFAULT_ASIC_EFFICIENCY_J_PER_TH,
    DEFAULT_PUE,
    DEFAULT_USD_PER_EUR,
    DEFAULT_EUR2026_TO_EUR2020_FACTOR,
    DEFAULT_NON_ELECTRIC_OPEX_EUR2020_PER_MWH,
    add_btc_mining_link,
    btc_hashprice_eur2020_per_th_day,
    btc_electricity_value_eur2020_per_mwh,
)


# =============================================================================
# General helpers
# =============================================================================

def add_carrier_if_missing(
    n,
    carrier,
):
    if carrier not in n.carriers.index:
        n.add(
            "Carrier",
            carrier,
        )


def ensure_free(
    n,
    component,
    name,
):
    table = {
        "Generator": n.generators.index,
        "Bus": n.buses.index,
        "Link": n.links.index,
        "Store": n.stores.index,
    }

    if name in table[component]:
        raise RuntimeError(
            f"{component} '{name}' already exists."
        )


# =============================================================================
# Hydrogen production
# =============================================================================

def add_hydrogen(
    n,
    bus,
    p_nom_mw,
    annual_h2_kt,
    efficiency_lhv,
    variable_cost_eur_per_mwh_el,
):
    """
    Add PEM electrolysis with an annual H2-production target.

    Electricity:
        AC bus
          |
          v
        PEM Link
          |
          v
        H2 product bus
          |
          v
        cumulative-product Store

    The Store is NOT physical H2 storage.

    It is a bookkeeping accumulator used to impose:

        sum_t H2_output_t = annual target

    on an LHV basis.
    """

    h2_bus = "THESIS H2 product bus"
    link = "THESIS PEM electrolyser"
    accumulator = "THESIS H2 annual product accumulator"

    carrier = "thesis H2 product"

    ensure_free(
        n,
        "Bus",
        h2_bus,
    )

    ensure_free(
        n,
        "Link",
        link,
    )

    ensure_free(
        n,
        "Store",
        accumulator,
    )

    add_carrier_if_missing(
        n,
        carrier,
    )

    # -------------------------------------------------------------------------
    # Annual H2 target
    # -------------------------------------------------------------------------

    # 1 kt = 1,000,000 kg
    h2_kg = (
        annual_h2_kt
        * 1_000_000
    )

    # Thesis LHV convention: 33.33 kWh/kg
    target_mwh_h2 = (
        h2_kg
        * 33.33
        / 1000
    )

    # Corresponding electricity requirement
    target_mwh_el = (
        target_mwh_h2
        / efficiency_lhv
    )

    # -------------------------------------------------------------------------
    # Feasibility check
    # -------------------------------------------------------------------------

    if isinstance(
        n.snapshot_weightings,
        pd.DataFrame,
    ):
        if (
            "generators"
            in n.snapshot_weightings.columns
        ):
            weighted_hours = float(
                n.snapshot_weightings[
                    "generators"
                ].sum()
            )
        else:
            weighted_hours = float(
                n.snapshot_weightings
                .iloc[:, 0]
                .sum()
            )
    else:
        weighted_hours = float(
            n.snapshot_weightings.sum()
        )

    max_mwh_el = (
        p_nom_mw
        * weighted_hours
    )

    if (
        target_mwh_el
        > max_mwh_el + 1e-6
    ):
        raise ValueError(
            "H2 target is infeasible with selected PEM capacity.\n"
            f"Required electricity : {target_mwh_el:.3f} MWh\n"
            f"Maximum electricity  : {max_mwh_el:.3f} MWh"
        )

    # -------------------------------------------------------------------------
    # H2 product bus
    # -------------------------------------------------------------------------

    n.add(
        "Bus",
        h2_bus,
        carrier=carrier,
    )

    # -------------------------------------------------------------------------
    # PEM electrolyser
    # -------------------------------------------------------------------------

    n.add(
        "Link",
        link,

        bus0=bus,
        bus1=h2_bus,

        carrier=carrier,

        p_nom=p_nom_mw,
        p_nom_extendable=False,

        # Simplified hourly flexibility
        p_min_pu=0.0,
        p_max_pu=1.0,

        efficiency=efficiency_lhv,

        # Stack-throughput / non-electric variable cost
        marginal_cost=(
            variable_cost_eur_per_mwh_el
        ),

        capital_cost=0.0,
    )

    # -------------------------------------------------------------------------
    # Annual H2 product accumulator
    # -------------------------------------------------------------------------

    e_min = pd.Series(
        0.0,
        index=n.snapshots,
    )

    e_max = pd.Series(
        1.0,
        index=n.snapshots,
    )

    # Force final cumulative product exactly to target.
    e_min.iloc[-1] = 1.0
    e_max.iloc[-1] = 1.0

    n.add(
        "Store",
        accumulator,

        bus=h2_bus,
        carrier=carrier,

        e_nom=target_mwh_h2,
        e_nom_extendable=False,

        e_initial=0.0,
        e_cyclic=False,

        e_min_pu=e_min,
        e_max_pu=e_max,

        standing_loss=0.0,

        capital_cost=0.0,
        marginal_cost=0.0,
    )

    return {
        "type": "H2",
        "link": link,
        "store": accumulator,
        "target_h2_mwh": target_mwh_h2,
        "target_el_mwh": target_mwh_el,
        "target_h2_kg": h2_kg,
    }


# =============================================================================
# CLI
# =============================================================================

parser = argparse.ArgumentParser(
    description=(
        "Add thesis BTC mining or PEM hydrogen production "
        "to a prepared PyPSA-Earth network."
    )
)

parser.add_argument(
    "--base-network",
    required=True,
)

parser.add_argument(
    "--output-network",
    required=True,
)

parser.add_argument(
    "--consumer",
    choices=[
        "btc",
        "h2",
    ],
    required=True,
)

parser.add_argument(
    "--bus",
    default="KZ0 0",
)

parser.add_argument(
    "--p-nom-mw",
    type=float,
    required=True,
)


# =============================================================================
# BTC source-chain parameters
#
# Defaults are the frozen thesis-central assumptions from thesis_btc.py.
# They are exposed here for transparent sensitivity/reproduction work.
#
# There is intentionally NO direct --btc-value-eur-per-mwh argument.
# The BTC electricity value must be derived from the underlying assumptions.
# =============================================================================

parser.add_argument(
    "--btc-hashprice-usd2026-per-ph-day",
    type=float,
    default=(
        DEFAULT_HASHPRICE_USD2026_PER_PH_DAY
    ),
)

parser.add_argument(
    "--btc-asic-efficiency-j-per-th",
    type=float,
    default=(
        DEFAULT_ASIC_EFFICIENCY_J_PER_TH
    ),
)

parser.add_argument(
    "--btc-pue",
    type=float,
    default=DEFAULT_PUE,
)

parser.add_argument(
    "--btc-usd-per-eur",
    type=float,
    default=DEFAULT_USD_PER_EUR,
)

parser.add_argument(
    "--btc-eur2026-to-eur2020-factor",
    type=float,
    default=(
        DEFAULT_EUR2026_TO_EUR2020_FACTOR
    ),
)

parser.add_argument(
    "--btc-non-electric-opex-eur2020-per-mwh",
    type=float,
    default=(
        DEFAULT_NON_ELECTRIC_OPEX_EUR2020_PER_MWH
    ),
)


# =============================================================================
# H2-specific parameters
# =============================================================================

parser.add_argument(
    "--annual-h2-kt",
    type=float,
    default=None,
)

parser.add_argument(
    "--h2-efficiency-lhv",
    type=float,
    default=0.640,
)

parser.add_argument(
    "--h2-variable-cost-eur-per-mwh-el",
    type=float,
    default=1.884615,
)


args = parser.parse_args()


# =============================================================================
# Input validation
# =============================================================================

base = Path(
    args.base_network
)

output = Path(
    args.output_network
)

if not base.exists():
    raise FileNotFoundError(
        base
    )

if args.p_nom_mw <= 0:
    raise ValueError(
        "p_nom_mw must be > 0."
    )


n = pypsa.Network(
    base
)


if args.bus not in n.buses.index:
    raise ValueError(
        f"Bus '{args.bus}' does not exist."
    )

if str(
    n.buses.at[
        args.bus,
        "carrier",
    ]
) != "AC":
    raise ValueError(
        f"Bus '{args.bus}' is not an AC bus."
    )


# =============================================================================
# Add selected consumer
# =============================================================================

if args.consumer == "btc":

    gross_value = (
        btc_electricity_value_eur2020_per_mwh(
            hashprice_usd_per_ph_day=(
                args.btc_hashprice_usd2026_per_ph_day
            ),
            asic_efficiency_j_per_th=(
                args.btc_asic_efficiency_j_per_th
            ),
            pue=args.btc_pue,
            usd_per_eur=(
                args.btc_usd_per_eur
            ),
            eur2026_to_eur2020_factor=(
                args.btc_eur2026_to_eur2020_factor
            ),
        )
    )

    hashprice_eur2020_per_th_day = (
        btc_hashprice_eur2020_per_th_day(
            hashprice_usd_per_ph_day=(
                args.btc_hashprice_usd2026_per_ph_day
            ),
            usd_per_eur=(
                args.btc_usd_per_eur
            ),
            eur2026_to_eur2020_factor=(
                args.btc_eur2026_to_eur2020_factor
            ),
        )
    )

    net_value = (
        gross_value
        - args.btc_non_electric_opex_eur2020_per_mwh
    )

    add_btc_mining_link(
        n=n,
        electricity_bus=args.bus,
        p_nom_mw=args.p_nom_mw,
        hashprice_usd_per_ph_day=(
            args.btc_hashprice_usd2026_per_ph_day
        ),
        asic_efficiency_j_per_th=(
            args.btc_asic_efficiency_j_per_th
        ),
        pue=args.btc_pue,
        usd_per_eur=(
            args.btc_usd_per_eur
        ),
        eur2026_to_eur2020_factor=(
            args.btc_eur2026_to_eur2020_factor
        ),
        non_electric_opex_eur2020_per_mwh=(
            args.btc_non_electric_opex_eur2020_per_mwh
        ),
    )

    result = {
        "type": "BTC",
        "link": BTC_LINK_NAME,
        "service_bus": BTC_BUS_NAME,
        "store": BTC_STORE_NAME,
        "hashprice_eur2020_per_th_day":
            hashprice_eur2020_per_th_day,
        "gross_value_eur2020_per_mwh":
            gross_value,
        "net_value_eur2020_per_mwh":
            net_value,
    }


elif args.consumer == "h2":

    if args.annual_h2_kt is None:
        raise ValueError(
            "--annual-h2-kt is required "
            "for H2 scenarios."
        )

    if not (
        0
        < args.h2_efficiency_lhv
        <= 1
    ):
        raise ValueError(
            "H2 efficiency must be in (0,1]."
        )

    result = add_hydrogen(
        n=n,
        bus=args.bus,
        p_nom_mw=args.p_nom_mw,
        annual_h2_kt=args.annual_h2_kt,
        efficiency_lhv=(
            args.h2_efficiency_lhv
        ),
        variable_cost_eur_per_mwh_el=(
            args.h2_variable_cost_eur_per_mwh_el
        ),
    )


# =============================================================================
# Metadata
# =============================================================================

n.meta = dict(
    getattr(
        n,
        "meta",
        {},
    )
)

n.meta.update({
    "thesis_flexible_consumer":
        args.consumer,

    "thesis_connection_bus":
        args.bus,

    "thesis_consumer_p_nom_mw":
        float(args.p_nom_mw),

    "thesis_model_note":
        (
            "Flexible-consumer extension added to frozen "
            "PyPSA-Earth prepared network."
        ),
})


if args.consumer == "btc":

    n.meta.update({
        "thesis_btc_representation":
            "Link-service-bus-accumulator",

        "thesis_btc_hashprice_usd2026_per_ph_day":
            float(
                args.btc_hashprice_usd2026_per_ph_day
            ),

        "thesis_btc_hashprice_eur2020_per_th_day":
            float(
                result[
                    "hashprice_eur2020_per_th_day"
                ]
            ),

        "thesis_btc_asic_efficiency_j_per_th":
            float(
                args.btc_asic_efficiency_j_per_th
            ),

        "thesis_btc_pue":
            float(
                args.btc_pue
            ),

        "thesis_btc_usd_per_eur":
            float(
                args.btc_usd_per_eur
            ),

        "thesis_btc_eur2026_to_eur2020_factor":
            float(
                args.btc_eur2026_to_eur2020_factor
            ),

        "thesis_btc_non_electric_opex_eur2020_per_mwh":
            float(
                args.btc_non_electric_opex_eur2020_per_mwh
            ),

        "thesis_btc_gross_value_eur2020_per_mwh":
            float(
                result[
                    "gross_value_eur2020_per_mwh"
                ]
            ),

        "thesis_btc_net_value_eur2020_per_mwh":
            float(
                result[
                    "net_value_eur2020_per_mwh"
                ]
            ),
    })


if args.consumer == "h2":

    n.meta.update({
        "thesis_h2_target_kt_per_year":
            float(
                args.annual_h2_kt
            ),

        "thesis_h2_efficiency_lhv":
            float(
                args.h2_efficiency_lhv
            ),

        "thesis_h2_variable_cost_eur_per_mwh_el":
            float(
                args.h2_variable_cost_eur_per_mwh_el
            ),
    })


# =============================================================================
# Export
# =============================================================================

output.parent.mkdir(
    parents=True,
    exist_ok=True,
)

n.export_to_netcdf(
    output
)


# =============================================================================
# Report
# =============================================================================

print()
print("=" * 100)
print(
    "FLEXIBLE-CONSUMER NETWORK CREATED"
)
print("=" * 100)

print(
    f"Consumer       : {result['type']}"
)

print(
    f"Connection bus : {args.bus}"
)

print(
    f"Capacity       : "
    f"{args.p_nom_mw:.3f} MW"
)

print(
    f"Base network   : {base}"
)

print(
    f"Output network : {output}"
)


if args.consumer == "btc":

    print()
    print(
        "BTC representation      : "
        "Link -> service bus -> accumulator Store"
    )

    print(
        "BTC hashprice           : "
        f"{args.btc_hashprice_usd2026_per_ph_day:.6f} "
        "USD2026/(PH/s)/day"
    )

    print(
        "BTC converted hashprice : "
        f"{result['hashprice_eur2020_per_th_day']:.12f} "
        "EUR2020/(TH/s)/day"
    )

    print(
        "ASIC efficiency         : "
        f"{args.btc_asic_efficiency_j_per_th:.6f} J/TH"
    )

    print(
        "PUE                     : "
        f"{args.btc_pue:.6f}"
    )

    print(
        "BTC gross electricity value : "
        f"{result['gross_value_eur2020_per_mwh']:.12f} "
        "EUR2020/MWh"
    )

    print(
        "BTC net electricity value   : "
        f"{result['net_value_eur2020_per_mwh']:.12f} "
        "EUR2020/MWh"
    )

    print(
        "BTC annual electricity  : ENDOGENOUS"
    )

    print(
        "BTC flexibility         : "
        "0-100%, no deferred-energy requirement"
    )


else:

    print()
    print(
        "Annual H2 target       : "
        f"{args.annual_h2_kt:.6f} kt/a"
    )

    print(
        "H2 target energy       : "
        f"{result['target_h2_mwh']/1000:.3f} "
        "GWh_H2 LHV/a"
    )

    print(
        "Required electricity   : "
        f"{result['target_el_mwh']/1000:.3f} "
        "GWh_el/a"
    )

    print(
        "PEM efficiency         : "
        f"{args.h2_efficiency_lhv:.4f} LHV"
    )

    print(
        "PEM hourly dispatch    : 0-100%"
    )


print()
print(
    "The network has NOT been solved."
)

print("=" * 100)
