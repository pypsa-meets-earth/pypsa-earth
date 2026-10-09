from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from thesis_flexible_consumers import (
    BTC_LINK_NAME,
    LEGACY_BTC_GENERATOR_NAME,
    btc_capacity_factor_pct,
    btc_capacity_mw,
    btc_consumption_mwh,
    btc_consumption_series_mw,
    btc_objective_contribution_eur,
    btc_representation,
    upstream_power_system_cost_eur,
)

# ============================================================
# Scenario files
# ============================================================

SCENARIO_CANDIDATES = {
    "S0": [
        "results/scenarios/final_S0_S6_Link_20260922/" "S0_KZ_2045_Reference.nc",
    ],
    "S1": [
        "results/scenarios/final_S0_S6_Link_20260922/" "S1_KZ_2045_ZeroDirectCO2.nc",
    ],
    "S2": [
        "results/scenarios/final_S0_S6_Link_20260922/" "S2_KZ_2045_Reference_BTC.nc",
    ],
    "S3": [
        "results/scenarios/final_S0_S6_Link_20260922/" "S3_KZ_2045_Reference_H2.nc",
    ],
    "S4": [
        "results/scenarios/final_S0_S6_Link_20260922/"
        "S4_KZ_2045_ZeroDirectCO2_BTC.nc",
    ],
    "S5": [
        "results/scenarios/final_S0_S6_Link_20260922/" "S5_KZ_2045_ZeroDirectCO2_H2.nc",
    ],
    "S6": [
        "results/scenarios/final_S0_S6_Link_20260922/"
        "S6_KZ_2045_ZeroDirectCO2_BTC_H2.nc",
    ],
}


BTC = LEGACY_BTC_GENERATOR_NAME
PEM = "THESIS PEM electrolyser"
H2STORE = "THESIS H2 annual product accumulator"


# ============================================================
# Carrier aliases
# ============================================================

ALIASES = {
    "solar": {
        "solar",
    },
    "wind": {
        "onwind",
        "offwind-ac",
        "offwind-dc",
    },
    "coal": {
        "coal",
    },
    "ccgt": {
        "ccgt",
    },
    "ocgt": {
        "ocgt",
    },
    "lignite": {
        "lignite",
    },
    "oil": {
        "oil",
    },
    "ror": {
        "ror",
        "run of river",
        "run-of-river",
    },
}


# ============================================================
# General helpers
# ============================================================


def find_file(candidates):
    for candidate in candidates:
        p = Path(candidate)

        if p.exists():
            return p

    raise FileNotFoundError(
        "None of these candidate files exists:\n" + "\n".join(candidates)
    )


def get_weights(n):
    if isinstance(n.snapshot_weightings, pd.DataFrame):
        if "generators" in n.snapshot_weightings.columns:
            return n.snapshot_weightings["generators"]

        return n.snapshot_weightings.iloc[:, 0]

    return n.snapshot_weightings


def optimal_values(
    table,
    opt_col,
    nominal_col,
):
    if opt_col in table.columns:
        s = table[opt_col].copy()

        if nominal_col in table.columns:
            s = s.fillna(table[nominal_col])

        return s

    return table[nominal_col].copy()


def carrier_mask(table, aliases):
    carriers = table.carrier.astype(str).str.lower()

    aliases = {x.lower() for x in aliases}

    return carriers.isin(aliases)


def dense_static_or_time(
    n,
    component,
    attr,
    names,
):
    if component == "Generator":
        static = n.generators
        dynamic = getattr(n.generators_t, attr)

    elif component == "Load":
        static = n.loads
        dynamic = getattr(n.loads_t, attr)

    else:
        raise ValueError(component)

    out = pd.DataFrame(
        index=n.snapshots,
        columns=names,
        dtype=float,
    )

    for name in names:

        if name in dynamic.columns:
            out[name] = dynamic[name]

        else:
            out[name] = float(static.at[name, attr])

    return out


# ============================================================
# Demand
# ============================================================


def annual_base_demand_twh(
    n,
    w,
):
    names = list(n.loads.index)

    if not names:
        return 0.0

    p = dense_static_or_time(
        n,
        "Load",
        "p_set",
        names,
    )

    return p.sum(axis=1).mul(w).sum() / 1e6


# ============================================================
# Generation capacity
# ============================================================


def generator_capacity_gw(
    n,
    aliases,
):
    mask = carrier_mask(n.generators, aliases)

    if not mask.any():
        return 0.0

    cap = optimal_values(
        n.generators,
        "p_nom_opt",
        "p_nom",
    )

    return cap.loc[mask].sum() / 1000


def reservoir_capacity_gw(n):
    if n.storage_units.empty:
        return 0.0

    c = n.storage_units.carrier.astype(str).str.lower()

    mask = c.str.contains("hydro") | c.str.contains("reservoir")

    if not mask.any():
        return 0.0

    cap = optimal_values(
        n.storage_units,
        "p_nom_opt",
        "p_nom",
    )

    return cap.loc[mask].sum() / 1000


# ============================================================
# Generation
# ============================================================


def generator_energy_twh(
    n,
    w,
    aliases,
):
    mask = carrier_mask(n.generators, aliases)

    idx = n.generators.index[mask]

    if len(idx) == 0:
        return 0.0

    p = n.generators_t.p[idx].clip(lower=0)

    return p.mul(w, axis=0).sum().sum() / 1e6


def reservoir_generation_twh(
    n,
    w,
):
    if n.storage_units.empty:
        return 0.0

    c = n.storage_units.carrier.astype(str).str.lower()

    mask = c.str.contains("hydro") | c.str.contains("reservoir")

    idx = n.storage_units.index[mask]

    if len(idx) == 0:
        return 0.0

    if hasattr(n.storage_units_t, "p_dispatch"):
        p = n.storage_units_t.p_dispatch[idx]
    else:
        p = n.storage_units_t.p[idx].clip(lower=0)

    return p.mul(w, axis=0).sum().sum() / 1e6


# ============================================================
# CO2
# ============================================================


def co2_for_generator_group_mt(
    n,
    w,
    aliases,
):
    mask = carrier_mask(n.generators, aliases)

    idx = n.generators.index[mask]

    if len(idx) == 0:
        return 0.0

    energy = n.generators_t.p[idx].clip(lower=0).mul(w, axis=0).sum()

    total_t = 0.0

    for name in idx:

        carrier = str(n.generators.at[name, "carrier"])

        efficiency = float(n.generators.at[name, "efficiency"])

        if efficiency <= 0:
            continue

        if carrier not in n.carriers.index:
            continue

        if "co2_emissions" not in n.carriers.columns:
            continue

        intensity = float(n.carriers.at[carrier, "co2_emissions"])

        total_t += energy[name] / efficiency * intensity

    return total_t / 1e6


# ============================================================
# Battery
# ============================================================


def battery_metrics(n):

    if n.stores.empty:
        energy_gwh = 0.0

    else:
        c = n.stores.carrier.astype(str).str.lower()

        names = n.stores.index.astype(str).str.lower()

        mask = c.str.contains("battery") | names.str.contains("battery")

        e = optimal_values(
            n.stores,
            "e_nom_opt",
            "e_nom",
        )

        energy_gwh = e.loc[mask].sum() / 1000

    if n.links.empty:
        return (
            energy_gwh,
            0.0,
            0.0,
        )

    carriers = n.links.carrier.astype(str).str.lower()

    names = n.links.index.astype(str).str.lower()

    charge_mask = carriers.str.contains("battery charger") | names.str.contains(
        "battery charger"
    )

    discharge_mask = carriers.str.contains("battery discharger") | names.str.contains(
        "battery discharger"
    )

    p = optimal_values(
        n.links,
        "p_nom_opt",
        "p_nom",
    )

    charge_gw = p.loc[charge_mask].sum() / 1000

    discharge_gw = p.loc[discharge_mask].sum() / 1000

    return (
        energy_gwh,
        charge_gw,
        discharge_gw,
    )


# ============================================================
# Curtailment
# ============================================================


def renewable_curtailment(
    n,
    w,
    aliases,
):
    mask = carrier_mask(n.generators, aliases)

    idx = n.generators.index[mask]

    if len(idx) == 0:
        return (
            0.0,
            np.nan,
        )

    p_nom = optimal_values(
        n.generators,
        "p_nom_opt",
        "p_nom",
    ).loc[idx]

    p_max_pu = dense_static_or_time(
        n,
        "Generator",
        "p_max_pu",
        list(idx),
    )

    available = p_max_pu.mul(p_nom, axis=1)

    actual = n.generators_t.p[idx].clip(lower=0)

    available_mwh = available.mul(w, axis=0).sum().sum()

    actual_mwh = actual.mul(w, axis=0).sum().sum()

    curtailed_mwh = max(
        available_mwh - actual_mwh,
        0.0,
    )

    if available_mwh > 0:
        rate = 100 * curtailed_mwh / available_mwh

    else:
        rate = np.nan

    return (
        curtailed_mwh / 1e6,
        rate,
    )


# ============================================================
# Prices
# ============================================================


def price_metrics(
    n,
    w,
):
    ac = n.buses.index[n.buses.carrier.astype(str).eq("AC")]

    prices = n.buses_t.marginal_price[ac]

    unweighted = float(prices.mean().mean())

    load_names = list(n.loads.index)

    if not load_names:
        return (
            unweighted,
            np.nan,
        )

    loads = dense_static_or_time(
        n,
        "Load",
        "p_set",
        load_names,
    )

    load_by_bus = pd.DataFrame(
        0.0,
        index=n.snapshots,
        columns=ac,
    )

    for load in load_names:

        bus = n.loads.at[load, "bus"]

        if bus in load_by_bus.columns:
            load_by_bus[bus] += loads[load]

    denominator = load_by_bus.sum(axis=1).mul(w).sum()

    numerator = prices.mul(load_by_bus).sum(axis=1).mul(w).sum()

    if denominator > 0:
        load_weighted = numerator / denominator

    else:
        load_weighted = np.nan

    return (
        unweighted,
        float(load_weighted),
    )


def consumer_weighted_price(
    n,
    w,
    bus,
    consumption,
):
    if consumption.sum() <= 0:
        return np.nan

    price = n.buses_t.marginal_price[bus]

    denominator = consumption.mul(w).sum()

    numerator = price.mul(consumption).mul(w).sum()

    if denominator <= 0:
        return np.nan

    return float(numerator / denominator)


# ============================================================
# Transmission
# ============================================================


def ac_line_volume(n):
    if n.lines.empty:
        return 0.0

    capacity = optimal_values(
        n.lines,
        "s_nom_opt",
        "s_nom",
    )

    return float((capacity * n.lines.length).sum() / 1e6)


# ============================================================
# Load shedding
# ============================================================


def load_shedding_twh(
    n,
    w,
):
    carriers = n.generators.carrier.astype(str).str.lower()

    names = n.generators.index.astype(str).str.lower()

    mask = carriers.str.contains("load shedding") | names.str.contains("load shedding")

    idx = n.generators.index[mask]

    if len(idx) == 0:
        return 0.0

    return n.generators_t.p[idx].clip(lower=0).mul(w, axis=0).sum().sum() / 1e6


# ============================================================
# Flexible consumers
# ============================================================


def btc_metrics(
    n,
    w,
):
    """
    Extract BTC-mining KPIs.

    Final thesis scenarios use the validated BTC Link
    representation. The legacy negative-Generator formulation
    remains readable only for historical regression networks.
    """

    result = {
        "BTC capacity [GW]": 0.0,
        "BTC consumption [TWh]": 0.0,
        "BTC FLH [h/a]": 0.0,
        "BTC capacity factor [%]": 0.0,
        "BTC operating hours [h]": 0.0,
        "BTC full-power hours [h]": 0.0,
        "BTC modeled value [EUR bn]": 0.0,
        "BTC weighted electricity price [EUR/MWh]": np.nan,
    }

    representation = btc_representation(n)

    if representation == "none":
        return result

    load = btc_consumption_series_mw(n)

    p_nom = btc_capacity_mw(n)

    mwh = btc_consumption_mwh(n)

    flh = mwh / p_nom if p_nom > 0 else np.nan

    cf = btc_capacity_factor_pct(n)

    # Magnitude of the actual BTC objective contribution.
    #
    # btc_objective_contribution_eur() is negative because
    # BTC provides utility/revenue to the optimization.
    # Using the solved coefficient is essential because the
    # PyPSA-Earth noisy_costs preparation slightly perturbs
    # marginal costs before optimization.
    modeled_value_bn = -btc_objective_contribution_eur(n) / 1e9

    if representation == "link":

        bus = n.links.at[
            BTC_LINK_NAME,
            "bus0",
        ]

    elif representation == "legacy_generator":

        bus = n.generators.at[
            LEGACY_BTC_GENERATOR_NAME,
            "bus",
        ]

    else:
        raise RuntimeError("Unexpected BTC representation: " f"{representation}")

    result.update(
        {
            "BTC capacity [GW]": p_nom / 1000,
            "BTC consumption [TWh]": mwh / 1e6,
            "BTC FLH [h/a]": flh,
            "BTC capacity factor [%]": cf,
            "BTC operating hours [h]": float((load > 1e-3).sum()),
            "BTC full-power hours [h]": float((load >= 0.999 * p_nom).sum()),
            "BTC modeled value [EUR bn]": modeled_value_bn,
            "BTC weighted electricity price [EUR/MWh]": consumer_weighted_price(
                n,
                w,
                bus,
                load,
            ),
        }
    )

    return result


def h2_metrics(
    n,
    w,
):
    result = {
        "PEM capacity [GW]": 0.0,
        "PEM electricity [TWh]": 0.0,
        "PEM FLH [h/a]": 0.0,
        "PEM capacity factor [%]": 0.0,
        "H2 output [TWh_LHV]": 0.0,
        "H2 production [kt/a]": 0.0,
        "PEM variable cost [EUR bn]": 0.0,
        "PEM weighted electricity price [EUR/MWh]": np.nan,
    }

    if PEM not in n.links.index:
        return result

    electricity = (n.links_t.p0[PEM]).clip(lower=0)

    output = (-n.links_t.p1[PEM]).clip(lower=0)

    p_nom = float(n.links.at[PEM, "p_nom"])

    variable_cost = float(n.links.at[PEM, "marginal_cost"])

    el_mwh = float(electricity.mul(w).sum())

    output_mwh = float(output.mul(w).sum())

    flh = el_mwh / p_nom if p_nom > 0 else np.nan

    cf = 100 * flh / float(w.sum()) if p_nom > 0 else np.nan

    h2_kg = output_mwh * 1000 / 33.33

    bus = n.links.at[PEM, "bus0"]

    result.update(
        {
            "PEM capacity [GW]": p_nom / 1000,
            "PEM electricity [TWh]": el_mwh / 1e6,
            "PEM FLH [h/a]": flh,
            "PEM capacity factor [%]": cf,
            "H2 output [TWh_LHV]": output_mwh / 1e6,
            "H2 production [kt/a]": h2_kg / 1e6,
            "PEM variable cost [EUR bn]": (el_mwh * variable_cost / 1e9),
            "PEM weighted electricity price [EUR/MWh]": consumer_weighted_price(
                n,
                w,
                bus,
                electricity,
            ),
        }
    )

    return result


# ============================================================
# One complete scenario
# ============================================================


def extract_scenario(
    scenario,
    path,
):
    print()
    print("=" * 110)
    print(f"Reading {scenario}: {path}")
    print("=" * 110)

    n = pypsa.Network(path)

    w = get_weights(n)

    # --------------------------------------------------------
    # Generation
    # --------------------------------------------------------

    solar_gen = generator_energy_twh(
        n,
        w,
        ALIASES["solar"],
    )

    wind_gen = generator_energy_twh(
        n,
        w,
        ALIASES["wind"],
    )

    ror_gen = generator_energy_twh(
        n,
        w,
        ALIASES["ror"],
    )

    reservoir_gen = reservoir_generation_twh(
        n,
        w,
    )

    coal_gen = generator_energy_twh(
        n,
        w,
        ALIASES["coal"],
    )

    ccgt_gen = generator_energy_twh(
        n,
        w,
        ALIASES["ccgt"],
    )

    ocgt_gen = generator_energy_twh(
        n,
        w,
        ALIASES["ocgt"],
    )

    lignite_gen = generator_energy_twh(
        n,
        w,
        ALIASES["lignite"],
    )

    oil_gen = generator_energy_twh(
        n,
        w,
        ALIASES["oil"],
    )

    renewable_generation = solar_gen + wind_gen + ror_gen + reservoir_gen

    fossil_generation = coal_gen + ccgt_gen + ocgt_gen + lignite_gen + oil_gen

    primary_generation = renewable_generation + fossil_generation

    if primary_generation > 0:
        renewable_share = 100 * renewable_generation / primary_generation
    else:
        renewable_share = np.nan

    # --------------------------------------------------------
    # Emissions
    # --------------------------------------------------------

    coal_co2 = co2_for_generator_group_mt(
        n,
        w,
        ALIASES["coal"],
    )

    ccgt_co2 = co2_for_generator_group_mt(
        n,
        w,
        ALIASES["ccgt"],
    )

    ocgt_co2 = co2_for_generator_group_mt(
        n,
        w,
        ALIASES["ocgt"],
    )

    lignite_co2 = co2_for_generator_group_mt(
        n,
        w,
        ALIASES["lignite"],
    )

    oil_co2 = co2_for_generator_group_mt(
        n,
        w,
        ALIASES["oil"],
    )

    total_co2 = coal_co2 + ccgt_co2 + ocgt_co2 + lignite_co2 + oil_co2

    # --------------------------------------------------------
    # Battery
    # --------------------------------------------------------

    (
        battery_e,
        battery_charge,
        battery_discharge,
    ) = battery_metrics(n)

    # --------------------------------------------------------
    # Curtailment
    # --------------------------------------------------------

    (
        solar_curt_twh,
        solar_curt_pct,
    ) = renewable_curtailment(
        n,
        w,
        ALIASES["solar"],
    )

    (
        wind_curt_twh,
        wind_curt_pct,
    ) = renewable_curtailment(
        n,
        w,
        ALIASES["wind"],
    )

    (
        ror_curt_twh,
        ror_curt_pct,
    ) = renewable_curtailment(
        n,
        w,
        ALIASES["ror"],
    )

    # --------------------------------------------------------
    # Prices
    # --------------------------------------------------------

    (
        unweighted_price,
        load_weighted_price,
    ) = price_metrics(
        n,
        w,
    )

    # --------------------------------------------------------
    # Flexible consumers
    # --------------------------------------------------------

    btc = btc_metrics(
        n,
        w,
    )

    h2 = h2_metrics(
        n,
        w,
    )

    flex_twh = btc["BTC consumption [TWh]"] + h2["PEM electricity [TWh]"]

    # --------------------------------------------------------
    # Corrected power-system cost
    # --------------------------------------------------------

    raw_objective_bn = float(n.objective) / 1e9

    objective_constant_bn = float(getattr(n, "objective_constant", np.nan)) / 1e9

    # Remove the BTC utility/revenue contribution using the
    # actual solved marginal-cost coefficient, then remove the
    # PEM-specific variable cost exactly as in the established
    # thesis system boundary.
    power_system_cost_bn = (
        upstream_power_system_cost_eur(n) / 1e9 - h2["PEM variable cost [EUR bn]"]
    )

    # --------------------------------------------------------
    # CO2 constraint
    # --------------------------------------------------------

    co2_limit = np.nan
    co2_mu = np.nan

    if not n.global_constraints.empty:

        gc = n.global_constraints

        candidate = gc[
            (gc.index.astype(str).str.lower().str.contains("co2"))
            | (gc["type"].astype(str).str.lower().str.contains("primary_energy"))
        ]

        if len(candidate):

            row = candidate.iloc[0]

            co2_limit = float(row["constant"]) / 1e6

            if pd.notna(row.get("mu", np.nan)):
                co2_mu = float(row["mu"])

    # --------------------------------------------------------
    # Final row
    # --------------------------------------------------------

    row = {
        "Scenario": scenario,
        "Snapshots [-]": len(n.snapshots),
        "Annual base demand [TWh]": annual_base_demand_twh(
            n,
            w,
        ),
        "Flexible consumption [TWh]": flex_twh,
        "Total electricity demand incl. flexible loads [TWh]": (
            annual_base_demand_twh(
                n,
                w,
            )
            + flex_twh
        ),
        "Raw optimization objective [EUR bn]": raw_objective_bn,
        "Objective constant [EUR bn]": objective_constant_bn,
        "Power-system cost proxy [EUR bn]": power_system_cost_bn,
        "CO2 emissions [MtCO2/a]": total_co2,
        "CO2 limit [MtCO2/a]": co2_limit,
        "CO2 dual mu [EUR/tCO2]": co2_mu,
        "CO2 shadow-price magnitude [EUR/tCO2]": (
            abs(co2_mu) if pd.notna(co2_mu) else np.nan
        ),
        # Capacity
        "Solar capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["solar"],
        ),
        "Wind capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["wind"],
        ),
        "Coal capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["coal"],
        ),
        "CCGT capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["ccgt"],
        ),
        "OCGT capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["ocgt"],
        ),
        "RoR capacity [GW]": generator_capacity_gw(
            n,
            ALIASES["ror"],
        ),
        "Reservoir hydro capacity [GW]": reservoir_capacity_gw(n),
        "Battery energy capacity [GWh]": battery_e,
        "Battery charging power [GW]": battery_charge,
        "Battery discharging power [GW]": battery_discharge,
        # Generation
        "Solar generation [TWh]": solar_gen,
        "Wind generation [TWh]": wind_gen,
        "RoR generation [TWh]": ror_gen,
        "Reservoir hydro generation [TWh]": reservoir_gen,
        "Coal generation [TWh]": coal_gen,
        "CCGT generation [TWh]": ccgt_gen,
        "OCGT generation [TWh]": ocgt_gen,
        "Lignite generation [TWh]": lignite_gen,
        "Oil generation [TWh]": oil_gen,
        "Renewable generation [TWh]": renewable_generation,
        "Fossil generation [TWh]": fossil_generation,
        "Renewable share of gross primary generation [%]": renewable_share,
        # Curtailment
        "Solar curtailment [TWh]": solar_curt_twh,
        "Solar curtailment [%]": solar_curt_pct,
        "Wind curtailment [TWh]": wind_curt_twh,
        "Wind curtailment [%]": wind_curt_pct,
        "RoR curtailment [TWh]": ror_curt_twh,
        "RoR curtailment [%]": ror_curt_pct,
        # Prices/grid/reliability
        "Load-weighted marginal price [EUR/MWh]": load_weighted_price,
        "Unweighted marginal price [EUR/MWh]": unweighted_price,
        "AC line volume [million MWkm]": ac_line_volume(n),
        "Load shedding [TWh]": load_shedding_twh(
            n,
            w,
        ),
    }

    row.update(btc)
    row.update(h2)

    return row


# ============================================================
# Extract all scenarios
# ============================================================

rows = []

resolved_paths = {}

for scenario in [
    "S0",
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
]:

    path = find_file(SCENARIO_CANDIDATES[scenario])

    resolved_paths[scenario] = str(path)

    rows.append(
        extract_scenario(
            scenario,
            path,
        )
    )


df = pd.DataFrame(rows).set_index("Scenario")


# ============================================================
# Final flexible-consumer representation invariants
# ============================================================

# S2/S4/S6 are the production BTC scenarios and MUST use the
# validated Link representation. Re-open only these three
# networks here so that a future accidental path change cannot
# silently reintroduce the legacy Generator formulation.
for scenario in [
    "S2",
    "S4",
    "S6",
]:
    check_n = pypsa.Network(resolved_paths[scenario])

    representation = btc_representation(check_n)

    if representation != "link":
        raise AssertionError(
            f"{scenario}: expected final BTC Link "
            f"representation, found "
            f"'{representation}'."
        )

    btc_capacity = float(
        df.loc[
            scenario,
            "BTC capacity [GW]",
        ]
    )

    if abs(btc_capacity - 1.0) > 1e-9:
        raise AssertionError(
            f"{scenario}: BTC capacity is " f"{btc_capacity} GW, expected 1.0 GW."
        )


# Scenarios without BTC must remain BTC-free.
for scenario in [
    "S0",
    "S1",
    "S3",
    "S5",
]:
    btc_energy = float(
        df.loc[
            scenario,
            "BTC consumption [TWh]",
        ]
    )

    if abs(btc_energy) > 1e-9:
        raise AssertionError(
            f"{scenario}: unexpected BTC consumption " f"{btc_energy} TWh/a."
        )


# S3/S5/S6 must retain the fixed annual 100 kt H2 service.
for scenario in [
    "S3",
    "S5",
    "S6",
]:
    h2_kt = float(
        df.loc[
            scenario,
            "H2 production [kt/a]",
        ]
    )

    if abs(h2_kt - 100.0) > 1e-6:
        raise AssertionError(
            f"{scenario}: H2 production is " f"{h2_kt} kt/a, expected 100 kt/a."
        )


print("\\nFinal flexible-consumer representation " "invariants: PASSED")


# ============================================================
# Baseline-relative effects
# ============================================================

BASELINES = {
    "S2": "S0",
    "S3": "S0",
    "S4": "S1",
    "S5": "S1",
    "S6": "S1",
}

df["Incremental power-system cost vs baseline [EUR bn]"] = np.nan

df["Incremental power-system cost per flexible MWh [EUR/MWh]"] = np.nan

df["Incremental CO2 vs baseline [MtCO2/a]"] = np.nan

for scenario, baseline in BASELINES.items():

    delta_cost_bn = (
        df.loc[scenario, "Power-system cost proxy [EUR bn]"]
        - df.loc[baseline, "Power-system cost proxy [EUR bn]"]
    )

    df.loc[scenario, "Incremental power-system cost vs baseline [EUR bn]"] = (
        delta_cost_bn
    )

    flex_twh = df.loc[scenario, "Flexible consumption [TWh]"]

    if flex_twh > 0:

        df.loc[scenario, "Incremental power-system cost per flexible MWh [EUR/MWh]"] = (
            delta_cost_bn * 1e9 / (flex_twh * 1e6)
        )

    df.loc[scenario, "Incremental CO2 vs baseline [MtCO2/a]"] = (
        df.loc[scenario, "CO2 emissions [MtCO2/a]"]
        - df.loc[baseline, "CO2 emissions [MtCO2/a]"]
    )


# ============================================================
# Pairwise thesis comparisons
# ============================================================

PAIRS = {
    "S2-S0 Reference BTC effect": ("S2", "S0"),
    "S3-S0 Reference H2 effect": ("S3", "S0"),
    "S4-S1 NetZero BTC effect": ("S4", "S1"),
    "S5-S1 NetZero H2 effect": ("S5", "S1"),
    "S6-S1 NetZero combined effect": ("S6", "S1"),
}


PAIR_KPIS = [
    "Power-system cost proxy [EUR bn]",
    "CO2 emissions [MtCO2/a]",
    "Solar capacity [GW]",
    "Wind capacity [GW]",
    "Battery energy capacity [GWh]",
    "Battery charging power [GW]",
    "Battery discharging power [GW]",
    "Renewable generation [TWh]",
    "Solar curtailment [TWh]",
    "Wind curtailment [TWh]",
    "RoR curtailment [TWh]",
    "Load-weighted marginal price [EUR/MWh]",
    "AC line volume [million MWkm]",
    "BTC consumption [TWh]",
    "BTC capacity factor [%]",
    "PEM electricity [TWh]",
    "H2 production [kt/a]",
]


pair_rows = []

for comparison, (
    scenario,
    baseline,
) in PAIRS.items():

    for kpi in PAIR_KPIS:

        scenario_value = df.loc[scenario, kpi]

        baseline_value = df.loc[baseline, kpi]

        absolute_change = scenario_value - baseline_value

        if pd.notna(baseline_value) and baseline_value != 0:
            relative_change = 100 * absolute_change / abs(baseline_value)
        else:
            relative_change = np.nan

        pair_rows.append(
            {
                "Comparison": comparison,
                "KPI": kpi,
                "Baseline": baseline,
                "Scenario": scenario,
                "Baseline value": baseline_value,
                "Scenario value": scenario_value,
                "Absolute change": absolute_change,
                "Relative change [%]": relative_change,
            }
        )


pair_df = pd.DataFrame(pair_rows)


# ============================================================
# Formal S6 interaction terms
# ============================================================

INTERACTION_KPIS = [
    "Power-system cost proxy [EUR bn]",
    "Solar capacity [GW]",
    "Wind capacity [GW]",
    "Battery energy capacity [GWh]",
    "Battery charging power [GW]",
    "Battery discharging power [GW]",
    "Renewable generation [TWh]",
    "Solar curtailment [TWh]",
    "Wind curtailment [TWh]",
    "RoR curtailment [TWh]",
    "Load-weighted marginal price [EUR/MWh]",
    "AC line volume [million MWkm]",
]


interaction_rows = []

for kpi in INTERACTION_KPIS:

    s1 = df.loc["S1", kpi]

    s4 = df.loc["S4", kpi]

    s5 = df.loc["S5", kpi]

    s6 = df.loc["S6", kpi]

    btc_effect = s4 - s1

    h2_effect = s5 - s1

    additive_expected = s1 + btc_effect + h2_effect

    interaction = s6 - additive_expected

    interaction_rows.append(
        {
            "KPI": kpi,
            "S1 baseline": s1,
            "S4 BTC": s4,
            "S5 H2": s5,
            "S6 BTC+H2": s6,
            "BTC effect S4-S1": btc_effect,
            "H2 effect S5-S1": h2_effect,
            "Additive expectation": additive_expected,
            "Interaction S6-(S4+S5-S1)": interaction,
        }
    )


interaction_df = pd.DataFrame(interaction_rows).set_index("KPI")


# ============================================================
# Direct BTC-H2 competition indicators
# ============================================================

competition = pd.DataFrame(
    {
        "Metric": [
            "BTC consumption S4 [TWh]",
            "BTC consumption S6 [TWh]",
            "BTC consumption change S6-S4 [TWh]",
            "BTC capacity factor S4 [%]",
            "BTC capacity factor S6 [%]",
            "BTC CF change S6-S4 [percentage points]",
            "PEM weighted price S5 [EUR/MWh]",
            "PEM weighted price S6 [EUR/MWh]",
            "PEM weighted price change S6-S5 [EUR/MWh]",
        ],
        "Value": [
            df.loc["S4", "BTC consumption [TWh]"],
            df.loc["S6", "BTC consumption [TWh]"],
            (
                df.loc["S6", "BTC consumption [TWh]"]
                - df.loc["S4", "BTC consumption [TWh]"]
            ),
            df.loc["S4", "BTC capacity factor [%]"],
            df.loc["S6", "BTC capacity factor [%]"],
            (
                df.loc["S6", "BTC capacity factor [%]"]
                - df.loc["S4", "BTC capacity factor [%]"]
            ),
            df.loc["S5", "PEM weighted electricity price [EUR/MWh]"],
            df.loc["S6", "PEM weighted electricity price [EUR/MWh]"],
            (
                df.loc["S6", "PEM weighted electricity price [EUR/MWh]"]
                - df.loc["S5", "PEM weighted electricity price [EUR/MWh]"]
            ),
        ],
    }
)


# ============================================================
# Save
# ============================================================

OUT = Path("results/scenarios")

OUT.mkdir(parents=True, exist_ok=True)


kpi_file = OUT / "final_S0_S6_flexible_kpis.csv"

pair_file = OUT / "final_S0_S6_pairwise_changes.csv"

interaction_file = OUT / "final_S0_S6_interaction_terms.csv"

competition_file = OUT / "final_S4_S6_BTC_H2_competition.csv"


df.to_csv(kpi_file)

pair_df.to_csv(
    pair_file,
    index=False,
)

interaction_df.to_csv(interaction_file)

competition.to_csv(
    competition_file,
    index=False,
)


# ============================================================
# Print headline table
# ============================================================

DISPLAY = [
    "Annual base demand [TWh]",
    "Flexible consumption [TWh]",
    "Power-system cost proxy [EUR bn]",
    "Incremental power-system cost vs baseline [EUR bn]",
    "CO2 emissions [MtCO2/a]",
    "Solar capacity [GW]",
    "Wind capacity [GW]",
    "Battery energy capacity [GWh]",
    "Solar curtailment [TWh]",
    "Wind curtailment [TWh]",
    "Load-weighted marginal price [EUR/MWh]",
    "AC line volume [million MWkm]",
    "BTC consumption [TWh]",
    "BTC capacity factor [%]",
    "PEM electricity [TWh]",
    "H2 production [kt/a]",
]


print()
print("=" * 150)
print("FINAL S0-S6 THESIS KPI TABLE")
print("=" * 150)

print(df[DISPLAY].T.to_string(float_format=lambda x: f"{x:.6f}"))

print()
print("=" * 150)
print("S6 FORMAL INTERACTION TERMS")
print("=" * 150)

print(interaction_df.to_string(float_format=lambda x: f"{x:.6f}"))

print()
print("=" * 150)
print("DIRECT BTC-H2 COMPETITION INDICATORS")
print("=" * 150)

print(competition.to_string(index=False, float_format=lambda x: f"{x:.6f}"))

print()
print("Resolved network files:")

for scenario, path in resolved_paths.items():
    print(f"{scenario}: {path}")

print()
print("Written:")
print(kpi_file)
print(pair_file)
print(interaction_file)
print(competition_file)

print("=" * 150)
