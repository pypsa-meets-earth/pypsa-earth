from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

SCENARIOS = {
    "S0_2045_Reference": Path("results/scenarios/S0_2045_1h_final/" "S0_KZ_2045_1h.nc"),
    "S1_2045_NetZero": Path(
        "results/scenarios/S1_2045_NetZero_1h_final/" "S1_KZ_2045_NetZero_1h.nc"
    ),
}

OUTPUT = Path("results/scenarios/final_scenario_kpis.csv")


# ============================================================
# Helpers
# ============================================================


def snapshot_weights(n):
    if isinstance(n.snapshot_weightings, pd.DataFrame):
        return n.snapshot_weightings["generators"]

    return n.snapshot_weightings


def optimal_capacity(df, opt_col, nominal_col):
    """
    Use optimized nominal capacity where available,
    otherwise fall back to nominal capacity.
    """
    if opt_col in df.columns:
        x = df[opt_col].copy()
        return x.where(x.notna(), df[nominal_col])

    return df[nominal_col].copy()


def get_generator_capacity(n, carrier):
    g = n.generators[n.generators.carrier == carrier]

    if g.empty:
        return 0.0

    cap = optimal_capacity(g, "p_nom_opt", "p_nom")

    return cap.sum() / 1000


def get_generation(n, carrier, w):
    idx = n.generators.index[n.generators.carrier == carrier]

    if len(idx) == 0:
        return 0.0

    return n.generators_t.p[idx].mul(w, axis=0).sum().sum() / 1e6


def get_curtailment(n, carrier, w):
    idx = n.generators.index[n.generators.carrier == carrier]

    if len(idx) == 0:
        return np.nan, 0.0, 0.0

    cap_all = optimal_capacity(n.generators, "p_nom_opt", "p_nom")

    cap = cap_all.loc[idx]

    available = (
        n.generators_t.p_max_pu[idx].mul(cap, axis=1).mul(w, axis=0).sum().sum() / 1e6
    )

    used = n.generators_t.p[idx].mul(w, axis=0).sum().sum() / 1e6

    curtailed = max(available - used, 0.0)

    pct = 100 * curtailed / available if available > 0 else np.nan

    return pct, curtailed, available


def get_store_energy(n, carrier):
    x = n.stores[n.stores.carrier == carrier]

    if x.empty:
        return 0.0

    cap = optimal_capacity(x, "e_nom_opt", "e_nom")

    return cap.sum() / 1000


def get_link_capacity(n, carrier):
    x = n.links[n.links.carrier == carrier]

    if x.empty:
        return 0.0

    cap = optimal_capacity(x, "p_nom_opt", "p_nom")

    return cap.sum() / 1000


def get_reservoir_hydro(n, w):
    h = n.storage_units[n.storage_units.carrier == "hydro"]

    if h.empty:
        return 0.0, 0.0, 0.0

    capacity_gw = h.p_nom.sum() / 1000

    energy_gwh = (h.p_nom * h.max_hours).sum() / 1000

    if hasattr(n.storage_units_t, "p_dispatch"):
        dispatch = n.storage_units_t.p_dispatch[h.index]
    else:
        dispatch = n.storage_units_t.p[h.index].clip(lower=0)

    generation_twh = dispatch.mul(w, axis=0).sum().sum() / 1e6

    return (
        capacity_gw,
        energy_gwh,
        generation_twh,
    )


def get_co2_emissions(n, w):
    total = 0.0
    by_carrier = {}

    for carrier, idx in n.generators.groupby("carrier").groups.items():
        if carrier not in n.carriers.index:
            continue

        co2 = n.carriers.at[carrier, "co2_emissions"]

        if pd.isna(co2) or co2 == 0:
            continue

        idx = list(idx)

        dispatch = n.generators_t.p[idx]

        efficiency = n.generators.loc[idx, "efficiency"]

        fuel_input = dispatch.div(efficiency, axis=1)

        emissions = fuel_input.mul(w, axis=0).sum().sum() * co2 / 1e6

        by_carrier[carrier] = emissions
        total += emissions

    return total, by_carrier


def get_co2_constraint(n):
    """
    Returns:
        limit [MtCO2/a]
        utilisation [%]
        dual mu
    """
    if n.global_constraints.empty:
        return np.nan, np.nan, np.nan

    candidates = n.global_constraints[
        n.global_constraints["type"]
        .astype(str)
        .str.contains(
            "co2|primary_energy",
            case=False,
            regex=True,
        )
    ]

    if candidates.empty:
        return np.nan, np.nan, np.nan

    row = candidates.iloc[0]

    limit_mt = float(row["constant"]) / 1e6

    mu = (
        float(row["mu"])
        if "mu" in candidates.columns and pd.notna(row["mu"])
        else np.nan
    )

    return limit_mt, np.nan, mu


def get_prices(n, w):
    price = n.buses_t.marginal_price

    load_by_bus = n.loads_t.p_set.T.groupby(n.loads.bus).sum().T

    common = price.columns.intersection(load_by_bus.columns)

    if len(common) == 0:
        return np.nan, np.nan

    p = price[common]
    l = load_by_bus[common]

    denominator = l.mul(w, axis=0).sum().sum()

    weighted = p.mul(l).mul(w, axis=0).sum().sum() / denominator

    active_buses = l.sum(axis=0) > 0

    unweighted = p.loc[:, active_buses].mean().mean()

    return weighted, unweighted


def get_ac_line_volume(n):
    line_cap = optimal_capacity(n.lines, "s_nom_opt", "s_nom")

    return (line_cap * n.lines.length).sum() / 1e6


# ============================================================
# Main scenario extraction
# ============================================================


def extract_scenario(name, path):
    print()
    print("=" * 100)
    print(f"Reading {name}")
    print(path)
    print("=" * 100)

    if not path.exists():
        raise FileNotFoundError(path)

    n = pypsa.Network(path)

    w = snapshot_weights(n)

    # --------------------------------------------------------
    # Demand
    # --------------------------------------------------------

    total_load = n.loads_t.p_set.sum(axis=1)

    demand_twh = total_load.mul(w).sum() / 1e6

    # --------------------------------------------------------
    # Hydro
    # --------------------------------------------------------

    (
        hydro_capacity_gw,
        hydro_energy_gwh,
        hydro_generation_twh,
    ) = get_reservoir_hydro(n, w)

    # --------------------------------------------------------
    # Generation
    # --------------------------------------------------------

    coal_gen = get_generation(n, "coal", w)
    ccgt_gen = get_generation(n, "CCGT", w)
    ocgt_gen = get_generation(n, "OCGT", w)
    lignite_gen = get_generation(n, "lignite", w)
    oil_gen = get_generation(n, "oil", w)
    solar_gen = get_generation(n, "solar", w)
    wind_gen = get_generation(n, "onwind", w)
    ror_gen = get_generation(n, "ror", w)
    load_shedding = get_generation(n, "load shedding", w)

    renewable_generation = solar_gen + wind_gen + ror_gen + hydro_generation_twh

    fossil_generation = coal_gen + ccgt_gen + ocgt_gen + lignite_gen + oil_gen

    total_primary_generation = renewable_generation + fossil_generation

    # Renewable share of electricity generation.
    #
    # This is preferable to renewable generation / demand,
    # because gross renewable output can exceed final demand
    # due to storage losses.
    renewable_share = (
        100 * renewable_generation / total_primary_generation
        if total_primary_generation > 0
        else np.nan
    )

    renewable_generation_to_demand = (
        100 * renewable_generation / demand_twh if demand_twh > 0 else np.nan
    )

    # --------------------------------------------------------
    # CO2
    # --------------------------------------------------------

    co2_total, co2_by_carrier = get_co2_emissions(n, w)

    (
        co2_limit,
        _,
        co2_mu,
    ) = get_co2_constraint(n)

    if np.isfinite(co2_limit):
        if abs(co2_limit) > 1e-12:
            co2_utilisation = 100 * co2_total / co2_limit
        else:
            co2_utilisation = np.nan
    else:
        co2_utilisation = np.nan

    # --------------------------------------------------------
    # Curtailment
    # --------------------------------------------------------

    (
        solar_curt_pct,
        solar_curt_twh,
        solar_available,
    ) = get_curtailment(
        n,
        "solar",
        w,
    )

    (
        wind_curt_pct,
        wind_curt_twh,
        wind_available,
    ) = get_curtailment(
        n,
        "onwind",
        w,
    )

    (
        ror_curt_pct,
        ror_curt_twh,
        ror_available,
    ) = get_curtailment(
        n,
        "ror",
        w,
    )

    # --------------------------------------------------------
    # Prices
    # --------------------------------------------------------

    (
        weighted_price,
        unweighted_price,
    ) = get_prices(n, w)

    # --------------------------------------------------------
    # Result row
    # --------------------------------------------------------

    objective = float(n.objective)

    objective_constant = float(
        getattr(
            n,
            "objective_constant",
            np.nan,
        )
    )

    row = {
        "Scenario": name,
        # General
        "Snapshots [-]": len(n.snapshots),
        "Annual demand [TWh]": demand_twh,
        # Objective
        "Objective [EUR bn]": objective / 1e9,
        "Objective constant [EUR bn]": (
            objective_constant / 1e9 if np.isfinite(objective_constant) else np.nan
        ),
        "Decision-dependent objective [EUR bn]": (
            (objective - objective_constant) / 1e9
            if np.isfinite(objective_constant)
            else np.nan
        ),
        # CO2
        "CO2 emissions [MtCO2/a]": co2_total,
        "CO2 limit [MtCO2/a]": co2_limit,
        "CO2 limit utilisation [%]": co2_utilisation,
        "CO2 constraint dual mu [EUR/tCO2]": co2_mu,
        "CO2 shadow-price magnitude [EUR/tCO2]": (
            abs(co2_mu) if np.isfinite(co2_mu) else np.nan
        ),
        "Coal CO2 [MtCO2/a]": co2_by_carrier.get("coal", 0.0),
        "CCGT CO2 [MtCO2/a]": co2_by_carrier.get("CCGT", 0.0),
        "OCGT CO2 [MtCO2/a]": co2_by_carrier.get("OCGT", 0.0),
        # Capacity
        "Solar capacity [GW]": get_generator_capacity(n, "solar"),
        "Wind capacity [GW]": get_generator_capacity(n, "onwind"),
        "Coal capacity [GW]": get_generator_capacity(n, "coal"),
        "CCGT capacity [GW]": get_generator_capacity(n, "CCGT"),
        "OCGT capacity [GW]": get_generator_capacity(n, "OCGT"),
        "Reservoir hydro capacity [GW]": hydro_capacity_gw,
        "Reservoir hydro energy [GWh]": hydro_energy_gwh,
        "RoR capacity [GW]": get_generator_capacity(n, "ror"),
        "Battery energy [GWh]": get_store_energy(n, "battery"),
        "Battery charger [GW]": get_link_capacity(n, "battery charger"),
        "Battery discharger [GW]": get_link_capacity(n, "battery discharger"),
        "H2 store [GWh]": get_store_energy(n, "H2"),
        "H2 electrolysis [GW]": get_link_capacity(n, "H2 electrolysis"),
        "H2 fuel cell [GW]": get_link_capacity(n, "H2 fuel cell"),
        # Generation
        "Coal generation [TWh]": coal_gen,
        "CCGT generation [TWh]": ccgt_gen,
        "OCGT generation [TWh]": ocgt_gen,
        "Lignite generation [TWh]": lignite_gen,
        "Oil generation [TWh]": oil_gen,
        "Solar generation [TWh]": solar_gen,
        "Wind generation [TWh]": wind_gen,
        "Reservoir hydro generation [TWh]": hydro_generation_twh,
        "RoR generation [TWh]": ror_gen,
        "Renewable generation [TWh]": renewable_generation,
        "Fossil generation [TWh]": fossil_generation,
        "Renewable share of gross generation [%]": renewable_share,
        "Renewable generation / demand [%]": renewable_generation_to_demand,
        "Load shedding [TWh]": load_shedding,
        # Curtailment
        "Solar curtailment [%]": solar_curt_pct,
        "Solar curtailed [TWh]": solar_curt_twh,
        "Wind curtailment [%]": wind_curt_pct,
        "Wind curtailed [TWh]": wind_curt_twh,
        "RoR curtailment [%]": ror_curt_pct,
        "RoR curtailed [TWh]": ror_curt_twh,
        # Prices
        "Load-weighted marginal price [EUR/MWh]": weighted_price,
        "Unweighted mean marginal price [EUR/MWh]": unweighted_price,
        # Grid
        "AC line volume [million MWkm]": get_ac_line_volume(n),
    }

    return row


# ============================================================
# Run
# ============================================================

rows = []

for scenario, path in SCENARIOS.items():
    rows.append(extract_scenario(scenario, path))

df = pd.DataFrame(rows)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(OUTPUT, index=False)

# ============================================================
# Scenario comparison
# ============================================================

COMPARISON_OUTPUT = Path("results/scenarios/final_scenario_comparison.csv")

s0 = df.loc[df["Scenario"] == "S0_2045_Reference"].iloc[0]

s1 = df.loc[df["Scenario"] == "S1_2045_NetZero"].iloc[0]

comparison_rows = []

for column in df.columns:

    if column == "Scenario":
        continue

    value_s0 = s0[column]
    value_s1 = s1[column]

    # Determine whether KPI itself is already expressed in percent
is_percentage_kpi = "[%]" in column

if is_percentage_kpi:
    percentage_point_change = abs_change
    relative_change = np.nan
else:
    percentage_point_change = np.nan

    if pd.notna(value_s0) and pd.notna(value_s1) and abs(value_s0) > 1e-4:
        relative_change = 100 * (value_s1 - value_s0) / abs(value_s0)
    else:
        relative_change = np.nan

comparison_rows.append(
    {
        "KPI": column,
        "S0_2045_Reference": value_s0,
        "S1_2045_NetZero": value_s1,
        "Absolute change S1-S0": abs_change,
        "Percentage-point change [pp]": percentage_point_change,
        "Relative change [%]": relative_change,
    }
)

comparison = pd.DataFrame(comparison_rows)

comparison.to_csv(COMPARISON_OUTPUT, index=False)

print()
print("=" * 100)
print("S0 vs S1 COMPARISON")
print("=" * 100)

print(comparison.set_index("KPI").round(6).to_string())

print()
print(f"Saved comparison to: " f"{COMPARISON_OUTPUT}")

print()
print("=" * 100)
print("FINAL KPI TABLE")
print("=" * 100)

print(df.set_index("Scenario").T.round(6).to_string())

print()
print("=" * 100)
print(f"Saved to: {OUTPUT}")
print("=" * 100)
