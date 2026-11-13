from pathlib import Path

import numpy as np
import pandas as pd
import pypsa

# ============================================================
# FINAL THESIS SCENARIOS
# ============================================================

BASE = Path("results/thesis_2045/networks")
OUT = Path("results/thesis_2045/final_analysis")
OUT.mkdir(parents=True, exist_ok=True)

SCENARIOS = {
    "S0": {
        "file": BASE / "elec_s_10_ec_lcopt_S0-1h.nc",
        "opts": "S0-1h",
        "co2": False,
        "btc": False,
        "pem": False,
        "regime": "reference",
    },
    "S1": {
        "file": BASE / "elec_s_10_ec_lcopt_S1-Co2L-1h.nc",
        "opts": "S1-Co2L-1h",
        "co2": True,
        "btc": False,
        "pem": False,
        "regime": "zero-direct-CO2",
    },
    "S2": {
        "file": BASE / "elec_s_10_ec_lcopt_S2-BTC-1h.nc",
        "opts": "S2-BTC-1h",
        "co2": False,
        "btc": True,
        "pem": False,
        "regime": "reference",
    },
    "S3": {
        "file": BASE / "elec_s_10_ec_lcopt_S3-PEM-1h.nc",
        "opts": "S3-PEM-1h",
        "co2": False,
        "btc": False,
        "pem": True,
        "regime": "reference",
    },
    "S4": {
        "file": BASE / "elec_s_10_ec_lcopt_S4-Co2L-BTC-1h.nc",
        "opts": "S4-Co2L-BTC-1h",
        "co2": True,
        "btc": True,
        "pem": False,
        "regime": "zero-direct-CO2",
    },
    "S5": {
        "file": BASE / "elec_s_10_ec_lcopt_S5-Co2L-PEM-1h.nc",
        "opts": "S5-Co2L-PEM-1h",
        "co2": True,
        "btc": False,
        "pem": True,
        "regime": "zero-direct-CO2",
    },
    "S6": {
        "file": BASE / "elec_s_10_ec_lcopt_S6-Co2L-BTC-PEM-1h.nc",
        "opts": "S6-Co2L-BTC-PEM-1h",
        "co2": True,
        "btc": True,
        "pem": True,
        "regime": "zero-direct-CO2",
    },
}


BTC_BUS = "KZ0 0 bitcoin"
BTC_LINK = "KZ0 0 Bitcoin mining"
BTC_SERVICE = "KZ0 0 Bitcoin service"

PEM_BUS = "THESIS H2 product bus"
PEM_LINK = "THESIS PEM electrolyser"
PEM_STORE = "THESIS H2 annual product accumulator"

CONNECTION_BUS = "KZ0 0"

EXPECTED_BASE_DEMAND_TWH = 187.0
EXPECTED_BTC_MW = 1000.0
EXPECTED_BTC_VALUE = 62.7967

EXPECTED_PEM_MW = 1000.0
EXPECTED_PEM_EFF = 0.64
EXPECTED_PEM_MC = 1.884615
EXPECTED_H2_MWH = 3_333_000.0
EXPECTED_PEM_EL_MWH = 5_207_812.5
EXPECTED_H2_KT = 100.0


# ============================================================
# HELPERS
# ============================================================


def nominal_capacity(df, base_col, opt_col):
    out = df[base_col].astype(float).copy()

    if opt_col in df.columns:
        mask = df[opt_col].notna()
        out.loc[mask] = df.loc[mask, opt_col].astype(float)

    return out


def weighted_sum(df_or_series, weights):
    return df_or_series.mul(weights, axis=0).sum()


def get_co2_results(n, weights):
    """
    Reconstruct the operational generator contribution to direct CO2
    and inspect Store / StorageUnit terms relevant to PyPSA's
    primary_energy global constraint.
    """

    if "co2_emissions" not in n.carriers.columns:
        co2 = pd.Series(0.0, index=n.carriers.index)
    else:
        co2 = n.carriers["co2_emissions"].fillna(0.0).astype(float)

    # --------------------------------------------------------
    # Generator contribution
    # --------------------------------------------------------

    gen_factor = n.generators["carrier"].map(co2).fillna(0.0).astype(float)

    emitting_gens = n.generators.index[gen_factor != 0]

    generator_co2 = 0.0
    emitting_generation_mwh = 0.0
    max_emitting_dispatch_mw = 0.0

    if len(emitting_gens):
        p = n.generators_t.p[emitting_gens].astype(float)

        # Positive generator dispatch only.
        p_positive = p.clip(lower=0.0)

        efficiency = n.get_switchable_as_dense("Generator", "efficiency")[
            emitting_gens
        ].astype(float)

        factors = gen_factor.loc[emitting_gens]

        emissions = p_positive.div(efficiency).mul(factors, axis=1).mul(weights, axis=0)

        generator_co2 = float(emissions.sum().sum())

        positive_factor_gens = gen_factor[gen_factor > 0].index.intersection(
            n.generators.index
        )

        if len(positive_factor_gens):
            p_emit = (
                n.generators_t.p[positive_factor_gens].astype(float).clip(lower=0.0)
            )

            emitting_generation_mwh = float(p_emit.mul(weights, axis=0).sum().sum())

            max_emitting_dispatch_mw = float(p_emit.max().max())

    # --------------------------------------------------------
    # Store contribution
    # --------------------------------------------------------

    store_co2 = 0.0
    emitting_stores = []

    if len(n.stores):
        store_factor = n.stores["carrier"].map(co2).fillna(0.0).astype(float)

        for name in n.stores.index[store_factor != 0]:
            emitting_stores.append(name)

            cyclic = bool(n.stores.at[name, "e_cyclic"])

            if cyclic:
                continue

            initial = float(
                n.stores.at[name, "e_initial"]
                if pd.notna(n.stores.at[name, "e_initial"])
                else 0.0
            )

            final = float(n.stores_t.e[name].iloc[-1])

            store_co2 += (initial - final) * float(store_factor.loc[name])

    # --------------------------------------------------------
    # StorageUnit contribution
    # --------------------------------------------------------

    storage_unit_co2 = 0.0
    emitting_storage_units = []

    if len(n.storage_units):
        su_factor = n.storage_units["carrier"].map(co2).fillna(0.0).astype(float)

        for name in n.storage_units.index[su_factor != 0]:
            emitting_storage_units.append(name)

            cyclic = bool(n.storage_units.at[name, "cyclic_state_of_charge"])

            if cyclic:
                continue

            initial = n.storage_units.at[name, "state_of_charge_initial"]

            initial = float(initial) if pd.notna(initial) else 0.0

            final = float(n.storage_units_t.state_of_charge[name].iloc[-1])

            storage_unit_co2 += (initial - final) * float(su_factor.loc[name])

    lhs = generator_co2 + store_co2 + storage_unit_co2

    return {
        "direct_generator_co2_t": generator_co2,
        "store_co2_term_t": store_co2,
        "storage_unit_co2_term_t": storage_unit_co2,
        "primary_energy_co2_lhs_t": lhs,
        "emitting_generation_mwh": emitting_generation_mwh,
        "max_emitting_dispatch_mw": max_emitting_dispatch_mw,
        "emitting_stores": "; ".join(emitting_stores),
        "emitting_storage_units": "; ".join(emitting_storage_units),
    }


def generation_and_capacity(n, weights, scenario):
    """
    Extract conventional/electricity-supply Generator capacities
    and annual generation by carrier.

    sign < 0 generators such as BTC service are intentionally
    excluded.
    """

    supply_gens = n.generators.index[
        (n.generators["sign"].astype(float) >= 0)
        & (n.generators["carrier"].astype(str).str.lower() != "load shedding")
    ]

    cap = nominal_capacity(
        n.generators,
        "p_nom",
        "p_nom_opt",
    ).loc[supply_gens]

    gen = n.generators_t.p[supply_gens].astype(float).mul(weights, axis=0).sum()

    carrier = n.generators.loc[supply_gens, "carrier"].astype(str)

    tmp = pd.DataFrame(
        {
            "carrier": carrier,
            "capacity_mw": cap,
            "generation_mwh": gen,
        }
    )

    capacity = tmp.groupby("carrier")["capacity_mw"].sum().reset_index()

    capacity.insert(0, "scenario", scenario)

    generation = tmp.groupby("carrier")["generation_mwh"].sum().reset_index()

    generation.insert(0, "scenario", scenario)

    return capacity, generation


def renewable_curtailment(n, weights):
    results = {}

    cap = nominal_capacity(
        n.generators,
        "p_nom",
        "p_nom_opt",
    )

    p_max_pu = n.get_switchable_as_dense(
        "Generator",
        "p_max_pu",
    ).astype(float)

    for carrier in ["solar", "onwind"]:
        gens = n.generators.index[
            (n.generators["carrier"] == carrier)
            & (n.generators["sign"].astype(float) >= 0)
        ]

        if not len(gens):
            results[f"{carrier}_potential_twh"] = 0.0
            results[f"{carrier}_generation_twh"] = 0.0
            results[f"{carrier}_curtailment_twh"] = 0.0
            results[f"{carrier}_curtailment_pct"] = np.nan
            continue

        potential = p_max_pu[gens].mul(cap.loc[gens], axis=1)

        actual = n.generators_t.p[gens].astype(float).clip(lower=0.0)

        potential_mwh = float(potential.mul(weights, axis=0).sum().sum())

        actual_mwh = float(actual.mul(weights, axis=0).sum().sum())

        curtailed_mwh = max(
            potential_mwh - actual_mwh,
            0.0,
        )

        results[f"{carrier}_potential_twh"] = potential_mwh / 1e6

        results[f"{carrier}_generation_twh"] = actual_mwh / 1e6

        results[f"{carrier}_curtailment_twh"] = curtailed_mwh / 1e6

        results[f"{carrier}_curtailment_pct"] = (
            100.0 * curtailed_mwh / potential_mwh if potential_mwh > 0 else np.nan
        )

    return results


def storage_results(n, weights, scenario):
    rows = []

    if not len(n.storage_units):
        return pd.DataFrame(
            columns=[
                "scenario",
                "carrier",
                "power_mw",
                "energy_mwh",
                "charge_mwh",
                "discharge_mwh",
            ]
        )

    cap = nominal_capacity(
        n.storage_units,
        "p_nom",
        "p_nom_opt",
    )

    max_hours = n.storage_units["max_hours"].fillna(0.0).astype(float)

    energy = cap * max_hours

    p = n.storage_units_t.p.astype(float)

    charge = (-p.clip(upper=0.0)).mul(weights, axis=0).sum()

    discharge = p.clip(lower=0.0).mul(weights, axis=0).sum()

    tmp = pd.DataFrame(
        {
            "carrier": n.storage_units["carrier"].astype(str),
            "power_mw": cap,
            "energy_mwh": energy,
            "charge_mwh": charge,
            "discharge_mwh": discharge,
        }
    )

    grouped = tmp.groupby("carrier").sum(numeric_only=True).reset_index()

    grouped.insert(0, "scenario", scenario)

    return grouped


def battery_results(n, weights, scenario):
    """
    Extract PyPSA-Earth battery results.

    Batteries are represented as:
      - Store: battery energy capacity
      - Link: battery charger
      - Link: battery discharger

    They are NOT StorageUnits in this model.
    """

    store_mask = n.stores["carrier"].astype(str).str.lower() == "battery"

    charger_mask = n.links["carrier"].astype(str).str.lower() == "battery charger"

    discharger_mask = n.links["carrier"].astype(str).str.lower() == "battery discharger"

    stores = n.stores.index[store_mask]
    chargers = n.links.index[charger_mask]
    dischargers = n.links.index[discharger_mask]

    if not len(stores):
        return pd.DataFrame(
            [
                {
                    "scenario": scenario,
                    "battery_store_count": 0,
                    "battery_energy_mwh": 0.0,
                    "battery_charge_power_mw": 0.0,
                    "battery_discharge_link_p_nom_mw": 0.0,
                    "battery_discharge_power_mw": 0.0,
                    "battery_charge_mwh": 0.0,
                    "battery_discharge_mwh": 0.0,
                }
            ]
        )

    # --------------------------------------------------------
    # Energy capacity
    # --------------------------------------------------------

    store_energy = nominal_capacity(
        n.stores,
        "e_nom",
        "e_nom_opt",
    )

    battery_energy_mwh = float(store_energy.loc[stores].sum())

    # --------------------------------------------------------
    # Power capacities
    # --------------------------------------------------------

    link_power = nominal_capacity(
        n.links,
        "p_nom",
        "p_nom_opt",
    )

    battery_charge_power_mw = float(link_power.loc[chargers].sum())

    battery_discharge_link_p_nom_mw = float(link_power.loc[dischargers].sum())

    # Discharger p_nom is defined at bus0 (battery side).
    # Convert to electricity delivered at bus1.
    discharge_efficiency = n.links.loc[dischargers, "efficiency"].astype(float)

    battery_discharge_power_mw = float(
        (link_power.loc[dischargers] * discharge_efficiency).sum()
    )

    # --------------------------------------------------------
    # Annual charging / discharging
    # --------------------------------------------------------

    if len(chargers):
        charge = n.links_t.p0[chargers].astype(float).clip(lower=0.0)

        battery_charge_mwh = float(charge.mul(weights, axis=0).sum().sum())
    else:
        battery_charge_mwh = 0.0

    if len(dischargers):
        # Electricity delivered to the AC buses is -p1.
        discharge = (-n.links_t.p1[dischargers].astype(float)).clip(lower=0.0)

        battery_discharge_mwh = float(discharge.mul(weights, axis=0).sum().sum())
    else:
        battery_discharge_mwh = 0.0

    return pd.DataFrame(
        [
            {
                "scenario": scenario,
                "battery_store_count": len(stores),
                "battery_energy_mwh": battery_energy_mwh,
                "battery_charge_power_mw": battery_charge_power_mw,
                "battery_discharge_link_p_nom_mw": battery_discharge_link_p_nom_mw,
                "battery_discharge_power_mw": battery_discharge_power_mw,
                "battery_charge_mwh": battery_charge_mwh,
                "battery_discharge_mwh": battery_discharge_mwh,
            }
        ]
    )


def transmission_results(n):
    if not len(n.lines):
        return {
            "ac_line_volume_gwkm": 0.0,
            "ac_line_expansion_gwkm": 0.0,
        }

    cap = nominal_capacity(
        n.lines,
        "s_nom",
        "s_nom_opt",
    )

    base = n.lines["s_nom"].astype(float)

    length = n.lines["length"].fillna(0.0).astype(float)

    volume = float((cap * length).sum()) / 1000.0

    expansion = float(((cap - base).clip(lower=0.0) * length).sum()) / 1000.0

    return {
        "ac_line_volume_gwkm": volume,
        "ac_line_expansion_gwkm": expansion,
    }


# ============================================================
# MAIN AUDIT / EXTRACTION
# ============================================================

headline_rows = []
capacity_tables = []
generation_tables = []
storage_tables = []
battery_tables = []
audit_rows = []

print("=" * 78)
print("FINAL S0-S6 THESIS AUDIT")
print("=" * 78)

for scenario, cfg in SCENARIOS.items():

    fn = cfg["file"]

    print(f"\n{'=' * 30} {scenario} {'=' * 30}")
    print(fn)

    assert fn.exists(), f"Missing result network: {fn}"

    n = pypsa.Network(fn)
    meta = n.meta or {}

    weights = n.snapshot_weightings["generators"].reindex(n.snapshots).astype(float)

    # --------------------------------------------------------
    # General invariants
    # --------------------------------------------------------

    assert len(n.snapshots) == 8760

    opts = meta.get("wildcards", {}).get("opts")

    assert opts == cfg["opts"], (
        f"{scenario}: opts mismatch: " f"{opts} != {cfg['opts']}"
    )

    assert meta.get("solving", {}).get("options", {}).get("noisy_costs") is False

    assert meta.get("sector", {}).get("chp") is False

    p_set = n.get_switchable_as_dense(
        "Load",
        "p_set",
    ).astype(float)

    base_demand_mwh = float(p_set.mul(weights, axis=0).sum().sum())

    base_demand_twh = base_demand_mwh / 1e6

    assert abs(base_demand_twh - EXPECTED_BASE_DEMAND_TWH) < 1e-6, (
        f"{scenario}: wrong base demand " f"{base_demand_twh} TWh"
    )

    # --------------------------------------------------------
    # Scenario isolation
    # --------------------------------------------------------

    btc_present = (
        BTC_BUS in n.buses.index
        and BTC_LINK in n.links.index
        and BTC_SERVICE in n.generators.index
    )

    pem_present = (
        PEM_BUS in n.buses.index
        and PEM_LINK in n.links.index
        and PEM_STORE in n.stores.index
    )

    assert btc_present == cfg["btc"], f"{scenario}: BTC presence mismatch"

    assert pem_present == cfg["pem"], f"{scenario}: PEM presence mismatch"

    # --------------------------------------------------------
    # CO2 constraint
    # --------------------------------------------------------

    co2_constraint_present = "CO2Limit" in n.global_constraints.index

    assert co2_constraint_present == cfg["co2"], f"{scenario}: CO2 constraint mismatch"

    if cfg["co2"]:
        gc = n.global_constraints.loc["CO2Limit"]

        assert gc["sense"] == "<="
        assert gc["type"] == "primary_energy"
        assert gc["carrier_attribute"] == "co2_emissions"
        assert abs(float(gc["constant"])) < 1e-12

    # --------------------------------------------------------
    # BTC
    # --------------------------------------------------------

    btc_mwh = 0.0
    btc_twh = 0.0
    btc_flh = 0.0
    btc_cf = 0.0
    btc_service_value = 0.0

    if cfg["btc"]:

        l = n.links.loc[BTC_LINK]
        g = n.generators.loc[BTC_SERVICE]

        assert l.bus0 == CONNECTION_BUS
        assert l.bus1 == BTC_BUS

        assert abs(float(l.p_nom) - EXPECTED_BTC_MW) < 1e-9

        assert not bool(l.p_nom_extendable)

        assert abs(float(l.efficiency) - 1.0) < 1e-12

        assert abs(float(g.sign) + 1.0) < 1e-12

        assert abs(float(g.marginal_cost) + EXPECTED_BTC_VALUE) < 1e-9

        p_btc = n.links_t.p0[BTC_LINK].astype(float)

        p_service = n.generators_t.p[BTC_SERVICE].astype(float)

        assert np.allclose(
            p_btc.values,
            p_service.values,
            atol=1e-4,
        )

        assert p_btc.min() >= -1e-4
        assert p_btc.max() <= 1000.0 + 1e-4

        btc_mwh = float((p_btc * weights).sum())

        btc_twh = btc_mwh / 1e6

        btc_flh = btc_mwh / EXPECTED_BTC_MW

        btc_cf = btc_flh / float(weights.sum())

        btc_service_value = btc_mwh * EXPECTED_BTC_VALUE

    # --------------------------------------------------------
    # PEM
    # --------------------------------------------------------

    pem_el_mwh = 0.0
    pem_el_twh = 0.0
    h2_mwh = 0.0
    h2_kt = 0.0
    pem_flh = 0.0
    pem_cf = 0.0
    pem_variable_cost = 0.0

    if cfg["pem"]:

        l = n.links.loc[PEM_LINK]
        s = n.stores.loc[PEM_STORE]

        assert l.bus0 == CONNECTION_BUS
        assert l.bus1 == PEM_BUS

        assert abs(float(l.p_nom) - EXPECTED_PEM_MW) < 1e-9

        assert not bool(l.p_nom_extendable)

        assert abs(float(l.efficiency) - EXPECTED_PEM_EFF) < 1e-12

        assert abs(float(l.marginal_cost) - EXPECTED_PEM_MC) < 1e-9

        assert abs(float(s.e_nom) - EXPECTED_H2_MWH) < 1e-6

        p_pem = n.links_t.p0[PEM_LINK].astype(float)

        pem_el_mwh = float((p_pem * weights).sum())

        pem_el_twh = pem_el_mwh / 1e6

        h2_mwh = pem_el_mwh * EXPECTED_PEM_EFF

        h2_kt = h2_mwh * 1000.0 / 33.33 / 1e6

        pem_flh = pem_el_mwh / EXPECTED_PEM_MW

        pem_cf = pem_flh / float(weights.sum())

        pem_variable_cost = pem_el_mwh * EXPECTED_PEM_MC

        final_h2 = float(n.stores_t.e[PEM_STORE].iloc[-1])

        assert abs(pem_el_mwh - EXPECTED_PEM_EL_MWH) < 1.0

        assert abs(h2_mwh - EXPECTED_H2_MWH) < 1.0

        assert abs(h2_kt - EXPECTED_H2_KT) < 1e-6

        assert abs(final_h2 - EXPECTED_H2_MWH) < 1.0

    # --------------------------------------------------------
    # CO2 reconstruction
    # --------------------------------------------------------

    co2_result = get_co2_results(
        n,
        weights,
    )

    if cfg["co2"]:

        assert (
            co2_result["direct_generator_co2_t"] < 1.0
        ), f"{scenario}: generator CO2 is not zero"

        assert (
            co2_result["emitting_generation_mwh"] < 1.0
        ), f"{scenario}: emitting generation is not zero"

        assert co2_result["primary_energy_co2_lhs_t"] <= 1.0, (
            f"{scenario}: reconstructed " f"primary-energy CO2 LHS > tolerance"
        )

    # --------------------------------------------------------
    # Corrected system-cost accounting
    # --------------------------------------------------------

    raw_objective = float(n.objective)

    upstream_system_cost = raw_objective + btc_service_value - pem_variable_cost

    total_electricity_withdrawal_twh = base_demand_twh + btc_twh + pem_el_twh

    # --------------------------------------------------------
    # Detailed electricity results
    # --------------------------------------------------------

    cap_table, gen_table = generation_and_capacity(
        n,
        weights,
        scenario,
    )

    capacity_tables.append(cap_table)
    generation_tables.append(gen_table)

    storage_table = storage_results(
        n,
        weights,
        scenario,
    )

    storage_tables.append(storage_table)

    curtailment = renewable_curtailment(
        n,
        weights,
    )

    transmission = transmission_results(n)

    # --------------------------------------------------------
    # Battery results: Store + charger/discharger Links
    # --------------------------------------------------------

    battery_table = battery_results(
        n,
        weights,
        scenario,
    )

    battery_tables.append(battery_table)

    battery = battery_table.iloc[0]

    battery_energy_mwh = float(battery["battery_energy_mwh"])

    battery_charge_power_mw = float(battery["battery_charge_power_mw"])

    battery_discharge_power_mw = float(battery["battery_discharge_power_mw"])

    battery_charge_mwh = float(battery["battery_charge_mwh"])

    battery_discharge_mwh = float(battery["battery_discharge_mwh"])

    # --------------------------------------------------------
    # Electricity price at common connection bus
    # --------------------------------------------------------

    price = n.buses_t.marginal_price[CONNECTION_BUS].astype(float)

    mean_price = float((price * weights).sum() / weights.sum())

    median_price = float(price.median())

    # --------------------------------------------------------
    # Store headline row
    # --------------------------------------------------------

    headline_rows.append(
        {
            "scenario": scenario,
            "regime": cfg["regime"],
            "base_demand_twh": base_demand_twh,
            "btc_twh": btc_twh,
            "btc_flh_h": btc_flh,
            "btc_cf_pct": btc_cf * 100.0,
            "pem_el_twh": pem_el_twh,
            "h2_twh_lhv": h2_mwh / 1e6,
            "h2_kt": h2_kt,
            "pem_flh_h": pem_flh,
            "pem_cf_pct": pem_cf * 100.0,
            "total_electricity_withdrawal_twh": total_electricity_withdrawal_twh,
            "raw_objective_eur": raw_objective,
            "btc_service_value_eur": btc_service_value,
            "pem_variable_cost_eur": pem_variable_cost,
            "upstream_system_cost_eur": upstream_system_cost,
            "direct_generator_co2_t": co2_result["direct_generator_co2_t"],
            "primary_energy_co2_lhs_t": co2_result["primary_energy_co2_lhs_t"],
            "emitting_generation_twh": co2_result["emitting_generation_mwh"] / 1e6,
            "solar_generation_twh": curtailment["solar_generation_twh"],
            "solar_curtailment_twh": curtailment["solar_curtailment_twh"],
            "solar_curtailment_pct": curtailment["solar_curtailment_pct"],
            "onwind_generation_twh": curtailment["onwind_generation_twh"],
            "onwind_curtailment_twh": curtailment["onwind_curtailment_twh"],
            "onwind_curtailment_pct": curtailment["onwind_curtailment_pct"],
            "battery_charge_power_gw": battery_charge_power_mw / 1000.0,
            "battery_discharge_power_gw": battery_discharge_power_mw / 1000.0,
            "battery_energy_gwh": battery_energy_mwh / 1000.0,
            "battery_charge_twh": battery_charge_mwh / 1e6,
            "battery_discharge_twh": battery_discharge_mwh / 1e6,
            "ac_line_volume_gwkm": transmission["ac_line_volume_gwkm"],
            "ac_line_expansion_gwkm": transmission["ac_line_expansion_gwkm"],
            "mean_price_kz0_eur_per_mwh": mean_price,
            "median_price_kz0_eur_per_mwh": median_price,
        }
    )

    audit_rows.append(
        {
            "scenario": scenario,
            "opts": opts,
            "snapshots": len(n.snapshots),
            "base_demand_twh": base_demand_twh,
            "co2_constraint": co2_constraint_present,
            "btc_present": btc_present,
            "pem_present": pem_present,
            "noisy_costs": False,
            "sector_chp": False,
            "status": "PASS",
        }
    )

    print(
        f"PASS | demand={base_demand_twh:.3f} TWh "
        f"| BTC={btc_twh:.6f} TWh "
        f"| PEM={pem_el_twh:.6f} TWh "
        f"| H2={h2_kt:.6f} kt "
        f"| CO2={co2_result['direct_generator_co2_t']:.3f} t "
        f"| Csys={upstream_system_cost / 1e9:.6f} bn EUR/a"
    )


# ============================================================
# COMBINE TABLES
# ============================================================

headline = pd.DataFrame(headline_rows)
audit = pd.DataFrame(audit_rows)

capacity = pd.concat(
    capacity_tables,
    ignore_index=True,
)

generation = pd.concat(
    generation_tables,
    ignore_index=True,
)

storage = pd.concat(
    storage_tables,
    ignore_index=True,
)

battery = pd.concat(
    battery_tables,
    ignore_index=True,
)


# ============================================================
# SCENARIO COMPARISONS
# ============================================================

comparison_pairs = [
    ("S2", "S0"),
    ("S3", "S0"),
    ("S4", "S1"),
    ("S5", "S1"),
    ("S6", "S1"),
    ("S6", "S4"),
]

numeric_cols = headline.select_dtypes(include=[np.number]).columns

comparison_rows = []

h = headline.set_index("scenario")

for case, baseline in comparison_pairs:
    for metric in numeric_cols:

        case_value = float(h.loc[case, metric])

        baseline_value = float(h.loc[baseline, metric])

        delta = case_value - baseline_value

        if np.isfinite(baseline_value) and abs(baseline_value) > 1e-12:
            delta_pct = 100.0 * delta / baseline_value
        else:
            delta_pct = np.nan

        comparison_rows.append(
            {
                "comparison": f"{case}-{baseline}",
                "scenario": case,
                "baseline": baseline,
                "metric": metric,
                "scenario_value": case_value,
                "baseline_value": baseline_value,
                "delta": delta,
                "delta_pct": delta_pct,
            }
        )

comparisons = pd.DataFrame(comparison_rows)


# ============================================================
# S6 INTERACTION EFFECTS
# I_X = X_S6 - X_S4 - X_S5 + X_S1
# ============================================================

interaction_rows = []

for metric in numeric_cols:

    values = [h.loc[s, metric] for s in ["S1", "S4", "S5", "S6"]]

    if not all(np.isfinite(values)):
        interaction = np.nan
    else:
        interaction = (
            float(h.loc["S6", metric])
            - float(h.loc["S4", metric])
            - float(h.loc["S5", metric])
            + float(h.loc["S1", metric])
        )

    interaction_rows.append(
        {
            "metric": metric,
            "interaction_S6_minus_S4_minus_S5_plus_S1": interaction,
        }
    )

interactions = pd.DataFrame(interaction_rows)


# ============================================================
# WRITE FINAL CSV FILES
# ============================================================

audit.to_csv(
    OUT / "table_00_final_scenario_audit.csv",
    index=False,
)

headline.to_csv(
    OUT / "table_01_final_headline_kpis.csv",
    index=False,
)

capacity.to_csv(
    OUT / "table_02_capacity_by_carrier.csv",
    index=False,
)

generation.to_csv(
    OUT / "table_03_generation_by_carrier.csv",
    index=False,
)

storage.to_csv(
    OUT / "table_04_storage_by_carrier.csv",
    index=False,
)

battery.to_csv(
    OUT / "table_04b_battery_results.csv",
    index=False,
)

comparisons.to_csv(
    OUT / "table_05_scenario_comparisons.csv",
    index=False,
)

interactions.to_csv(
    OUT / "table_06_interaction_effects.csv",
    index=False,
)


# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

print("\n" + "=" * 78)
print("FINAL AUDIT SUMMARY")
print("=" * 78)

print(
    audit[
        [
            "scenario",
            "base_demand_twh",
            "co2_constraint",
            "btc_present",
            "pem_present",
            "status",
        ]
    ].to_string(index=False)
)

print("\n" + "=" * 78)
print("HEADLINE KPI SUMMARY")
print("=" * 78)

display_cols = [
    "scenario",
    "upstream_system_cost_eur",
    "direct_generator_co2_t",
    "btc_twh",
    "btc_cf_pct",
    "pem_el_twh",
    "h2_kt",
    "pem_cf_pct",
    "solar_curtailment_twh",
    "onwind_curtailment_twh",
    "battery_charge_power_gw",
    "battery_discharge_power_gw",
    "battery_energy_gwh",
    "ac_line_expansion_gwkm",
    "mean_price_kz0_eur_per_mwh",
]

print(
    headline[display_cols].to_string(
        index=False,
        float_format=lambda x: f"{x:,.6f}",
    )
)

print("\n" + "=" * 78)
print("FILES WRITTEN")
print("=" * 78)

for f in sorted(OUT.glob("*.csv")):
    print(f)

print("\nPASS: all final S0-S6 thesis scenarios audited and extracted.")
