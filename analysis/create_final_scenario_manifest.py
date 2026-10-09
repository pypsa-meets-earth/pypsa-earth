import hashlib
from pathlib import Path

import pandas as pd
import pypsa
from thesis_flexible_consumers import (
    btc_capacity_mw,
    btc_consumption_twh,
    btc_representation,
    h2_product_kt,
    pem_electricity_twh,
)

ROOT = Path("results/scenarios/final_S0_S6_Link_20260922")

SCENARIOS = {
    "S0": ("S0_KZ_2045_Reference.nc", "reference", False, False),
    "S1": ("S1_KZ_2045_ZeroDirectCO2.nc", "zero-direct-CO2", False, False),
    "S2": ("S2_KZ_2045_Reference_BTC.nc", "reference", True, False),
    "S3": ("S3_KZ_2045_Reference_H2.nc", "reference", False, True),
    "S4": ("S4_KZ_2045_ZeroDirectCO2_BTC.nc", "zero-direct-CO2", True, False),
    "S5": ("S5_KZ_2045_ZeroDirectCO2_H2.nc", "zero-direct-CO2", False, True),
    "S6": ("S6_KZ_2045_ZeroDirectCO2_BTC_H2.nc", "zero-direct-CO2", True, True),
}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def weights(n):
    if isinstance(n.snapshot_weightings, pd.DataFrame):
        return n.snapshot_weightings["generators"]
    return n.snapshot_weightings


def annual_base_demand_twh(n):
    w = weights(n)
    return float(n.loads_t.p_set.sum(axis=1).mul(w).sum() / 1e6)


def direct_co2_mt(n):
    w = weights(n)
    total_t = 0.0

    for name, row in n.generators.iterrows():
        carrier = str(row.carrier)

        if carrier not in n.carriers.index:
            continue

        if "co2_emissions" not in n.carriers.columns:
            continue

        intensity = float(n.carriers.at[carrier, "co2_emissions"])

        efficiency = float(row.efficiency)

        if intensity == 0 or efficiency <= 0:
            continue

        generation_mwh = float(n.generators_t.p[name].clip(lower=0).mul(w).sum())

        total_t += generation_mwh / efficiency * intensity

    return total_t / 1e6


rows = []

for scenario, (
    filename,
    regime,
    btc_expected,
    pem_expected,
) in SCENARIOS.items():

    path = ROOT / filename

    if not path.exists():
        raise FileNotFoundError(path)

    print(f"Checking {scenario}: {path}")

    n = pypsa.Network(path)

    rows.append(
        {
            "scenario": scenario,
            "system_regime": regime,
            "network_file": str(path),
            "sha256": sha256(path),
            "snapshots": len(n.snapshots),
            "weighted_hours": float(weights(n).sum()),
            "annual_base_demand_twh": annual_base_demand_twh(n),
            "direct_co2_mt_per_a": direct_co2_mt(n),
            "btc_expected": btc_expected,
            "btc_representation": btc_representation(n),
            "btc_capacity_mw": btc_capacity_mw(n),
            "btc_electricity_twh": btc_consumption_twh(n),
            "pem_expected": pem_expected,
            "pem_electricity_twh": pem_electricity_twh(n),
            "h2_production_kt_per_a": h2_product_kt(n),
            "raw_objective_eur": float(n.objective),
        }
    )


df = pd.DataFrame(rows).set_index("scenario")


# ------------------------------------------------------------------
# Global invariants
# ------------------------------------------------------------------

for scenario in df.index:

    if int(df.loc[scenario, "snapshots"]) != 8760:
        raise AssertionError(f"{scenario}: snapshots != 8760")

    if abs(float(df.loc[scenario, "weighted_hours"]) - 8760.0) > 1e-6:
        raise AssertionError(f"{scenario}: weighted hours != 8760")

    if abs(float(df.loc[scenario, "annual_base_demand_twh"]) - 187.0) > 1e-6:
        raise AssertionError(f"{scenario}: base demand != 187 TWh/a")


# ------------------------------------------------------------------
# BTC invariants
# ------------------------------------------------------------------

for scenario in ["S2", "S4", "S6"]:

    if df.loc[scenario, "btc_representation"] != "link":
        raise AssertionError(f"{scenario}: BTC is not represented by Link")

    if abs(float(df.loc[scenario, "btc_capacity_mw"]) - 1000.0) > 1e-9:
        raise AssertionError(f"{scenario}: BTC capacity != 1000 MW")


for scenario in ["S0", "S1", "S3", "S5"]:

    if df.loc[scenario, "btc_representation"] != "none":
        raise AssertionError(f"{scenario}: unexpected BTC component")


# ------------------------------------------------------------------
# H2 invariants
# ------------------------------------------------------------------

for scenario in ["S3", "S5", "S6"]:

    if abs(float(df.loc[scenario, "h2_production_kt_per_a"]) - 100.0) > 1e-6:
        raise AssertionError(f"{scenario}: H2 production != 100 kt/a")


# ------------------------------------------------------------------
# Zero-direct-CO2 invariants
# ------------------------------------------------------------------

for scenario in ["S1", "S4", "S5", "S6"]:

    if abs(float(df.loc[scenario, "direct_co2_mt_per_a"])) > 1e-8:
        raise AssertionError(f"{scenario}: direct CO2 is not zero")


output = ROOT / "final_scenario_manifest.csv"

df.reset_index().to_csv(
    output,
    index=False,
)


print()
print("=" * 110)
print("FINAL S0-S6 ARCHIVE MANIFEST")
print("=" * 110)

print(
    df[
        [
            "system_regime",
            "snapshots",
            "annual_base_demand_twh",
            "direct_co2_mt_per_a",
            "btc_representation",
            "btc_capacity_mw",
            "btc_electricity_twh",
            "pem_electricity_twh",
            "h2_production_kt_per_a",
        ]
    ].to_string(float_format=lambda x: f"{x:.9f}")
)

print()
print("ALL FINAL ARCHIVE INVARIANTS PASSED")
print()
print(f"Saved: {output}")
