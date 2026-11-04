from pathlib import Path

import pandas as pd
import pypsa
from thesis_flexible_consumers import summary

SCENARIOS = {
    "S0": ("results/scenarios/final_S0_S6_Link_20260922/" "S0_KZ_2045_Reference.nc"),
    "S1": (
        "results/scenarios/final_S0_S6_Link_20260922/" "S1_KZ_2045_ZeroDirectCO2.nc"
    ),
    "S2": (
        "results/scenarios/final_S0_S6_Link_20260922/" "S2_KZ_2045_Reference_BTC.nc"
    ),
    "S3": ("results/scenarios/final_S0_S6_Link_20260922/" "S3_KZ_2045_Reference_H2.nc"),
    "S4": (
        "results/scenarios/final_S0_S6_Link_20260922/" "S4_KZ_2045_ZeroDirectCO2_BTC.nc"
    ),
    "S5": (
        "results/scenarios/final_S0_S6_Link_20260922/" "S5_KZ_2045_ZeroDirectCO2_H2.nc"
    ),
    "S6": (
        "results/scenarios/final_S0_S6_Link_20260922/"
        "S6_KZ_2045_ZeroDirectCO2_BTC_H2.nc"
    ),
}


rows = []

for scenario, path_string in SCENARIOS.items():

    path = Path(path_string)

    if not path.exists():
        raise FileNotFoundError(f"{scenario}: {path}")

    print(f"Loading {scenario}: " f"{path}")

    n = pypsa.Network(path)

    row = {
        "scenario": scenario,
        "network": str(path),
    }

    row.update(summary(n))

    rows.append(row)


df = pd.DataFrame(rows)

columns = [
    "scenario",
    "btc_representation",
    "btc_capacity_mw",
    "btc_electricity_twh",
    "btc_capacity_factor_pct",
    "pem_present",
    "pem_electricity_twh",
    "h2_product_twh_lhv",
    "h2_product_kt",
    "flexible_electricity_twh",
    "raw_objective_eur",
    "btc_objective_contribution_eur",
    "upstream_power_system_cost_eur",
    "network",
]

df = df[columns]


print("\n" "============================================================")

print("FINAL FLEXIBLE-CONSUMER EXTRACTION CHECK")

print("============================================================")

print(
    df.to_string(
        index=False,
        float_format=lambda x: f"{x:.9f}",
    )
)


output = Path("results/scenarios/" "final_flexible_consumer_extraction_check.csv")

output.parent.mkdir(
    parents=True,
    exist_ok=True,
)

df.to_csv(
    output,
    index=False,
)


# =============================================================================
# Thesis invariants
# =============================================================================

by_scenario = df.set_index("scenario")


# BTC presence
assert (
    by_scenario.loc[
        ["S0", "S1", "S3", "S5"],
        "btc_capacity_mw",
    ]
    == 0.0
).all()

assert (
    by_scenario.loc[
        ["S2", "S4", "S6"],
        "btc_capacity_mw",
    ]
    == 1000.0
).all()


# Final BTC scenarios must use Link representation.
assert (
    by_scenario.loc[
        ["S2", "S4", "S6"],
        "btc_representation",
    ]
    == "link"
).all()


# PEM presence
assert (
    ~by_scenario.loc[
        ["S0", "S1", "S2", "S4"],
        "pem_present",
    ]
).all()

assert (
    by_scenario.loc[
        ["S3", "S5", "S6"],
        "pem_present",
    ]
).all()


# H2 target
for scenario in [
    "S3",
    "S5",
    "S6",
]:
    h2_kt = float(
        by_scenario.loc[
            scenario,
            "h2_product_kt",
        ]
    )

    if abs(h2_kt - 100.0) > 1e-6:

        raise AssertionError(
            f"{scenario}: " f"H2 target is {h2_kt} kt/a, " "expected 100 kt/a."
        )


print("\n" "============================================================")

print("ALL FLEXIBLE-CONSUMER EXTRACTION " "INVARIANTS PASSED")

print("============================================================")

print("\nSaved:")

print(output)
