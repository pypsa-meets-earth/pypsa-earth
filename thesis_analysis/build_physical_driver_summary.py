from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("results/thesis_2045/final_analysis")

headline = pd.read_csv(BASE / "table_01_final_headline_kpis.csv").set_index("scenario")

capacity = pd.read_csv(BASE / "table_02_capacity_by_carrier.csv")

generation = pd.read_csv(BASE / "table_03_generation_by_carrier.csv")


# ============================================================
# PIVOTS
# ============================================================

cap = (
    capacity.pivot(
        index="scenario",
        columns="carrier",
        values="capacity_mw",
    ).fillna(0.0)
    / 1000.0
)

gen = (
    generation.pivot(
        index="scenario",
        columns="carrier",
        values="generation_mwh",
    ).fillna(0.0)
    / 1e6
)


def get(df, scenario, carrier):
    if carrier not in df.columns:
        return 0.0
    return float(df.loc[scenario, carrier])


# ============================================================
# COMPARISONS
# ============================================================

pairs = [
    ("S2", "S0"),
    ("S3", "S0"),
    ("S4", "S1"),
    ("S5", "S1"),
    ("S6", "S1"),
]

rows = []

for case, baseline in pairs:

    added_load = (
        headline.loc[case, "total_electricity_withdrawal_twh"]
        - headline.loc[baseline, "total_electricity_withdrawal_twh"]
    )

    base_curtailment = (
        headline.loc[baseline, "solar_curtailment_twh"]
        + headline.loc[baseline, "onwind_curtailment_twh"]
    )

    case_curtailment = (
        headline.loc[case, "solar_curtailment_twh"]
        + headline.loc[case, "onwind_curtailment_twh"]
    )

    curtailment_reduction = base_curtailment - case_curtailment

    if abs(added_load) > 1e-12:
        offset_ratio = 100.0 * curtailment_reduction / added_load
    else:
        offset_ratio = np.nan

    fossil_carriers = [
        "coal",
        "CCGT",
        "OCGT",
    ]

    fossil_case = sum(get(gen, case, c) for c in fossil_carriers)

    fossil_base = sum(get(gen, baseline, c) for c in fossil_carriers)

    rows.append(
        {
            "comparison": f"{case}-{baseline}",
            "added_flexible_load_twh": added_load,
            "delta_system_cost_m_eur": (
                headline.loc[case, "upstream_system_cost_eur"]
                - headline.loc[baseline, "upstream_system_cost_eur"]
            )
            / 1e6,
            "delta_direct_co2_mt": (
                headline.loc[case, "direct_generator_co2_t"]
                - headline.loc[baseline, "direct_generator_co2_t"]
            )
            / 1e6,
            "delta_solar_capacity_gw": get(cap, case, "solar")
            - get(cap, baseline, "solar"),
            "delta_onwind_capacity_gw": get(cap, case, "onwind")
            - get(cap, baseline, "onwind"),
            "delta_solar_generation_twh": get(gen, case, "solar")
            - get(gen, baseline, "solar"),
            "delta_onwind_generation_twh": get(gen, case, "onwind")
            - get(gen, baseline, "onwind"),
            "delta_fossil_generation_twh": fossil_case - fossil_base,
            "curtailment_reduction_twh": curtailment_reduction,
            "curtailment_offset_share_pct": offset_ratio,
            "delta_battery_energy_gwh": (
                headline.loc[case, "battery_energy_gwh"]
                - headline.loc[baseline, "battery_energy_gwh"]
            ),
            "delta_battery_power_gw": (
                headline.loc[case, "battery_charge_power_gw"]
                - headline.loc[baseline, "battery_charge_power_gw"]
            ),
            "delta_transmission_expansion_gwkm": (
                headline.loc[case, "ac_line_expansion_gwkm"]
                - headline.loc[baseline, "ac_line_expansion_gwkm"]
            ),
            "delta_mean_price_eur_per_mwh": (
                headline.loc[case, "mean_price_kz0_eur_per_mwh"]
                - headline.loc[baseline, "mean_price_kz0_eur_per_mwh"]
            ),
        }
    )


out = pd.DataFrame(rows)

fn = BASE / "table_07_physical_driver_summary.csv"

out.to_csv(
    fn,
    index=False,
)

print("=" * 130)
print("PHYSICAL DRIVER SUMMARY")
print("=" * 130)

print(
    out.to_string(
        index=False,
        float_format=lambda x: f"{x:,.4f}",
    )
)

print("\nWritten:")
print(fn)
