#!/usr/bin/env python3
"""Build thesis economic tables from existing, optimized PyPSA networks.

Read-only input: results/thesis_economic_2045/networks/*.nc
Outputs only: results/thesis_economic_2045/final_analysis/
No Snakemake or optimization is invoked.

Reproducibility reference: ba-pypsa-earth, commit a89178b9.
"""

import math
from pathlib import Path

import pandas as pd
import pypsa
import yaml

ROOT = Path(__file__).resolve().parents[1]
NETWORK_DIR = ROOT / "results/thesis_economic_2045/networks"
OUT_DIR = ROOT / "results/thesis_economic_2045/final_analysis"
CONFIG = ROOT / "configs/thesis_economic_scenarios.yaml"
BASE_CONFIG = ROOT / "config.yaml"
PREFIX = "elec_s_10_ec_lcopt_"
BTC_NAME = "KZ0 0 Bitcoin mining optimized"
PEM_NAME = "THESIS PEM electrolyser"
H2_STORE = "THESIS H2 annual product accumulator"

ORDER = [
    "E0",
    "E1",
    "E2",
    "E3",
    "E4R35",
    "E4R38",
    "E4R40",
    "E4R100",
    "E5C75",
    "E5C125",
    "E6CAP2000",
]
LABELS = {
    "E0": "Reference",
    "E1": "BTC optimized",
    "E2": "PEM optimized",
    "E3": "BTC + PEM optimized",
    "E4R35": "BTC revenue 35 EUR/MWh",
    "E4R38": "BTC revenue 38 EUR/MWh",
    "E4R40": "BTC revenue 40 EUR/MWh",
    "E4R100": "BTC revenue 100 EUR/MWh",
    "E5C75": "BTC CAPEX x0.75",
    "E5C125": "BTC CAPEX x1.25",
    "E6CAP2000": "BTC cap 2000 MW",
}


def weighted_mwh(series, weights):
    """Integrate a power series [MW] using snapshot duration [h]."""
    return float(series.mul(weights, axis=0).sum())


def read_scenario(key, option, h2_lhv_kwh_per_kg):
    path = NETWORK_DIR / f"{PREFIX}{option}.nc"
    if not path.is_file():
        raise FileNotFoundError(f"Missing optimized network: {path}")

    n = pypsa.Network(path)
    weights = n.snapshot_weightings.generators
    objective = float(n.objective)
    if not math.isfinite(objective):
        raise ValueError(f"Invalid objective in {path}")

    row = {
        "scenario": key,
        "description": LABELS[key],
        "option": option,
        "btc_capacity_MW": 0.0,
        "btc_electricity_TWh": 0.0,
        "btc_full_load_hours": 0.0,
        "btc_revenue_MEUR_per_a": 0.0,
        "btc_annualized_capex_MEUR_per_a": 0.0,
        "btc_revenue_assumption_EUR_per_MWh": "",
        "btc_capacity_bound_MW": "",
        "pem_capacity_MW": 0.0,
        "pem_electricity_TWh": 0.0,
        "pem_full_load_hours": 0.0,
        "h2_production_kt_per_a": 0.0,
        "pem_efficiency_LHV": "",
        "system_objective_MEUR_per_a": objective / 1e6,
    }

    if BTC_NAME in n.generators.index:
        g = n.generators.loc[BTC_NAME]
        capacity = float(g.p_nom_opt)
        electricity_mwh = weighted_mwh(n.generators_t.p[BTC_NAME], weights)
        revenue_rate = -float(g.marginal_cost)
        row.update(
            {
                "btc_capacity_MW": capacity,
                "btc_electricity_TWh": electricity_mwh / 1e6,
                "btc_full_load_hours": (
                    electricity_mwh / capacity if capacity > 0 else 0.0
                ),
                "btc_revenue_MEUR_per_a": electricity_mwh * revenue_rate / 1e6,
                "btc_annualized_capex_MEUR_per_a": capacity
                * float(g.capital_cost)
                / 1e6,
                "btc_revenue_assumption_EUR_per_MWh": revenue_rate,
                "btc_capacity_bound_MW": float(g.p_nom_max),
            }
        )

    if PEM_NAME in n.links.index:
        if H2_STORE not in n.stores.index:
            raise ValueError(f"PEM present but H2 accumulator missing: {key}")
        pem = n.links.loc[PEM_NAME]
        capacity = (
            float(pem.p_nom_opt) if bool(pem.p_nom_extendable) else float(pem.p_nom)
        )
        electricity_mwh = weighted_mwh(n.links_t.p0[PEM_NAME], weights)
        # This dedicated Store represents cumulative H2 product, not a real tank.
        h2_mwh = float(n.stores_t.e[H2_STORE].iloc[-1])
        row.update(
            {
                "pem_capacity_MW": capacity,
                "pem_electricity_TWh": electricity_mwh / 1e6,
                "pem_full_load_hours": (
                    electricity_mwh / capacity if capacity > 0 else 0.0
                ),
                "h2_production_kt_per_a": h2_mwh / h2_lhv_kwh_per_kg / 1000.0,
                "pem_efficiency_LHV": float(pem.efficiency),
            }
        )
    return row


def write_latex(df, path):
    lines = [
        "% Auto-generated from optimized PyPSA networks; do not edit manually.",
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Economic scenario results for the Kazakhstan electricity system (2045).}",
        r"\label{tab:economic_scenario_results}",
        r"\small",
        r"\begin{tabular}{lrrrrrr}",
        r"\hline",
        r"Scenario & BTC (MW) & PEM (MW) & BTC (TWh/a) & H$_2$ (kt/a) & Objective (MEUR/a) & $\Delta$E0 (MEUR/a) \\",
        r"\hline",
    ]
    for _, r in df.iterrows():
        vals = [
            str(r["scenario"]),
            f'{r["btc_capacity_MW"]:,.1f}'.replace(",", "\\,"),
            f'{r["pem_capacity_MW"]:,.1f}'.replace(",", "\\,"),
            f'{r["btc_electricity_TWh"]:.3f}',
            f'{r["h2_production_kt_per_a"]:.2f}',
            f'{r["system_objective_MEUR_per_a"]:,.2f}'.replace(",", "\\,"),
            f'{r["delta_objective_vs_E0_MEUR_per_a"]:+.2f}',
        ]
        lines.append(" & ".join(vals) + r" \\")
    lines += [
        r"\hline",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    with CONFIG.open(encoding="utf-8") as fh:
        config = yaml.safe_load(fh)
    with BASE_CONFIG.open(encoding="utf-8") as fh:
        base = yaml.safe_load(fh)
    h2_lhv = float(
        base["electricity"]["flexible_consumers"]["pem"]["h2_lhv_kwh_per_kg"]
    )
    options = {opt.split("-", 1)[0]: opt for opt in config["scenario"]["opts"]}
    missing = sorted(set(ORDER) - set(options))
    if missing:
        raise ValueError(f"Missing options in economic scenario config: {missing}")
    if len(options) != len(ORDER):
        raise ValueError("Unexpected or duplicate scenarios in configuration.")

    df = pd.DataFrame([read_scenario(key, options[key], h2_lhv) for key in ORDER])
    obj0 = float(df.loc[df.scenario == "E0", "system_objective_MEUR_per_a"].iloc[0])
    obj1 = float(df.loc[df.scenario == "E1", "system_objective_MEUR_per_a"].iloc[0])
    df["delta_objective_vs_E0_MEUR_per_a"] = df["system_objective_MEUR_per_a"] - obj0
    df["delta_objective_vs_E1_MEUR_per_a"] = df["system_objective_MEUR_per_a"] - obj1

    effects = df[
        [
            "scenario",
            "description",
            "btc_revenue_assumption_EUR_per_MWh",
            "btc_capacity_bound_MW",
            "btc_capacity_MW",
            "btc_electricity_TWh",
            "btc_revenue_MEUR_per_a",
            "btc_annualized_capex_MEUR_per_a",
            "system_objective_MEUR_per_a",
            "delta_objective_vs_E0_MEUR_per_a",
            "delta_objective_vs_E1_MEUR_per_a",
        ]
    ].copy()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    file1 = OUT_DIR / "table_01_economic_scenario_kpis.csv"
    file2 = OUT_DIR / "table_02_economic_sensitivity_effects.csv"
    file3 = OUT_DIR / "table_economic_scenario_results.tex"
    df.to_csv(file1, index=False, float_format="%.8f")
    effects.to_csv(file2, index=False, float_format="%.8f")
    write_latex(df, file3)

    by_scenario = df.set_index("scenario")
    interaction = (
        by_scenario.at["E3", "system_objective_MEUR_per_a"]
        - by_scenario.at["E1", "system_objective_MEUR_per_a"]
        - by_scenario.at["E2", "system_objective_MEUR_per_a"]
        + by_scenario.at["E0", "system_objective_MEUR_per_a"]
    )
    manifest = OUT_DIR / "ECONOMIC_RESULTS_MANIFEST.txt"
    manifest.write_text(
        "THESIS ECONOMIC PYPSA-EARTH RESULTS\n"
        "Code reference: ba-pypsa-earth @ a89178b9\n"
        "Source: optimized network NetCDF files in results/thesis_economic_2045/networks/\n"
        "Generation: python scripts/build_thesis_economic_results.py\n"
        "Scenarios: " + ", ".join(ORDER) + "\n"
        "Files: table_01_economic_scenario_kpis.csv, "
        "table_02_economic_sensitivity_effects.csv, "
        "table_economic_scenario_results.tex\n"
        "Objective: annualized net system objective including modeled BTC revenues; "
        "not gross electricity supply cost.\n"
        "H2 annual volume: final energy of bookkeeping Store, converted "
        f"using configured H2 LHV {h2_lhv:g} kWh/kg.\n"
        "BTC power: negative-sign Generator dispatch; positive consumption in results.\n"
        f"E3 interaction E3-E1-E2+E0: {interaction:.6f} MEUR/a\n",
        encoding="utf-8",
    )
    print(f"Completed {len(df)} economic scenarios; no optimization performed.")
    print(
        df[
            [
                "scenario",
                "btc_capacity_MW",
                "pem_capacity_MW",
                "h2_production_kt_per_a",
                "system_objective_MEUR_per_a",
            ]
        ].to_string(index=False, float_format=lambda x: f"{x:.3f}")
    )
    print(f"Economic interaction E3-E1-E2+E0: {interaction:+.4f} MEUR/a")
    print(f"Files written to: {OUT_DIR}")
    for path in (file1, file2, file3, manifest):
        print(" -", path.name)


if __name__ == "__main__":
    main()
