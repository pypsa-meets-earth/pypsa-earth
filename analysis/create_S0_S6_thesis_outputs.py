from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# Paths
# ============================================================

ROOT = Path("results/scenarios")

KPI_FILE = ROOT / "final_S0_S6_flexible_kpis.csv"
INTERACTION_FILE = ROOT / "final_S0_S6_interaction_terms.csv"
COMPETITION_FILE = ROOT / "final_S4_S6_BTC_H2_competition.csv"

FIG_DIR = ROOT / "thesis_figures"
TABLE_DIR = ROOT / "thesis_tables"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Load final archived results
# ============================================================

df = pd.read_csv(
    KPI_FILE,
    index_col=0,
)

interaction = pd.read_csv(
    INTERACTION_FILE,
    index_col=0,
)

competition = pd.read_csv(
    COMPETITION_FILE,
)


SCENARIOS = [
    "S0",
    "S1",
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
]

LABELS = {
    "S0": "S0\nReference",
    "S1": "S1\nZero direct CO$_2$",
    "S2": "S2\nRef. + BTC",
    "S3": "S3\nRef. + H$_2$",
    "S4": "S4\nZero CO$_2$ + BTC",
    "S5": "S5\nZero CO$_2$ + H$_2$",
    "S6": "S6\nZero CO$_2$ + BTC + H$_2$",
}


# ============================================================
# Plot style
# ============================================================

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.dpi": 120,
        "savefig.dpi": 350,
    }
)


def clean_axis(ax):
    ax.grid(
        axis="y",
        alpha=0.25,
        linewidth=0.7,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig, stem):
    png = FIG_DIR / f"{stem}.png"
    pdf = FIG_DIR / f"{stem}.pdf"

    fig.tight_layout()

    fig.savefig(
        png,
        bbox_inches="tight",
        dpi=350,
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Written: {png}")
    print(f"Written: {pdf}")


def annotate_bars(
    ax,
    bars,
    fmt="{:.2f}",
    rotation=0,
):
    for bar in bars:

        height = bar.get_height()

        if not np.isfinite(height):
            continue

        ax.annotate(
            fmt.format(height),
            xy=(
                bar.get_x() + bar.get_width() / 2,
                height,
            ),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=rotation,
        )


# ============================================================
# TABLE 00 — Scenario matrix
# ============================================================

scenario_matrix = pd.DataFrame(
    {
        "Scenario": SCENARIOS,
        "Electricity system": [
            "2045 Reference",
            "2045 zero-direct-CO2",
            "2045 Reference",
            "2045 Reference",
            "2045 zero-direct-CO2",
            "2045 zero-direct-CO2",
            "2045 zero-direct-CO2",
        ],
        "Bitcoin mining": [
            "None",
            "None",
            "1 GW flexible BTC",
            "None",
            "1 GW flexible BTC",
            "None",
            "1 GW flexible BTC",
        ],
        "Hydrogen production": [
            "None",
            "None",
            "None",
            "1 GW PEM; 100 kt H2/a",
            "None",
            "1 GW PEM; 100 kt H2/a",
            "1 GW PEM; 100 kt H2/a",
        ],
    }
)

scenario_matrix.to_csv(
    TABLE_DIR / "table_00_scenario_matrix.csv",
    index=False,
)

scenario_matrix.to_latex(
    TABLE_DIR / "table_00_scenario_matrix.tex",
    index=False,
    escape=True,
    caption=("Definition of the principal PyPSA-Earth scenarios."),
    label="tab:scenario_matrix",
    position="htbp",
)


# ============================================================
# TABLE 01 — Headline KPIs
# ============================================================

headline = df.loc[
    SCENARIOS,
    [
        "Annual base demand [TWh]",
        "Flexible consumption [TWh]",
        "Power-system cost proxy [EUR bn]",
        "CO2 emissions [MtCO2/a]",
        "Solar capacity [GW]",
        "Wind capacity [GW]",
        "Battery energy capacity [GWh]",
        "Solar curtailment [TWh]",
        "Wind curtailment [TWh]",
        "Load-weighted marginal price [EUR/MWh]",
        "BTC consumption [TWh]",
        "BTC capacity factor [%]",
        "PEM electricity [TWh]",
        "PEM capacity factor [%]",
        "H2 production [kt/a]",
    ],
].copy()

headline.columns = [
    "Base demand [TWh]",
    "Flexible demand [TWh]",
    "Upstream system cost [EUR bn]",
    "CO2 [Mt/a]",
    "Solar [GW]",
    "Wind [GW]",
    "Battery [GWh]",
    "Solar curtail. [TWh]",
    "Wind curtail. [TWh]",
    "Load-weighted price [EUR/MWh]",
    "BTC electricity [TWh]",
    "BTC CF [%]",
    "PEM electricity [TWh]",
    "PEM CF [%]",
    "H2 [kt/a]",
]

headline.index.name = "Scenario"

headline.round(3).to_csv(TABLE_DIR / "table_01_headline_kpis.csv")

headline.round(3).to_latex(
    TABLE_DIR / "table_01_headline_kpis.tex",
    escape=True,
    na_rep="--",
    caption=("Headline results of the principal " "PyPSA-Earth scenarios."),
    label="tab:headline_kpis",
    position="htbp",
)


# ============================================================
# TABLE 02 — Baseline-relative system effects
# ============================================================

BASELINES = {
    "S2": "S0",
    "S3": "S0",
    "S4": "S1",
    "S5": "S1",
    "S6": "S1",
}

effect_rows = []

for scenario, baseline in BASELINES.items():

    flex = df.loc[scenario, "Flexible consumption [TWh]"]

    cost_delta = (
        df.loc[scenario, "Power-system cost proxy [EUR bn]"]
        - df.loc[baseline, "Power-system cost proxy [EUR bn]"]
    )

    co2_delta = (
        df.loc[scenario, "CO2 emissions [MtCO2/a]"]
        - df.loc[baseline, "CO2 emissions [MtCO2/a]"]
    )

    solar_curt_delta = (
        df.loc[scenario, "Solar curtailment [TWh]"]
        - df.loc[baseline, "Solar curtailment [TWh]"]
    )

    wind_curt_delta = (
        df.loc[scenario, "Wind curtailment [TWh]"]
        - df.loc[baseline, "Wind curtailment [TWh]"]
    )

    total_curt_delta = solar_curt_delta + wind_curt_delta

    if flex > 0:
        cost_per_flex_mwh = cost_delta * 1e9 / (flex * 1e6)

        co2_per_flex_mwh = co2_delta * 1e6 / (flex * 1e6)

    else:
        cost_per_flex_mwh = np.nan
        co2_per_flex_mwh = np.nan

    effect_rows.append(
        {
            "Scenario": scenario,
            "Baseline": baseline,
            "Flexible electricity [TWh]": flex,
            "Incremental upstream cost [EUR million]": cost_delta * 1000,
            "Incremental cost [EUR/MWh_flex]": cost_per_flex_mwh,
            "Incremental CO2 [Mt/a]": co2_delta,
            "Incremental CO2 [t/MWh_flex]": co2_per_flex_mwh,
            "Change in solar+wind curtailment [TWh]": total_curt_delta,
            "Change in solar capacity [GW]": (
                df.loc[scenario, "Solar capacity [GW]"]
                - df.loc[baseline, "Solar capacity [GW]"]
            ),
            "Change in wind capacity [GW]": (
                df.loc[scenario, "Wind capacity [GW]"]
                - df.loc[baseline, "Wind capacity [GW]"]
            ),
            "Change in battery energy [GWh]": (
                df.loc[scenario, "Battery energy capacity [GWh]"]
                - df.loc[baseline, "Battery energy capacity [GWh]"]
            ),
        }
    )


effects = pd.DataFrame(effect_rows).set_index("Scenario")

effects.round(4).to_csv(TABLE_DIR / "table_02_baseline_relative_effects.csv")

effects.round(4).to_latex(
    TABLE_DIR / "table_02_baseline_relative_effects.tex",
    escape=True,
    na_rep="--",
    caption=(
        "Incremental system effects of Bitcoin mining "
        "and hydrogen production relative to the "
        "corresponding electricity-system baseline."
    ),
    label="tab:flexible_load_effects",
    position="htbp",
)


# ============================================================
# TABLE 03 — Synergy / conflict summary
# ============================================================

interaction_col = "Interaction S6-(S4+S5-S1)"

cost_interaction_bn = float(
    interaction.loc[
        "Power-system cost proxy [EUR bn]",
        interaction_col,
    ]
)

solar_curt_interaction = float(
    interaction.loc[
        "Solar curtailment [TWh]",
        interaction_col,
    ]
)

wind_curt_interaction = float(
    interaction.loc[
        "Wind curtailment [TWh]",
        interaction_col,
    ]
)

ror_curt_interaction = float(
    interaction.loc[
        "RoR curtailment [TWh]",
        interaction_col,
    ]
)

total_curt_interaction = (
    solar_curt_interaction + wind_curt_interaction + ror_curt_interaction
)


def competition_value(metric):
    row = competition.loc[competition["Metric"] == metric]

    if row.empty:
        return np.nan

    return float(row.iloc[0]["Value"])


btc_consumption_change = competition_value("BTC consumption change S6-S4 [TWh]")

btc_cf_change = competition_value("BTC CF change S6-S4 [percentage points]")

pem_price_change = competition_value("PEM weighted price change S6-S5 [EUR/MWh]")


synergy = pd.DataFrame(
    [
        {
            "Indicator": "Combined cost non-additivity",
            "Value": cost_interaction_bn * 1e6,
            "Unit": "EUR thousand/a",
            "Interpretation": (
                "Effectively neutral; very small " "positive interaction."
            ),
        },
        {
            "Indicator": "Combined renewable-curtailment non-additivity",
            "Value": total_curt_interaction,
            "Unit": "TWh/a",
            "Interpretation": (
                "Effectively neutral; small additional " "curtailment reduction."
            ),
        },
        {
            "Indicator": "BTC electricity change when H2 is added",
            "Value": btc_consumption_change,
            "Unit": "TWh/a",
            "Interpretation": ("Negligible direct competition for " "electricity."),
        },
        {
            "Indicator": "BTC capacity-factor change when H2 is added",
            "Value": btc_cf_change,
            "Unit": "percentage points",
            "Interpretation": ("Negligible operational competition."),
        },
        {
            "Indicator": "PEM weighted-price change when BTC is added",
            "Value": pem_price_change,
            "Unit": "EUR/MWh",
            "Interpretation": ("Negligible increase in marginal " "electricity value."),
        },
    ]
)

synergy.round(6).to_csv(
    TABLE_DIR / "table_03_synergy_conflict_summary.csv",
    index=False,
)

synergy.round(6).to_latex(
    TABLE_DIR / "table_03_synergy_conflict_summary.tex",
    index=False,
    escape=True,
    caption=(
        "Indicators of non-additivity and direct competition "
        "between Bitcoin mining and hydrogen production "
        "in the net-zero scenarios."
    ),
    label="tab:synergy_conflict",
    position="htbp",
)


# ============================================================
# FIGURE 01 — Total upstream power-system cost
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.7))

values = df.loc[SCENARIOS, "Power-system cost proxy [EUR bn]"]

bars = ax.bar(
    range(len(SCENARIOS)),
    values,
)

ax.set_xticks(range(len(SCENARIOS)))

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Upstream power-system cost [EUR bn/a]")

clean_axis(ax)

annotate_bars(
    ax,
    bars,
    fmt="{:.2f}",
)

save_figure(
    fig,
    "figure_01_upstream_power_system_cost",
)


# ============================================================
# FIGURE 02 — Incremental cost per flexible MWh
# ============================================================

flex_scenarios = [
    "S2",
    "S3",
    "S4",
    "S5",
    "S6",
]

fig, ax = plt.subplots(figsize=(7.5, 4.7))

values = df.loc[
    flex_scenarios, "Incremental power-system cost per flexible MWh [EUR/MWh]"
]

bars = ax.bar(
    range(len(flex_scenarios)),
    values,
)

ax.set_xticks(range(len(flex_scenarios)))

ax.set_xticklabels([LABELS[s] for s in flex_scenarios])

ax.set_ylabel("Incremental upstream cost [EUR/MWh$_{flex}$]")

clean_axis(ax)

annotate_bars(
    ax,
    bars,
    fmt="{:.3f}",
)

save_figure(
    fig,
    "figure_02_incremental_cost_per_flexible_MWh",
)


# ============================================================
# FIGURE 03 — CO2 emissions
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.7))

values = df.loc[SCENARIOS, "CO2 emissions [MtCO2/a]"]

bars = ax.bar(
    range(len(SCENARIOS)),
    values,
)

ax.set_xticks(range(len(SCENARIOS)))

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Direct fossil CO$_2$ emissions [MtCO$_2$/a]")

clean_axis(ax)

annotate_bars(
    ax,
    bars,
    fmt="{:.2f}",
)

save_figure(
    fig,
    "figure_03_CO2_emissions",
)


# ============================================================
# FIGURE 04 — Renewable curtailment
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.8))

x = np.arange(len(SCENARIOS))

solar = df.loc[SCENARIOS, "Solar curtailment [TWh]"].to_numpy()

wind = df.loc[SCENARIOS, "Wind curtailment [TWh]"].to_numpy()

ax.bar(
    x,
    solar,
    label="Solar",
)

ax.bar(
    x,
    wind,
    bottom=solar,
    label="Wind",
)

ax.set_xticks(x)

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Renewable curtailment [TWh/a]")

ax.legend(frameon=False)

clean_axis(ax)

save_figure(
    fig,
    "figure_04_solar_wind_curtailment",
)


# ============================================================
# FIGURE 05 — Solar and wind capacity
# ============================================================

fig, ax = plt.subplots(figsize=(8.2, 4.8))

x = np.arange(len(SCENARIOS))

width = 0.38

solar = df.loc[SCENARIOS, "Solar capacity [GW]"].to_numpy()

wind = df.loc[SCENARIOS, "Wind capacity [GW]"].to_numpy()

ax.bar(
    x - width / 2,
    solar,
    width=width,
    label="Solar",
)

ax.bar(
    x + width / 2,
    wind,
    width=width,
    label="Wind",
)

ax.set_xticks(x)

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Installed generation capacity [GW]")

ax.legend(frameon=False)

clean_axis(ax)

save_figure(
    fig,
    "figure_05_solar_wind_capacity",
)


# ============================================================
# FIGURE 06 — Battery energy capacity
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.7))

values = df.loc[SCENARIOS, "Battery energy capacity [GWh]"]

bars = ax.bar(
    range(len(SCENARIOS)),
    values,
)

ax.set_xticks(range(len(SCENARIOS)))

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Battery energy capacity [GWh]")

clean_axis(ax)

annotate_bars(
    ax,
    bars,
    fmt="{:.1f}",
)

save_figure(
    fig,
    "figure_06_battery_energy_capacity",
)


# ============================================================
# FIGURE 07 — Flexible-consumer utilization
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.8))

x = np.arange(len(SCENARIOS))

width = 0.38

btc_cf = df.loc[SCENARIOS, "BTC capacity factor [%]"].to_numpy()

pem_cf = df.loc[SCENARIOS, "PEM capacity factor [%]"].to_numpy()

ax.bar(
    x - width / 2,
    btc_cf,
    width=width,
    label="Bitcoin mining",
)

ax.bar(
    x + width / 2,
    pem_cf,
    width=width,
    label="PEM electrolysis",
)

ax.set_xticks(x)

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Capacity factor [%]")

ax.set_ylim(
    0,
    105,
)

ax.legend(frameon=False)

clean_axis(ax)

save_figure(
    fig,
    "figure_07_flexible_consumer_capacity_factor",
)


# ============================================================
# FIGURE 08 — Flexible electricity consumption
# ============================================================

fig, ax = plt.subplots(figsize=(8.0, 4.8))

btc_e = df.loc[SCENARIOS, "BTC consumption [TWh]"].to_numpy()

pem_e = df.loc[SCENARIOS, "PEM electricity [TWh]"].to_numpy()

ax.bar(
    x,
    btc_e,
    label="Bitcoin mining",
)

ax.bar(
    x,
    pem_e,
    bottom=btc_e,
    label="PEM electrolysis",
)

ax.set_xticks(x)

ax.set_xticklabels([LABELS[s] for s in SCENARIOS])

ax.set_ylabel("Flexible electricity consumption [TWh/a]")

ax.legend(frameon=False)

clean_axis(ax)

save_figure(
    fig,
    "figure_08_flexible_electricity_consumption",
)


# ============================================================
# Summary
# ============================================================

print()
print("=" * 100)
print("THESIS OUTPUT GENERATION COMPLETE")
print("=" * 100)

print()
print(f"Figures: {FIG_DIR}")
print(f"Tables : {TABLE_DIR}")

print()
print("Recommended core figures for Chapter 4:")
print("  figure_02_incremental_cost_per_flexible_MWh")
print("  figure_03_CO2_emissions")
print("  figure_04_solar_wind_curtailment")
print("  figure_05_solar_wind_capacity")
print("  figure_06_battery_energy_capacity")
print("  figure_07_flexible_consumer_capacity_factor")

print()
print("Recommended core tables:")
print("  table_00_scenario_matrix")
print("  table_01_headline_kpis")
print("  table_02_baseline_relative_effects")
print("  table_03_synergy_conflict_summary")

print("=" * 100)
