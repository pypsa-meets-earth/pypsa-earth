from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# PATHS
# ============================================================

BASE = Path("results/thesis_2045/final_analysis")

FIG = BASE / "figures_thesis"
TEX = BASE / "tables_thesis"

FIG.mkdir(parents=True, exist_ok=True)
TEX.mkdir(parents=True, exist_ok=True)


# ============================================================
# INPUTS
# ============================================================

headline = pd.read_csv(BASE / "table_01_final_headline_kpis.csv").set_index("scenario")

capacity = pd.read_csv(BASE / "table_02_capacity_by_carrier.csv")

generation = pd.read_csv(BASE / "table_03_generation_by_carrier.csv")

storage = pd.read_csv(BASE / "table_04_storage_by_carrier.csv")

battery = pd.read_csv(BASE / "table_04b_battery_results.csv").set_index("scenario")

drivers = pd.read_csv(BASE / "table_07_physical_driver_summary.csv").set_index(
    "comparison"
)


SCENARIOS = ["S0", "S1", "S2", "S3", "S4", "S5", "S6"]

SCENARIO_NAMES = {
    "S0": "Reference",
    "S1": "Zero-direct-CO$_2$",
    "S2": "Reference + BTC",
    "S3": "Reference + PEM",
    "S4": "Zero-direct-CO$_2$ + BTC",
    "S5": "Zero-direct-CO$_2$ + PEM",
    "S6": "Zero-direct-CO$_2$ + BTC + PEM",
}

COMPARISONS = [
    "S2-S0",
    "S3-S0",
    "S4-S1",
    "S5-S1",
    "S6-S1",
]


# ============================================================
# FIGURE STYLE
# ============================================================

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
    }
)

CARRIER_COLORS = {
    "coal": "#4d4d4d",
    "CCGT": "#d95f02",
    "OCGT": "#fdb863",
    "ror": "#5ab4ac",
    "hydro": "#2b8cbe",
    "onwind": "#41ab5d",
    "solar": "#fdd049",
}

CURTAILMENT_COLORS = {
    "solar": "#fdd049",
    "onwind": "#41ab5d",
}


def finish_figure(fig, name):
    fig.savefig(
        FIG / f"{name}.pdf",
        bbox_inches="tight",
    )
    fig.savefig(
        FIG / f"{name}.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(fig)


def add_values(ax, values, fmt="{:.1f}", offset=0.02):
    ymax = max(
        max(values),
        1.0,
    )

    for i, value in enumerate(values):
        ax.text(
            i,
            value + ymax * offset,
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=7,
        )


# ============================================================
# PREPARE CAPACITY AND GENERATION
# ============================================================

cap = (
    capacity.pivot(
        index="scenario",
        columns="carrier",
        values="capacity_mw",
    )
    .reindex(SCENARIOS)
    .fillna(0.0)
    / 1000.0
)

gen = (
    generation.pivot(
        index="scenario",
        columns="carrier",
        values="generation_mwh",
    )
    .reindex(SCENARIOS)
    .fillna(0.0)
    / 1e6
)

# Add reservoir hydro, represented as StorageUnit rather than Generator.
hydro = storage[storage["carrier"].astype(str).str.lower() == "hydro"].copy()

hydro_power = (
    hydro.set_index("scenario")["power_mw"].reindex(SCENARIOS).fillna(0.0) / 1000.0
)

hydro_generation = (
    hydro.set_index("scenario")["discharge_mwh"].reindex(SCENARIOS).fillna(0.0) / 1e6
)

cap["hydro"] = hydro_power
gen["hydro"] = hydro_generation

carrier_order = [
    "coal",
    "CCGT",
    "OCGT",
    "ror",
    "hydro",
    "onwind",
    "solar",
]

carrier_order = [c for c in carrier_order if c in cap.columns or c in gen.columns]


# ============================================================
# FIGURE 1 — INSTALLED POWER CAPACITY
# ============================================================

fig, ax = plt.subplots(figsize=(7.2, 4.3))

bottom = np.zeros(len(SCENARIOS))

for carrier in carrier_order:
    values = cap[carrier].values if carrier in cap.columns else np.zeros(len(SCENARIOS))

    ax.bar(
        SCENARIOS,
        values,
        bottom=bottom,
        label=carrier,
        color=CARRIER_COLORS.get(carrier),
        edgecolor="white",
        linewidth=0.4,
    )

    bottom += values

ax.set_ylabel("Installed capacity [GW]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

ax.legend(
    ncol=4,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.16),
    frameon=False,
)

finish_figure(
    fig,
    "fig_01_installed_capacity",
)


# ============================================================
# FIGURE 2 — ANNUAL ELECTRICITY GENERATION
# ============================================================

fig, ax = plt.subplots(figsize=(7.2, 4.3))

bottom = np.zeros(len(SCENARIOS))

for carrier in carrier_order:
    values = gen[carrier].values if carrier in gen.columns else np.zeros(len(SCENARIOS))

    ax.bar(
        SCENARIOS,
        values,
        bottom=bottom,
        label=carrier,
        color=CARRIER_COLORS.get(carrier),
        edgecolor="white",
        linewidth=0.4,
    )

    bottom += values

ax.set_ylabel("Annual electricity generation [TWh/a]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

ax.legend(
    ncol=4,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.16),
    frameon=False,
)

finish_figure(
    fig,
    "fig_02_generation_mix",
)


# ============================================================
# FIGURE 3 — RENEWABLE CURTAILMENT
# ============================================================

solar_curt = headline.reindex(SCENARIOS)["solar_curtailment_twh"].astype(float)

wind_curt = headline.reindex(SCENARIOS)["onwind_curtailment_twh"].astype(float)

fig, ax = plt.subplots(figsize=(7.2, 4.0))

ax.bar(
    SCENARIOS,
    solar_curt,
    label="Solar",
    color=CURTAILMENT_COLORS["solar"],
    edgecolor="white",
)

ax.bar(
    SCENARIOS,
    wind_curt,
    bottom=solar_curt,
    label="Onshore wind",
    color=CURTAILMENT_COLORS["onwind"],
    edgecolor="white",
)

ax.set_ylabel("Renewable curtailment [TWh/a]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

ax.legend(
    loc="upper left",
    frameon=False,
)

finish_figure(
    fig,
    "fig_03_renewable_curtailment",
)


# ============================================================
# FIGURE 4 — BATTERY ENERGY CAPACITY
# ============================================================

battery_energy = headline.reindex(SCENARIOS)["battery_energy_gwh"].astype(float)

fig, ax = plt.subplots(figsize=(7.2, 3.8))

ax.bar(
    SCENARIOS,
    battery_energy.values,
    edgecolor="black",
    linewidth=0.4,
)

add_values(
    ax,
    battery_energy.values,
    "{:.1f}",
)

ax.set_ylabel("Battery energy capacity [GWh]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

finish_figure(
    fig,
    "fig_04_battery_energy_capacity",
)


# ============================================================
# FIGURE 5 — TRANSMISSION EXPANSION
# ============================================================

transmission = headline.reindex(SCENARIOS)["ac_line_expansion_gwkm"].astype(float)

fig, ax = plt.subplots(figsize=(7.2, 3.8))

ax.bar(
    SCENARIOS,
    transmission.values,
    edgecolor="black",
    linewidth=0.4,
)

add_values(
    ax,
    transmission.values,
    "{:,.0f}",
)

ax.set_ylabel("AC transmission expansion [GW km]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

finish_figure(
    fig,
    "fig_05_transmission_expansion",
)


# ============================================================
# FIGURE 6 — UPSTREAM SYSTEM COST
# ============================================================

cost = headline.reindex(SCENARIOS)["upstream_system_cost_eur"].astype(float) / 1e9

fig, ax = plt.subplots(figsize=(7.2, 3.8))

ax.bar(
    SCENARIOS,
    cost.values,
    edgecolor="black",
    linewidth=0.4,
)

add_values(
    ax,
    cost.values,
    "{:.2f}",
    offset=0.015,
)

ax.set_ylabel("Upstream electricity-system cost [bn EUR/a]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

finish_figure(
    fig,
    "fig_06_upstream_system_cost",
)


# ============================================================
# FIGURE 7 — DIRECT CO2 EMISSIONS
# ============================================================

co2 = headline.reindex(SCENARIOS)["direct_generator_co2_t"].astype(float) / 1e6

fig, ax = plt.subplots(figsize=(7.2, 3.8))

ax.bar(
    SCENARIOS,
    co2.values,
    edgecolor="black",
    linewidth=0.4,
)

add_values(
    ax,
    co2.values,
    "{:.1f}",
)

ax.set_ylabel("Direct CO$_2$ emissions [MtCO$_2$/a]")
ax.set_xlabel("Scenario")
ax.grid(axis="y", alpha=0.2)

finish_figure(
    fig,
    "fig_07_direct_co2_emissions",
)


# ============================================================
# FIGURE 8 — FLEXIBLE LOAD VS CURTAILMENT REDUCTION
# ============================================================

d = drivers.reindex(COMPARISONS)

added = d["added_flexible_load_twh"].astype(float)

reduced = d["curtailment_reduction_twh"].astype(float)

x = np.arange(len(COMPARISONS))
width = 0.36

fig, ax = plt.subplots(figsize=(7.4, 4.2))

ax.bar(
    x - width / 2,
    added.values,
    width,
    label="Added flexible electricity demand",
)

ax.bar(
    x + width / 2,
    reduced.values,
    width,
    label="Reduction in wind + solar curtailment",
)

ax.axhline(
    0.0,
    linewidth=0.8,
    color="black",
)

ax.set_xticks(x)
ax.set_xticklabels(COMPARISONS)

ax.set_ylabel("Annual electricity [TWh/a]")
ax.set_xlabel("Scenario comparison")
ax.grid(axis="y", alpha=0.2)

ax.legend(
    frameon=False,
    loc="upper left",
)

finish_figure(
    fig,
    "fig_08_flexible_load_vs_curtailment",
)


# ============================================================
# LATEX TABLE 1 — SCENARIO MATRIX
# ============================================================

scenario_tex = r"""
\begin{table}[H]
    \centering
    \caption{Final national PyPSA-Earth scenario matrix.}
    \label{tab:national_scenario_matrix}
    \small
    \renewcommand{\arraystretch}{1.15}
    \begin{tabular}{llll}
        \toprule
        Scenario & Electricity-system regime & BTC mining & PEM electrolysis \\
        \midrule
        S0 & Reference & -- & -- \\
        S1 & Zero-direct-\ce{CO2} & -- & -- \\
        S2 & Reference & 1 GW & -- \\
        S3 & Reference & -- & 1 GW; 100 kt \ce{H2}/a \\
        S4 & Zero-direct-\ce{CO2} & 1 GW & -- \\
        S5 & Zero-direct-\ce{CO2} & -- & 1 GW; 100 kt \ce{H2}/a \\
        S6 & Zero-direct-\ce{CO2} & 1 GW & 1 GW; 100 kt \ce{H2}/a \\
        \bottomrule
    \end{tabular}
\end{table}
""".strip()

(TEX / "table_national_scenario_matrix.tex").write_text(scenario_tex)


# ============================================================
# LATEX TABLE 2 — HEADLINE KPIs
# ============================================================

rows = []

for s in SCENARIOS:
    total_curt = float(
        headline.loc[
            s,
            "solar_curtailment_twh",
        ]
    ) + float(
        headline.loc[
            s,
            "onwind_curtailment_twh",
        ]
    )

    rows.append(
        (
            s,
            headline.loc[
                s,
                "upstream_system_cost_eur",
            ]
            / 1e9,
            headline.loc[
                s,
                "direct_generator_co2_t",
            ]
            / 1e6,
            headline.loc[
                s,
                "btc_twh",
            ],
            headline.loc[
                s,
                "pem_el_twh",
            ],
            headline.loc[
                s,
                "h2_kt",
            ],
            total_curt,
            headline.loc[
                s,
                "battery_energy_gwh",
            ],
            headline.loc[
                s,
                "ac_line_expansion_gwkm",
            ],
        )
    )

lines = []

for row in rows:
    lines.append(
        f"{row[0]} & "
        f"{row[1]:.3f} & "
        f"{row[2]:.3f} & "
        f"{row[3]:.3f} & "
        f"{row[4]:.3f} & "
        f"{row[5]:.1f} & "
        f"{row[6]:.2f} & "
        f"{row[7]:.1f} & "
        f"{row[8]:,.0f} \\\\"
    )

headline_tex = (
    r"""
\begin{table}[H]
    \centering
    \caption{Headline results of the national PyPSA-Earth scenarios.}
    \label{tab:national_headline_results}
    \scriptsize
    \renewcommand{\arraystretch}{1.15}
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lrrrrrrrr}
        \toprule
        Scenario
        & System cost
        & Direct \ce{CO2}
        & BTC electricity
        & PEM electricity
        & \ce{H2}
        & RE curtailment
        & Battery energy
        & AC expansion \\
        & [bn EUR/a]
        & [Mt/a]
        & [TWh/a]
        & [TWh/a]
        & [kt/a]
        & [TWh/a]
        & [GWh]
        & [GW km] \\
        \midrule
"""
    + "\n".join(lines)
    + r"""
        \bottomrule
    \end{tabular}
    }
\end{table}
"""
)

(TEX / "table_national_headline_results.tex").write_text(headline_tex.strip())


# ============================================================
# LATEX TABLE 3 — PHYSICAL DRIVER COMPARISONS
# ============================================================

driver_lines = []

for comparison in COMPARISONS:
    r = drivers.loc[comparison]

    driver_lines.append(
        f"{comparison} & "
        f"{r['added_flexible_load_twh']:.3f} & "
        f"{r['delta_system_cost_m_eur']:.3f} & "
        f"{r['delta_direct_co2_mt']:.3f} & "
        f"{r['delta_solar_capacity_gw']:.3f} & "
        f"{r['delta_onwind_capacity_gw']:.3f} & "
        f"{r['curtailment_reduction_twh']:.3f} & "
        f"{r['delta_battery_energy_gwh']:.2f} & "
        f"{r['delta_transmission_expansion_gwkm']:.1f} \\\\"
    )

drivers_tex = (
    r"""
\begin{table}[H]
    \centering
    \caption{System-level effects of the flexible-consumer scenarios relative to their corresponding reference cases.}
    \label{tab:national_physical_drivers}
    \scriptsize
    \renewcommand{\arraystretch}{1.15}
    \resizebox{\textwidth}{!}{
    \begin{tabular}{lrrrrrrrr}
        \toprule
        Comparison
        & Added load
        & $\Delta C_{\mathrm{sys}}$
        & $\Delta$\ce{CO2}
        & $\Delta$ solar
        & $\Delta$ wind
        & Curtailment reduction
        & $\Delta E_{\mathrm{bat}}$
        & $\Delta$ AC expansion \\
        & [TWh/a]
        & [M EUR/a]
        & [Mt/a]
        & [GW]
        & [GW]
        & [TWh/a]
        & [GWh]
        & [GW km] \\
        \midrule
"""
    + "\n".join(driver_lines)
    + r"""
        \bottomrule
    \end{tabular}
    }
\end{table}
"""
)

(TEX / "table_national_physical_drivers.tex").write_text(drivers_tex.strip())


# ============================================================
# MANIFEST
# ============================================================

manifest = """
FINAL NATIONAL PYPSA-EARTH THESIS RESULTS PACKAGE

FIGURES
-------
fig_01_installed_capacity
    Installed electricity-supply capacity by carrier.

fig_02_generation_mix
    Annual electricity generation by carrier.

fig_03_renewable_curtailment
    Solar and onshore-wind curtailment.

fig_04_battery_energy_capacity
    Optimized battery energy capacity.

fig_05_transmission_expansion
    AC grid expansion.

fig_06_upstream_system_cost
    Corrected upstream electricity-system cost.

fig_07_direct_co2_emissions
    Direct generator CO2 emissions.

fig_08_flexible_load_vs_curtailment
    Added flexible electricity demand versus reduction in renewable
    curtailment for the principal scenario comparisons.

TABLES
------
table_national_scenario_matrix.tex
table_national_headline_results.tex
table_national_physical_drivers.tex
""".strip()

(BASE / "THESIS_RESULTS_MANIFEST.txt").write_text(manifest)


print("=" * 78)
print("THESIS RESULTS PACKAGE CREATED")
print("=" * 78)

print("\nFigures:")
for f in sorted(FIG.glob("*.pdf")):
    print(f)

print("\nLaTeX tables:")
for f in sorted(TEX.glob("*.tex")):
    print(f)

print("\nManifest:")
print(BASE / "THESIS_RESULTS_MANIFEST.txt")
