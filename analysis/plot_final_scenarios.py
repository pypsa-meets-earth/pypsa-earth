from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ============================================================
# Paths
# ============================================================

INPUT = Path("results/scenarios/final_scenario_kpis.csv")

OUTPUT_DIR = Path("results/scenarios/figures")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(INPUT)

df = df.set_index("Scenario")

scenario_names = {
    "S0_2045_Reference": "S0 Reference",
    "S1_2045_NetZero": "S1 Zero direct CO$_2$",
}

scenarios = list(scenario_names.keys())
labels = [scenario_names[s] for s in scenarios]


# ============================================================
# Helper
# ============================================================


def save_figure(filename):
    path = OUTPUT_DIR / filename

    plt.tight_layout()

    plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.close()

    print(f"Saved: {path}")


# ============================================================
# Figure 1
# Installed generation capacity
# ============================================================

capacity_kpis = {
    "Solar": "Solar capacity [GW]",
    "Wind": "Wind capacity [GW]",
    "Coal": "Coal capacity [GW]",
    "CCGT": "CCGT capacity [GW]",
    "OCGT": "OCGT capacity [GW]",
    "Reservoir hydro": "Reservoir hydro capacity [GW]",
    "Run-of-river": "RoR capacity [GW]",
}

x = np.arange(len(scenarios))
width = 0.11

plt.figure(figsize=(11, 6))

offsets = (np.arange(len(capacity_kpis)) - (len(capacity_kpis) - 1) / 2) * width

for offset, (technology, column) in zip(offsets, capacity_kpis.items()):
    values = [df.loc[s, column] for s in scenarios]

    plt.bar(
        x + offset,
        values,
        width,
        label=technology,
    )

plt.xticks(x, labels)

plt.ylabel("Installed capacity [GW]")

plt.title("Installed Generation Capacity – Kazakhstan 2045")

plt.legend(ncol=2)

plt.grid(axis="y", alpha=0.3)

save_figure("01_installed_generation_capacity.png")


# ============================================================
# Figure 2
# Electricity generation
# ============================================================

generation_kpis = {
    "Solar": "Solar generation [TWh]",
    "Wind": "Wind generation [TWh]",
    "Reservoir hydro": "Reservoir hydro generation [TWh]",
    "Run-of-river": "RoR generation [TWh]",
    "Coal": "Coal generation [TWh]",
    "CCGT": "CCGT generation [TWh]",
    "OCGT": "OCGT generation [TWh]",
}

plt.figure(figsize=(9, 6))

bottom = np.zeros(len(scenarios))

for technology, column in generation_kpis.items():

    values = np.array([df.loc[s, column] for s in scenarios])

    plt.bar(
        labels,
        values,
        bottom=bottom,
        label=technology,
    )

    bottom += values

plt.axhline(
    df.loc["S0_2045_Reference", "Annual demand [TWh]"],
    linestyle="--",
    linewidth=1.3,
    label="Annual demand",
)

plt.ylabel("Electricity generation [TWh]")

plt.title("Electricity Generation Mix – Kazakhstan 2045")

plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

plt.grid(axis="y", alpha=0.3)

save_figure("02_generation_mix.png")


# ============================================================
# Figure 3
# Battery deployment
# ============================================================

battery_kpis = {
    "Battery energy [GWh]": "Battery energy [GWh]",
    "Battery charger [GW]": "Battery charger [GW]",
    "Battery discharger [GW]": "Battery discharger [GW]",
}

for title, column in battery_kpis.items():

    plt.figure(figsize=(7, 5))

    values = [df.loc[s, column] for s in scenarios]

    plt.bar(labels, values)

    plt.ylabel(title)

    plt.title(f"{title} – Kazakhstan 2045")

    plt.grid(axis="y", alpha=0.3)

    filename = (
        title.lower()
        .replace(" ", "_")
        .replace("[", "")
        .replace("]", "")
        .replace("/", "_")
        + ".png"
    )

    save_figure(filename)


# ============================================================
# Figure 4
# Renewable curtailment
# ============================================================

curtailment = {
    "Solar": "Solar curtailment [%]",
    "Wind": "Wind curtailment [%]",
    "Run-of-river": "RoR curtailment [%]",
}

x = np.arange(len(scenarios))

width = 0.24

plt.figure(figsize=(9, 5.5))

for i, (technology, column) in enumerate(curtailment.items()):

    values = [df.loc[s, column] for s in scenarios]

    plt.bar(
        x + (i - 1) * width,
        values,
        width,
        label=technology,
    )

plt.xticks(x, labels)

plt.ylabel("Curtailment [%]")

plt.title("Renewable Curtailment – Kazakhstan 2045")

plt.legend()

plt.grid(axis="y", alpha=0.3)

save_figure("04_renewable_curtailment.png")


# ============================================================
# Figure 5
# CO2 emissions
# ============================================================

plt.figure(figsize=(7, 5))

values = [df.loc[s, "CO2 emissions [MtCO2/a]"] for s in scenarios]

plt.bar(labels, values)

plt.ylabel("CO$_2$ emissions [MtCO$_2$/a]")

plt.title("Electricity-Sector CO$_2$ Emissions – Kazakhstan 2045")

plt.grid(axis="y", alpha=0.3)

save_figure("05_co2_emissions.png")


# ============================================================
# Figure 6
# Marginal electricity prices
# ============================================================

price_kpis = {
    "Load-weighted": "Load-weighted marginal price [EUR/MWh]",
    "Unweighted": "Unweighted mean marginal price [EUR/MWh]",
}

x = np.arange(len(scenarios))

width = 0.32

plt.figure(figsize=(8, 5))

for i, (name, column) in enumerate(price_kpis.items()):

    values = [df.loc[s, column] for s in scenarios]

    offset = (i - 0.5) * width

    plt.bar(
        x + offset,
        values,
        width,
        label=name,
    )

plt.xticks(x, labels)

plt.ylabel("Marginal electricity price [EUR/MWh]")

plt.title("Mean Marginal Electricity Prices – Kazakhstan 2045")

plt.legend()

plt.grid(axis="y", alpha=0.3)

save_figure("06_marginal_prices.png")


# ============================================================
# Figure 7
# Objective
# ============================================================

plt.figure(figsize=(7, 5))

values = [df.loc[s, "Objective [EUR bn]"] for s in scenarios]

plt.bar(labels, values)

plt.ylabel("Objective [EUR billion]")

plt.title("Optimized System Objective – Kazakhstan 2045")

plt.grid(axis="y", alpha=0.3)

save_figure("07_system_objective.png")


# ============================================================
# Figure 8
# AC transmission line volume
# ============================================================

plt.figure(figsize=(7, 5))

values = [df.loc[s, "AC line volume [million MWkm]"] for s in scenarios]

plt.bar(labels, values)

plt.ylabel("AC line volume [million MWkm]")

plt.title("AC Transmission Line Volume – Kazakhstan 2045")

plt.grid(axis="y", alpha=0.3)

save_figure("08_ac_line_volume.png")


print()
print("=" * 80)
print("All final scenario figures created successfully.")
print(f"Output directory: {OUTPUT_DIR}")
print("=" * 80)
