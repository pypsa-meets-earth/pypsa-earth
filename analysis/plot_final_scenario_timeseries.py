from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa

# ============================================================
# Paths
# ============================================================

SCENARIOS = {
    "S0_2045_Reference": Path("results/scenarios/S0_2045_1h_final/" "S0_KZ_2045_1h.nc"),
    "S1_2045_NetZero": Path(
        "results/scenarios/S1_2045_NetZero_1h_final/" "S1_KZ_2045_NetZero_1h.nc"
    ),
}

SCENARIO_LABELS = {
    "S0_2045_Reference": "S0 Reference",
    "S1_2045_NetZero": "S1 Zero direct CO$_2$",
}

OUTPUT_DIR = Path("results/scenarios/figures")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Helpers
# ============================================================


def load_network(path):
    if not path.exists():
        raise FileNotFoundError(path)

    return pypsa.Network(path)


def total_load(n):
    """
    National hourly electricity demand [MW].
    """
    return n.loads_t.p_set.sum(axis=1)


def generator_dispatch(n, carrier):
    """
    National hourly generator dispatch by carrier [MW].
    """
    idx = n.generators.index[n.generators.carrier == carrier]

    if len(idx) == 0:
        return pd.Series(0.0, index=n.snapshots)

    return n.generators_t.p[idx].sum(axis=1)


def reservoir_dispatch(n):
    """
    National hourly reservoir-hydro electricity output [MW].
    """
    idx = n.storage_units.index[n.storage_units.carrier == "hydro"]

    if len(idx) == 0:
        return pd.Series(0.0, index=n.snapshots)

    if hasattr(n.storage_units_t, "p_dispatch"):
        return n.storage_units_t.p_dispatch[idx].sum(axis=1)

    return n.storage_units_t.p[idx].clip(lower=0).sum(axis=1)


def battery_state_of_charge(n):
    """
    Aggregate battery Store energy [GWh].
    """
    idx = n.stores.index[n.stores.carrier == "battery"]

    if len(idx) == 0:
        return pd.Series(0.0, index=n.snapshots)

    return n.stores_t.e[idx].sum(axis=1) / 1000


def battery_charge_discharge(n):
    """
    Electricity-side battery charging/discharging [GW].

    charger:
        p0 > 0 means electricity consumed for charging

    discharger:
        -p1 > 0 means electricity delivered back to AC buses
    """

    charger_idx = n.links.index[n.links.carrier == "battery charger"]

    discharger_idx = n.links.index[n.links.carrier == "battery discharger"]

    if len(charger_idx):
        charge = n.links_t.p0[charger_idx].sum(axis=1).clip(lower=0) / 1000
    else:
        charge = pd.Series(0.0, index=n.snapshots)

    if len(discharger_idx):
        discharge = -n.links_t.p1[discharger_idx].sum(axis=1).clip(lower=0) / 1000
    else:
        discharge = pd.Series(0.0, index=n.snapshots)

    return charge, discharge


def hourly_load_weighted_price(n):
    """
    Hourly national load-weighted marginal electricity price [EUR/MWh].
    """
    price = n.buses_t.marginal_price

    load_by_bus = n.loads_t.p_set.T.groupby(n.loads.bus).sum().T

    common = price.columns.intersection(load_by_bus.columns)

    p = price[common]
    l = load_by_bus[common]

    denominator = l.sum(axis=1)

    numerator = (p * l).sum(axis=1)

    result = numerator.div(denominator.replace(0, np.nan))

    return result


def build_operation_dataframe(n):
    """
    Assemble national hourly operation.
    Units:
      generation/load/residual = GW
      battery state of charge = GWh
      price = EUR/MWh
    """

    df = pd.DataFrame(index=n.snapshots)

    df["Demand"] = total_load(n) / 1000

    df["Solar"] = generator_dispatch(n, "solar") / 1000

    df["Wind"] = generator_dispatch(n, "onwind") / 1000

    df["Run-of-river"] = generator_dispatch(n, "ror") / 1000

    df["Reservoir hydro"] = reservoir_dispatch(n) / 1000

    df["Coal"] = generator_dispatch(n, "coal") / 1000

    df["CCGT"] = generator_dispatch(n, "CCGT") / 1000

    df["OCGT"] = generator_dispatch(n, "OCGT") / 1000

    df["Load shedding"] = generator_dispatch(n, "load shedding") / 1000

    charge, discharge = battery_charge_discharge(n)

    df["Battery charge"] = charge
    df["Battery discharge"] = discharge

    df["Battery energy"] = battery_state_of_charge(n)

    df["Marginal price"] = hourly_load_weighted_price(n)

    renewable_direct = (
        df["Solar"] + df["Wind"] + df["Run-of-river"] + df["Reservoir hydro"]
    )

    df["Residual demand"] = df["Demand"] - renewable_direct

    return df


def save_figure(filename):
    path = OUTPUT_DIR / filename

    plt.tight_layout()

    plt.savefig(path, dpi=300, bbox_inches="tight")

    plt.close()

    print(f"Saved: {path}")


# ============================================================
# Load both final networks
# ============================================================

networks = {name: load_network(path) for name, path in SCENARIOS.items()}

operations = {name: build_operation_dataframe(n) for name, n in networks.items()}


# ============================================================
# Automatically identify the S1 stress week
# ============================================================

s1 = operations["S1_2045_NetZero"]

ROLLING_HOURS = 168

rolling_residual = (
    s1["Residual demand"].rolling(ROLLING_HOURS, min_periods=ROLLING_HOURS).mean()
)

stress_end = rolling_residual.idxmax()

stress_end_pos = s1.index.get_loc(stress_end)

stress_start_pos = stress_end_pos - ROLLING_HOURS + 1

stress_index = s1.index[stress_start_pos : stress_end_pos + 1]

stress_start = stress_index[0]
stress_end = stress_index[-1]

print()
print("=" * 90)
print("AUTOMATICALLY SELECTED S1 STRESS WEEK")
print("=" * 90)
print(f"Start : {stress_start}")
print(f"End   : {stress_end}")
print(
    "Criterion: highest 168-hour mean residual demand "
    "after solar, wind, RoR and reservoir hydro."
)
print("=" * 90)


# ============================================================
# Figure 09a / 09b
# Representative stress-week generation dispatch
# ============================================================

generation_columns = [
    "Solar",
    "Wind",
    "Run-of-river",
    "Reservoir hydro",
    "Coal",
    "CCGT",
    "OCGT",
]

for scenario, df in operations.items():

    week = df.loc[stress_index]

    plt.figure(figsize=(13, 6))

    bottom = np.zeros(len(week))

    for column in generation_columns:

        values = week[column].to_numpy()

        plt.fill_between(
            week.index,
            bottom,
            bottom + values,
            label=column,
            alpha=0.8,
        )

        bottom += values

    plt.plot(
        week.index,
        week["Demand"],
        linewidth=1.8,
        label="Demand",
    )

    plt.ylabel("Power [GW]")

    plt.xlabel("Time")

    plt.title(f"Hourly Electricity Dispatch – " f"{SCENARIO_LABELS[scenario]}")

    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.grid(axis="y", alpha=0.25)

    filename = (
        "09a_dispatch_stress_week_S0.png"
        if scenario == "S0_2045_Reference"
        else "09b_dispatch_stress_week_S1.png"
    )

    save_figure(filename)


# ============================================================
# Figure 10
# Battery operation in stress week
# ============================================================

plt.figure(figsize=(13, 6))

for scenario, df in operations.items():

    week = df.loc[stress_index]

    plt.plot(
        week.index,
        week["Battery energy"],
        linewidth=1.6,
        label=SCENARIO_LABELS[scenario],
    )

plt.ylabel("Battery state of charge [GWh]")

plt.xlabel("Time")

plt.title("Battery State of Charge During Stress Week")

plt.legend()

plt.grid(alpha=0.25)

save_figure("10_battery_state_of_charge_stress_week.png")


# ============================================================
# Figure 11a / 11b
# Battery charging and discharging
# ============================================================

for scenario, df in operations.items():

    week = df.loc[stress_index]

    plt.figure(figsize=(13, 5.5))

    plt.plot(
        week.index,
        week["Battery discharge"],
        label="Battery discharge",
        linewidth=1.4,
    )

    plt.plot(
        week.index,
        -week["Battery charge"],
        label="Battery charge",
        linewidth=1.4,
    )

    plt.axhline(0, linewidth=0.8)

    plt.ylabel("Battery power [GW]")

    plt.xlabel("Time")

    plt.title(f"Battery Charging and Discharging – " f"{SCENARIO_LABELS[scenario]}")

    plt.legend()

    plt.grid(alpha=0.25)

    filename = (
        "11a_battery_operation_S0.png"
        if scenario == "S0_2045_Reference"
        else "11b_battery_operation_S1.png"
    )

    save_figure(filename)


# ============================================================
# Figure 12
# Residual demand comparison
# ============================================================

plt.figure(figsize=(13, 5.5))

for scenario, df in operations.items():

    week = df.loc[stress_index]

    plt.plot(
        week.index,
        week["Residual demand"],
        linewidth=1.5,
        label=SCENARIO_LABELS[scenario],
    )

plt.axhline(0, linewidth=0.8)

plt.ylabel("Residual demand [GW]")

plt.xlabel("Time")

plt.title("Residual Demand During S1 Stress Week")

plt.legend()

plt.grid(alpha=0.25)

save_figure("12_residual_demand_stress_week.png")


# ============================================================
# Figure 13
# Marginal-price duration curve
# ============================================================

plt.figure(figsize=(9, 5.5))

for scenario, df in operations.items():

    values = df["Marginal price"].dropna().sort_values(ascending=False).to_numpy()

    duration = np.arange(1, len(values) + 1) / len(values) * 100

    plt.plot(
        duration,
        values,
        linewidth=1.5,
        label=SCENARIO_LABELS[scenario],
    )

plt.xlabel("Share of hours exceeded [%]")

plt.ylabel("Load-weighted marginal price [EUR/MWh]")

plt.title("Marginal Price Duration Curve – Kazakhstan 2045")

plt.legend()

plt.grid(alpha=0.25)

save_figure("13_marginal_price_duration_curve.png")


# ============================================================
# Figure 14
# Residual-demand duration curve
# ============================================================

plt.figure(figsize=(9, 5.5))

for scenario, df in operations.items():

    values = df["Residual demand"].sort_values(ascending=False).to_numpy()

    duration = np.arange(1, len(values) + 1) / len(values) * 100

    plt.plot(
        duration,
        values,
        linewidth=1.5,
        label=SCENARIO_LABELS[scenario],
    )

plt.axhline(0, linewidth=0.8)

plt.xlabel("Share of hours exceeded [%]")

plt.ylabel("Residual demand [GW]")

plt.title("Residual-Demand Duration Curve – Kazakhstan 2045")

plt.legend()

plt.grid(alpha=0.25)

save_figure("14_residual_demand_duration_curve.png")


print()
print("=" * 90)
print("Hourly scenario analysis completed successfully.")
print(f"Stress week: {stress_start} to {stress_end}")
print(f"Figures saved to: {OUTPUT_DIR}")
print("=" * 90)
