from pathlib import Path

import numpy as np
import pandas as pd
import pypsa


# ============================================================
# Paths
# ============================================================

SCENARIOS = {
    "S0_2045_Reference": Path(
        "results/scenarios/S0_2045_1h_final/"
        "S0_KZ_2045_1h.nc"
    ),
    "S1_2045_NetZero": Path(
        "results/scenarios/S1_2045_NetZero_1h_final/"
        "S1_KZ_2045_NetZero_1h.nc"
    ),
}

OUTPUT = Path(
    "results/scenarios/final_operational_metrics.csv"
)


# ============================================================
# Helpers
# ============================================================

def weights(n):
    if isinstance(n.snapshot_weightings, pd.DataFrame):
        return n.snapshot_weightings["generators"]

    return n.snapshot_weightings


def generator_dispatch(n, carrier):
    idx = n.generators.index[
        n.generators.carrier == carrier
    ]

    if len(idx) == 0:
        return pd.Series(
            0.0,
            index=n.snapshots
        )

    return (
        n.generators_t.p[idx]
        .sum(axis=1)
    )


def reservoir_dispatch(n):
    idx = n.storage_units.index[
        n.storage_units.carrier == "hydro"
    ]

    if len(idx) == 0:
        return pd.Series(
            0.0,
            index=n.snapshots
        )

    if hasattr(
        n.storage_units_t,
        "p_dispatch"
    ):
        return (
            n.storage_units_t
            .p_dispatch[idx]
            .sum(axis=1)
        )

    return (
        n.storage_units_t
        .p[idx]
        .clip(lower=0)
        .sum(axis=1)
    )


def battery_energy(n):
    idx = n.stores.index[
        n.stores.carrier == "battery"
    ]

    if len(idx) == 0:
        return pd.Series(
            0.0,
            index=n.snapshots
        )

    return (
        n.stores_t.e[idx]
        .sum(axis=1)
        / 1000
    )


def battery_capacity(n):
    x = n.stores[
        n.stores.carrier == "battery"
    ]

    if x.empty:
        return 0.0

    if "e_nom_opt" in x.columns:
        cap = x.e_nom_opt.where(
            x.e_nom_opt.notna(),
            x.e_nom
        )
    else:
        cap = x.e_nom

    return cap.sum() / 1000


def battery_charge_discharge(n):
    charger = n.links.index[
        n.links.carrier == "battery charger"
    ]

    discharger = n.links.index[
        n.links.carrier == "battery discharger"
    ]

    if len(charger):
        charge = (
            n.links_t.p0[charger]
            .sum(axis=1)
            .clip(lower=0)
            / 1000
        )
    else:
        charge = pd.Series(
            0.0,
            index=n.snapshots
        )

    if len(discharger):
        discharge = (
            -n.links_t.p1[discharger]
            .sum(axis=1)
            .clip(lower=0)
            / 1000
        )
    else:
        discharge = pd.Series(
            0.0,
            index=n.snapshots
        )

    return charge, discharge


def load_weighted_price(n):
    price = n.buses_t.marginal_price

    load_by_bus = (
        n.loads_t.p_set
        .T.groupby(n.loads.bus)
        .sum()
        .T
    )

    common = price.columns.intersection(
        load_by_bus.columns
    )

    p = price[common]
    l = load_by_bus[common]

    numerator = (
        p * l
    ).sum(axis=1)

    denominator = (
        l.sum(axis=1)
        .replace(0, np.nan)
    )

    return numerator / denominator


def operation_dataframe(n):
    df = pd.DataFrame(
        index=n.snapshots
    )

    df["Demand"] = (
        n.loads_t.p_set.sum(axis=1)
        / 1000
    )

    df["Solar"] = (
        generator_dispatch(
            n,
            "solar"
        )
        / 1000
    )

    df["Wind"] = (
        generator_dispatch(
            n,
            "onwind"
        )
        / 1000
    )

    df["RoR"] = (
        generator_dispatch(
            n,
            "ror"
        )
        / 1000
    )

    df["Reservoir"] = (
        reservoir_dispatch(n)
        / 1000
    )

    df["Coal"] = (
        generator_dispatch(
            n,
            "coal"
        )
        / 1000
    )

    df["CCGT"] = (
        generator_dispatch(
            n,
            "CCGT"
        )
        / 1000
    )

    df["OCGT"] = (
        generator_dispatch(
            n,
            "OCGT"
        )
        / 1000
    )

    df["Load shedding"] = (
        generator_dispatch(
            n,
            "load shedding"
        )
        / 1000
    )

    df["Renewable direct"] = (
        df["Solar"]
        + df["Wind"]
        + df["RoR"]
        + df["Reservoir"]
    )

    df["Residual demand"] = (
        df["Demand"]
        - df["Renewable direct"]
    )

    df["Renewable surplus"] = (
        -df["Residual demand"]
    ).clip(lower=0)

    charge, discharge = (
        battery_charge_discharge(n)
    )

    df["Battery charge"] = charge
    df["Battery discharge"] = discharge
    df["Battery energy"] = battery_energy(n)
    df["Price"] = load_weighted_price(n)

    return df


# ============================================================
# Load networks
# ============================================================

networks = {
    name: pypsa.Network(path)
    for name, path in SCENARIOS.items()
}

operations = {
    name: operation_dataframe(n)
    for name, n in networks.items()
}


# ============================================================
# Identify same S1 stress week as plotting script
# ============================================================

s1 = operations[
    "S1_2045_NetZero"
]

ROLLING_HOURS = 168

rolling = (
    s1["Residual demand"]
    .rolling(
        ROLLING_HOURS,
        min_periods=ROLLING_HOURS
    )
    .mean()
)

stress_end = rolling.idxmax()

end_pos = s1.index.get_loc(
    stress_end
)

start_pos = (
    end_pos
    - ROLLING_HOURS
    + 1
)

stress_index = s1.index[
    start_pos:
    end_pos + 1
]

stress_start = stress_index[0]
stress_end = stress_index[-1]


# ============================================================
# Extract metrics
# ============================================================

rows = []

for scenario, n in networks.items():

    df = operations[scenario]
    w = weights(n)

    # Align weights
    w = w.reindex(df.index)

    week = df.loc[
        stress_index
    ]

    w_week = w.loc[
        stress_index
    ]

    battery_cap = (
        battery_capacity(n)
    )

    annual_charge_twh = (
        df["Battery charge"]
        .mul(w)
        .sum()
        / 1000
    )

    annual_discharge_twh = (
        df["Battery discharge"]
        .mul(w)
        .sum()
        / 1000
    )

    week_charge_twh = (
        week["Battery charge"]
        .mul(w_week)
        .sum()
        / 1000
    )

    week_discharge_twh = (
        week["Battery discharge"]
        .mul(w_week)
        .sum()
        / 1000
    )

    equivalent_cycles = (
        annual_discharge_twh
        * 1000
        / battery_cap
        if battery_cap > 1e-9
        else np.nan
    )

    load_shed_twh = (
        df["Load shedding"]
        .mul(w)
        .sum()
        / 1000
    )

    row = {
        "Scenario": scenario,

        # Stress-week definition
        "Stress week start":
            str(stress_start),

        "Stress week end":
            str(stress_end),

        # Demand
        "Annual peak demand [GW]":
            df["Demand"].max(),

        "Stress-week peak demand [GW]":
            week["Demand"].max(),

        "Stress-week mean demand [GW]":
            week["Demand"].mean(),

        # Residual demand
        "Annual maximum residual demand [GW]":
            df["Residual demand"].max(),

        "Annual minimum residual demand [GW]":
            df["Residual demand"].min(),

        "Stress-week max residual demand [GW]":
            week["Residual demand"].max(),

        "Stress-week mean residual demand [GW]":
            week["Residual demand"].mean(),

        "Hours residual demand > 0 [h]":
            int(
                (
                    df["Residual demand"]
                    > 0
                ).sum()
            ),

        "Hours renewable surplus [h]":
            int(
                (
                    df["Residual demand"]
                    < 0
                ).sum()
            ),

        "Maximum renewable surplus [GW]":
            df["Renewable surplus"].max(),

        # Battery
        "Battery energy capacity [GWh]":
            battery_cap,

        "Battery maximum state of charge [GWh]":
            df["Battery energy"].max(),

        "Battery minimum state of charge [GWh]":
            df["Battery energy"].min(),

        "Annual battery charging [TWh]":
            annual_charge_twh,

        "Annual battery discharge [TWh]":
            annual_discharge_twh,

        "Stress-week battery charging [TWh]":
            week_charge_twh,

        "Stress-week battery discharge [TWh]":
            week_discharge_twh,

        "Approx. battery equivalent full cycles [1/a]":
            equivalent_cycles,

        # Load shedding
        "Annual load shedding [TWh]":
            load_shed_twh,

        "Maximum hourly load shedding [GW]":
            df["Load shedding"].max(),

        # Price
        "Minimum load-weighted price [EUR/MWh]":
            df["Price"].min(),

        "Median load-weighted price [EUR/MWh]":
            df["Price"].median(),

        "95th percentile price [EUR/MWh]":
            df["Price"].quantile(0.95),

        "99th percentile price [EUR/MWh]":
            df["Price"].quantile(0.99),

        "Maximum load-weighted price [EUR/MWh]":
            df["Price"].max(),

        "Hours with negative price [h]":
            int(
                (
                    df["Price"]
                    < 0
                ).sum()
            ),

        "Hours price > 100 EUR/MWh [h]":
            int(
                (
                    df["Price"]
                    > 100
                ).sum()
            ),
    }

    rows.append(row)


# ============================================================
# Save
# ============================================================

result = pd.DataFrame(
    rows
)

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True
)

result.to_csv(
    OUTPUT,
    index=False
)


# ============================================================
# Terminal output
# ============================================================

print()
print("=" * 110)
print("FINAL OPERATIONAL METRICS")
print("=" * 110)

print(
    result
    .set_index("Scenario")
    .T
    .to_string()
)

print()
print("=" * 110)
print(
    f"Stress week: "
    f"{stress_start} to {stress_end}"
)
print(
    f"Saved to: {OUTPUT}"
)
print("=" * 110)
