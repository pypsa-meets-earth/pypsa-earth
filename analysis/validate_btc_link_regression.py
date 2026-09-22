from pathlib import Path

import numpy as np
import pandas as pd
import pypsa


OLD_SOLVED = Path(
    "results/networks/"
    "elec_s_10_ec_lcopt_1h-S2BTC.nc"
)

NEW_SOLVED = Path(
    "results/networks/"
    "elec_s_10_ec_lcopt_1h-S2BTCLINKV2.nc"
)

OLD_INPUT = Path(
    "networks/"
    "elec_s_10_ec_lcopt_1h-S2BTC.nc"
)

NEW_INPUT = Path(
    "networks/"
    "elec_s_10_ec_lcopt_1h-S2BTCLINKV2.nc"
)

OUTPUT = Path(
    "results/scenarios/"
    "btc_link_regression_S2.csv"
)

OLD_BTC = "THESIS BTC flexible mining"
NEW_BTC = "THESIS BTC mining link"
BTC_STORE = "THESIS BTC service accumulator"


def weights(n, column):
    return (
        n.snapshot_weightings[column]
        .reindex(n.snapshots)
        .astype(float)
    )


def annual_sum(series, n, column="generators"):
    w = weights(n, column)

    return float(
        series.reindex(n.snapshots)
        .mul(w)
        .sum()
    )


def effective_generator_capacity(n):
    capacity = n.generators["p_nom"].astype(float).copy()

    if "p_nom_opt" in n.generators.columns:
        extendable = (
            n.generators["p_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        capacity.loc[extendable] = (
            n.generators.loc[
                extendable,
                "p_nom_opt",
            ].astype(float)
        )

    return capacity


def effective_store_capacity(n):
    capacity = n.stores["e_nom"].astype(float).copy()

    if "e_nom_opt" in n.stores.columns:
        extendable = (
            n.stores["e_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        capacity.loc[extendable] = (
            n.stores.loc[
                extendable,
                "e_nom_opt",
            ].astype(float)
        )

    return capacity


def effective_line_capacity(n):
    capacity = n.lines["s_nom"].astype(float).copy()

    if "s_nom_opt" in n.lines.columns:
        extendable = (
            n.lines["s_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        capacity.loc[extendable] = (
            n.lines.loc[
                extendable,
                "s_nom_opt",
            ].astype(float)
        )

    return capacity


def dense_generator_attribute(n, attribute):
    static = pd.DataFrame(
        np.tile(
            n.generators[attribute]
            .astype(float)
            .to_numpy(),
            (len(n.snapshots), 1),
        ),
        index=n.snapshots,
        columns=n.generators.index,
    )

    dynamic = getattr(
        n.generators_t,
        attribute,
    )

    if not dynamic.empty:
        cols = dynamic.columns.intersection(
            n.generators.index
        )

        static.loc[:, cols] = (
            dynamic.loc[
                n.snapshots,
                cols,
            ]
        )

    return static


def dense_load_p_set(n):
    static = pd.DataFrame(
        np.tile(
            n.loads["p_set"]
            .astype(float)
            .to_numpy(),
            (len(n.snapshots), 1),
        ),
        index=n.snapshots,
        columns=n.loads.index,
    )

    dynamic = n.loads_t.p_set

    if not dynamic.empty:
        cols = dynamic.columns.intersection(
            n.loads.index
        )

        static.loc[:, cols] = (
            dynamic.loc[
                n.snapshots,
                cols,
            ]
        )

    return static


def generator_capacity_gw(n, carrier):
    cap = effective_generator_capacity(n)

    mask = n.generators.carrier.eq(carrier)

    return float(
        cap.loc[mask].sum()
        / 1000.0
    )


def store_capacity_gwh(n, carrier):
    cap = effective_store_capacity(n)

    mask = n.stores.carrier.eq(carrier)

    return float(
        cap.loc[mask].sum()
        / 1000.0
    )


def generation_twh(n, carrier):
    names = n.generators.index[
        n.generators.carrier.eq(carrier)
    ]

    if len(names) == 0:
        return 0.0

    dispatch = (
        n.generators_t.p[names]
        .sum(axis=1)
    )

    return (
        annual_sum(dispatch, n)
        / 1e6
    )


def curtailment_twh(n, carrier):
    names = n.generators.index[
        n.generators.carrier.eq(carrier)
    ]

    if len(names) == 0:
        return 0.0

    cap = effective_generator_capacity(n)

    p_max_pu = dense_generator_attribute(
        n,
        "p_max_pu",
    )

    potential = (
        p_max_pu[names]
        .mul(
            cap.loc[names],
            axis=1,
        )
    )

    actual = n.generators_t.p[names]

    curtailed = (
        potential
        - actual
    ).clip(lower=0.0)

    hourly_total = curtailed.sum(axis=1)

    return (
        annual_sum(hourly_total, n)
        / 1e6
    )


def direct_co2_mt(n):
    if "co2_emissions" not in n.carriers.columns:
        return float("nan")

    carrier_co2 = (
        n.carriers["co2_emissions"]
        .astype(float)
    )

    co2_intensity = (
        n.generators["carrier"]
        .map(carrier_co2)
        .fillna(0.0)
    )

    efficiency = (
        n.generators["efficiency"]
        .astype(float)
        .replace(0.0, np.nan)
    )

    co2_per_mwh_el = (
        co2_intensity
        / efficiency
    ).fillna(0.0)

    dispatch = (
        n.generators_t.p[
            n.generators.index
        ]
        .clip(lower=0.0)
    )

    hourly = (
        dispatch
        .mul(
            co2_per_mwh_el,
            axis=1,
        )
        .sum(axis=1)
    )

    return (
        annual_sum(hourly, n)
        / 1e6
    )


def load_weighted_price(n):
    load = dense_load_p_set(n)

    bus_load = pd.DataFrame(
        0.0,
        index=n.snapshots,
        columns=n.buses.index,
    )

    groups = n.loads.groupby("bus").groups

    for bus, load_names in groups.items():
        if bus not in bus_load.columns:
            continue

        bus_load[bus] = (
            load[list(load_names)]
            .sum(axis=1)
        )

    prices = (
        n.buses_t.marginal_price
        .reindex(
            index=n.snapshots,
            columns=n.buses.index,
        )
    )

    w = weights(n, "generators")

    numerator = (
        (
            prices
            * bus_load
        )
        .sum(axis=1)
        .mul(w)
        .sum()
    )

    denominator = (
        bus_load
        .sum(axis=1)
        .mul(w)
        .sum()
    )

    return float(
        numerator
        / denominator
    )


def line_volume_million_mwkm(n):
    cap = effective_line_capacity(n)

    volume = (
        cap
        * n.lines["length"].astype(float)
    ).sum()

    return float(
        volume
        / 1e6
    )


def btc_metrics(old, new):
    old_dispatch_signed = (
        old.generators_t.p[OLD_BTC]
    )

    old_consumption = (
        -old_dispatch_signed
    )

    new_consumption = (
        new.links_t.p0[NEW_BTC]
        .clip(lower=0.0)
    )

    old_energy = annual_sum(
        old_consumption,
        old,
    )

    new_energy = annual_sum(
        new_consumption,
        new,
    )

    old_hours = float(
        weights(
            old,
            "generators",
        ).sum()
    )

    new_hours = float(
        weights(
            new,
            "generators",
        ).sum()
    )

    old_cf = (
        100.0
        * old_energy
        / (
            old.generators.at[
                OLD_BTC,
                "p_nom",
            ]
            * old_hours
        )
    )

    new_cf = (
        100.0
        * new_energy
        / (
            new.links.at[
                NEW_BTC,
                "p_nom",
            ]
            * new_hours
        )
    )

    return (
        old_energy,
        new_energy,
        old_cf,
        new_cf,
    )


def btc_objective_contribution_old(n):
    p = n.generators_t.p[OLD_BTC]

    mc = float(
        n.generators.at[
            OLD_BTC,
            "marginal_cost",
        ]
    )

    w = weights(n, "objective")

    return float(
        p.mul(w).sum()
        * mc
    )


def btc_objective_contribution_new(n):
    p0 = n.links_t.p0[NEW_BTC]

    mc = float(
        n.links.at[
            NEW_BTC,
            "marginal_cost",
        ]
    )

    w = weights(n, "objective")

    return float(
        p0.mul(w).sum()
        * mc
    )


for path in [
    OLD_SOLVED,
    NEW_SOLVED,
    OLD_INPUT,
    NEW_INPUT,
]:
    if not path.exists():
        raise FileNotFoundError(path)


print("Loading networks...")

old = pypsa.Network(OLD_SOLVED)
new = pypsa.Network(NEW_SOLVED)

old_input = pypsa.Network(OLD_INPUT)
new_input = pypsa.Network(NEW_INPUT)


(
    old_btc_energy,
    new_btc_energy,
    old_btc_cf,
    new_btc_cf,
) = btc_metrics(
    old,
    new,
)


old_btc_obj = (
    btc_objective_contribution_old(old)
)

new_btc_obj = (
    btc_objective_contribution_new(new)
)


old_system_cost = (
    float(old.objective)
    - old_btc_obj
)

new_system_cost = (
    float(new.objective)
    - new_btc_obj
)


metrics = {
    "BTC electricity [TWh/a]": (
        old_btc_energy / 1e6,
        new_btc_energy / 1e6,
    ),

    "BTC capacity factor [%]": (
        old_btc_cf,
        new_btc_cf,
    ),

    "Solar capacity [GW]": (
        generator_capacity_gw(
            old,
            "solar",
        ),
        generator_capacity_gw(
            new,
            "solar",
        ),
    ),

    "Wind capacity [GW]": (
        generator_capacity_gw(
            old,
            "onwind",
        ),
        generator_capacity_gw(
            new,
            "onwind",
        ),
    ),

    "Battery energy [GWh]": (
        store_capacity_gwh(
            old,
            "battery",
        ),
        store_capacity_gwh(
            new,
            "battery",
        ),
    ),

    "Solar generation [TWh/a]": (
        generation_twh(
            old,
            "solar",
        ),
        generation_twh(
            new,
            "solar",
        ),
    ),

    "Wind generation [TWh/a]": (
        generation_twh(
            old,
            "onwind",
        ),
        generation_twh(
            new,
            "onwind",
        ),
    ),

    "Coal generation [TWh/a]": (
        generation_twh(
            old,
            "coal",
        ),
        generation_twh(
            new,
            "coal",
        ),
    ),

    "CCGT generation [TWh/a]": (
        generation_twh(
            old,
            "CCGT",
        ),
        generation_twh(
            new,
            "CCGT",
        ),
    ),

    "OCGT generation [TWh/a]": (
        generation_twh(
            old,
            "OCGT",
        ),
        generation_twh(
            new,
            "OCGT",
        ),
    ),

    "Solar curtailment [TWh/a]": (
        curtailment_twh(
            old,
            "solar",
        ),
        curtailment_twh(
            new,
            "solar",
        ),
    ),

    "Wind curtailment [TWh/a]": (
        curtailment_twh(
            old,
            "onwind",
        ),
        curtailment_twh(
            new,
            "onwind",
        ),
    ),

    "Direct CO2 [Mt/a]": (
        direct_co2_mt(old),
        direct_co2_mt(new),
    ),

    "Load-weighted price [EUR/MWh]": (
        load_weighted_price(old),
        load_weighted_price(new),
    ),

    "AC line volume [million MWkm]": (
        line_volume_million_mwkm(old),
        line_volume_million_mwkm(new),
    ),

    "Raw objective [EUR bn/a]": (
        float(old.objective) / 1e9,
        float(new.objective) / 1e9,
    ),

    "Upstream system-cost proxy [EUR bn/a]": (
        old_system_cost / 1e9,
        new_system_cost / 1e9,
    ),
}


rows = []

for metric, (old_value, new_value) in metrics.items():
    delta = new_value - old_value

    if abs(old_value) > 1e-12:
        relative = (
            100.0
            * delta
            / abs(old_value)
        )
    else:
        relative = np.nan

    rows.append(
        {
            "Metric": metric,
            "Old Generator S2": old_value,
            "New Link S2": new_value,
            "Delta": delta,
            "Delta [%]": relative,
        }
    )


df = pd.DataFrame(rows)

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

df.to_csv(
    OUTPUT,
    index=False,
)


print(
    "\n"
    "============================================================"
)

print(
    "S2 BTC GENERATOR -> LINK FULL REGRESSION"
)

print(
    "============================================================"
)

print(
    df.to_string(
        index=False,
        float_format=lambda x: f"{x:.9f}",
    )
)


print(
    "\n"
    "============================================================"
)

print(
    "ECONOMIC INPUT CHECK"
)

print(
    "============================================================"
)

print(
    "Old unsolved Generator marginal cost:",
    old_input.generators.at[
        OLD_BTC,
        "marginal_cost",
    ],
)

print(
    "New unsolved Link marginal cost:",
    new_input.links.at[
        NEW_BTC,
        "marginal_cost",
    ],
)

print(
    "Old solved Generator marginal cost:",
    old.generators.at[
        OLD_BTC,
        "marginal_cost",
    ],
)

print(
    "New solved Link marginal cost:",
    new.links.at[
        NEW_BTC,
        "marginal_cost",
    ],
)


print(
    "\n"
    "============================================================"
)

print(
    "BTC BOOKKEEPING STORE CHECK"
)

print(
    "============================================================"
)

store_p = (
    new.stores_t.p[BTC_STORE]
)

link_p0 = (
    new.links_t.p0[NEW_BTC]
)

eff = float(
    new.links.at[
        NEW_BTC,
        "efficiency",
    ]
)

service_balance_error = (
    (-store_p)
    - link_p0 * eff
).abs().max()

store_e = (
    new.stores_t.e[BTC_STORE]
)

print(
    "Maximum Store discharge [MW]:",
    float(store_p.max()),
)

print(
    "Maximum Store charging [MW]:",
    float(-store_p.min()),
)

print(
    "Maximum service-balance error [MW]:",
    float(service_balance_error),
)

print(
    "Final accumulated BTC service [MWh]:",
    float(store_e.iloc[-1]),
)

print(
    "Store energy capacity [MWh]:",
    float(
        new.stores.at[
            BTC_STORE,
            "e_nom",
        ]
    ),
)

print(
    "\nSaved regression table:"
)

print(OUTPUT)
