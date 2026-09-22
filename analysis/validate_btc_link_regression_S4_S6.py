from pathlib import Path

import numpy as np
import pandas as pd
import pypsa


# =============================================================================
# Files
# =============================================================================

S1_PATH = Path(
    "results/networks/"
    "elec_s_10_ec_lcopt_Co2L-1h.nc"
)

S5_PATH = Path(
    "results/networks/"
    "elec_s_10_ec_lcopt_Co2L-1h-S5PEM.nc"
)

SCENARIOS = {
    "S4": {
        "old": Path(
            "results/networks/"
            "elec_s_10_ec_lcopt_Co2L-1h-S4BTC.nc"
        ),
        "new": Path(
            "results/networks/"
            "elec_s_10_ec_lcopt_Co2L-1h-S4BTCLINKV2.nc"
        ),
    },
    "S6": {
        "old": Path(
            "results/networks/"
            "elec_s_10_ec_lcopt_Co2L-1h-S6BOTH.nc"
        ),
        "new": Path(
            "results/networks/"
            "elec_s_10_ec_lcopt_Co2L-1h-S6BTCLINKV2.nc"
        ),
    },
}

OLD_BTC = "THESIS BTC flexible mining"
NEW_BTC = "THESIS BTC mining link"
BTC_STORE = "THESIS BTC service accumulator"

OUTPUT = Path(
    "results/scenarios/"
    "btc_link_regression_S4_S6.csv"
)


# =============================================================================
# Generic helpers
# =============================================================================

def weights(n, column="generators"):
    return (
        n.snapshot_weightings[column]
        .reindex(n.snapshots)
        .astype(float)
    )


def annual_sum(series, n, column="generators"):
    return float(
        series.reindex(n.snapshots)
        .mul(weights(n, column))
        .sum()
    )


def effective_generator_capacity(n):
    cap = n.generators["p_nom"].astype(float).copy()

    if "p_nom_opt" in n.generators.columns:
        ext = (
            n.generators["p_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        cap.loc[ext] = (
            n.generators.loc[
                ext,
                "p_nom_opt",
            ].astype(float)
        )

    return cap


def effective_store_capacity(n):
    cap = n.stores["e_nom"].astype(float).copy()

    if "e_nom_opt" in n.stores.columns:
        ext = (
            n.stores["e_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        cap.loc[ext] = (
            n.stores.loc[
                ext,
                "e_nom_opt",
            ].astype(float)
        )

    return cap


def effective_line_capacity(n):
    cap = n.lines["s_nom"].astype(float).copy()

    if "s_nom_opt" in n.lines.columns:
        ext = (
            n.lines["s_nom_extendable"]
            .fillna(False)
            .astype(bool)
        )

        cap.loc[ext] = (
            n.lines.loc[
                ext,
                "s_nom_opt",
            ].astype(float)
        )

    return cap


def dense_generator_attribute(n, attr):
    df = pd.DataFrame(
        np.tile(
            n.generators[attr]
            .astype(float)
            .to_numpy(),
            (len(n.snapshots), 1),
        ),
        index=n.snapshots,
        columns=n.generators.index,
    )

    dynamic = getattr(
        n.generators_t,
        attr,
    )

    if not dynamic.empty:
        cols = dynamic.columns.intersection(
            n.generators.index
        )

        df.loc[:, cols] = (
            dynamic.loc[
                n.snapshots,
                cols,
            ]
        )

    return df


def dense_load_p_set(n):
    df = pd.DataFrame(
        np.tile(
            n.loads["p_set"]
            .astype(float)
            .to_numpy(),
            (len(n.snapshots), 1),
        ),
        index=n.snapshots,
        columns=n.loads.index,
    )

    if not n.loads_t.p_set.empty:
        cols = (
            n.loads_t.p_set.columns
            .intersection(n.loads.index)
        )

        df.loc[:, cols] = (
            n.loads_t.p_set.loc[
                n.snapshots,
                cols,
            ]
        )

    return df


# =============================================================================
# System KPIs
# =============================================================================

def generator_capacity_gw(n, carrier):
    mask = n.generators.carrier.eq(carrier)

    return float(
        effective_generator_capacity(n)
        .loc[mask]
        .sum()
        / 1000.0
    )


def store_capacity_gwh(n, carrier):
    mask = n.stores.carrier.eq(carrier)

    return float(
        effective_store_capacity(n)
        .loc[mask]
        .sum()
        / 1000.0
    )


def generation_twh(n, carrier):
    names = n.generators.index[
        n.generators.carrier.eq(carrier)
    ]

    if len(names) == 0:
        return 0.0

    p = (
        n.generators_t.p[names]
        .sum(axis=1)
    )

    return annual_sum(p, n) / 1e6


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

    return (
        annual_sum(
            curtailed.sum(axis=1),
            n,
        )
        / 1e6
    )


def direct_co2_mt(n):
    carrier_co2 = (
        n.carriers
        .get(
            "co2_emissions",
            pd.Series(
                0.0,
                index=n.carriers.index,
            ),
        )
        .astype(float)
    )

    intensity = (
        n.generators["carrier"]
        .map(carrier_co2)
        .fillna(0.0)
    )

    efficiency = (
        n.generators["efficiency"]
        .astype(float)
        .replace(0.0, np.nan)
    )

    specific = (
        intensity
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
            specific,
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

    for bus, names in (
        n.loads.groupby("bus").groups.items()
    ):
        if bus in bus_load.columns:
            bus_load[bus] = (
                load[list(names)]
                .sum(axis=1)
            )

    prices = (
        n.buses_t.marginal_price
        .reindex(
            index=n.snapshots,
            columns=n.buses.index,
        )
    )

    w = weights(n)

    numerator = (
        (prices * bus_load)
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
    volume = (
        effective_line_capacity(n)
        * n.lines["length"].astype(float)
    ).sum()

    return float(volume / 1e6)


# =============================================================================
# BTC
# =============================================================================

def btc_energy_mwh(n, representation):
    if representation == "old":
        p = (
            -n.generators_t.p[OLD_BTC]
        ).clip(lower=0.0)

    else:
        p = (
            n.links_t.p0[NEW_BTC]
        ).clip(lower=0.0)

    return annual_sum(p, n)


def btc_capacity_mw(n, representation):
    if representation == "old":
        return float(
            n.generators.at[
                OLD_BTC,
                "p_nom",
            ]
        )

    return float(
        n.links.at[
            NEW_BTC,
            "p_nom",
        ]
    )


def btc_capacity_factor(n, representation):
    energy = btc_energy_mwh(
        n,
        representation,
    )

    hours = float(weights(n).sum())

    return (
        100.0
        * energy
        / (
            btc_capacity_mw(
                n,
                representation,
            )
            * hours
        )
    )


def btc_objective_term(n, representation):
    w = weights(n, "objective")

    if representation == "old":
        p = n.generators_t.p[OLD_BTC]

        mc = float(
            n.generators.at[
                OLD_BTC,
                "marginal_cost",
            ]
        )

    else:
        p = n.links_t.p0[NEW_BTC]

        mc = float(
            n.links.at[
                NEW_BTC,
                "marginal_cost",
            ]
        )

    return float(
        p.mul(w).sum()
        * mc
    )


def upstream_system_cost(n, representation):
    return (
        float(n.objective)
        - btc_objective_term(
            n,
            representation,
        )
    )


# =============================================================================
# H2 identification from S5 relative to S1
# =============================================================================

def identify_pem_components(s1, s5):
    new_links = (
        s5.links.index
        .difference(s1.links.index)
    )

    new_stores = (
        s5.stores.index
        .difference(s1.stores.index)
    )

    if len(new_links) == 0:
        raise RuntimeError(
            "No S5-specific Link found."
        )

    # Prefer explicit PEM/electrolyser naming.
    candidates = [
        x for x in new_links
        if (
            "pem" in str(x).lower()
            or "electro" in str(x).lower()
            or "h2" in str(x).lower()
        )
    ]

    if len(candidates) == 1:
        pem_link = candidates[0]

    elif len(new_links) == 1:
        pem_link = new_links[0]

    else:
        raise RuntimeError(
            "Could not uniquely identify PEM Link. "
            f"S5-specific Links: {list(new_links)}"
        )

    # Same logic for bookkeeping product Store.
    store_candidates = [
        x for x in new_stores
        if (
            "pem" in str(x).lower()
            or "hydrogen" in str(x).lower()
            or "h2" in str(x).lower()
            or "product" in str(x).lower()
        )
    ]

    if len(store_candidates) == 1:
        pem_store = store_candidates[0]

    elif len(new_stores) == 1:
        pem_store = new_stores[0]

    else:
        pem_store = None

    return (
        pem_link,
        pem_store,
        list(new_links),
        list(new_stores),
    )


def pem_electricity_twh(n, pem_link):
    return (
        annual_sum(
            n.links_t.p0[pem_link]
            .clip(lower=0.0),
            n,
        )
        / 1e6
    )


def pem_h2_output_twh(n, pem_link):
    # PyPSA Link convention:
    # p1 is negative when energy is injected into bus1.
    h2 = (
        -n.links_t.p1[pem_link]
    ).clip(lower=0.0)

    return (
        annual_sum(h2, n)
        / 1e6
    )


# =============================================================================
# Load networks
# =============================================================================

required = [
    S1_PATH,
    S5_PATH,
]

for s in SCENARIOS.values():
    required += [
        s["old"],
        s["new"],
    ]

for path in required:
    if not path.exists():
        raise FileNotFoundError(path)


print("Loading S1 and S5...")

s1 = pypsa.Network(S1_PATH)
s5 = pypsa.Network(S5_PATH)

(
    PEM_LINK,
    PEM_STORE,
    S5_NEW_LINKS,
    S5_NEW_STORES,
) = identify_pem_components(
    s1,
    s5,
)


print(
    "\n"
    "============================================================"
)
print("PEM COMPONENT IDENTIFICATION")
print(
    "============================================================"
)

print("S5-specific Links:")
for x in S5_NEW_LINKS:
    print(" ", x)

print("S5-specific Stores:")
for x in S5_NEW_STORES:
    print(" ", x)

print("\nIdentified PEM Link:")
print(" ", PEM_LINK)

print("Identified PEM Store:")
print(" ", PEM_STORE)


# =============================================================================
# Regression
# =============================================================================

rows = []


def add_row(
    scenario,
    metric,
    old_value,
    new_value,
):
    delta = new_value - old_value

    relative = (
        100.0
        * delta
        / abs(old_value)
        if abs(old_value) > 1e-12
        else np.nan
    )

    rows.append(
        {
            "Scenario": scenario,
            "Metric": metric,
            "Old Generator": old_value,
            "New Link": new_value,
            "Delta": delta,
            "Delta [%]": relative,
        }
    )


for scenario, paths in SCENARIOS.items():

    print(
        "\n"
        "============================================================"
    )
    print(f"LOADING {scenario}")
    print(
        "============================================================"
    )

    old = pypsa.Network(
        paths["old"]
    )

    new = pypsa.Network(
        paths["new"]
    )

    metrics = {
        "BTC electricity [TWh/a]": (
            btc_energy_mwh(
                old,
                "old",
            )
            / 1e6,

            btc_energy_mwh(
                new,
                "new",
            )
            / 1e6,
        ),

        "BTC capacity factor [%]": (
            btc_capacity_factor(
                old,
                "old",
            ),

            btc_capacity_factor(
                new,
                "new",
            ),
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
            upstream_system_cost(
                old,
                "old",
            )
            / 1e9,

            upstream_system_cost(
                new,
                "new",
            )
            / 1e9,
        ),
    }

    if scenario == "S6":

        metrics[
            "PEM electricity [TWh/a]"
        ] = (
            pem_electricity_twh(
                old,
                PEM_LINK,
            ),
            pem_electricity_twh(
                new,
                PEM_LINK,
            ),
        )

        metrics[
            "H2 product [TWh_LHV/a]"
        ] = (
            pem_h2_output_twh(
                old,
                PEM_LINK,
            ),
            pem_h2_output_twh(
                new,
                PEM_LINK,
            ),
        )

        if (
            PEM_STORE is not None
            and PEM_STORE in old.stores.index
            and PEM_STORE in new.stores.index
        ):
            metrics[
                "Final H2 product Store [TWh_LHV]"
            ] = (
                float(
                    old.stores_t.e[
                        PEM_STORE
                    ].iloc[-1]
                )
                / 1e6,

                float(
                    new.stores_t.e[
                        PEM_STORE
                    ].iloc[-1]
                )
                / 1e6,
            )

    for metric, (
        old_value,
        new_value,
    ) in metrics.items():

        add_row(
            scenario,
            metric,
            old_value,
            new_value,
        )

    # ---------------------------------------------------------
    # BTC bookkeeping validation
    # ---------------------------------------------------------

    store_p = (
        new.stores_t.p[
            BTC_STORE
        ]
    )

    service_error = (
        (-store_p)
        - new.links_t.p0[
            NEW_BTC
        ]
        * float(
            new.links.at[
                NEW_BTC,
                "efficiency",
            ]
        )
    ).abs().max()

    print(
        f"\n{scenario} BTC bookkeeping:"
    )

    print(
        "  max Store discharge [MW]:",
        float(store_p.max()),
    )

    print(
        "  max Store charge [MW]:",
        float(-store_p.min()),
    )

    print(
        "  max service balance error [MW]:",
        float(service_error),
    )

    print(
        "  final BTC service [TWh]:",
        float(
            new.stores_t.e[
                BTC_STORE
            ].iloc[-1]
        )
        / 1e6,
    )

    # ---------------------------------------------------------
    # Zero-direct-CO2 constraint check
    # ---------------------------------------------------------

    print(
        f"\n{scenario} global constraints:"
    )

    print(
        new.global_constraints.to_string()
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
    "S4/S6 GENERATOR -> LINK REGRESSION"
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


# =============================================================================
# Automatic validation summary
# =============================================================================

print(
    "\n"
    "============================================================"
)
print("VALIDATION SUMMARY")
print(
    "============================================================"
)

# Physical-system metrics should remain extremely close.
physical = df[
    ~df["Metric"].str.contains(
        "Raw objective",
        regex=False,
    )
]

max_abs_pct = (
    physical["Delta [%]"]
    .abs()
    .replace([np.inf], np.nan)
    .max()
)

print(
    "Largest absolute relative difference "
    "among compared physical KPIs [%]:",
    max_abs_pct,
)

for scenario in ["S4", "S6"]:
    co2 = df[
        (df["Scenario"] == scenario)
        & (
            df["Metric"]
            == "Direct CO2 [Mt/a]"
        )
    ]["New Link"].iloc[0]

    print(
        f"{scenario} new direct CO2 [Mt/a]:",
        co2,
    )

if "S6" in SCENARIOS:
    pem_el = df[
        (df["Scenario"] == "S6")
        & (
            df["Metric"]
            == "PEM electricity [TWh/a]"
        )
    ]["New Link"].iloc[0]

    h2 = df[
        (df["Scenario"] == "S6")
        & (
            df["Metric"]
            == "H2 product [TWh_LHV/a]"
        )
    ]["New Link"].iloc[0]

    print(
        "S6 PEM electricity [TWh/a]:",
        pem_el,
    )

    print(
        "S6 H2 product [TWh_LHV/a]:",
        h2,
    )

print(
    "\nSaved:"
)
print(OUTPUT)
