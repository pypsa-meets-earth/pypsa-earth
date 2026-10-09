from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path("results/scenarios")

KPI_FILE = ROOT / "final_S0_S6_flexible_kpis.csv"

INTERACTION_FILE = ROOT / "final_S0_S6_interaction_terms.csv"

FIG_DIR = ROOT / "thesis_figures"

TABLE_DIR = ROOT / "thesis_tables"

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


if not KPI_FILE.exists():
    raise FileNotFoundError(KPI_FILE)

if not INTERACTION_FILE.exists():
    raise FileNotFoundError(INTERACTION_FILE)


df = pd.read_csv(
    KPI_FILE,
    index_col=0,
)

interaction = pd.read_csv(
    INTERACTION_FILE,
    index_col=0,
)


# =============================================================================
# Plot style
# =============================================================================

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
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


def set_signed_limits(
    ax,
    values,
    lower_pad=0.18,
    upper_pad=0.18,
):
    """
    Explicitly reserve space above and below bars for labels.
    """

    values = np.asarray(
        values,
        dtype=float,
    )

    values = values[np.isfinite(values)]

    if len(values) == 0:
        return

    vmin = min(
        0.0,
        float(values.min()),
    )

    vmax = max(
        0.0,
        float(values.max()),
    )

    span = vmax - vmin

    if span == 0.0:
        span = max(
            abs(vmin),
            abs(vmax),
            1.0,
        )

    if vmin < 0.0 and vmax > 0.0:
        ymin = vmin - lower_pad * span

        ymax = vmax + upper_pad * span

    elif vmax <= 0.0:

        ymin = vmin - 0.22 * span

        ymax = 0.12 * span

    else:

        ymin = -0.10 * span

        ymax = vmax + upper_pad * span

    ax.set_ylim(
        ymin,
        ymax,
    )


def annotate_bars(
    ax,
    bars,
    decimals=3,
    small_threshold=None,
):
    """
    Place all value labels clearly outside the bars.

    If small_threshold is given, non-zero values with an
    absolute magnitude below that threshold are shown as
    <threshold or >-threshold instead of appearing as zero.
    """

    ymin, ymax = ax.get_ylim()

    span = ymax - ymin

    tiny = 0.015 * span

    for bar in bars:

        value = float(bar.get_height())

        if not np.isfinite(value):
            continue

        x = bar.get_x() + bar.get_width() / 2

        if (
            small_threshold is not None
            and value != 0.0
            and abs(value) < small_threshold
        ):
            if value > 0:
                label = f"<{small_threshold:.3f}"
            else:
                label = f">-{small_threshold:.3f}"
        else:
            label = f"{value:.{decimals}f}"

        if abs(value) < tiny:

            if value < 0:

                y = -0.045 * span

                va = "top"

            else:

                y = 0.045 * span

                va = "bottom"

        elif value > 0:

            y = value + 0.025 * span

            va = "bottom"

        else:

            y = value - 0.025 * span

            va = "top"

        ax.text(
            x,
            y,
            label,
            ha="center",
            va=va,
            fontsize=8,
            clip_on=False,
            zorder=20,
        )


def save_figure(
    fig,
    stem,
    rect,
):
    png = FIG_DIR / f"{stem}.png"

    pdf = FIG_DIR / f"{stem}.pdf"

    fig.tight_layout(rect=rect)

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


# =============================================================================
# FIGURE 09
# Baseline-relative system response
# =============================================================================

BASELINES = {
    "S2": "S0",
    "S3": "S0",
    "S4": "S1",
    "S5": "S1",
    "S6": "S1",
}


comparison_labels = {
    "S2": "S2-S0",
    "S3": "S3-S0",
    "S4": "S4-S1",
    "S5": "S5-S1",
    "S6": "S6-S1",
}


response_rows = []


for scenario, baseline in BASELINES.items():

    solar_delta = (
        df.loc[
            scenario,
            "Solar capacity [GW]",
        ]
        - df.loc[
            baseline,
            "Solar capacity [GW]",
        ]
    )

    wind_delta = (
        df.loc[
            scenario,
            "Wind capacity [GW]",
        ]
        - df.loc[
            baseline,
            "Wind capacity [GW]",
        ]
    )

    battery_delta = (
        df.loc[
            scenario,
            "Battery energy capacity [GWh]",
        ]
        - df.loc[
            baseline,
            "Battery energy capacity [GWh]",
        ]
    )

    curtailment_delta = (
        df.loc[
            scenario,
            "Solar curtailment [TWh]",
        ]
        + df.loc[
            scenario,
            "Wind curtailment [TWh]",
        ]
        - df.loc[
            baseline,
            "Solar curtailment [TWh]",
        ]
        - df.loc[
            baseline,
            "Wind curtailment [TWh]",
        ]
    )

    response_rows.append(
        {
            "Scenario": scenario,
            "Baseline": baseline,
            "Comparison": comparison_labels[scenario],
            "Solar capacity change [GW]": solar_delta,
            "Wind capacity change [GW]": wind_delta,
            "Battery energy capacity change [GWh]": battery_delta,
            "Solar+wind curtailment change [TWh/a]": curtailment_delta,
        }
    )


response = pd.DataFrame(response_rows).set_index("Scenario")


response.to_csv(TABLE_DIR / "figure_09_baseline_relative_system_response_data.csv")


xlabels = [
    response.loc[
        scenario,
        "Comparison",
    ]
    for scenario in BASELINES
]


panels = [
    (
        "Solar capacity change [GW]",
        r"$\Delta$ solar capacity [GW]",
        "Solar capacity",
        3,
    ),
    (
        "Wind capacity change [GW]",
        r"$\Delta$ wind capacity [GW]",
        "Wind capacity",
        3,
    ),
    (
        "Battery energy capacity change [GWh]",
        r"$\Delta$ battery energy capacity [GWh]",
        "Battery energy capacity",
        3,
    ),
    (
        "Solar+wind curtailment change [TWh/a]",
        r"$\Delta$ solar + wind curtailment [TWh/a]",
        "Variable-renewable curtailment",
        3,
    ),
]


fig, axes = plt.subplots(
    2,
    2,
    figsize=(
        11,
        7.8,
    ),
)

axes = axes.ravel()


for ax, (
    column,
    ylabel,
    title,
    decimals,
) in zip(
    axes,
    panels,
):

    values = [
        float(
            response.loc[
                scenario,
                column,
            ]
        )
        for scenario in BASELINES
    ]

    bars = ax.bar(
        xlabels,
        values,
    )

    ax.axhline(
        0.0,
        linewidth=0.9,
    )

    ax.set_ylabel(ylabel)

    ax.set_title(
        title,
        pad=8,
    )

    clean_axis(ax)

    set_signed_limits(
        ax,
        values,
    )

    annotate_bars(
        ax,
        bars,
        decimals=decimals,
        small_threshold=0.001,
    )


fig.suptitle(
    (
        "System response to flexible electricity demand "
        "relative to the corresponding baseline"
    ),
    fontsize=11,
    y=0.985,
)


fig.text(
    0.5,
    0.018,
    (
        "S2 and S3 are relative to S0 (reference); "
        "S4, S5 and S6 are relative to S1 "
        "(zero-direct-CO$_2$)."
    ),
    ha="center",
    fontsize=9,
)


save_figure(
    fig,
    "figure_09_baseline_relative_system_response",
    rect=(
        0.02,
        0.055,
        0.99,
        0.955,
    ),
)


# =============================================================================
# FIGURE 10
# S6 combined response vs additive expectation
# =============================================================================

required_columns = [
    "S1 baseline",
    "S6 BTC+H2",
    "Additive expectation",
    "Interaction S6-(S4+S5-S1)",
]


for column in required_columns:

    if column not in interaction.columns:

        raise KeyError(f"Missing interaction column: {column}")


def actual_effect(kpi):

    return (
        interaction.loc[
            kpi,
            "S6 BTC+H2",
        ]
        - interaction.loc[
            kpi,
            "S1 baseline",
        ]
    )


def additive_effect(kpi):

    return (
        interaction.loc[
            kpi,
            "Additive expectation",
        ]
        - interaction.loc[
            kpi,
            "S1 baseline",
        ]
    )


def interaction_effect(kpi):

    return interaction.loc[
        kpi,
        "Interaction S6-(S4+S5-S1)",
    ]


rows = []


cost_kpi = "Power-system cost proxy [EUR bn]"


rows.append(
    {
        "Metric": "Power-system cost",
        "Unit": "million EUR/a",
        "Actual combined effect": actual_effect(cost_kpi) * 1000.0,
        "Additive expectation": additive_effect(cost_kpi) * 1000.0,
        "Interaction": interaction_effect(cost_kpi) * 1000.0,
    }
)


for (
    metric_name,
    kpi,
    unit,
) in [
    (
        "Solar capacity",
        "Solar capacity [GW]",
        "GW",
    ),
    (
        "Wind capacity",
        "Wind capacity [GW]",
        "GW",
    ),
    (
        "Battery energy capacity",
        "Battery energy capacity [GWh]",
        "GWh",
    ),
]:

    rows.append(
        {
            "Metric": metric_name,
            "Unit": unit,
            "Actual combined effect": actual_effect(kpi),
            "Additive expectation": additive_effect(kpi),
            "Interaction": interaction_effect(kpi),
        }
    )


solar_curt = "Solar curtailment [TWh]"

wind_curt = "Wind curtailment [TWh]"


rows.append(
    {
        "Metric": "Solar + wind curtailment",
        "Unit": "TWh/a",
        "Actual combined effect": actual_effect(solar_curt) + actual_effect(wind_curt),
        "Additive expectation": additive_effect(solar_curt)
        + additive_effect(wind_curt),
        "Interaction": interaction_effect(solar_curt) + interaction_effect(wind_curt),
    }
)


interaction_plot_data = pd.DataFrame(rows).set_index("Metric")


interaction_plot_data.to_csv(
    TABLE_DIR / "figure_10_S6_interaction_nonadditivity_data.csv"
)


metrics = [
    "Power-system cost",
    "Solar capacity",
    "Wind capacity",
    "Battery energy capacity",
    "Solar + wind curtailment",
]


decimals = {
    "Power-system cost": 3,
    "Solar capacity": 4,
    "Wind capacity": 4,
    "Battery energy capacity": 4,
    "Solar + wind curtailment": 4,
}


fig, axes = plt.subplots(
    2,
    3,
    figsize=(
        12,
        8.2,
    ),
)

axes = axes.ravel()


for ax, metric in zip(
    axes[:5],
    metrics,
):

    row = interaction_plot_data.loc[metric]

    values = [
        float(row["Actual combined effect"]),
        float(row["Additive expectation"]),
    ]

    bars = ax.bar(
        [
            "Actual S6\nvs S1",
            "Additive\nS4 + S5",
        ],
        values,
    )

    ax.axhline(
        0.0,
        linewidth=0.9,
    )

    ax.set_title(
        metric,
        pad=8,
    )

    ax.set_ylabel(row["Unit"])

    clean_axis(ax)

    set_signed_limits(
        ax,
        values,
        lower_pad=0.20,
        upper_pad=0.22,
    )

    annotate_bars(
        ax,
        bars,
        decimals=decimals[metric],
    )


# =============================================================================
# Sixth panel: interaction summary
# =============================================================================

axes[-1].axis("off")


summary_lines = [
    (r"Interaction term: " r"$I=(S6-S1)-[(S4-S1)+(S5-S1)]$"),
    "",
    (
        "Power-system cost: "
        f"{interaction_plot_data.loc['Power-system cost', 'Interaction']:+.4f} "
        "million EUR/a"
    ),
    (
        "Solar capacity: "
        f"{interaction_plot_data.loc['Solar capacity', 'Interaction']:+.4f} GW"
    ),
    (
        "Wind capacity: "
        f"{interaction_plot_data.loc['Wind capacity', 'Interaction']:+.4f} GW"
    ),
    (
        "Battery energy capacity: "
        f"{interaction_plot_data.loc['Battery energy capacity', 'Interaction']:+.4f} "
        "GWh"
    ),
    (
        "Solar + wind curtailment: "
        f"{interaction_plot_data.loc['Solar + wind curtailment', 'Interaction']:+.4f} "
        "TWh/a"
    ),
]


axes[-1].text(
    0.02,
    0.94,
    "Non-additive interaction",
    transform=axes[-1].transAxes,
    ha="left",
    va="top",
    fontsize=10.5,
    fontweight="bold",
)

axes[-1].text(
    0.02,
    0.86,
    "\n".join(summary_lines),
    transform=axes[-1].transAxes,
    ha="left",
    va="top",
    fontsize=9,
    linespacing=1.45,
)


fig.suptitle(
    ("Combined BTC-H$_2$ response in S6: " "actual effect versus additive expectation"),
    fontsize=11,
    y=0.985,
)


fig.text(
    0.5,
    0.018,
    ("Actual combined effect: S6-S1. " "Additive expectation: (S4-S1) + (S5-S1)."),
    ha="center",
    fontsize=9,
)


save_figure(
    fig,
    "figure_10_S6_interaction_nonadditivity",
    rect=(
        0.02,
        0.055,
        0.99,
        0.955,
    ),
)


print()
print("=" * 100)
print("ADDITIONAL THESIS FIGURES REGENERATED")
print("=" * 100)
