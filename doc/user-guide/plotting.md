<!-- SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Plotting & Summary Visualization

Visualization is a central part of the PyPSA-Earth workflow. After solving a
network, built-in scripts let you generate geographic maps of the power system,
produce aggregated cost and energy summaries, and explore results
interactively via Jupyter notebooks. This page documents those tools, explains
how to invoke them via Snakemake, and highlights the relevant configuration
options.

## Overview of Visualization Scripts

PyPSA-Earth currently ships three main scripts for plotting and result
analysis, located in the `scripts/` directory:

| Script | Snakemake Rule | Purpose |
|---|---|---|
| `plot_network.py` | `plot_network` | Renders geographic maps of solved networks |
| `make_summary.py` | `make_summary` | Aggregates optimization results into `.csv` summary tables |
| `plot_summary.py` | `plot_summary` | Generates summary charts from the `costs` and `energy` tables |

All three scripts are driven by Snakemake wildcards and the `plotting:`
section of `config.yaml`, so no manual Python invocation is needed for
standard use.

## Network Maps: `plot_network`

The `plot_network` rule creates static geographic maps of a solved network. Bus
sizes represent installed capacity and line widths represent transfer capacity.

### Running it

```bash
# Plot a solved network using the supported attr wildcard p_nom
snakemake -j 1 "results/plots/elec_s_10_ec_lcopt_Co2L-3H_p_nom.pdf"

# Plot the same map as a PNG
snakemake -j 1 "results/plots/elec_s_10_ec_lcopt_Co2L-3H_p_nom.png"
```

The `{attr}` wildcard is currently required and only supports `p_nom`. The
`{ext}` wildcard controls the output format. Supported formats depend on your
matplotlib backend, typically including `pdf`, `png`, and `svg`.

### What is plotted

- Buses as proportionally sized circles, split into pie segments by carrier.
- Lines and links with widths proportional to their transfer capacity.
- A geographic base map whose colors are controlled in the plotting config.

### Relevant config options

```yaml
plotting:
  map:
    figsize: [7, 7]
    boundaries: [-10.2, 29, 35, 72]
    p_nom:
      bus_size_factor: 5.e+4
      linewidth_factor: 3.e+3
    color_geomap:
      ocean: white
      land: whitesmoke
```

!!! tip
    Adjust `boundaries` to zoom into your region of interest. If bus circles or
    line widths look too large or too small for your clustering level, tune
    `bus_size_factor` and `linewidth_factor`.

## Summary Generation: `make_summary`

The `make_summary` rule aggregates solved networks into `.csv` files stored
under `results/summaries/`. These files can be consumed by `plot_summary.py`
or read directly in a notebook or data analysis workflow.

### Running it

```bash
# Summarize all solved networks defined in config.yaml
snakemake -j 1 make_summary

# Summarize a specific scenario for one country (for example Nigeria)
snakemake -j 1 "results/summaries/elec_s_all_ec_lall_Co2L-3H_NG"
```

The `{country}` wildcard narrows the summary to buses belonging to a given
two-letter ISO country code. Use `all` to include every country in the model.

### Output files

`make_summary.py` writes several `.csv` files per scenario combination:

| File | Description |
|---|---|
| `costs.csv` | Annualized system costs by technology |
| `curtailment.csv` | Curtailed energy by carrier |
| `energy.csv` | Annual generation and dispatch by carrier |
| `capacity.csv` | Installed capacity by carrier |
| `supply.csv` | Supply grouped by bus carrier and component |
| `supply_energy.csv` | Supply energy grouped by bus carrier and component |
| `prices.csv` | Mean marginal prices by bus carrier |
| `weighted_prices.csv` | Load-weighted prices by carrier |
| `metrics.csv` | High-level system metrics such as line volume and CO2 shadow prices |

### Relevant config options

```yaml
costs:
  default_exchange_rate: 0.7532
  discountrate: 0.07

electricity:
  max_hours:
    battery: 6
    H2: 168
```

## Summary Plots: `plot_summary`

Once `make_summary` has run, the `plot_summary` rule reads the summary tables
and produces charts for the `costs` and `energy` summaries.

### Running it

```bash
# Generate all configured summary plots
snakemake -j 1 plot_summary

# Plot the cost summary for a specific country as PDF
snakemake -j 1 "results/plots/summary_costs_elec_s_all_ec_lall_Co2L-3H_NG.pdf"
```

### What is plotted

- Cost summary as a stacked bar chart of annualized system cost by technology.
- Energy summary as a stacked chart of annual generation and storage dispatch.

Small contributors can be grouped into `Other` using the configured thresholds:

```yaml
plotting:
  costs_max: 10
  costs_threshold: 0.2
  energy_max: 20000
  energy_min: -20000
  energy_threshold: 15
```

## Technology Color Palette

All plotting scripts share a color palette defined under
`plotting: tech_colors:` in `config.default.yaml`. These colors are applied
consistently across network maps and summary charts.

Key colors from the default configuration:

| Carrier | Config key | Default color |
|---|---|---|
| Onshore Wind | `onwind` | `#235ebc` |
| Offshore Wind (AC) | `offwind-ac` | `#6895dd` |
| Offshore Wind (DC) | `offwind-dc` | `#74c6f2` |
| Solar PV | `solar` | `#f9d002` |
| Hydro | `hydro` | `#08ad97` |
| Pumped Hydro Storage | `PHS` | `#08ad97` |
| Run of River | `ror` | `#4adbc8` |
| OCGT | `OCGT` | `#d35050` |
| CCGT | `CCGT` | `#b80404` |
| nuclear | `nuclear` | `#ff9000` |
| Coal | `coal` | `#707070` |
| Oil | `oil` | `#262626` |
| Battery | `battery` | `slategray` |
| Hydrogen | `H2` | `#ea048a` |

To override a color, add or edit the entry under `plotting: tech_colors:` in
your `config.yaml`:

```yaml
plotting:
  tech_colors:
    solar: "#e8c800"
    onwind: "#1a3e8a"
```

The helper `_helpers.sanitize_carriers()` in `scripts/_helpers.py` reads these
values and assigns them to the PyPSA network's `carriers` component, making
them automatically available to `n.plot()`.

## Technology Groupings

Plotting scripts sort carriers into named groups that control ordering and
labelling of stacked charts. The default groupings are defined under
`plotting:` in `config.default.yaml` and can be overridden in `config.yaml`.

```yaml
plotting:
  vre_techs:
    - onwind
    - offwind-ac
    - offwind-dc
    - solar
    - ror
  conv_techs:
    - OCGT
    - CCGT
    - nuclear
    - Nuclear
    - coal
    - oil
  storage_techs:
    - hydro+PHS
    - battery
    - H2
  renewable_storage_techs:
    - PHS
    - hydro
```

## Jupyter Notebook Examples

The [pypsa-meets-earth/documentation](https://github.com/pypsa-meets-earth/documentation)
repository contains Jupyter notebooks for exploring PyPSA-Earth outputs. These
examples are useful for validation workflows and custom post-processing beyond
the built-in Snakemake targets.

```bash
git clone https://github.com/pypsa-meets-earth/documentation.git
cd documentation
jupyter lab
```

## Custom Plotting with PyPSA

For analyses beyond the built-in scripts, load a solved network directly in
Python and use PyPSA's `n.plot()` method.

### Load a solved network

```python
import pypsa

n = pypsa.Network("results/networks/elec_s_10_ec_lcopt_Co2L-3H.nc")
```

### Static network map

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 8))

bus_sizes = n.generators.groupby(["bus", "carrier"])["p_nom_opt"].sum() / 5e4

n.plot(
    ax=ax,
    bus_sizes=bus_sizes,
    bus_colors=n.carriers["color"],
    line_widths=n.lines["s_nom"] / 3e3,
    geomap=True,
)
plt.tight_layout()
plt.savefig("my_network_map.png", dpi=150)
```

### Generation mix bar chart

```python
import matplotlib.pyplot as plt

gen = (
    n.generators_t.p.multiply(n.snapshot_weightings.generators, axis=0)
    .sum()
    .groupby(n.generators.carrier)
    .sum()
    / 1e6
)

colors = n.carriers.loc[gen.index, "color"]

gen.plot.bar(color=colors, figsize=(10, 5))
plt.ylabel("Annual generation (TWh)")
plt.title("Generation mix")
plt.tight_layout()
plt.savefig("generation_mix.png", dpi=150)
```

### Line loading map

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 8))

loading = n.lines_t.p0.abs().divide(n.lines.s_nom, axis=1).mean()

n.plot(
    ax=ax,
    line_colors=loading,
    line_cmap="YlOrRd",
    line_widths=1.5,
    geomap=True,
)
plt.savefig("line_loading.png", dpi=150)
```

!!! note
    For a full reference of `n.plot()` parameters, see the
    [PyPSA plotting documentation](https://pypsa.readthedocs.io/en/latest/user-guide/plotting/).

## Quick Reference

| Task | Command |
|---|---|
| Plot a solved network map | `snakemake -j 1 "results/plots/elec_s_10_ec_lcopt_Co2L-3H_p_nom.pdf"` |
| Generate summary CSVs | `snakemake -j 1 make_summary` |
| Country-specific summary | `snakemake -j 1 "results/summaries/elec_s_all_ec_lall_Co2L-3H_NG"` |
| Plot summary charts | `snakemake -j 1 plot_summary` |
| Plot a specific summary figure | `snakemake -j 1 "results/plots/summary_costs_elec_s_all_ec_lall_Co2L-3H_NG.pdf"` |

Explore valid wildcard combinations for your configuration with a dry-run:

```bash
snakemake -j 1 plot_summary -n
```

For more on the wildcards used in these rules, see [Wildcards](wildcards.md).
