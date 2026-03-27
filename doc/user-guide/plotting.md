<!-- SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

<!-- SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors -->
<!-- SPDX-License-Identifier: CC-BY-4.0 -->

# Plotting & Summary Visualization

Visualization is a central part of the PyPSA-Earth workflow. After solving a network, several built-in scripts let you generate geographic maps of the power system, produce aggregated cost and energy summaries, and explore results interactively via Jupyter notebooks. This page documents those tools, explains how to invoke them via Snakemake, and provides guidance on customizing outputs for your own analysis.

---

## Overview of Visualization Scripts

PyPSA-Earth ships four main scripts for plotting and result analysis, located in the `scripts/` directory:

| Script | Snakemake Rule | Purpose |
|---|---|---|
| `plot_network.py` | `plot_network` | Renders geographic maps of solved networks (capacity, loading, etc.) |
| `make_summary.py` | `make_summary` | Aggregates optimization results into `.csv` summary tables |
| `plot_summary.py` | `plot_summary` | Generates bar/area charts from the summary `.csv` files |
| `plot_p_nom_max.py` | `plot_p_nom_max` | Maps regionally disaggregated renewable potentials by technology |

All scripts are driven by Snakemake wildcards and the `plotting:` section of `config.yaml`, so no manual Python invocation is needed for standard use.

---

## Network Maps: `plot_network`

The `plot_network` rule creates **static geographic maps** of a solved network. Bus sizes represent installed capacity (`p_nom_opt`) and line widths represent transfer capacity.

### Running it

```bash
# Plot a solved network (PDF output)
snakemake -j 1 results/plots/elec_s_10_ec_lcopt_Co2L-3H.pdf

# Plot with PNG output
snakemake -j 1 results/plots/elec_s_10_ec_lcopt_Co2L-3H.png
```

The `{ext}` wildcard controls the output format. Supported formats depend on your matplotlib backend — typically `pdf`, `png`, and `svg`.

### What is plotted

The map shows:

- **Buses** as proportionally sized circles, split into pie segments by carrier (generation mix at each node). Size is scaled by `bus_size_factor` in the config.
- **Lines and links** with widths proportional to their transfer capacity, scaled by `linewidth_factor`.
- **Geographic base map** with ocean and land colors set by `color_geomap`.

### Relevant config options

```yaml
plotting:
  map:
    figsize: [7, 7]
    boundaries: [-10.2, 29, 35, 72]   # [lon_min, lon_max, lat_min, lat_max]
    p_nom:
      bus_size_factor: 5.e+4
      linewidth_factor: 3.e+3
    color_geomap:
      ocean: white
      land: whitesmoke
```

!!! tip
    Adjust `boundaries` to zoom into your region of interest. Adjust `bus_size_factor` and `linewidth_factor` if bus circles or line widths appear too large or too small for your network resolution.

---

## Summary Generation: `make_summary`

The `make_summary` rule aggregates key results from one or more solved networks into **`.csv` files** stored under `results/summaries/`. These files contain quantities such as:

- Total system cost broken down by technology
- Installed capacity (optimal `p_nom_opt`) per carrier
- Annual energy generation and curtailment per carrier
- CO₂ emissions

### Running it

```bash
# Summarize all solved networks defined in config.yaml
snakemake -j 1 make_summary

# Summarize only buses for a specific country (e.g. Nigeria)
snakemake -j 1 results/summaries/elec_s_all_lall_Co2L-3H_NG
```

The `{country}` wildcard narrows the summary to buses belonging to the specified two-letter ISO country code. Use `all` to include every country in the model.

### Output files

`make_summary.py` writes several `.csv` files per scenario combination:

| File | Description |
|---|---|
| `costs.csv` | Annualized system costs in EUR/year by technology |
| `capacities.csv` | Installed capacity in GW by carrier |
| `energy.csv` | Annual generation in TWh by carrier |
| `curtailment.csv` | Curtailed energy in TWh by carrier |

These files feed directly into `plot_summary.py` but can also be read into a Jupyter notebook or any data analysis tool.

### Relevant config options

```yaml
costs:
  USD2013_to_EUR2013: 0.7532
  discountrate: 0.07

electricity:
  max_hours:
    battery: 6
    H2: 168
```

---

## Summary Plots: `plot_summary`

Once `make_summary` has run, the `plot_summary` rule reads the `.csv` files and produces **bar and area charts** of system costs and energy mixes.

### Running it

```bash
# Generate all summary plots
snakemake -j 1 plot_summary

# Plot for a specific country
snakemake -j 1 results/plots/summary_elec_s_all_lall_Co2L-3H_NG.pdf
```

### What is plotted

- **Cost summary** — stacked bar chart of annualized system cost in EUR/year, broken down by technology category (VRE, conventional, storage, transmission).
- **Energy summary** — stacked bar/area chart of annual electricity generation and storage dispatch in TWh.

Thresholds control which small contributors are shown explicitly vs. grouped into "Other":

```yaml
plotting:
  costs_max: 10          # upper limit of cost axis (billion EUR)
  costs_threshold: 0.2   # technologies below this share are grouped as "Other"
  energy_max: 20000      # upper limit of energy axis (TWh)
  energy_min: -20000     # lower limit (negative = storage charging / curtailment)
  energy_threshold: 15   # technologies below this TWh value are grouped as "Other"
```

---

## Renewable Potential Maps: `plot_p_nom_max`

This rule maps the **maximum installable capacity** for a given renewable technology across all network regions.

### Running it

```bash
# Plot solar potential map
snakemake -j 1 results/plots/p_nom_max_solar.pdf

# Plot onshore wind potential map
snakemake -j 1 results/plots/p_nom_max_onwind.pdf
```

The `{technology}` wildcard accepts `onwind`, `offwind-ac`, `offwind-dc`, and `solar`.

---

## Technology Color Palette

All plotting scripts share a color palette defined under `plotting: tech_colors:` in `config.default.yaml`. Colors are specified as hex codes and applied consistently across maps and summary charts.

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
| OCGT | `OCGT` | `#e0986c` |
| CCGT | `CCGT` | `#e8d04b` |
| Nuclear | `nuclear` | `#ff8c00` |
| Coal | `coal` | `#545454` |
| Oil | `oil` | `#808080` |
| Battery | `battery` | `#b8ea04` |
| Hydrogen | `H2` | `#ea048a` |

To override a color, add or edit the entry under `plotting: tech_colors:` in your `config.yaml`:

```yaml
plotting:
  tech_colors:
    solar: "#e8c800"    # custom yellow for solar
    onwind: "#1a3e8a"   # custom dark blue for onshore wind
```

The helper `_helpers.sanitize_carriers()` in `scripts/_helpers.py` reads these values and assigns them to the PyPSA network's `carriers` component, making them automatically available to `n.plot()`.

---

## Technology Groupings

Plotting scripts sort carriers into named groups that control the ordering and labelling of stacked charts. Edit these lists in `config.yaml` to match your carrier set:

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

---

## Jupyter Notebook Examples

The [pypsa-meets-earth/documentation](https://github.com/pypsa-meets-earth/documentation) repository contains hands-on Jupyter notebooks for exploring PyPSA-Earth outputs. Relevant notebooks include:

- **Network validation** — compares model outputs to public power system statistics (e.g. Nigeria case study).
- **Capacity validation** — inspects installed capacity against reference datasets.
- **Demand validation** — checks power demand values and load profiles.
- **Renewable potential analysis** — visualizes wind and solar potential time series.

```bash
git clone https://github.com/pypsa-meets-earth/documentation.git
cd documentation
jupyter lab
```

---

## Custom Plotting with PyPSA

For analyses beyond the built-in scripts, load any solved network directly in Python and use PyPSA's built-in `n.plot()` method.

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
    / 1e6  # MWh → TWh
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
    For a full reference of `n.plot()` parameters (bus sizes, colormaps, flow arrows, etc.) see the [PyPSA plotting documentation](https://pypsa.readthedocs.io/en/latest/user-guide/plotting/).

---

## Quick Reference

| Task | Command |
|---|---|
| Plot a solved network map | `snakemake -j 1 results/plots/<network_name>.pdf` |
| Generate summary CSVs | `snakemake -j 1 make_summary` |
| Country-specific summary | `snakemake -j 1 results/summaries/<network_name>_<CC>` |
| Plot summary charts | `snakemake -j 1 plot_summary` |
| Plot renewable potential | `snakemake -j 1 results/plots/p_nom_max_<technology>.pdf` |

Explore all valid wildcard combinations for your `config.yaml` with a dry-run:

```bash
snakemake -j 1 plot_summary -n
```

For more on the wildcards used in these rules, see [Wildcards](wildcards.md).
