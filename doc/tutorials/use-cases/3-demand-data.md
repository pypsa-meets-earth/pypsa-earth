<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 3: Calibrate Electricity Demand

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) and [Part 2](2-analyze-results.md). You should still have `config.KZ.yaml` in the project root and a baseline solved network for comparison.

## Introduction

In [Part 2](2-analyze-results.md) we compared the baseline model against **2020** national statistics. Total electricity demand looked roughly plausible (~107 TWh), but we had not chosen that number on purpose — it came from PyPSA-Earth defaults. In this tutorial we take control of demand: where the hourly profiles come from, how to adjust them with a multiplier, and how to align the annual total with a reference year before we calibrate generation and capacities in later parts.

This series builds an **electricity-only (power system) model**. We do not switch on sector coupling (heat, transport, industry, hydrogen, and so on). Everything in this part lives under `load_options` in the config and the `build_demand_profiles` rule — not under `sector`.

---

## Where demand enters the workflow

Electricity demand is not read directly from national statistics inside the optimiser. Two rules handle it — both before the network is simplified or solved:

```
build_demand_profiles  →  resources/KZ/demand_profiles.csv
        ↓
add_electricity        →  networks/KZ/elec.nc
        ↓
… simplify, cluster, prepare, solve …
```

| Rule | What it does |
|---|---|
| **`build_demand_profiles`** | Loads GEGIS or DemandCast, applies `scale`, and writes an hourly profile per bus → `demand_profiles.csv` |
| **`add_electricity`** | Attaches those values as PyPSA **Load** components on the network |

All demand settings live under **`load_options`** in the config. They affect **`build_demand_profiles` only**.

```yaml
--8<-- "configtables/snippets/load_options.yaml"
```

For the full parameter list see [load_options in the Configuration reference](../../user-guide/configuration.md#load_options).

---

## Step 1: Choose a demand source

PyPSA-Earth builds hourly load from a **global demand dataset** — pre-calculated country-level electricity consumption, disaggregated using weather patterns — not from national statistics directly. Two datasets are available via `load_options.source`:

| Source | Config value | Weather years | Data location |
|---|---|---|---|
| **[GEGIS](http://dx.doi.org/10.1016/j.esr.2020.100606)** | `gegis` | 2011, 2013, 2018 | `data/ssp2-2.6/{prediction_year}/era5_{weather_year}/` (regional `.nc` / `.csv`, e.g. `Asia.nc` for KZ) |
| **[DemandCast](https://github.com/open-energy-transition/demandcast)** | `demcast` | 2000–2024 (must match `snapshots` year) | `data/demand/forecasts_on_historical_period.parquet` (via databundle) |

**GEGIS** is the default. Profiles are pre-calculated from population and GDP using the [synde](https://github.com/euronion/synde) workflow; `prediction_year` selects the SSP pathway folder (e.g. `2030` under `ssp2-2.6`).

**DemandCast** ([Steijn et al., 2025](https://arxiv.org/abs/2510.08000)) adds wider **time** coverage (25 years vs three for GEGIS) and better **spatial** coverage for countries where GEGIS had gaps or zero demand. For most countries with good reference data, hourly patterns correlate strongly with GEGIS; annual totals are still adjusted with `scale` either way — see the comparison and Ember validation in [#1724](https://github.com/pypsa-meets-earth/pypsa-earth/issues/1724).

For Kazakhstan we stay with the default **`gegis`**. The country is part of the **Asia** regional file, which Snakemake selects automatically from your `countries` list.

Add the source to `config.KZ.yaml`:

```yaml
load_options:
  source: gegis
```

---

## Step 2: Set weather and prediction year (GEGIS)

These two settings apply when `source: gegis`. They pick **which pre-calculated file** Snakemake loads from:

```
data/ssp2-2.6/{prediction_year}/era5_{weather_year}/Asia.nc
```

| Key | Meaning |
|---|---|
| `weather_year` | **Weather year** — which calendar year's temperature patterns shape the hourly load curve. GEGIS provides 2011, 2013, and 2018. This must match the **year your model simulates**. With the default config that is **2013** (`start: "2013-01-01"`, `end: "2014-01-01"` near the top of `config.default.yaml`), so `weather_year: 2013` is the right choice. |
| `prediction_year` | **Pathway year** — which GEGIS future-demand scenario to use. `prediction_year: 2030` loads demand built for **2030** under SSP2-2.6 (population and GDP assumptions for that year). It sets the socioeconomic level in the profile; it is **not** the calendar year you validate against in Part 2. |

In short: **`weather_year`** drives *when* the hours look like (hot summers, cold winters); **`prediction_year`** drives *how much* demand GEGIS expects in that SSP future (e.g. 2030).

```yaml
load_options:
  source: gegis
  weather_year: 2013
  prediction_year: 2030
```

!!! note "Why `2030`, not `2020`?"
    GEGIS does **not** ship a 2020 pathway folder. Under `data/ssp2-2.6/` the available `prediction_year` values are **2030, 2040, 2050, and 2100** — future SSP scenario years, not historical calendar years.

    We keep **`prediction_year: 2030`** because it is the PyPSA-Earth default, matches `planning_horizons: [2030]` in `config.default.yaml`, and is the folder included in the standard databundle. That is why the Part 1 baseline was already near **107 TWh** without any tuning.

    We still validate against **2020** statistics via **`scale`** (Step 3), not by changing these years. DemandCast can use **`weather_year: 2020`**, but only if you also move **`snapshots` `start` / `end`** to 2020 — see below.

!!! tip
    Validating against **2020** statistics does not mean setting `weather_year: 2020` or `prediction_year: 2020` under GEGIS — GEGIS offers neither. Keep `weather_year: 2013` and `prediction_year: 2030`, then adjust the **annual total** with `scale` (next step). If you change GEGIS `weather_year` to 2011 or 2018, update `snapshots` `start` / `end` to that same calendar year.

If you experiment with **DemandCast** (`source: demcast`), set `weather_year` to the **calendar year you simulate** — the same year as `snapshots` `start` / `end`. The parquet holds 2000–2024, but `build_demand_profiles` keeps only rows inside the snapshot window; a mismatch (e.g. `weather_year: 2020` with default 2013 snapshots) produces an **empty** `demand_profiles.csv`.

---

## Step 3: Calibrate annual demand with `scale`

Steps 1–2 fixed the **hourly shape** of demand. Step 3 sets the **annual total** to match **KEGOC 2020 (107.3 TWh)**.

**`load_options.scale`** multiplies every hour after the profile loads. The curve shape is unchanged; only the annual total moves. The same `scale` logic applies whether you use GEGIS or DemandCast from Step 1.

You can pass `scale` in **two forms**:

| Form | When to use | Example |
|---|---|---|
| **Single float** | One country in `countries` (this tutorial) | `scale: 1.005` |
| **Per-country dictionary** | Multi-country run — different multiplier per country | `scale: { DEFAULT: 1.0, KZ: 1.005 }` |

With the dictionary form, each key is an ISO country code. **`DEFAULT`** is the fallback multiplier for any country in `countries` that does not have its own entry; if you omit it, PyPSA-Earth uses `1.0`.

**1. Pick target and baseline** — both from [Part 2](2-analyze-results.md):

| Role | Source | TWh |
|---|---|---|
| **Target** (numerator) | [KEGOC 2020](https://ar2020.kegoc.kz/eng/index.html) | **107.3** |
| **Baseline** (denominator) | Model (Part 2) | **106.8** |

Other targets: [EIA](https://www.eia.gov/international/data/country/KAZ) 103.4 TWh, [Ember](https://ember-climate.org/data/data-explorer/) 107.9 TWh.

!!! note "Where does the baseline come from?"
    The GEGIS input file (`data/ssp2-2.6/2030/era5_2013/Asia.csv`) contains about **108.0 TWh** for Kazakhstan. After simplification, load on electrically isolated sub-networks is dropped from the main grid, leaving about **106.8 TWh** in the solved network from Part 2 — a gap of roughly **1.2 TWh**. The `scale` factor in this step absorbs that loss together with any GEGIS-vs-actual discrepancy in one multiplier. Later in this series, when you improve network topology to reconnect isolated regions, that topology gap largely disappears and `scale` mainly corrects the dataset-vs-statistics difference.

**2. Compute the multiplier:**

$$
\text{scale} = \frac{107.3}{106.8} \approx 1.005
$$

```python
scale = 107.3 / 106.8
print(f"scale: {scale:.4f}")  # → 1.005
```

**3. Add `scale` to your `load_options` block**

For Kazakhstan alone, a single float is enough:

```yaml
load_options:
  source: gegis
  weather_year: 2013
  prediction_year: 2030
  scale: 1.005
```

You can [download the file](snippets/config.KZ.demand.yaml){: download="config.KZ.yaml"} and merge it with your existing `config.KZ.yaml`, or add the `load_options` block by hand.

---

## Step 4: Re-run the workflow

You updated `load_options` in Steps 1–3. Run the same target as [Part 1](1-baseline-model.md):

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Snakemake compares your config with the last run and **rebuilds only what changed**. Because `load_options` feeds `build_demand_profiles`, expect that rule to run again, followed by `add_electricity` and everything downstream through `solve_network`. Cutouts, OSM data, and the base network stay cached from Part 1.

!!! tip "Optional: disable cutout download on re-runs"
    If Part 1 completed successfully, add to `config.KZ.yaml`:

    ```yaml
    enable:
      retrieve_cutout: false
    ```

    Demand changes do not require a new ERA5 cutout. This tells Snakemake to use the cached file under `cutouts/KZ/` and avoids another download attempt. See the [FAQ](../../community/faq.md#cutout-download-failed-retrieve_cutout) if cutout retrieval caused trouble in Part 1.

When the run finishes, the updated solved network overwrites the same file you analysed in Part 2:

```
results/KZ/networks/elec_s_10_ec_lcopt_6h.nc
```

**Expected runtime:** much faster than Part 1 — on the order of **7-10 minutes** with HiGHS, because most upstream rules are skipped.

---

## Step 5: Verify in your notebook

Reopen the notebook from Part 2 and reload the network. Total demand should match your KEGOC target:

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")
weights = n.snapshot_weightings.generators
total_TWh = n.loads_t.p_set.multiply(weights, axis=0).sum().sum() / 1e6
print(f"Total annual demand: {total_TWh:.1f} TWh")
```

Expected output:

```
Total annual demand: 107.3 TWh
```

The hourly profile shape is unchanged from Part 2 — only the annual total moved.

---

## Recap

| Step | Config key | Value | Role |
|---|---|---|---|
| 1 | `load_options.source` | `gegis` | Demand dataset |
| 2 | `load_options.weather_year` | `2013` | Hourly shape (match simulation year) |
| 2 | `load_options.prediction_year` | `2030` | GEGIS SSP level (no 2020 folder) |
| 3 | `load_options.scale` | `1.005` or `{ DEFAULT: 1.0, KZ: 1.005 }` | Match 2020 KEGOC total (float for single country; dict for multi-country) |

Demand is now anchored to **2020** consumption statistics. Generation and installed capacity are calibrated in **[Part 4](4-generation-data.md)** — lock the 2020 fleet, filter powerplantmatching, and compare installed capacities to the same validation tables.
