<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 4: Calibrate the Generation Fleet

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 3](3-demand-data.md). Demand should already be calibrated (`load_options.scale: 1.005`) and `config.KZ.yaml` should include the Part 3 settings.

## Introduction

In [Part 2](2-analyze-results.md) the **installed capacity** table showed a large gap for solar and wind: the model reported **~20 GW solar** and **~6 GW wind** against **~1 GW** each in 2020 statistics. Coal and gas were closer to the national report by KEGOC, but still imperfect.

PyPSA-Earth loads the existing fleet from **[powerplantmatching](https://github.com/FRESNA/powerplantmatching)** — the same global database Part 2 already reflected for coal and gas. Solar and wind farms are in there too. What inflated the baseline was not missing data, but default settings that **add** IRENA capacity on top and **allow the optimiser to build more** (IRENA year **2023**, renewables in `extendable_carriers`).

In this tutorial we lock the model to the **2020 fleet**:

1. See how Snakemake builds the plant list and inspect it.
2. Filter plants to units operating in **2020** (`powerplants_filter`).
3. Turn off new build and set IRENA totals for **2020** (`extendable_carriers`, `estimate_renewable_capacities`).
4. Optionally **replace** powerplantmatching with a custom plant list where the global database is incomplete (often gas and smaller thermal units).

After re-running, compare installed capacity against the same Part 2 validation tables. With no extendable carriers, **`p_nom`** and **`p_nom_opt`** should match — unlike the Part 2 baseline, where optimal capacity (~20 GW solar, ~6 GW wind) was far above installed (~1–2 GW each).

!!! note "Part 2 vs Part 4"
    Part 2's validation table used **`p_nom_opt`**, which counted **new build on top of** existing solar and wind. Those already had **`p_nom` ~1–2 GW** from powerplantmatching and IRENA top-up — not zero — but defaults also marked them **extendable**, so the optimiser added much more. After locking the fleet here, installed and optimal capacity should match.

This series stays in the **electricity workflow** (`solve_all_networks`). All settings below live under **`electricity`** in the config.

---

## Where generation enters the workflow

Generation is assembled in two Snakemake rules — both before the network is simplified or solved. **`build_powerplants`** pulls data from powerplantmatching (and optionally your own CSV), then writes one row per plant to disk:

```
build_powerplants  →  resources/KZ/powerplants.csv   ← one row per plant
        ↓
add_electricity    →  networks/KZ/elec.nc
        ↓
… simplify, cluster, prepare, solve …
```

| Rule | What it does |
|---|---|
| **`build_powerplants`** | Pulls powerplantmatching data for your `countries`, applies `powerplants_filter`, optionally **`merge`** or **`replace`** with **`data/custom_powerplants.csv`** → `powerplants.csv` |
| **`add_electricity`** | Reads `powerplants.csv`, attaches demand and generators → `networks/KZ/elec.nc` |

---

## Step 1: Inspect the plant list

The rule **`build_powerplants`** writes **`resources/KZ/powerplants.csv`** — one row per plant. A few lines from a Part 1 run (columns trimmed for readability):

```csv
Name,Fueltype,Technology,Capacity,DateIn,DateOut,Country
Ekibastuz,Hard Coal,Steam Turbine,4000,1980,2025,KZ
Shulbi,Hydro,Reservoir,702,1987,2087,KZ
Almaty,CCGT,CCGT,120,2017,2057,KZ
Zhangiz Tobe Solar Farm,Solar,PV,30,2019,2044,KZ
Aktogay Wind Farm,Wind,Onshore,150,2022,2047,KZ
Taraz Dzhambul Kazakhstan,Oil,Steam Turbine,185,1996,2036,KZ
```

The full file has more columns (`lat`, `lon`, `bus`, …) and on the order of **100 rows** for Kazakhstan. Totals by fuel type (GW):

| Fueltype | GW |
|---|---|
| Hard Coal | 11.72 |
| CCGT | 2.43 |
| Hydro | 3.12 |
| Solar | 1.27 |
| Wind | 1.21 |
| Oil | 1.14 |

Solar and wind are already listed — not only conventional plants. Your file should look similar after Part 1 simulation; open `resources/KZ/powerplants.csv` locally if you want the complete list.

Compare with the **2020 installed capacity** table in [Part 2](2-analyze-results.md). Coal and hydro are roughly the right order of magnitude; **gas** is often under-represented in global databases. The Part 2 **model** solar/wind totals were much higher than KEGOC because of **IRENA gap-fill and extendable build**, not because powerplantmatching had no renewables.

!!! tip
    The log file `logs/KZ/build_powerplants.log` shows how many plants survived the filter. After you change `powerplants_filter` in Step 3, check this log to confirm rows were dropped as expected.

---

## Step 2: Understand the defaults we are changing

PyPSA-Earth defaults in `config.default.yaml` target a **forward-looking** model (plants active around 2023, room for new build). For **2020 validation** we override three keys:

| Config key | Default | Part 4 (2020 validation) | Effect |
|---|---|---|---|
| `electricity.powerplants_filter` | plants active in **2022–2023** | plants operating in **2020** | Filters **every row** in `powerplants.csv` at build time |
| `electricity.extendable_carriers.Generator` | `[solar, onwind, offwind-ac, offwind-dc, OCGT]` | `[]` | No new solar, wind, or gas build |
| `electricity.estimate_renewable_capacities.year` | `2023` | `2020` | IRENA solar/wind totals match the validation year |

Kazakhstan is landlocked — offshore wind profiles are empty and irrelevant here. An empty `Generator` list under `extendable_carriers` is appropriate for locking the full 2020 fleet.

---

## Step 3: Filter plants for 2020

`powerplants_filter` is a [pandas.query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) expression applied in **`build_powerplants`** after loading powerplantmatching data. It runs on **all fuel types** — coal, gas, hydro, solar, and wind.

Keep plants that were **operating in 2020**:

- **`DateIn <= 2020`** — commissioned on or before 2020 (or unknown commission date).
- **`DateOut >= 2020`** — still operating in 2020 (or unknown decommission date).

The `DateOut != DateOut` / `DateIn != DateIn` clauses treat **missing** dates (NaN) as unknown — those plants are kept.

Add to `config.KZ.yaml`:

```yaml
electricity:
  powerplants_filter: (DateOut >= 2020 or DateOut != DateOut) and (DateIn <= 2020 or DateIn != DateIn)
```

---

## Step 4: Lock renewable capacities

Two settings work together: stop the optimiser from expanding, and pin solar/wind to IRENA **2020** totals.

**4a. Remove carriers from `extendable_carriers`**

Default includes solar, wind, and **OCGT**. For a fixed 2020 fleet, leave empty lists:

```yaml
electricity:
  extendable_carriers:
    Generator: []
    StorageUnit: []
    Store: []
    Link: []
```

**4b. Point IRENA stats at 2020**

PyPSA-Earth sums existing solar and wind from `powerplants.csv`, then fills any national shortfall to match IRENA using `estimate_renewable_capacities`:

```yaml
electricity:
  estimate_renewable_capacities:
    stats: irena
    year: 2020
```

- **`stats: irena`** — top up solar and wind to IRENASTAT national installed capacity.
- **`year: 2020`** — match the validation year.

---

## Step 5: Complete your config

Merge the generation settings with your Part 3 demand block (`load_options`, `enable.retrieve_cutout: false`, etc.):

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.generation.yaml"
```

Or [download the snippet](snippets/config.KZ.generation.yaml){: download="config.KZ.yaml"}.

Your file should contain **both** `load_options` (Part 3) and `electricity` (Part 4). Snakemake merges this file on top of `config.default.yaml` — you only need the keys you are overriding.

---

## Step 6: Re-run the workflow

Generation settings affect **`build_powerplants`** and **`add_electricity`** and everything downstream:

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Snakemake rebuilds `powerplants.csv`, re-runs `add_electricity`, and re-solves. Expect a similar runtime to [Part 3](3-demand-data.md) (~7–10 minutes with cached cutouts).

---

## Step 7: Verify installed capacity

Reopen the notebook from Part 2 and reload the solved network. Use the same **`statistics()`** call as in [Part 2](2-analyze-results.md) — with a locked fleet, **Installed** and **Optimal** capacity should match for generators:

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")

caps = n.statistics()["Installed Capacity"].dropna() / 1e3  # GW
print(caps.sort_values(ascending=False).to_string())
```

After Part 4 you should see something like (all values in **GW** — from the `/ 1e3` above):

```
Line         AC                    31.19
Generator    Load shedding         17.70
             Coal                  11.60
             Combined-Cycle Gas     2.47
StorageUnit  Reservoir & Dam        2.34
Generator    Oil                    1.14
             Solar                  0.91
             Onshore Wind           0.43
             Run of River           0.13
```

Ignore **Load shedding** and **AC** line capacity — they are not part of the generation fleet. **Solar (~0.9 GW)** and **onshore wind (~0.4 GW)** should now sit near KEGOC/IRENA 2020 levels, not the ~20 / ~6 GW from the Part 2 baseline.

**Comparison vs. KEGOC 2020 (GW):**

| Carrier | KEGOC 2020 | Baseline (Part 2) | After Part 4 |
|---|---|---|---|
| Coal | 13.41 | 11.72 | 11.60 |
| Gas | 6.01 | 2.43 | 2.47 |
| Hydro | 2.95 | 2.46 | 2.46 |
| Solar | 0.96 | 20.06 | 0.91 |
| Wind | 0.51 | 6.19 | 0.43 |

Solar and wind capacities are now aligned with statistics — that was the main gap. **Coal and gas** installed capacity are still below KEGOC; powerplantmatching under-reports gas in particular. The **Advanced** section below uses the custom fleet to replace powerplantmatching (or merge extra plants).

### Generation mix (still approximate)

Locking **capacities** does not fix **dispatch**. An optional check — same `statistics()` call, **Supply** column in **TWh**:

```python
supply = n.statistics()["Supply"].dropna() / 1e6  # TWh
print(supply.sort_values(ascending=False).to_string())
```

After Part 4 you might see (in GW):

```
Line         AC                    110.29
Generator    Coal                   90.85
             Load shedding           7.99
StorageUnit  Reservoir & Dam         2.83
Generator    Combined-Cycle Gas      2.40
             Onshore Wind            1.28
             Solar                   1.16
             Run of River            0.78
             Oil                     0.00
```

**Solar (~1.2 TWh)** and **wind (~1.3 TWh)** are now in the right ballpark vs KEGOC 2020 generation (~1.3 / ~1.1 TWh). **Coal (~91 TWh)** is too high and **gas (~2.4 TWh)** too low — the model still lacks gas capacity and over-relies on coal. **Hydro (~3.6 TWh)** remains under-generated. Non-zero **load shedding (~8 TWh)** means the solve still struggles to meet demand in some hours.

Improving the **generation** table in [Part 2](2-analyze-results.md) is a follow-up step (custom plants, hydro inflow, marginal costs, and so on).

---

## Advanced: Custom powerplants

[Step 7](#step-7-verify-installed-capacity) left **gas** short: **~2.5 GW** CCGT installed and **~2.4 TWh** annual generation vs **~6 GW** gas in KEGOC's 2020 statistics, and almost no **OCGT** peakers in powerplantmatching. For Kazakhstan, the global database misses many industrial gas turbines and smaller thermal units that show up in national statistics.

This tutorial ships a copy of [`custom_powerplants.csv`](https://github.com/pypsa-meets-earth/pypsa-kz-data/blob/main/data/custom_powerplants.csv) from the [**pypsa-kz-data**](https://github.com/pypsa-meets-earth/pypsa-kz-data) repository — **144 rows** in powerplantmatching format — as a **full Kazakhstan fleet** instead of patching gaps by hand.

### What is in the file

[Download the file](snippets/custom_powerplants.KZ.csv){: download="custom_powerplants.csv"} and save it as **`data/custom_powerplants.csv`**, replacing the existing placeholder. Snakemake always reads that path in **`build_powerplants`**. With the default `custom_powerplants: false` it is ignored until you set **`replace`** (use this list instead of powerplantmatching) or **`merge`** (append to powerplantmatching) below — we recommend **`replace`** for this full fleet.

Installed capacity in the custom list (MW):

| Fueltype | MW |
|---|---|
| Hard Coal | 13,171 |
| CCGT | 3,580 |
| Hydro | 2,751 |
| OCGT | 1,626 |
| Solar | 822 |
| Wind | 649 |

The main gain over powerplantmatching is **OCGT** (**~1.6 GW**) and fuller **CCGT** coverage. **`replace`** makes this list the sole plant source (no double-counting against overlapping powerplantmatching rows). **`merge`** is available if you only want to add missing units.

Example rows (industrial peakers absent from powerplantmatching):

```csv
Name,Fueltype,Technology,Set,Country,Capacity,DateIn,DateOut,lat,lon
ES AFP TNK Kazchrome,OCGT,OCGT,PP,KZ,135,1996,2146,50.34808,57.13376
ZHGTS 56 JSC CNPC-Aktobe,OCGT,OCGT,PP,KZ,118,1999,2149,48.385282,57.434115
Aktobe CHP,CCGT,CCGT,CHP,KZ,88,1962,2112,50.33618,57.14072
```

The same **`powerplants_filter`** from Step 3 still applies — plants with `DateIn > 2020` are dropped.

### Enable custom powerplants

Add to `config.KZ.yaml` (alongside the Part 4 `electricity` block):

```yaml
electricity:
  custom_powerplants: replace   # "false" | "merge" | "replace"
```

| Value | Behaviour |
|---|---|
| **`false`** | powerplantmatching only (Steps 1–7) — `data/custom_powerplants.csv` is ignored |
| **`merge`** | Append custom rows to powerplantmatching |
| **`replace`** | Use only your CSV — **recommended** for this list |

### Re-run and verify

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Snakemake rebuilds **`build_powerplants`** when `custom_powerplants` or the CSV changes. Re-check installed capacity (Step 7): **CCGT** and **OCGT** should move toward KEGOC 2020, and annual **Combined-Cycle Gas** / **Open-Cycle Gas** supply should rise relative to the coal-heavy mix from Step 7.

After **`replace`**, installed capacity (GW) tracks KEGOC 2020 much more closely:

```
Line         AC                    31.81
Generator    Load shedding         17.70
             Coal                  13.17
             Combined-Cycle Gas     3.58
StorageUnit  Reservoir & Dam        2.31
Generator    Open-Cycle Gas         1.60
             Solar                  0.80
             Onshore Wind           0.46
             Run of River           0.06
```

**CCGT (~3.6 GW)** and **OCGT (~1.6 GW)** together reach **~5.2 GW** — near KEGOC **~6 GW** gas. **Coal (~13.2 GW)** matches KEGOC **~13.4 GW**.

Annual **Supply** (TWh) is still far from KEGOC 2020:

```
Line         AC                    109.37
Generator    Coal                   92.35
             Load shedding           7.67
StorageUnit  Reservoir & Dam         3.93
Generator    Onshore Wind            1.34
             Solar                   0.99
             Combined-Cycle Gas      0.63
             Run of River            0.40
             Open-Cycle Gas          0.00
Load         -                       0.00
```

**Solar (~1.0 TWh)** and **wind (~1.3 TWh)** look reasonable vs KEGOC (~1.3 / ~1.1 TWh). **Coal (~92 TWh)** is far too high and **gas (~0.6 TWh from CCGT; OCGT idle)** far too low vs KEGOC (~75 / ~22 TWh). **Hydro (~4.3 TWh combined)** is still under-generated vs KEGOC (~9.5 TWh).

**Load shedding (~7.7 TWh)** is still high — roughly unchanged from Step 7 even after capacity improves. Several factors may contribute, alone or together:

- **Network topology** — a bus may be poorly connected (missing line or no path from local generation to load).
- **Transmission limits** — line capacities may block power from reaching demand in some hours or regions.
- **Spatial mismatch** — demand and generation may be on different **10-cluster** buses; national totals can look fine while some regions are short.
- **Dispatch economics** — cheap coal is used first; gas, hydro, and renewables may not fill remaining gaps.

Capacity alignment does not automatically fix the **generation mix** either. The optimiser still dispatches the cheapest available energy. With default settings that often means:

- **No CO₂ or fuel constraints** — nothing pushes the model away from cheap coal.
- **Marginal costs** — coal is typically cheaper than gas in the default cost tables, so CCGT/OCGT sit idle even when capacity is there.
- **Spatial resolution** — how load and plants are distributed across clusters may not match real geography.
- **Transmission** — line limits may be loose enough that northern coal serves most of the country, or tight enough to cause regional shortfalls.

---

## Recap

| Step | Config key | Role |
|---|---|---|
| 3 | `electricity.powerplants_filter` | Keep plants operating in **2020** |
| 4a | `electricity.extendable_carriers` | Empty `Generator` — no new solar/wind/gas build |
| 4b | `electricity.estimate_renewable_capacities` | IRENA **2020** solar/wind totals |
| Adv. | `electricity.custom_powerplants: replace` | Use **`data/custom_powerplants.csv`** as the full fleet |

Demand is calibrated (Part 3); the **2020 generation fleet** is now locked for capacity comparisons. The fleet matches KEGOC, yet the model still sheds **~8 TWh** of load — a network problem, not a generation one. In **[Part 5](5-network-topology.md)** we diagnose the electrically isolated sub-networks behind that shedding and fix them through simplification settings.
