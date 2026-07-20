<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 7: Adapt Fuel and Generation Costs

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 6](6-transmission-network.md). Demand should be near KEGOC (`load_options.scale: 0.994`), the 2020 fleet locked with `custom_powerplants: replace`, load shedding at ~0 TWh after Parts 5–6, and `config.KZ.yaml` should include those settings.

## Introduction

By the end of [Part 6](6-transmission-network.md) the model has the right **demand**, **fleet size**, and **network connectivity**. Installed coal and gas capacity sit near KEGOC 2020. Dispatch does **not**.

After the custom fleet (and without RoR rows trapping hydro inflow), annual **Supply** still looks roughly like:

| Carrier | Model (TWh) | KEGOC 2020 (TWh) |
|---|---|---|
| Coal | ~89 | ~75 |
| Gas (CCGT + OCGT) | ~0.6 | ~22 |
| Hydro | ~7 | ~9.5 |
| Solar / wind | ~1 / ~1.3 | ~1.3 / ~1.1 |

Capacity is available — **~5.2 GW** gas sits in the network — but the optimiser barely uses it. With a fixed fleet, the LP is a **merit-order** problem: generators with lower **marginal cost** (EUR/MWh) are dispatched first. Default **technology-data** fuel prices make **coal cheaper than gas**, so coal fills most hours.

In this tutorial we:

1. Inspect the cost database and the marginal costs attached to generators.
2. Override **fuel prices** for Kazakhstan (preferred).
3. Optionally set **marginal costs** directly as a shortcut.
4. Re-solve and check whether gas displaces some coal.

Everything lives under **`costs`** in the config. Changing it re-runs **`process_cost_data`** and everything that consumes `resources/KZ/costs_*_elec.csv` — including `add_electricity` and the solve.

---

## Where costs enter the workflow

```
retrieve_cost_data   →  resources/KZ/costs_{year}.csv   (technology-data)
        ↓
process_cost_data    →  resources/KZ/costs_{year}_elec.csv
                        capital_cost, marginal_cost, co2_emissions
                        ← costs.fuel / costs.VOM / costs.marginal_cost overrides
        ↓
add_electricity      →  attaches generators with those marginal costs
        ↓
… simplify, cluster, prepare, solve …
```

| Step | What happens |
|---|---|
| **`retrieve_cost_data`** | Downloads [PyPSA/technology-data](https://github.com/PyPSA/technology-data) for `costs.year` (default **2030**) |
| **`process_cost_data`** | Builds `capital_cost` and `marginal_cost`; applies config overrides |
| **`add_electricity`** | Copies those values onto generators / storage units |

**Formula used for thermal plants:**

\[
\text{marginal\_cost} = \text{VOM} + \frac{\text{fuel}}{\text{efficiency}}
\]

Fuel is in **EUR/MWh<sub>th</sub>**; efficiency converts it to electrical MWh. CCGT and OCGT inherit the **`gas`** fuel price; coal uses the **`coal`** fuel price.

See the [costs user guide](../../user-guide/costs.md) and [configuration reference](../../user-guide/configuration.md#costs).

---

## Step 1: Diagnose the default merit order

Reload the Part 6 solved network and print average variable costs by carrier:

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")

print(n.generators.groupby("carrier")[["marginal_cost", "capital_cost"]].mean())
```

With European defaults you should see something like (**EUR/MWh** for `marginal_cost`):

```
               marginal_cost   capital_cost
carrier
CCGT               35.01      104788
OCGT               64.69       47719
coal               30.11      337208
onwind              0.02      ...
solar               0.02      ...
```

**Coal (~30)** is cheaper than **CCGT (~35)** and much cheaper than **OCGT (~65)**. Renewables are near zero. With no CO₂ price and no heat demand for CHPs, the optimiser runs coal first and leaves gas idle — exactly the Part 4/6 generation gap.

You can also inspect the processed cost file:

```python
import pandas as pd

costs = pd.read_csv("resources/KZ/costs_2030_elec.csv", index_col=0)
print(
    costs.loc[
        ["coal", "gas", "CCGT", "OCGT"], ["fuel", "VOM", "efficiency", "marginal_cost"]
    ]
)
```

Expect **gas fuel ≫ coal fuel** in the database — that is the root cause.

---

## Step 2: Choose the cost year

Defaults use **`costs.year: 2030`** (forward-looking technology-data). For a **2020 validation** study, pin the same year as your fleet and IRENA settings:

```yaml
costs:
  year: 2020
```

This switches which `costs_{year}.csv` is retrieved and processed. Absolute fuel numbers may shift slightly by vintage; the **relative** coal-vs-gas gap usually remains. You still need Step 3 to change fuel prices.

---

## Step 3: Override fuel prices (preferred)

Config keys under **`costs.fuel`** overwrite the **fuel** column **before** marginal costs are computed. CCGT and OCGT automatically pick up the new **`gas`** price.

Add a starting point that makes gas **cheaper on the margin** than coal (illustrative — tune to your region and data):

```yaml
costs:
  year: 2020
  fuel:                    # EUR/MWh_th
    coal: 15.0             # raise from technology-data (~8–10)
    gas: 12.0              # lower from technology-data (~20–25)
```

Rough check with typical efficiencies (coal ~0.46, CCGT ~0.50, VOM ~6 / ~4):

| Carrier | Approximate \( \mathrm{VOM} + \mathrm{fuel}/\eta \) |
|---|---|
| Coal | \(6 + 15/0.46 \approx 39\) EUR/MWh<sub>el</sub> |
| CCGT | \(4 + 12/0.50 = 28\) EUR/MWh<sub>el</sub> |

After this change, **CCGT should dispatch before coal** whenever transmission and capacity allow. That is what we want for moving toward KEGOC’s ~22 TWh gas.

!!! tip "How much to change"
    Start with a clear merit-order flip (as above), re-solve, then nudge **`coal`** / **`gas`** until annual supply is closer to KEGOC. Perfect match is not required in one shot — Part 8 (CO₂) can close remaining policy-driven gaps.

---

## Step 4: Optional shortcut — override `marginal_cost` directly

If you already know target **EUR/MWh<sub>el</sub>** values, write them under **`costs.marginal_cost`**. This runs **after** the fuel formula and **replaces** the computed marginal cost for that technology:

```yaml
costs:
  year: 2020
  marginal_cost:           # EUR/MWh_el — replaces VOM + fuel/η
    coal: 40.0
    CCGT: 28.0
    OCGT: 45.0
```

Use this for quick what-if tests. For a documented regional study, prefer **`costs.fuel`** so efficiencies and VOM stay consistent with technology-data.

!!! warning "Do not mix blindly"
    If you set both `fuel` and `marginal_cost` for the same technology, the **`marginal_cost`** overwrite wins. Keep one approach per carrier.

---

## Step 5: Complete your config

Merge with Parts 3–6 (`load_options`, `electricity`, `cluster_options`, `clean_osm_data_options`, `lines`, …):

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.costs.yaml"
```

Or [download the snippet](snippets/config.KZ.costs.yaml){: download="config.KZ.yaml"}.

The new block is only **`costs`**. Everything else should already be in your working `config.KZ.yaml`.

---

## Step 6: Re-run the workflow

Cost changes invalidate **`process_cost_data`** and every downstream rule that reads the elec cost table (including **`add_electricity`** and **`solve_network`**). Cutouts and OSM stay cached:

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Expect a runtime similar to Part 4 (~7–10 minutes with cached data), not a full Part 1 rebuild.

Confirm the override in the log and processed CSV:

```text
logs/KZ/build_cost_data_2020_elec.log
```

```python
costs = pd.read_csv("resources/KZ/costs_2020_elec.csv", index_col=0)
print(costs.loc[["coal", "gas", "CCGT", "OCGT"], ["fuel", "marginal_cost"]])
```

---

## Step 7: Verify the generation mix

Reload the solved network and compare **Supply** (TWh) and **marginal costs**:

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")

print(n.generators.groupby("carrier")["marginal_cost"].mean().sort_values())

supply = n.statistics()["Supply"].dropna() / 1e6  # TWh
print(supply.sort_values(ascending=False).to_string())
```

**What to look for**

| Check | Expectation |
|---|---|
| Merit order | `CCGT.marginal_cost` **&lt;** `coal.marginal_cost` (if you followed Step 3) |
| Gas supply | **Up** from ~0.6 TWh toward KEGOC **~22 TWh** |
| Coal supply | **Down** from ~89 TWh toward KEGOC **~75 TWh** |
| Load shedding | Still **~0 TWh** (costs do not re-isolate buses) |
| Hydro / VRE | Largely unchanged (near-zero MC; constrained by profiles) |

If gas barely moves, open the merit-order printout again — another carrier may still undercut CCGT, or **transmission** may limit western gas from serving northern load. If gas overshoots and coal collapses too far, raise `costs.fuel.gas` slightly or lower `costs.fuel.coal`.

Exact TWh numbers depend on your fuel values; treat Step 3 as a **calibration loop**, not a single magic number.

---

## Advanced: Edit the cost CSV

Config overrides cover the common cases (`fuel`, `VOM`, `efficiency`, `investment`, `FOM`, `lifetime`, `marginal_cost`, `capital_cost`). For deeper edits:

1. Set **`enable.retrieve_cost_data: false`** and edit **`data/costs.csv`**, **or**
2. Edit the retrieved file under **`resources/KZ/costs_{year}.csv`** and keep retrieve on only when you intend to refresh from GitHub.

Prefer config overrides for the tutorial so your changes stay visible in `config.KZ.yaml` and Git diffs.

---

## Recap

| Step | Config key | Role |
|---|---|---|
| 2 | `costs.year` | Use **2020** technology-data for the validation year |
| 3 | `costs.fuel.coal` / `costs.fuel.gas` | Preferred — sets thermal fuel prices; MC recomputed via VOM + fuel/η |
| 4 | `costs.marginal_cost.*` | Optional shortcut — replace EUR/MWh<sub>el</sub> directly |
| Adv. | `data/costs.csv` / retrieve flag | Full database edits when config is not enough |

Fleet and network were already calibrated; **dispatch economics** now push gas into the merit order. Remaining gaps (policy-driven coal shares, exact CO₂ totals) are a natural next step with an explicit **CO₂ limit or emission price** in a follow-up tutorial.
