<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 5: Adapt Fuel and Generation Costs

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 4](4-generation-data.md). Demand should be calibrated (`load_options.scale: 1.005`), the 2020 fleet locked with `custom_powerplants: replace`, and `config.KZ.yaml` should include those settings. The model will still shed **~7.7 TWh** of load on isolated buses at this point; that is expected, and it gets fixed in [Part 6](6-network-topology.md) and [Part 7](7-transmission-network.md).

## Introduction

By the end of [Part 4](4-generation-data.md) the model has the right **demand** and **fleet size**: installed coal and gas capacity sit near KEGOC 2020. Dispatch does **not**: the model runs almost entirely on coal. This tutorial checks the piece that is independent of network topology: **is the coal-over-gas dispatch actually explained by fuel and O&amp;M costs?** Rather than assume the answer, we replace the generic [technology-data](https://github.com/PyPSA/technology-data) defaults (European reference assumptions) with Kazakhstan-specific fuel and O&amp;M data and let the numbers speak.

After the custom fleet, annual **Supply** still looks roughly like:

| Carrier | Model (TWh) | KEGOC 2020 (TWh) |
|---|---|---|
| Coal | ~89 | ~75 |
| Gas (CCGT + OCGT) | ~0.6 | ~22 |
| Hydro | ~7 | ~9.5 |
| Solar / wind | ~1 / ~1.3 | ~1.3 / ~1.1 |

**~5.2 GW** of gas capacity sits in the network, but the optimiser barely uses it. With a fixed fleet, the underlying **linear program (LP)** is a **merit-order** problem: generators with lower **marginal cost** (EUR/MWh) are dispatched first. Default **technology-data** fuel prices make **coal cheaper than gas**, so coal fills most hours. That holds regardless of network topology, because cost ranking rather than connectivity decides dispatch order between connected buses, so it is worth checking against Kazakhstan-specific cost data now, before Parts 6–7 rebuild the grid.

In this tutorial we:

1. Inspect the default **marginal costs** on the generators in the solved network.
2. Compute country-specific **marginal costs** for Kazakhstan from sourced fuel and O&amp;M data.
3. Override **`costs.marginal_cost`** with computed marginal costs from step 2.
4. Re-solve and check whether country-specific costs actually explain the coal/gas split.

Everything lives under **`costs`** in the config. Changing it re-runs **`process_cost_data`** and everything that consumes `resources/KZ/costs_*_elec.csv`, including `add_electricity` and the solve.

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
\text{marginal_cost} = \text{VOM} + \frac{\text{fuel cost}}{\text{efficiency}}
\]


- **Fuel cost** (`fuel` in the cost tables) is in **EUR/MWh<sub>th</sub>**.
- **Efficiency** converts that thermal fuel cost to electrical MWh.
- **CCGT** and **OCGT** inherit the **`gas`** fuel price.
- **Coal** uses the **`coal`** fuel price.

See the [costs user guide](../../user-guide/costs.md) and [configuration reference](../../user-guide/configuration.md#costs).

---

## Step 1: Diagnose the default merit order

Reload the Part 4 solved network and print average variable costs by carrier:

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")

print(n.generators.groupby("carrier")["marginal_cost"].mean())
```

With European defaults you should see something like (**EUR/MWh**):

```
carrier
CCGT                46.81
OCGT                64.69
coal                30.11
load shedding   100000.01
offwind-ac           0.03
offwind-dc           0.02
onwind               0.02
solar                0.02
```

Generators also carry a **`capital_cost`** (annuitized investment), but with the Part 4 locked fleet it does not enter the objective; see the note in Step 4 for when it would matter.

**Coal (~30)** is cheaper than **CCGT (~47)** and much cheaper than **OCGT (~65)**. Renewables are near zero. **`load shedding`** sits at the default penalty price (**100,000 EUR/MWh**): it only dispatches when nothing else can meet demand, which is why it stays out of the normal merit order entirely. With no CO₂ price and no heat demand for CHPs, the optimiser runs coal first and leaves gas idle, exactly the Part 4 generation gap.

You can also inspect the processed cost file:

```python
import pandas as pd

costs = pd.read_csv("resources/KZ/costs_2030_elec.csv", index_col=0)
print(
    costs.loc[["coal", "CCGT", "OCGT"], ["fuel", "VOM", "efficiency", "marginal_cost"]]
)
```

You should see something like:

```
technology       fuel     VOM   efficiency   marginal_cost
coal            9.55    3.26        0.356           30.10
CCGT           24.57    4.44        0.580           46.80
OCGT           24.57    4.76        0.410           64.68
```

**Gas fuel cost (~24.6 EUR/MWh<sub>th</sub>) is about 2.6 times coal fuel cost (~9.6 EUR/MWh<sub>th</sub>)**: that gap is the root cause. CCGT and OCGT both pay the gas price, but OCGT's lower efficiency than CCGT's **(OCGT - 0.41 ; CCGT - 0.58)** stretches the same fuel cost over less electrical output, pushing its marginal cost well above CCGT's.

---

## Step 2: Choose the cost year

Defaults use **`costs.year: 2030`** (forward-looking technology-data). For a **2020 validation** study, pin the same year as your fleet and IRENA settings:

```yaml
costs:
  year: 2020
```

This switches which `costs_{year}.csv` is retrieved and processed. The fuel prices in the 2020 file may differ a little from the 2030 defaults, but coal usually stays cheaper than gas relative to each other. You still need Step 3 to bring in country-specific cost data.

---

## Step 3: Compute country-specific marginal costs for Kazakhstan

Generic technology-data defaults are a reasonable starting point, but country-specific data is better when you can find it. As an example, [Gasilov (2025), "Cost-optimal energy system development pathways for Kazakhstan," QazaqGreen](https://qazaqgreen.com/en/journal-qazaqgreen/industry-news/2994/) gives, for existing plants:

| Carrier | Fuel price | O&amp;M (OPEX) |
|---|---|---|
| Coal | 2.8 KZT/kWh | 6.0 KZT/kWh |
| CCGT | 24.5 KZT/kWh | 3.0 KZT/kWh |
| OCGT | 29.9 KZT/kWh | 4.0 KZT/kWh |

with 1 USD ≈ 500 KZT (the article's own exchange rate).

**Watch the units.** This fuel price already accounts for efficiency, unlike the fuel-cost term in the formula above (`VOM + fuel cost / η`), which is EUR/MWh<sub>th</sub> *before* efficiency. So rather than overriding `costs.fuel`, we compute `marginal_cost` ourselves and override that directly.

Since Fuel and O&amp;M are already on an electricity basis, no efficiency division is needed: just add them together to get the marginal cost:

| Carrier | Marginal cost | → USD/MWh | → EUR/MWh |
|---|---|---|---|
| Coal | 2.8 + 6.0 = 8.8 KZT/kWh | 17.6 | ≈ 16 |
| CCGT | 24.5 + 3.0 = 27.5 KZT/kWh | 55.0 | ≈ 51 |
| OCGT | 29.9 + 4.0 = 33.9 KZT/kWh | 67.8 | ≈ 62 |

!!! note "Unit conversion"
    Multiply KZT/kWh by **2** to get USD/MWh: 1 KZT/kWh = 1,000 KZT/MWh, and 1,000 ÷ 500 KZT/USD = 2 USD/MWh. The EUR column is an exemplary conversion (`output_currency` default is **EUR**); use your own current USD→EUR rate for your study.

**This does not flip the merit order.** Coal (~16 EUR/MWh) is still roughly **3×** cheaper than CCGT (~51). That gap is even wider than with technology-data defaults (~30 vs ~47 in Step 1), because the sourced Kazakhstan coal price is lower than the European default, while gas is not much higher. So Kazakhstan-specific fuel and O&amp;M data alone do not explain why KEGOC dispatches ~22 TWh/yr of gas.

---

## Step 4: Apply the country-specific costs with `costs.marginal_cost`

Because the sourced numbers are already electricity-basis EUR/MWh<sub>el</sub> values (not a thermal fuel price), write them under **`costs.marginal_cost`** rather than `costs.fuel`. This runs **after** the fuel formula and **replaces** the computed marginal cost for that technology. **`costs.capital_cost`** works the same way: a direct override of the annuitized investment cost, if you have a sourced CAPEX/annuity figure instead of raw `investment` / `lifetime` / `FOM`:

```yaml
costs:
  year: 2020
  marginal_cost:           # EUR/MWh_el, sourced Kazakhstan fuel + O&M data
    coal: 16.0
    CCGT: 51.0
    OCGT: 62.0
```

Use **`costs.marginal_cost`** whenever your source already reports cost per unit of electricity generated, as here. Use **`costs.fuel`** instead when your source gives a thermal-basis fuel price (EUR/MWh<sub>th</sub>, before efficiency); that way efficiency and VOM stay consistent with technology-data.

The same per-technology override pattern works for **`costs.VOM`**, **`costs.FOM`**, **`costs.efficiency`**, **`costs.investment`**, and **`costs.lifetime`**. That is useful if you find sourced values for those too.

!!! note "Why `capital_cost` would not change this solve"
    PyPSA only puts `capital_cost` into the objective for **extendable** capacity; fixed `p_nom` is excluded from that term entirely. Since Part 4 locked the fleet (`electricity.extendable_carriers` all empty), nothing here is extendable, so a `costs.capital_cost` override would not move this solution at all; it only matters once you allow new capacity, e.g. in a capacity-expansion study.

!!! warning "Do not mix blindly"
    If you set both `fuel` and `marginal_cost` for the same technology, the **`marginal_cost`** overwrite wins. Keep one approach per carrier.

---

## Step 5: Complete your config

Merge the `costs` changes above with your Parts 3–4 settings (`load_options`, `electricity`, and so on). The combined config looks like:

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.costs.yaml"
```

Or [download the snippet](snippets/config.KZ.costs.yaml){: download="config.KZ.yaml"}.

The new block is only **`costs`**. Everything else should already be in your working `config.KZ.yaml`.

---

## Step 6: Re-run the workflow

We now re-run the electricity workflow with the updated `config.KZ.yaml`. Because the `costs` options changed, Snakemake will re-execute **`process_cost_data`** and every downstream rule that reads the elec cost table (including **`add_electricity`** and **`solve_network`**). Cutouts and OSM stay cached:

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Expect a runtime similar to Part 4 (~7–10 minutes with cached data), not a full Part 1 rebuild.

While the full solve runs, do a quick sanity check that the override actually landed: `process_cost_data` writes the post-override table to `resources/KZ/costs_2020_elec.csv`, so it should already show your Step 4 values (16.0 / 51.0 / 62.0):

```python
costs = pd.read_csv("resources/KZ/costs_2020_elec.csv", index_col=0)
print(costs.loc[["coal", "CCGT", "OCGT"], ["marginal_cost"]])
```

If the numbers don't match, check for a typo in the carrier name or config indentation before waiting on the full solve.

---

## Step 7: Verify the generation mix

To check whether the country-specific costs changed the coal/gas dispatch, and to compare the resulting energy balance with KEGOC 2020, reload the solved network and inspect **marginal costs** and **Supply** (TWh):

```python
import pypsa

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")

print(n.generators.groupby("carrier")["marginal_cost"].mean().sort_values())

supply = n.statistics()["Supply"].dropna() / 1e6  # TWh
supply = supply.drop(["Line", "Load"], errors="ignore")  # dropping line and load values
print(supply.sort_values(ascending=False).to_string())
```

You should see something like:

```
carrier
solar                 0.02
offwind-dc            0.02
onwind                0.02
offwind-ac            0.03
coal                 16.01
CCGT                 51.01
OCGT                 62.01
load shedding    100000.01
Name: marginal_cost, dtype: float64
```

```
Generator    Coal                   89.89
             Load shedding           7.72
StorageUnit  Reservoir & Dam         7.15
Generator    Onshore Wind            1.34
             Solar                   0.99
             Combined-Cycle Gas      0.20
             Open-Cycle Gas          0.00
```

**What to look for**

| Check | Expectation |
|---|---|
| Merit order | `coal.marginal_cost` (~16 EUR/MWh) still **below** `CCGT`/`OCGT` (~51 / ~62 EUR/MWh); confirmed, country-specific data does **not** flip it |
| Gas supply | Stays low (**~0.2 TWh** above), even lower than Part 4's ~0.6 TWh; country-specific costs do not push gas up |
| Coal supply | Stays close to **~89 TWh** |
| Load shedding | Roughly **unchanged** from Part 4 (~7.7 TWh), a **network** issue (isolated buses), not a costs issue; fixed in [Part 6](6-network-topology.md) |
| Hydro / VRE | Largely unchanged (near-zero marginal cost; constrained by profiles) |

That is the honest result: sourced Kazakhstan fuel and O&amp;M data still make coal the cheaper option on the margin. This model does **not** reproduce KEGOC's real ~22 TWh/yr of gas generation through costs alone; if anything, gas supply drops slightly further. In practice much of Kazakhstan's gas fleet is combined-heat-and-power (CHP) serving district heat demand, or otherwise committed under contracts rather than pure economic merit order: a behaviour that marginal-cost override cannot capture. Closing this gap credibly would mean modelling that constraint directly (for example a must-run floor or heat-linked dispatch for gas CHP), not further fuel-price tuning.

A second confusion that is not cost-related sits in the network. By default, PyPSA-Earth lets the solver expand every transmission line at capital cost (`scenario.ll: ["copt"]`). Therefore, coal can reach any bus without a local gas plant ever needing to dispatch. This is further discussed in [Part 7](7-transmission-network.md#step-5-stop-lines-from-being-cost-optimized-past-their-rating) where this setting is turned off.

---

## Recap

| Step | Config key | Role |
|---|---|---|
| 2 | `costs.year` | Use **2020** technology-data for the validation year |
| 3 | *(sourced Kazakhstan fuel + O&amp;M data)* | Compute the country-specific marginal cost per technology (already EUR/MWh<sub>el</sub>) |
| 4 | `costs.marginal_cost.*` | Apply the computed values, replacing technology-data's generic figures |

Marginal costs now reflect sourced Kazakhstan fuel and O&amp;M data instead of generic technology-data defaults, but the coal/gas dispatch split barely moves, because country-specific costs still favor coal. The gap to KEGOC's observed gas generation is therefore **not primarily an economics problem** in this model. The model still sheds **~7.7 TWh** of load on electrically isolated buses; that is a separate **network** problem. **[Part 6](6-network-topology.md)** diagnoses and fixes those isolated sub-networks, and **[Part 7](7-transmission-network.md)** then rebuilds the transmission network itself with Kazakhstan-specific voltage levels and line ratings.
