<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 5: Fix Isolated Nodes and Load Shedding

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 4](4-generation-data.md). Your `config.KZ.yaml` should include the Part 3 demand settings and the Part 4 generation fleet, and you should have a solved network at `results/KZ/networks/elec_s_10_ec_lcopt_6h.nc`.

## Introduction

By the end of [Part 4](4-generation-data.md) the fleet matched KEGOC's 2020 capacities, yet the model still **shed about 8 TWh of load** — roughly 7–8% of Kazakhstan's annual demand. Load shedding means the optimiser could not serve demand at some buses in some hours, so it "dropped" that load at a very high penalty price. A validated fleet that still sheds load is a signal that the problem is **not generation** but the **network** the generation sits on.

In this tutorial we diagnose the cause — one or more **electrically isolated sub-networks** — and fix it by changing **how simplification handles islands**. This is the fastest lever: it re-runs in minutes and does not touch the OSM base network. A follow-up tutorial will go one level deeper and improve the base transmission network itself (voltage levels and line ratings) so fewer islands appear in the first place.

Everything in this part lives under **`cluster_options.simplify_network`** in the config and the **`simplify_network`** rule. It does not change the demand or the fleet.

---

## Where isolation comes from

PyPSA-Earth builds the transmission grid from **[OpenStreetMap (OSM)](https://www.openstreetmap.org/)** — volunteer-mapped substations and power lines, downloaded in [Part 1](1-baseline-model.md). The model topology is whatever OSM contains at download time, not an official KEGOC schematic.

Two facts about the workflow combine to create the problem:

```
build_demand_profiles  →  splits national demand across every substation bus
        ↓                  (weighted by population and GDP — no grid check)
simplify_network       →  merges/clusters buses, then handles leftover islands
        ↓
… cluster, prepare, solve …
```

1. **`build_demand_profiles`** takes the national annual total and distributes it across **all** substation buses using population and GDP weights. It does **not** check whether a bus is electrically connected to the rest of the grid.
2. If the OSM-derived network has a **gap** — a region whose lines never link to the national backbone — that region becomes its own **sub-network** (an electrical island). It still carries its share of demand, but the only generation available to serve it is whatever sits inside the island.

When an island's local generation cannot cover its local demand, the optimiser has no way to import power across the missing lines, so it **sheds** the unmet load. Nationally the capacity and demand totals look fine, but a some part of the country is starved.

!!! note "A PyPSA sub-network"
    After `n.determine_network_topology()`, every bus is labelled with a `sub_network`. Buses in the same `sub_network` are electrically connected; buses in different sub-networks cannot exchange power. A healthy country model has **one** dominant AC sub-network (the "backbone") carrying almost all load.

---

## Step 1: See the island on a map

Open your Part 2 notebook (`analyze_kz.ipynb`) and reload the solved network. First, colour buses by whether they belong to the main grid or to an island:

```python
import pypsa
import matplotlib.pyplot as plt

n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")
n.determine_network_topology()

# The backbone is the sub-network carrying the most load
load_by_sub = n.loads_t.p_set.mean().groupby(n.buses.sub_network).sum()
backbone = load_by_sub.idxmax()

is_isolated = n.buses.sub_network != backbone
bus_colors = is_isolated.map({True: "crimson", False: "seagreen"})

# Map extent: [lon_min, lon_max, lat_min, lat_max] — frame Kazakhstan
boundaries = [46, 88, 40, 56]

n.plot(
    bus_colors=bus_colors,
    bus_sizes=0.03,
    line_colors="lightgray",
    boundaries=boundaries,
    title="Green = main grid   |   Red = isolated",
)
plt.show()
```

Any **red** bus is electrically isolated from the green backbone — load there cannot import power from the rest of the grid. You may see **one red bus** rather than a whole western region; Step 3 explains why. After the fix in Step 4, re-plot — you should see no red buses.

![Kazakhstan network coloured by sub-network (before fix)](figures/kz_subnetworks.png)

*Green = main grid; red = isolated bus. Example from the Part 4 solved network with default simplification settings.*

---

## Step 2: Confirm load shedding on the island

The map shows *where* the problem is. Confirm that the **~8 TWh load shedding** from Part 4 sits on the **red** sub-network, not on the green backbone.

Load-shedding generators are named `<bus> load`:

```python
weights = n.snapshot_weightings.generators
shed = n.generators_t.p.filter(like="load").multiply(weights, axis=0).sum() / 1e6  # TWh
shed_by_sub = shed.groupby(
    n.generators.loc[shed.index, "bus"].map(n.buses.sub_network)
).sum()
print(shed_by_sub[shed_by_sub > 0].round(2))
```

Expected output — shedding on the **isolated** sub-network only (the backbone should not appear):

```
1    7.67
dtype: float64
```

The index is the **`sub_network`** id from Step 1 (the red bus). **~7.7 TWh** here matches Part 4's **~8 TWh** total load shedding: almost all unmet load sits on that one island, not on the green backbone.

!!! note "Why not measure island size to pick a threshold?"
    The sub-network shares you see in a **solved** network already went through simplification — including the default **`p_threshold_merge_isolated: 300`**, which collapses many small western pockets into **one** large isolated bus (~7% of national load). Measuring that share and trying to size `s_threshold_fetch_isolated` from it is misleading. The fix is to change simplification **first** (turn merge off, then fetch). Step 4 gives the values that work for Kazakhstan.

---

## Step 3: The three isolation thresholds

Simplification has three settings for isolated sub-networks, applied in this fixed order inside `simplify_network`:

| Order | Parameter | Unit | What it does |
|---|---|---|---|
| 1 | `p_threshold_drop_isolated` | **MW** (mean load) | **Deletes** islands whose total mean load is below the threshold. The load *disappears* from the model. |
| 2 | `p_threshold_merge_isolated` | **MW** (mean load) | Collapses small islands into **one isolated bus per country** — still disconnected from the backbone. |
| 3 | `s_threshold_fetch_isolated` | **share** of country load | Attaches islands whose share is below the threshold to the **nearest backbone bus**, creating a real electrical connection. |

The critical distinction:

- **`drop`** removes load (avoids shedding by throwing demand away).
- **`merge`** is the **default trap** (`p_threshold_merge_isolated: 300` in `config.default.yaml`): it collapses many small islands into **one isolated bus per country**. That bus still cannot import power, so it **still sheds** — and it can hold a **large** share of national load (e.g. ~7%), which makes later fixes harder to reason about.
- **`fetch`** is the **actual fix**: it wires stranded load onto the closest **connected** bus so the backbone's generation can serve it.

The relevant excerpt from `scripts/simplify_network.py`:

```python
if p_threshold_drop_isolated:
    n = drop_isolated_networks(n, threshold=p_threshold_drop_isolated)

if p_threshold_merge_isolated:
    n, merged_nodes_map = merge_isolated_networks(
        n,
        threshold=p_threshold_merge_isolated,
        aggregation_strategies=aggregation_strategies,
    )
    busmaps.append(merged_nodes_map)

if s_threshold_fetch_isolated:
    n, fetched_nodes_map = merge_into_network(
        n,
        threshold=s_threshold_fetch_isolated,
        aggregation_strategies=aggregation_strategies,
    )
    busmaps.append(fetched_nodes_map)
```

**Order matters in the config, not only in the script:** turn **`merge` off first**, then enable **`fetch`**. With merge disabled, western regions stay as **separate** small islands; `fetch` can attach each one whose share is below the threshold.

!!! note "This is a modelling simplification, not a real line"
    `fetch` does **not** build a physical line. It re-assigns the stranded load (and any generation) to the geographically nearest connected bus so the linear program can balance it. The load is served, but its exact electrical location is approximated. The physically faithful fix — adding the missing transmission — comes in a follow-up tutorial.

---

## Step 4: Settings for Kazakhstan

**1. `p_threshold_merge_isolated: false` — do this first.**

The PyPSA-Earth default merges islands below **300 MW** mean load onto one stranded bus. For KZ that stacks much of the western pocket into a single ~7% island. Disabling merge leaves smaller fragments that `fetch` can reconnect individually:

```yaml
p_threshold_merge_isolated: false
```

**2. `s_threshold_fetch_isolated: 0.05` — attach small islands to the backbone.**

`fetch` reconnects any sub-network whose load share is **below** the threshold (here **5%** of national load) to the nearest backbone bus. With merge off, the western fragments are each below 5%, so **`0.05` is enough** for Kazakhstan:

```yaml
s_threshold_fetch_isolated: 0.05
```

If you left the default merge on, one mega-island could exceed 5% and `fetch` would skip it — another reason to disable merge first.

**3. `p_threshold_drop_isolated: 10` — drop only tiny artefacts.**

Keep this low so real regional load is preserved for `fetch`, not deleted:

```yaml
p_threshold_drop_isolated: 10   # [MW] drop only negligible OSM islands
```

The default is **20 MW**; **10 MW** is slightly stricter on artefacts only. You can keep **20** if you prefer the default — it does not change the merge/fetch logic above.

## Step 5: Add the settings to `config.KZ.yaml`

Add a `cluster_options` block:

```yaml
cluster_options:
  simplify_network:
    p_threshold_merge_isolated: false # FIRST: do not stack islands on one stranded bus
    s_threshold_fetch_isolated: 0.05  # attach islands < 5% of national load to the nearest backbone bus
    p_threshold_drop_isolated: 10     # [MW] drop only tiny artefact islands
```

You can [download the file](snippets/config.KZ.topology.yaml){: download="config.KZ.yaml"} and merge it with your existing `config.KZ.yaml`, or add the `cluster_options` block by hand.

---

## Step 6: Re-run the workflow

Run the same target as before:

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

Because these settings feed **`simplify_network`**, Snakemake rebuilds from there — `simplify_network` → `cluster_network` → `add_extra_components` → `prepare_network` → `solve_network`. The OSM base network, cutouts, demand profiles, and powerplant list all stay cached from earlier parts.

**Expected runtime:** a few minutes — much faster than Part 1, since only the clustering-and-solve tail re-runs.

---

## Step 7: Verify the fix

Reload the solved network in your notebook and repeat the checks from Steps 1–2.

**Load shedding** — should collapse toward zero:

```python
n = pypsa.Network("results/KZ/networks/elec_s_10_ec_lcopt_6h.nc")
weights = n.snapshot_weightings.generators
shed_TWh = (
    n.generators_t.p.filter(like="load").multiply(weights, axis=0).sum().sum() / 1e6
)
print(f"Load shedding: {shed_TWh:.2f} TWh")
```

```
Load shedding: 0.0 TWh
```

Re-plotting with the Step 1 code should now show a single-colour map — no red pocket — because the previously isolated buses were fetched onto the backbone and can now import power.

!!! note "Total demand may tick up slightly"
    Reconnecting stranded load (instead of dropping it) restores a little of the demand that simplification used to discard — the gap discussed in [Part 3](3-demand-data.md#step-3-calibrate-annual-demand-with-scale). If you calibrated `scale` against the pre-fix baseline, re-check the annual total and adjust `scale` if needed.

---

## Recap

| Step | Config key | Value | Role |
|---|---|---|---|
| 4 | `cluster_options.simplify_network.p_threshold_merge_isolated` | `false` | **First:** do not collapse islands onto one stranded bus (default 300 MW) |
| 4 | `cluster_options.simplify_network.s_threshold_fetch_isolated` | `0.05` | Attach islands below 5% of national load to the nearest backbone bus |
| 4 | `cluster_options.simplify_network.p_threshold_drop_isolated` | `10` | Drop only tiny (<10 MW) artefact islands |

Load shedding caused by electrical islands is resolved: the stranded demand is now wired onto the main grid and served by national generation. This is a **simplification-level** fix — fast and effective, but it approximates *where* that load connects.

A follow-up tutorial (Part 6) will address the root cause in the base network: lowering the OSM voltage threshold to capture Kazakhstan's sub-transmission grid, adding the 110 kV level, and choosing line ratings that reflect the real, older KZ conductors — so the western region is genuinely connected rather than fetched.
