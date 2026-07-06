<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 6: Improve the Transmission Network

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 5](5-network-topology.md). Your `config.KZ.yaml` should include the Part 3 demand settings (`load_options.scale: 1.005`), the Part 4 generation fleet, and the Part 5 `cluster_options.simplify_network` settings.

## Introduction

[Part 5](5-network-topology.md) eliminated load shedding by **fetching** electrically isolated buses onto the main grid. That works, but it is a simplification shortcut — load is electrically re-assigned to the nearest connected bus, not carried over real missing lines.

In this tutorial we improve the **base transmission network** that PyPSA-Earth builds from OpenStreetMap (OSM):

1. **Kazakhstan's voltage hierarchy** — post-Soviet levels instead of European defaults.
2. **A lower OSM voltage floor (35 kV)** — keep regional sub-transmission links that the default 51 kV cutoff drops.
3. **KZ line types** — conductor ratings that reflect older, single-circuit lines at 110 kV and 220 kV.

These settings live upstream in the workflow (`clean_osm_data`, `base_network`, `electricity`, `lines`). Changing them triggers a **full rebuild** of the grid — slower than Part 5, but the grid topology and line capacities are closer to KZ reality before clustering and solve.

Keep the Part 5 **`cluster_options`** block. Even with a better OSM grid, some small islands may remain; merge-off and fetch still help.

---

## Where the grid enters the workflow

The transmission network is assembled long before `solve_network`:

```
download_osm_data  →  raw OSM extract
        ↓
clean_osm_data     →  filter by voltage, country, tags  ← threshold_voltage
        ↓
build_osm_network  →  buses, lines, transformers
        ↓
base_network       →  PyPSA components, line types, s_nom  ← min_voltage_rebase_voltage, ac_types
        ↓
add_electricity    →  attach demand and generators
        ↓
… simplify, cluster, prepare, solve …
```

| Config block | Snakemake rule | What it controls |
|---|---|---|
| **`clean_osm_data_options`** | `clean_osm_data` | Which OSM lines and substations survive the first filter |
| **`base_network`** | `base_network` | Minimum voltage kept when building `base.nc` |
| **`electricity`** | `simplify_network`, `cluster_network`, … | Nominal voltage levels the model simplifies toward |
| **`lines`** | `base_network` | Conductor type per voltage → thermal **`s_nom`** |

See the [configuration reference](../../user-guide/configuration.md) for every key.

---

## Step 1: See why KZ uses a different voltage scale

Kazakhstan follows the **post-Soviet/CIS** voltage scale, not the European one PyPSA-Earth uses by default:

| Role | Kazakhstan (typical) | PyPSA-Earth default |
|---|---|---|
| Sub-transmission / regional mesh | **35 kV**, **110 kV** | *(35 kV dropped; 132 kV)* |
| Transmission | **220 kV** | 220 kV, 300 kV, 380 kV |
| Extra-high-voltage backbone | **500 kV** | 500 kV, 750 kV |

**35 kV** is the critical layer for connectivity: in many regions it is the voltage that ties smaller substations into the wider grid. **110 kV** carries regional load between cities; **220 kV** and **500 kV** form the long-distance backbone.

PyPSA-Earth defaults assume a **European** stack (`electricity.voltages: [132, 220, 300, 380, 500, 750]`) and discard OSM assets below **51 kV** (`threshold_voltage: 51000`). For KZ that removes much of the sub-transmission mesh — one reason western pockets appeared as electrical islands in Parts 4–5.

**`electricity.voltages`** lists the nominal levels used while building **`base.nc`** from OSM (line-type matching on 35 / 110 / 220 / 500 kV assets). During **`simplify_network`**, every bus and line is mapped onto **`base_voltage`** (500 kV) and transformers are removed. The **solved** 10-cluster network therefore has **one nominal bus voltage** — not separate 110 / 220 / 500 kV buses. The low-voltage OSM layers still matter: they add regional lines and substation paths in `base.nc` **before** that collapse, which can change connectivity and line capacities even though you no longer see 35 kV or 110 kV in the final file.

---

## Step 2: Lower the OSM voltage floor to 35 kV

The only active voltage floor in the current workflow is **`clean_osm_data_options.threshold_voltage`** — applied in **`clean_osm_data`** to lines, cables, and substations:

| Key | Block | Rule |
|---|---|---|
| `threshold_voltage` | `clean_osm_data_options` | Filters raw OSM during `clean_osm_data` |

```yaml
clean_osm_data_options:
  threshold_voltage: 35000   # [V] keep 35 kV+ lines and substations from OSM
```

---

## Step 3: Set KZ nominal voltage levels

```yaml
electricity:
  base_voltage: 500
  voltages: [110, 220, 500]
```

- **`base_voltage`** — voltage level that **`simplify_network`** maps the entire grid onto (500 kV for KZ).
- **`voltages`** — nominal levels used when assigning **`lines.ac_types`** to OSM lines in **`base_network`** (each line matched to the closest kV in this list).

Replace the European default list for KZ — do not mix in 132 kV or 380 kV unless you have a reason to keep them.

---

## Step 4: Set KZ line types and ratings

Line thermal limits (**`s_nom`**) come from **`lines.ac_types`**: each voltage maps to a named conductor in PyPSA's line register. PyPSA-Earth picks the closest voltage match and computes capacity from the conductor's rated current and **`v_nom`**.

Defaults assume **modern European** bundles (e.g. two bundles at 220 kV). Kazakhstan's older corridors often run **single circuits** with lighter conductors. We encode that as a **modelling assumption** — lower **`s_nom`** at 220 kV and an explicit **110 kV** type (missing from the default table):

```yaml
lines:
  ac_types:
    110.: "149-AL1/24-ST1A 110.0"
    220.: "305-AL1/39-ST1A 110.0"
    500.: "490-AL1/64-ST1A 380.0"
```

| kV | KZ type | Effect |
|---|---|---|
| **110** | `149-AL1/24-ST1A 110.0` | Enables 110 kV line ratings |
| **220** | `305-AL1/39-ST1A 110.0` | Lighter conductor vs default **2-bundle** → lower **`s_nom`** |
| **500** | `490-AL1/64-ST1A 380.0` | Single-circuit EHV vs default **4-bundle** |

Tighter line limits can bind regional flows in the LP; looser defaults would under-state congestion on older 220 kV corridors.

---

## Step 5: Add the settings to `config.KZ.yaml`

Merge the blocks below with your existing Part 3–5 settings (`load_options`, `electricity` fleet keys, `cluster_options`, …):

```yaml
clean_osm_data_options:
  threshold_voltage: 35000

electricity:
  base_voltage: 500
  voltages: [110, 220, 500]

lines:
  ac_types:
    110.: "149-AL1/24-ST1A 110.0"
    220.: "305-AL1/39-ST1A 110.0"
    500.: "490-AL1/64-ST1A 380.0"
```

You can [download the file](snippets/config.KZ.transmission.yaml){: download="config.KZ.yaml"} and merge it with your existing `config.KZ.yaml`, or add the blocks by hand.

---

## Step 6: Re-run the workflow

Because these keys feed **`clean_osm_data`** and **`base_network`**, Snakemake rebuilds the grid from OSM onward — through demand, generation, simplify, cluster, and solve. Expect a **longer run** than Part 5 (similar to Part 1 in scope, but OSM and cutouts stay cached):

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

!!! tip "Use cached data"
    Keep from earlier parts:

    ```yaml
    enable:
      retrieve_databundle: false
      retrieve_cutout: false
    ```

The solved network still overwrites `results/KZ/networks/elec_s_10_ec_lcopt_6h.nc`.

---

## Step 7: Verify the solve

Reload the solved network in `analyze_kz.ipynb`. The OSM rebuild changes line capacities and topology — confirm [Part 5](5-network-topology.md#step-7-verify-the-fix) still holds and demand is near the KEGOC target:

```python
weights = n.snapshot_weightings.generators
total_TWh = n.loads_t.p_set.multiply(weights, axis=0).sum().sum() / 1e6
shed_TWh = (
    n.generators_t.p.filter(like="load").multiply(weights, axis=0).sum().sum() / 1e6
)

print(f"Total annual demand: {total_TWh:.2f} TWh")
print(f"Load shedding: {shed_TWh:.2f} TWh")
```

```
Total annual demand: 107.87 TWh
Load shedding: 0.00 TWh
```

- **Load shedding** at zero — no new islands after the grid rebuild.
- **Total demand 107.87 TWh** with **`scale: 1.005`** — slightly above KEGOC **107.34 TWh**. A one-time tweak in Step 8 brings it back.

`n.statistics()` should show the same story: **Load shedding** supply at **0 TWh**; generator **Installed Capacity** unchanged from Part 4. Dispatch stays coal-heavy (~96 TWh) — expected until marginal costs or CO₂ limits are tuned.

---

## Step 8: Final calibration of `scale`

The Part 3 **`scale: 1.005`** was set before the Part 5–6 grid changes. Rescale once from your Step 7 total:

```python
target_TWh = 107.34  # KEGOC 2020
measured_TWh = 107.87  # Step 7 (with scale 1.005)
new_scale = 1.005 * target_TWh / measured_TWh
print(f"scale: {new_scale:.4f}")  # → 1.0000
```

```yaml
load_options:
  scale: 1.0
```

`1.005 × 107.34 / 107.87` rounds to **1.0**.

In [Part 3](3-demand-data.md#step-3-calibrate-annual-demand-with-scale), **`scale: 1.005`** bundled two gaps in one multiplier: GEGIS vs KEGOC statistics, **and** demand lost when simplification dropped load on isolated buses (~1.2 TWh in the Part 2 grid). Parts 5–6 fixed the second problem — **fetch** reconnects islands, **35 kV** keeps regional lines — so that load stays in the model. With the grid settled, the extra **0.5%** scale is no longer needed: **`scale: 1.0`** lands on the KEGOC total without a topology fudge.

Re-run the workflow if you change `scale`. Load shedding should stay at **~0 TWh** — scaling demand does not re-isolate buses.

---

## Recap

| Step | Config key | Value | Role |
|---|---|---|---|
| 2 | `clean_osm_data_options.threshold_voltage` | `35000` | Keep 35 kV+ OSM assets |
| 3 | `electricity.base_voltage` | `500` | Single voltage after `simplify_network` |
| 3 | `electricity.voltages` | `[110, 220, 500]` | Line-type matching on multi-voltage `base.nc` |
| 4 | `lines.ac_types` | KZ conductor map | 110 kV type + derated 220 kV **`s_nom`** in build |
| 5 | *(keep Part 5)* | `cluster_options.simplify_network` | Merge off, fetch on — catch remaining islands |
| 8 | `load_options.scale` | `1.0` (from `1.005` × 107.34/107.87) | **Final** match to **107.34 TWh** after grid is settled |

The base grid now reflects Kazakhstan's post-Soviet voltage scale and more conservative line ratings. Part 5's fetch remains a useful safety net for OSM gaps that 35 kV alone cannot close.

Demand ends at **`scale: 1.0`** — the GEGIS profile at its native total, on a grid that no longer drops regional load. That is a cleaner baseline for the next tutorials (CO₂ limits, regional costs).
