<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 7: Improve the Transmission Network

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 6](6-network-topology.md). Your `config.KZ.yaml` should include the Part 3 demand settings (`load_options.scale: 1.005`), the Part 4 generation fleet, the Part 5 `costs` overrides, and the Part 6 `cluster_options.simplify_network` settings.

## Introduction

[Part 6](6-network-topology.md) eliminated load shedding by reconnecting electrically isolated buses to the main grid (**merge off**, **drop off**, **fetch on**). That fix works, but it is still a simplification shortcut: load is reassigned to the nearest connected bus, not carried over explicitly modelled transmission lines.

In this tutorial we improve the **base transmission network** that PyPSA-Earth builds from OpenStreetMap (OSM):

1. **Kazakhstan's voltage hierarchy** — post-Soviet levels instead of European defaults.
2. **A lower OSM voltage floor (35 kV)** — keep regional sub-transmission links that the default 51 kV cutoff drops.
3. **KZ line types** — conductor ratings that reflect older, single-circuit lines at 110 kV and 220 kV.

These settings live upstream in the workflow (`clean_osm_data`, `base_network`, `electricity`, `lines`). Changing them triggers a **full rebuild** of the grid — slower than Part 6, but it gives a base topology and line capacities that are closer to Kazakhstan's real network before clustering and solving.

Keep the Part 6 **`cluster_options`** block and the Part 5 **`costs`** block. Even with a better OSM grid, some small islands may remain; `merge` off, `fetch` on, and `drop` off still help preserve both served load and island generation. The `costs.marginal_cost` overrides from Part 5 are generator attributes, independent of topology, so they carry through the rebuild unchanged.

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

PyPSA-Earth defaults assume a **European** stack (`electricity.voltages: [132, 220, 300, 380, 500, 750]`) and discard OSM assets below **51 kV** (`threshold_voltage: 51000`). For KZ that removes much of the sub-transmission mesh — one reason western pockets appeared as electrical islands in Parts 4–6.

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

In practice, this means `electricity.voltages` controls which voltage layers are represented in **`base.nc`** (for KZ: 110 / 220 / 500 kV, with 35 kV retained by the OSM filter in Step 2), while `base_voltage` controls the later simplification collapse. During **`simplify_network`**, buses and lines are mapped to a single nominal level (500 kV) and transformers are removed. The solved clustered network therefore shows one nominal bus voltage, but the lower-voltage layers still matter upstream because they affect connectivity and line capacities before that collapse.

For KZ, use `[110, 220, 500]` as the default working set. Add other levels (for example 132 or 380 kV) only if your OSM data or study scope clearly requires them.

---

## Step 4: Set KZ line types and ratings

Line thermal limits (**`s_nom`**) come from **`lines.ac_types`**: each voltage maps to a named conductor in PyPSA's line register. PyPSA-Earth picks the closest voltage match and computes capacity from the conductor's rated current and **`v_nom`**.

Defaults assume **modern European** bundles (e.g. two bundles at 220 kV). Kazakhstan's older corridors often run **single circuits** with lighter conductors. We encode that as a **modelling assumption** — lower **`s_nom`** at 220 kV and an explicit **110 kV** type (missing from the default table):

```yaml
lines:
  ac_types:
    110.: "149-AL1/24-ST1A 110.0"
    220.: "243-AL1/39-ST1A 110.0"
    500.: "490-AL1/64-ST1A 380.0"
```

| kV | KZ type | Effect |
|---|---|---|
| **110** | `149-AL1/24-ST1A 110.0` | Enables 110 kV line ratings — the lightest conductor PyPSA-Earth's pinned PyPSA registers at this voltage |
| **220** | `243-AL1/39-ST1A 110.0` | Lighter conductor vs default **2-bundle** → lower **`s_nom`** |
| **500** | `490-AL1/64-ST1A 380.0` | Single-circuit EHV vs default **4-bundle** |

!!! note
    `ac_types` values must exist in your installed PyPSA's line register (`pypsa<=0.30.3` here — see the [version-pinned line types table](https://docs.pypsa.org/v0.30.3/user-guide/components.html#line-types), not the `latest` docs) — an unmatched name silently produces `NaN` **`s_nom`** instead of raising an error.

Tighter line limits can bind regional flows in the LP; looser defaults would under-state congestion on older 220 kV corridors.

**`lines.s_max_pu`** (default `0.7`) derates every line's usable capacity as a security margin. We lower it to `0.5` as a further conservative assumption for KZ's less-meshed grid:

```yaml
lines:
  s_max_pu: 0.5
```

---

## Step 5: Add the settings to `config.KZ.yaml`

Merge the blocks below with your existing Part 3–6 settings (`load_options`, `electricity` fleet keys, `costs`, `cluster_options`, …):

```yaml
clean_osm_data_options:
  threshold_voltage: 35000

electricity:
  base_voltage: 500
  voltages: [110, 220, 500]

lines:
  ac_types:
    110.: "149-AL1/24-ST1A 110.0"
    220.: "243-AL1/39-ST1A 110.0"
    500.: "490-AL1/64-ST1A 380.0"
  s_max_pu: 0.5
```

You can [download the file](snippets/config.KZ.transmission.yaml){: download="config.KZ.yaml"} and merge it with your existing `config.KZ.yaml`, or add the blocks by hand.

---

## Step 6: Re-run the workflow

Because these keys feed **`clean_osm_data`** and **`base_network`**, Snakemake rebuilds the grid from OSM onward — through demand, generation, simplify, cluster, and solve. Expect a **longer run** than Part 6 (similar to Part 1 in scope, but OSM and cutouts stay cached):

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

The solved network still overwrites `results/KZ/networks/elec_s_10_ec_lcopt_6h.nc`.

---

## Step 7: Verify the solve

Reload the solved network in `analyze_kz.ipynb`. The OSM rebuild changes line capacities and topology — confirm [Part 6](6-network-topology.md#step-7-verify-the-fix) still holds and demand is near the KEGOC target:

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
Total annual demand: 108.54 TWh
Load shedding: 0.00 TWh
```

- **Load shedding** at zero — no new islands after the grid rebuild.
- **Total demand 108.54 TWh** with **`scale: 1.005`** — above KEGOC **107.34 TWh**. A one-time tweak in Step 8 brings it back.

`n.statistics()` should show the same story: **Load shedding** supply at **0 TWh**. Installed capacity can differ from earlier runs if previously dropped island plants are now kept (`p_threshold_drop_isolated: false`), especially for hydro; compare against the Part 4/KEGOC reference table. Dispatch still reflects the **Part 5** `costs.marginal_cost` override — but, as Part 5 found, sourced Kazakhstan costs still favor coal, so don't expect coal share to drop much from that alone. Closing the remaining gap to KEGOC's ~22 TWh gas likely needs a must-run/CHP-style constraint or a CO₂ limit or emission price, both left as options for a follow-up tutorial.

---

## Step 8: Final calibration of `scale`

The Part 3 **`scale: 1.005`** was set before the Part 6–7 grid changes. Rescale once from your Step 7 total:

```python
target_TWh = 107.34  # KEGOC 2020
measured_TWh = 108.54  # Step 7 (with scale 1.005)
new_scale = 1.005 * target_TWh / measured_TWh
print(f"scale: {new_scale:.4f}")  # → 0.9939
```

```yaml
load_options:
  scale: 0.994
```

`1.005 × 107.34 / 108.54` gives **0.9939**, which rounds to **0.994**.

In [Part 3](3-demand-data.md#step-3-calibrate-annual-demand-with-scale), **`scale: 1.005`** bundled two gaps in one multiplier: GEGIS vs KEGOC statistics, **and** demand lost when simplification dropped load on isolated buses (~1.2 TWh in the Part 2 grid). Parts 6–7 fixed the second problem — **fetch** reconnects islands, **35 kV** keeps regional lines, and **drop off** preserves island plants — so that load and capacity stay in the model. With the grid settled, re-calibration gives a slightly lower final multiplier: **`scale: 0.994`**.

Re-run the workflow if you change `scale`. Load shedding should stay at **~0 TWh** — scaling demand does not re-isolate buses.

---

## Recap

| Step | Config key | Value | Role |
|---|---|---|---|
| 2 | `clean_osm_data_options.threshold_voltage` | `35000` | Keep 35 kV+ OSM assets |
| 3 | `electricity.base_voltage` | `500` | Single voltage after `simplify_network` |
| 3 | `electricity.voltages` | `[110, 220, 500]` | Line-type matching on multi-voltage `base.nc` |
| 4 | `lines.ac_types` | KZ conductor map | 110 kV type + derated 220 kV **`s_nom`** in build |
| 4 | `lines.s_max_pu` | `0.5` | Conservative N-1 margin for a less-meshed grid (default `0.7`) |
| 6 | *(keep Part 6)* | `cluster_options.simplify_network` | Merge off, fetch on, drop off — catch remaining islands without deleting plants |
| 8 | `load_options.scale` | `0.994` (from `1.005` × 107.34/108.54) | **Final** match to **107.34 TWh** after grid is settled |

The base grid now reflects Kazakhstan's post-Soviet voltage scale and more conservative line ratings. Part 6's simplification settings (merge off, fetch on, drop off) remain a useful safety net for OSM gaps that 35 kV alone cannot close.

Demand ends at **`scale: 0.994`** — a re-calibrated match to KEGOC on a grid that no longer drops regional load. Demand, fleet, real-world costs (Part 5), network topology (Part 6), and now transmission detail are all grounded in Kazakhstan-specific data and KEGOC 2020 evidence — though the coal/gas dispatch split, as Part 5 found, still needs more than costs alone to close. A natural next step from this validated baseline is exploring a must-run/CHP-style constraint or policy scenarios — an explicit CO₂ limit or emission price — in a follow-up tutorial.
