<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 7: Improve the Transmission Network

!!! note
    This tutorial assumes you have completed [Part 1](1-baseline-model.md) through [Part 6](6-network-topology.md). Your `config.KZ.yaml` should include the Part 4 generation fleet, the Part 5 `costs` overrides, and the Part 6 `clustering.simplify_network` settings and calibrated `load_options.scale: 0.994`.

## Introduction

[Part 6](6-network-topology.md) eliminated load shedding by reconnecting electrically isolated buses to the main grid (**merge off**, **drop off**, **fetch on**). That fix works, but it is still a simplification shortcut: load is reassigned to the nearest connected bus, not carried over explicitly modelled transmission lines.

In this tutorial we improve the **base transmission network** that PyPSA-Earth builds from OpenStreetMap (OSM):

1. **Kazakhstan's voltage hierarchy**: represent its post-Soviet/CIS voltage levels.
2. **A suitable OSM voltage floor (35 kV for KZ)**: retain regional sub-transmission links that a higher cutoff would drop.
3. **KZ line types**: conductor ratings that reflect older, single-circuit lines at 110 kV and 220 kV.
4. **Non-extendable lines**: stop the solver from cost-optimizing past the ratings from (3).

These settings live upstream in the workflow (`clean_osm_data`, `base_network`, `electricity`, `lines`, `scenario.ll`). Changing them triggers a **full rebuild** of the grid, slower than Part 6, but it gives a base topology and line capacities that are closer to Kazakhstan's physical network before clustering and solving.

Keep the Part 6 **`clustering`** block and the Part 5 **`costs`** block. Even with a better OSM grid, some small islands may remain; `merge` off, `fetch` on, and `drop` off still help preserve both served load and island generation. The `costs.marginal_cost` overrides from Part 5 are generator attributes, independent of topology, so they carry through the rebuild unchanged.

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
| **`osm.clean_osm_data`** | `clean_osm_data` | Which OSM lines and substations survive the first filter |
| **`base_network`** | `base_network` | Minimum voltage kept when building `base.nc` |
| **`electricity`** | `simplify_network`, `cluster_network`, … | Nominal voltage levels the model simplifies toward |
| **`lines`** | `base_network` | Conductor type per voltage → thermal **`s_nom`** |
| **`scenario.ll`** | `prepare_network` | How far lines may expand beyond that built **`s_nom`** |

See the [configuration reference](../../user-guide/configuration.md) for every key.

---

## Step 1: See why KZ uses a different voltage scale

Kazakhstan follows the **post-Soviet/CIS** voltage scale:

| Role | Kazakhstan (typical) | PyPSA-Earth default nominal levels |
|---|---|---|
| Sub-transmission / regional mesh | **35 kV**, **110 kV** | 132 kV |
| Transmission | **220 kV** | 220 kV, 300 kV, 380 kV |
| Extra-high-voltage backbone | **500 kV** | 500 kV, 750 kV |

**35 kV** is the critical layer for connectivity: in many regions it is the voltage that ties smaller substations into the wider grid. **110 kV** carries regional load between cities; **220 kV** and **500 kV** form the long-distance backbone.

The minimum voltage used to retain OSM lines and substations is configured with **`osm.clean_osm_data.threshold_voltage`**; check its current default in `config.default.yaml`. Cross-check this value against the voltage hierarchy and OSM data quality of the country being modelled. For Kazakhstan, **35 kV** provides better results because a higher threshold drops additional lines that connect some regions, removing parts of the sub-transmission mesh and contributing to the electrical islands observed in Parts 4–6.

---

## Step 2: Lower the OSM voltage floor to 35 kV

The only active voltage floor in the current workflow is **`osm.clean_osm_data.threshold_voltage`**, applied in **`clean_osm_data`** to lines, cables, and substations:

| Key | Block | Rule |
|---|---|---|
| `threshold_voltage` | `osm.clean_osm_data` | Filters raw OSM during `clean_osm_data` |

```yaml
osm:
  clean_osm_data:
    threshold_voltage: 35000   # [V] keep 35 kV+ lines and substations from OSM
```

---

## Step 3: Set KZ nominal voltage levels

```yaml
electricity:
  base_voltage: 500
  voltages: [110, 220, 500]
```

- **`base_voltage`**: voltage level that **`simplify_network`** maps the entire grid onto (500 kV for KZ).
- **`voltages`**: nominal levels used when assigning **`lines.ac_types`** to OSM lines in **`base_network`** (each line matched to the closest kV in this list).

In practice, this means `electricity.voltages` controls which voltage layers are represented in **`base.nc`** (for KZ: 110 / 220 / 500 kV, with 35 kV retained by the OSM filter in Step 2), while `base_voltage` controls the later simplification collapse. During **`simplify_network`**, buses and lines are mapped to a single nominal level (500 kV) and transformers are removed. The solved clustered network therefore shows one nominal bus voltage, but the lower-voltage layers still matter upstream because they affect connectivity and line capacities before that collapse.

For KZ, use `[110, 220, 500]` as the default working set. Add other levels (for example 132 or 380 kV) only if your OSM data or study scope clearly requires them.

---

## Step 4: Set KZ line types and ratings

Line thermal limits (**`s_nom`**) come from **`lines.ac_types`**: each voltage maps to a named conductor in PyPSA's line register. PyPSA-Earth picks the closest voltage match and computes capacity from the conductor's rated current and **`v_nom`**.

Defaults assume **modern European** bundles (e.g. two bundles at 220 kV). Kazakhstan's older corridors often run **single circuits** with lighter conductors. We encode that as a **modelling assumption**: lower **`s_nom`** at 220 kV and an explicit **110 kV** type (missing from the default table):

```yaml
lines:
  ac_types:
    110.: "149-AL1/24-ST1A 110.0"
    220.: "305-AL1/39-ST1A 110.0"
    500.: "490-AL1/64-ST1A 380.0"
```

| kV | KZ type | Effect |
|---|---|---|
| **110** | `149-AL1/24-ST1A 110.0` | Enables 110 kV line ratings, the lightest conductor PyPSA-Earth's pinned PyPSA registers at this voltage |
| **220** | `305-AL1/39-ST1A 110.0` | Single conductor vs default **2-bundle** (`i_nom` 0.74 vs 1.29 kA) → lower **`s_nom`** |
| **500** | `490-AL1/64-ST1A 380.0` | Single-circuit EHV vs default **4-bundle** |

!!! note
    `ac_types` values must exist in your installed PyPSA's line register (`pypsa<=0.30.3` here; see the [version-pinned line types table](https://docs.pypsa.org/v0.30.3/user-guide/components.html#line-types), not the `latest` docs). An unmatched name silently produces `NaN` **`s_nom`** instead of raising an error.

Tighter line limits can bind regional flows in the LP; looser defaults would under-state congestion on older 220 kV corridors.

---

## Step 5: Stop lines from being cost-optimized past their rating

**`ac_types`** sets what a line is *built* with, but by default PyPSA-Earth still lets the solver expand every line **beyond** that rating at capital cost. That happens in **`prepare_network`**, controlled by **`scenario.ll`** (the [`{ll}` wildcard](../../user-guide/wildcards.md#the-ll-wildcard) in every network filename, e.g. `..._lcopt_...`):

```python
if factor == "opt" or float(factor) > 1.0:
    n.lines["s_nom_extendable"] = True
```

The default **`ll: ["copt"]`** sets `factor = "opt"`, so every line, including the ones we just derated, becomes extendable. The solver can simply pay to build past our KZ conductor choices instead of respecting them, which defeats the point of Steps 1–4. Set `factor` to a fixed value ≤ 1.0 instead, so lines stay at their built capacity:

```yaml
scenario:
  ll: ["v1.0"]
```

The letter (`v`, `c`, or `l`) only decides *how* extra capacity would be allocated once `factor` is `"opt"` or `> 1`: `v` caps total build-out **volume** (MW·km) globally, `c` caps total **cost** globally, `l` caps each line individually at `factor` × its own capacity (full detail in the [wildcards reference](../../user-guide/wildcards.md#the-ll-wildcard)). At `factor = 1.0`, none of the three make any line extendable, so the letter is cosmetic here; `v1.0` reads most naturally as "no expansion of grid volume."

This changes the network filename from `..._lcopt_6h.nc` to `..._lv1.0_6h.nc`, a new file, not an overwrite of Part 6's.

---

## Step 6: Add the settings to `config.KZ.yaml`

Merge the blocks below with your existing Part 3–6 settings (`load_options`, `electricity` fleet keys, `costs`, `clustering`, …):

```yaml
scenario:
  ll: ["v1.0"]

osm:
  clean_osm_data:
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

## Step 7: Re-run the workflow

Because these keys feed **`clean_osm_data`**, **`base_network`**, and **`prepare_network`**, Snakemake rebuilds the grid from OSM onward, through demand, generation, simplify, cluster, and solve. Expect a **longer run** than Part 6 (similar to Part 1 in scope, but OSM and cutouts stay cached):

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

The solved network is written to `results/KZ/networks/elec_s_10_ec_lv1.0_6h.nc`.

---

## Step 8: Verify the solve

Reload the solved network in `analyze_kz.ipynb`. The OSM rebuild changes line capacities and topology, so check whether [Part 6](6-network-topology.md#step-7-verify-the-fix)'s zero-shedding result still holds now that lines carry country-specific, non-extendable ratings, and that demand lands on the KEGOC target. Since `scale: 0.994` already carried over from [Part 6](6-network-topology.md#step-8-final-calibration-of-scale), no further calibration is needed here:

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
Total annual demand: 107.35 TWh
Load shedding: 0.08 TWh
```

- **Total demand** lands almost exactly on KEGOC's **107.34 TWh** target, confirming `scale: 0.994` (carried over from Part 6) still holds once the transmission network is rebuilt.
- **Load shedding** stays small at **0.08 TWh**: non-extendable lines occasionally bind somewhere in the network, unlike Part 6's simplification shortcut, which reached exactly zero.

### Generation mix: the transmission fix, not costs, closes the gap

```python
supply = n.statistics()["Supply"].dropna() / 1e6  # TWh
supply = supply.drop(["Line", "Load"], errors="ignore")
print(supply.sort_values(ascending=False).to_string())
```

```
Generator    Coal                  75.02
             Combined-Cycle Gas    19.55
StorageUnit  Reservoir & Dam        7.90
Generator    Open-Cycle Gas         2.34
             Onshore Wind           1.37
             Solar                  1.09
             Load shedding          0.08
```

| Carrier | KEGOC 2020 | Part 5 (costs only) | Part 7 (+ country-specific transmission) |
|---|---|---|---|
| Coal | ~74.5 TWh | 89.89 TWh | **75.02 TWh** |
| Gas (CCGT + OCGT) | **~21.7 TWh** | 0.20 TWh | **21.89 TWh** |
| Hydro | ~9.5 TWh | 7.15 TWh | 7.90 TWh |
| Solar | ~1.25 TWh | 0.99 TWh | 1.09 TWh |
| Wind | ~1.1 TWh | 1.34 TWh | 1.37 TWh |
| Load shedding | n/a | 7.72 TWh | 0.08 TWh |
| **Total** | **~108.1 TWh** | ~107.3 TWh | **107.35 TWh** |

This is the payoff of Part 7: Part 5 already showed that sourced Kazakhstan costs alone do **not** explain KEGOC's gas dispatch; coal stayed cheaper on the margin regardless. Once lines carry country-specific KZ ratings (Step 4) and can no longer be cost-optimized around (Step 5), the picture flips: **gas lands at ~21.9 TWh**, matching KEGOC's ~21.7 TWh almost exactly, and coal drops by ~15 TWh. The gap was predominantly a **transmission-congestion** problem, not an economics problem: coal-heavy corridors simply can't deliver everywhere once the grid uses those ratings, so local gas plants pick up the difference. The residual **0.08 TWh** load shedding is the honest cost of that fidelity: a still-imperfect base network occasionally can't route around a binding line at all.

---

## Recap

| Step | Config key | Value | Role |
|---|---|---|---|
| 2 | `osm.clean_osm_data.threshold_voltage` | `35000` | Keep 35 kV+ OSM assets |
| 3 | `electricity.base_voltage` | `500` | Single voltage after `simplify_network` |
| 3 | `electricity.voltages` | `[110, 220, 500]` | Line-type matching on multi-voltage `base.nc` |
| 4 | `lines.ac_types` | KZ conductor map | 110 kV type + derated 220 kV **`s_nom`** in build |
| 5 | `scenario.ll` | `["v1.0"]` | Lines stay at built capacity, no cost-optimized expansion past Steps 1–4 (default `["copt"]`) |
| *(Part 6)* | *(keep Part 6)* | `clustering.simplify_network` | Merge off, fetch on, drop off; catch remaining islands without deleting plants |
| *(Part 6)* | *(keep Part 6)* | `load_options.scale` | `0.994`, already final; demand no longer changes here |

The base grid now reflects Kazakhstan's post-Soviet voltage scale and more conservative line ratings. Part 6's simplification settings (merge off, fetch on, drop off) remain a useful safety net for OSM gaps that 35 kV alone cannot close, and `scale: 0.994` carries through unchanged since Part 6 already finalized it.

Demand, fleet, country-specific costs (Part 5), network topology (Part 6), and now transmission detail are all grounded in Kazakhstan-specific data and KEGOC 2020 evidence. The coal/gas dispatch split, the open question since Part 4, is now resolved: Part 5 showed costs alone don't explain it (coal stays cheaper on the margin), but with country-specific, non-extendable transmission (Steps 4–5 here), gas lands close to KEGOC's ~22 TWh target. The dispatch gap was a **network** problem, not an economics problem. A small residual of load shedding remains as the honest cost of that fidelity.

---

## Conclusion

This completes Parts 1–7 of the Kazakhstan use-case series: from a baseline solve through demand and fleet calibration, country-specific costs, isolated-node fixes, and transmission detail. You now have a validated electricity-only model grounded in KEGOC 2020 evidence.

For other countries or further studies you may want to try options not covered here, for example:

- **CO₂ limits or prices**, applied via `scenario.opts` (e.g. `Co2L`) and the `co2` block
- **A custom busmap**, a user-supplied bus-to-cluster mapping instead of the algorithm-based clustering used so far, applied via `enable.custom_busmap`
- **Alternative clustering**, grouping buses by administrative regions (e.g. GADM) instead of the electrical-distance clustering used so far, applied via `clustering.alternative_clustering`
- **Capacity expansion**, re-enabled via `electricity.extendable_carriers` and `scenario.ll: ["copt"]`, once the baseline is trusted

See the [configuration](../../user-guide/configuration.md) and [wildcards](../../user-guide/wildcards.md) guides for details. Feel free to experiment with these on your own.
