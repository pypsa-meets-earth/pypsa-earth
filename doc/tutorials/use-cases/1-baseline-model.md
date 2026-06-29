<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Part 1: Build a Baseline Model

!!! note
    This tutorial assumes PyPSA-Earth is already installed and your environment is active.
    If not, see the [installation guide](../../home/installation.md) first.

## Introduction

Kazakhstan (KZ) is the world's ninth-largest country by area, stretching from the Caspian Sea in the west to the Altai Mountains in the east. Its power system is dominated by coal and gas, yet it holds vast untapped potential for wind and solar — making it a compelling case for energy transition modelling.

![Topographic map of Kazakhstan](https://upload.wikimedia.org/wikipedia/commons/thumb/3/3a/Kazakhstan_topographic_location_map_conic.png/960px-Kazakhstan_topographic_location_map_conic.png)

*Topographic map of Kazakhstan. Source: [NordNordWest](https://commons.wikimedia.org/wiki/User:NordNordWest) et al., [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Kazakhstan_topographic_location_map_conic.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0).*

In this first tutorial we will build and solve a baseline electricity model for Kazakhstan. We will create a minimal configuration file from scratch, adding only the settings that differ from the defaults. By the end you will have a solved network file you can inspect and a reproducible configuration to build on in the rest of the series.

---

## How configuration works

Before writing a single line, it is worth understanding how PyPSA-Earth reads configuration — because it explains why the file we are about to create is so short.

PyPSA-Earth loads configuration in layers. The Snakefile reads these files in order, with each subsequent file overriding only the keys it explicitly sets:

```
config.default.yaml   ← full set of defaults (do not edit this)
      ↓  merged with
config.yaml           ← base overrides (keep this minimal or empty)
      ↓  merged with
--configfile flag     ← your study-specific overrides (highest priority)
```

This means your study config can be very short — often just 10–20 lines. You never need to copy the entire default file. Pass your study config explicitly on the command line; Snakemake merges it on top of everything else.

!!! tip "Browsing the defaults"
    All default values live in `config.default.yaml` in the project root. Open it any time you want to know what a key does or what its default value is. For a structured, searchable version see the [Configuration reference](../../user-guide/configuration.md).

## Step 1: Create your configuration file

Time to get started. Create an empty file called `config.KZ.yaml` in the project root — this is the only file we will create in this tutorial. We will fill it in gradually over the next few steps, adding one setting at a time and explaining the reasoning behind each one.

=== "Linux / macOS"

    ```bash
    touch config.KZ.yaml
    ```

=== "Windows (PowerShell)"

    ```powershell
    New-Item config.KZ.yaml
    ```

=== "Windows (any)"

    Create the file manually: in VS Code open the project folder, right-click the root and choose **New File**, name it `config.KZ.yaml`, and save it empty.

---

## Step 2: Set the country

This is the most important line in any country study. Add it to `config.KZ.yaml`:

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.yaml:5:5"
```

That single line tells every rule in the workflow — shape building, OSM download, cutout extraction, demand estimation — to focus on Kazakhstan. Two-letter [ISO 3166-1 alpha-2](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2) codes are used throughout.

---

## Step 3: Give the run a name

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.yaml:7:9"
```

**Why does `run.name` matter?**
It is a prefix that PyPSA-Earth appends to every intermediate and result path. Without it all results land in a flat `results/` and `resources/` directory, so running two countries back-to-back would overwrite each other. With `name: "KZ"` all files are isolated under `results/KZ/` and `resources/KZ/`, making it straightforward to manage multiple studies in the same repository.

Setting `shared_cutouts: false` stores the weather cutout under `cutouts/KZ/` as well. Use `shared_cutouts: true` only when multiple runs cover the same geographic region and you want to reuse the same ERA5 file.

---

## Step 4: Set temporal resolution to 6-hourly

The temporal resolution of the optimisation is one of the settings controlled by the `opts` wildcard in `scenario`. The default `Co2L-3h` resamples the 8760-hour year to 3-hourly snapshots (2920 time steps). For a first run, 6-hourly resolution (`6h`) cuts the problem size in half and reduces solve time from roughly an hour to around ten minutes with the HiGHS solver.

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.yaml:11:12"
```

Once you are satisfied with the model setup you can increase resolution progressively: `6h` → `3h` → full hourly. The trade-off is always solve time against temporal accuracy.

---

## Step 5: Choose a solver

PyPSA-Earth supports several open-source and commercial LP solvers. For a baseline run HiGHS is recommended: it is free, ships with the conda environment, and handles networks of this size efficiently.

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.yaml:14:16"
```

**Overview of supported solvers:**

| Solver | License | Notes |
|---|---|---|
| **HiGHS** | Open source (MIT) | Recommended for most users; included in the conda environment |
| **Gurobi** | Commercial | Fastest for very large models; requires a licence |
| **CPLEX** | Commercial | IBM solver; requires a licence |
| **GLPK** | Open source (GPL) | Very slow; only for small test cases |
| **CBC** | Open source (EPL) | Reasonable for small/medium models |

Solver-specific options (threads, tolerances, barrier method) can be tuned under `solving.solver_options` in `config.default.yaml`. The HiGHS defaults already use the interior-point method with 4 threads, which is a good starting point.

---

## Your complete `config.KZ.yaml`

Take a moment to look at what you have built:

```yaml
--8<-- "tutorials/use-cases/snippets/config.KZ.yaml"
```

That is all — every other parameter (network voltages, renewable technology settings, cost assumptions, solver tolerances, …) is inherited unchanged from `config.default.yaml`. When you need details on a specific key, the [Configuration reference](../../user-guide/configuration.md) has the full list — later parts of this series will come back to individual sections as we calibrate the model.

Save this as `config.KZ.yaml` in the project root before running Snakemake, or [download the file](snippets/config.KZ.yaml){: download="config.KZ.yaml"} and copy it there.

---

## Step 6: Run the model

### Dry run first

Before executing anything, ask Snakemake to print the full execution plan without running a single rule:

```bash
conda activate pypsa-earth
snakemake --cores 1 solve_all_networks --configfile config.KZ.yaml --dryrun
```

**How to read the output:**

A successful dry run prints a numbered list of jobs, for example:

```
Job stats:
job                          count
-------------------------  -------
add_electricity                  1
add_extra_components             1
base_network                     1
build_bus_regions                1
...
solve_network                    1
total                           22
```

This tells you Snakemake has resolved the full dependency graph and knows exactly what to run. Check that:

- The rule list ends with `solve_network` — this confirms the target is reachable
- The `total` count looks reasonable (a fresh KZ run typically requires 15–25 jobs)
- No `MissingInputException` or `AmbiguousRuleException` errors appear

If you see errors, the most common causes are a missing `config.yaml` file (create an empty one with `touch config.yaml`) or a typo in `config.KZ.yaml`.

### Execute the workflow

Once the dry run looks correct, run for real:

```bash
snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
```

`--cores 4` allows Snakemake to run up to 4 independent rules in parallel and limits each rule to 4 CPU threads. Reduce to `--cores 1` if memory is limited.

**Expected runtime:** on a typical laptop or workstation, the full pipeline takes around **1 hour** with HIGHs solver, including the ERA5 cutout download (~20 GB). Note that this is a coarse 10-node model (`clusters: [10]` by default) — a realistic high-resolution model would take significantly longer.

!!! note "First run downloads data"
    On the first run, PyPSA-Earth downloads OSM network data, cost data, and a pre-compiled ERA5 weather cutout for 2013. Subsequent runs reuse the cached data.

### What a successful run looks like

When all rules complete, Snakemake prints a summary to the terminal:

```
22 of 22 steps (100%) done
Complete log: .snakemake/log/...
```

The solved network file is written to:

```
results/KZ/networks/elec_s_10_ec_lcopt_6h.nc
```

This is a PyPSA `Network` object stored as NetCDF. In **[Part 2](2-analyze-results.md)** we open it, inspect installed capacities and dispatch, and validate against national statistics.

## Recap

| Config key | Value set | Why |
|---|---|---|
| `countries` | `["KZ"]` | Scope all data to Kazakhstan |
| `run.name` | `"KZ"` | Isolate results in `results/KZ/` and `resources/KZ/` |
| `run.shared_cutouts` | `false` | Store cutout under `cutouts/KZ/` |
| `scenario.opts` | `[6h]` | 6-hourly resolution for a fast first solve |
| `solving.solver.name` | `highs` | Free, included, adequate for this model size |

That is everything needed to run a first Kazakhstan model. In **[Part 2](2-analyze-results.md)** we open the solved network and find out what the optimiser actually decided to build.
