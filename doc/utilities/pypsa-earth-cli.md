<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# PyPSA-Earth CLI

The PyPSA-Earth CLI is an interactive, menu-driven command line application
that helps new users navigate the workflow without having to hand-edit YAML
files or memorize `snakemake` commands. It is built with
[`typer`](https://typer.tiangolo.com/), [`InquirerPy`](https://inquirerpy.readthedocs.io/)
and [`rich`](https://rich.readthedocs.io/), and lives at
`scripts/non_workflow/pypsa_earth_cli.py`.

The CLI is organized as a main menu with four modules: a **Quiz Zone** built
around the use-case tutorials, a **Configuration Setup** wizard, a **Data
Retrieval** checker, and a **Run Model** launcher.

## Environment

The CLI's dependencies (`typer`, `inquirerpy`, `rich`) are not part of the
default `pypsa-earth` environment. They are bundled in the dedicated
`pypsa-earth-cli` [pixi](https://pixi.sh) environemnt. 

## Launching the CLI

The recommended way to launch the CLI is with pixi, using the
`pypsa-earth-cli` environment:

```bash
.../pypsa-earth % pixi run -e pypsa-earth-cli python scripts/non_workflow/pypsa_earth_cli.py
```

Running the script with no arguments prints a welcome message and opens the
interactive **main menu** as seen in the following image.

![PyPSA-Earth CLI main menu](pypsa_earth_cli_menu.png)

Each module is also exposed as a `typer`
subcommand, so it can be jumped to directly:

```bash
.../pypsa-earth % pixi run -e pypsa-earth-cli python scripts/non_workflow/pypsa_earth_cli.py quiz_zone
.../pypsa-earth % pixi run -e pypsa-earth-cli python scripts/non_workflow/pypsa_earth_cli.py config-setup
.../pypsa-earth % pixi run -e pypsa-earth-cli python scripts/non_workflow/pypsa_earth_cli.py run-model
```

> **Tip:** Append `--help` to list all available subcommands and options,
> e.g. `python scripts/non_workflow/pypsa_earth_cli.py --help`.

## Main menu

| # | Module | Description |
|---|--------|--------------|
| 1 | Quiz Zone | Quiz based on the use-case guide for PyPSA-Kazakhstan |
| 2 | Configuration Setup | Customize configuration parameters for a PyPSA-Earth model run |
| 3 | Data Retrieval | Fetch external data required for a PyPSA-Earth model run |
| 4 | Run Model | Run the PyPSA-Earth model |
| 5 | Exit | Close the application |

## Modules

### 1. Quiz Zone (`quiz_zone`)

An interactive quiz that mirrors the PyPSA-Kazakhstan use-case
documentation, split into eight use-cases: Baseline model, Analyze results,
Demand, Generation, Transmission, CO2 limits, Costs, and Return (to the main
menu). Questions and answers are defined in
`configs/pypsa_earth_cli/tutorial_questions.yaml`.

For each question, the user is prompted (free text or a multiple-choice
menu) until the correct answer is entered, and can optionally request a
hint if one is defined. Scoring per question starts at 1 point, with a
0.25-point penalty for each wrong attempt and a further 0.25-point penalty
if a hint is requested; the score is capped between 0 and 1 and summed
across the module.

For some of the quiz modules for e.g., the "Baseline model" use-case, the answers are also mapped onto config
parameters: once the questionnaire is completed, the CLI asks for a run
name and a solver choice and writes the result to `config.KZ.yaml`, then
offers to launch the model run immediately (module 4) with that file.

### 2. Configuration Setup (`config-setup`)

A guided editor for `config.default.yaml` that avoids exposing the full
config file at once. The workflow is:

1. Select a **user group**, defined in
   `configs/pypsa_earth_cli/user_groups.yaml`, which determines which
   top-level config sections are shown:

   | User group | Config sections exposed |
   |------------|--------------------------|
   | Sandbox | `countries`, `enable`, `run`, `scenario`, `snapshots` |
   | Specialist | `countries`, `enable`, `run`, `scenario`, `snapshots`, `foresight`, `cluster_options`, `build_shape_options`, `clean_osm_data_options` |
   | Strategist | `results_dir`, `summary_dir`, `countries`, `enable`, `custom_data` |

2. Pick a config option from that list. Nested options (e.g. `scenario`)
   open a submenu of their flattened child parameters (e.g.
   `scenario.simpl`, `scenario.ll`, `scenario.clusters`); flat options are
   edited directly.
3. Enter a new value for the selected parameter. Values for parameters that
   are lists in the original config are split on `,` into a list.
4. Repeat, or select **return to main menu** to finish.

The updated configuration is written to `config.cli_updated.yaml`, leaving
`config.default.yaml` untouched.

> **Note:** To add or change which parameters a user group can edit, modify
> `configs/pypsa_earth_cli/user_groups.yaml`.

### 3. Data Retrieval

Launches `scripts/non_workflow/databundle_cli.py` as a subprocess. It checks
the databundles listed in the active config against the files already
present on disk, prints a checklist table (files present are shown in
green, missing files in red), and lets the user interactively retrieve
missing data:

- `check` — refresh the checklist
- `all` — retrieve all missing databundles
- `rerun` — retrieve every databundle again
- `bundle_...` — retrieve one or more specific databundles by name
- **Enter** (empty input) — exit the retrieval loop

See [Retrieve Data Bundle](../user-guide/rules-reference/download-and-filter/retrieve-databundle-light.md)
for details on the underlying databundle configuration.

### 4. Run Model (`run-model`)

Runs the tutorial electricity model via `snakemake`. The user is prompted
for:

- a config file (any `config*.yaml` file found in the repository root, or
  the path passed to `run-model` programmatically),
- the number of cores to use,
- whether to run through `pixi` or `conda`.

The CLI always targets the `solve_all_networks` rule and builds the
equivalent command, for example:

```bash
pixi run snakemake -c 4 solve_all_networks --configfile config.KZ.yaml --rerun-incomplete
```

> **Note:** `run-model` is scoped to the elec-only tutorial workflow
> (`solve_all_networks`). For other targets, run `snakemake` directly as
> described in [Running Models](../user-guide/customization/running-models.md).

### 5. Exit

Prints a farewell message and closes the CLI.

## Related files

- CLI entry point: `scripts/non_workflow/pypsa_earth_cli.py`
- Data retrieval submodule: `scripts/non_workflow/databundle_cli.py`
- Quiz questions: `configs/pypsa_earth_cli/tutorial_questions.yaml`
- User group definitions: `configs/pypsa_earth_cli/user_groups.yaml`
- Environment definition: `pixi.toml` (`pypsa-earth-cli` feature/environment)
