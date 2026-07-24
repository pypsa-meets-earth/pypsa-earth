<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# FAQ

Common issues when running PyPSA-Earth, especially on a first workflow run. For general tooling tips, see also [Software Hints](software-hints.md).

---

## Cutout download failed (`retrieve_cutout`)

### Symptoms

- Snakemake stops at the `retrieve_cutout` rule.
- The log mentions a timeout, connection reset, SSL error, or incomplete download.
- The expected cutout file is missing under `cutouts/`.

On a first run this step downloads a **pre-built ERA5 weather cutout** (often **15–25 GB**, depending on region). A slow or unstable connection can interrupt the download.

### Where the file should end up

The path depends on your config:

| Setting | Expected cutout path (default name) |
|---|---|
| `run.shared_cutouts: true` | `cutouts/cutout-2013-era5.nc` |
| `run.shared_cutouts: false` | `cutouts/<run.name>/cutout-2013-era5.nc` |

For the [Kazakhstan use case](../tutorials/use-cases/1-baseline-model.md) (`run.name: "KZ"`, `shared_cutouts: false`), the file must be:

```
cutouts/KZ/cutout-2013-era5.nc
```

Check the log for the exact path:

```
logs/<run.name>/retrieve_cutout/cutout-2013-era5.log
```

### What to do

Download the cutout manually and place it where Snakemake expects it:

1. **Find your bundle** in [`configs/bundle_config.yaml`](https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/configs/bundle_config.yaml). Search for your country code under `countries:` in the `cutouts` bundles. For Kazakhstan, the bundle is `bundle_cutouts_northeurasia`.
2. **Download the archive** from the URL listed under that bundle's `urls:` (for example `gdrive` or `zenodo`).
3. **Unzip at the project root** (the folder that contains `Snakefile`). The archive creates `cutouts/cutout-2013-era5.nc`.
4. **Move the file if your run uses a per-country cutout folder** (`shared_cutouts: false`):

    ```bash
    mkdir -p cutouts/KZ
    mv cutouts/cutout-2013-era5.nc cutouts/KZ/
    ```

    Replace `KZ` with your `run.name` if different.

5. **Turn off automatic cutout retrieval** so Snakemake does not try to download again. Add to your config file (for example `config.KZ.yaml`):

    ```yaml
    enable:
      retrieve_cutout: false
    ```

6. **Re-run the workflow**:

    ```bash
    snakemake --cores 4 solve_all_networks --configfile config.KZ.yaml
    ```

If the download was interrupted but no partial file remains, you can also simply **retry** the Snakemake run — a stable connection often succeeds on a second attempt.

!!! note
    Do not set both `retrieve_cutout: true` and `build_cutout: true` while a cutout file already exists — Snakemake will refuse to overwrite it. After a successful manual download, keep `retrieve_cutout: false`.

---

## Other data download failed (`retrieve_databundle_light`)

If Snakemake fails while downloading **non-cutout** files (demand, OSM helpers, costs, etc.), try the databundle CLI from the project root:

```bash
python scripts/non_workflow/databundle_cli.py
```

It lists missing files, shows download links, and can retry the download. See also [Basic Setup](../user-guide/customization/basic-setup.md).

This helper applies to `retrieve_databundle_light` only — it does **not** download weather cutouts. Use the section above for cutout problems.
