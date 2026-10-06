<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Basic Setup

A good starting point to customize your model are settings of the default configuration file `config.default.yaml`. You may want to do a reserve copy of your current configuration file and then overwrite it by a default configuration:

```bash
cp config.default.yaml config.yaml
```

## Specify the country/region of interest

The model can be adapted to include any country, multiple countries (e.g. `Nigeria` and `Benin`) or full continents (currently whole regions, such as `Africa`, `Asia`, `Europe`, `Oceania`, `NorthAmerica`, and `SouthAmerica`, are available for simulation) using `countries` argument:

```yaml
countries: ["NG", "BJ"]
```

The PyPSA meets Earth initiative also supports dedicated regional models:

| Countries | Regional model |
|---|---|
| BN, KH, ID, LA, MY, MM, PH, SG, TH, TL, VN | [PyPSA-ASEAN](https://github.com/pypsa-meets-earth/pypsa-asean) |
| US | [PyPSA-NorthAmerica](https://github.com/pypsa-meets-earth/pypsa-northamerica) |

These community-maintained models are soft forks of PyPSA-Earth. Users are encouraged to test and contribute to them, and developers of other regional soft forks are welcome to have their models listed here.

## Configure `enable` section to download/build data

For a successful model run, ensure the download of essential open-source data, including databundle and cost data, is activated in the `enable` section:

```yaml
enable:
retrieve_databundle: true  #  Recommended 'true', for the first run. Otherwise data might be missing.
retrieve_cost_data: true  # true: retrieves cost data from technology data and saves in resources/costs.csv, false: uses cost data in data/costs.csv
download_osm_data: true  # If 'true', OpenStreetMap data will be downloaded for the above given countries
build_natura_raster: false # If True, than an exclusion raster will be build
build_cutout: false
# If "build_cutout" : true, then environmental data is extracted according to `snapshots` date range and `countries`
```

After the initial run, it is recommended to set the retrieval of databundle and cost data to `false` to prevent unnecessary redownloading of data.

When `build_natura_raster: false` is utilized, the exclusion raster for protected areas is sourced from the pre-compiled `data/natura.tiff` file downloaded with the databundle. Conversely, if `build_natura_raster` is set to true, the exclusion raster, delineating areas where renewables cannot be installed, is computed using the `build_natura_raster rule`.

When using the weather year 2013, it is recommended to use default `build_cutout: false` because pre-compiled cutouts are automatically downloaded with `retrieve_databundle: true`.
On contrary, when simulating a weather year other than 2013, it is crucial to set `build_cutout: true` in order to generate custom cutouts. However, it is essential to first configure the `Copernicus Climate Data Store (CDS) API`. Detailed instructions for setting up the `Copernicus API` can be found in [customization_copernicus](copernicus-data.md).
After initial run and successful generation of custom cutouts, `build_cutout` can be switched to false to avoid reconstructing the cutout.

!!! note
    No need to configure the `Copernicus API` if the weather year 2013 is used, as pre-compiled cutouts are automatically downloaded.

!!! tip
    Additionally, if you encounter issues with failed `retrieve_databundle`, you can use the following script to debug it through the command line interface (CLI):

    ```bash
    python scripts/non_workflow/databundle_cli.py
    ```

## Specify the weather year scope

Likewise, the example's temporal scope can be restricted (e.g. to 7 days):

```yaml
snapshots:
  start: "2013-03-01"
  end: "2013-03-07"
  inclusive: "left" # end is not inclusive
```

Ensure that the selected date range aligns with the dates available in the cutout dataset. If the weather data within the cutouts corresponds to the year 2013, then the range of snapshots should fall within that same year.

## Specify the demand year

By default, the demand weather year is inferred from `snapshots.start`:

```yaml
load_options:
  source: "gegis"
  weather_year: derive_from_snapshots
  prediction_year: 2030
  scale: 1
```

For GEGIS, `weather_year` selects the weather conditions used to generate the demand profile, while `prediction_year` selects the socioeconomic scenario year. Supported weather years are 2011, 2013, and 2018. DemandCast provides demand years from 2000 to 2024 and does not use `prediction_year`.

An explicit integer, such as `weather_year: 2013`, selects the demand year independently of snapshots and cutout selection. The configured snapshot range is never changed.

Demand is aligned to snapshots by month, day, and hour. Weekdays are not preserved when the years differ. If the demand year has no February 29, its February 28 profile is reused for February 29. Missing required timestamps or demand values raise an error.

When deriving a year from snapshots, the range must remain within one calendar year. January 1 of the following year is accepted as an exclusive end boundary.

## Configure `atlite` section

PyPSA-Earth processes historical weather data using [atlite](https://atlite.readthedocs.io/en/latest/). The selected cutout must cover the snapshot dates.

Set `atlite.default: derive_from_snapshots` to generate the cutout name from the snapshot year and the module configured under `atlite.cutouts.derive_from_snapshots`:

```yaml
snapshots:
  start: "2018-01-01"
  end: "2019-01-01"
  inclusive: "left"

load_options:
  weather_year: derive_from_snapshots

atlite:
  nprocesses: 4
  default: derive_from_snapshots
  cutouts:
    derive_from_snapshots:
      module: era5
      dx: 0.3
      dy: 0.3
```

This selects `cutout-2018-era5`. Tutorial runs append `-tutorial`. Automatic cutout naming uses the snapshot year independently of an explicitly selected demand year. The planning horizon remains unchanged.

For a custom cutout, specify its name and definition explicitly:

```yaml
atlite:
  default: my_fancy_cutout
  cutouts:
    my_fancy_cutout:
      module: era5
      dx: 0.1
      dy: 0.1
      x: [130, 145]
      y: [30, 45]
```

Explicit names and definitions are preserved, including additional named cutouts. Renewable technologies with `cutout: auto` use the resolved `atlite.default`; explicit technology-specific cutout names remain unchanged.

Enable `retrieve_cutout` to download an available pre-built cutout. Available names are read from the configured non-tutorial databundles in the `cutouts` category. To build locally, disable `retrieve_cutout` and enable `build_cutout`. Choosing to build a cutout that is also available for retrieval produces an informational log message. For an existing local cutout, both options can be disabled.

Demand-source year restrictions apply to demand profiles. The weather dataset used for the cutout must independently cover the configured snapshot dates.

To explore additional settings, refer to the [configuration](../configuration.md) page.
