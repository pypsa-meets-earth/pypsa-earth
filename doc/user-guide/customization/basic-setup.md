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

Year-related parameters are also being used when specifying `load_options`:

```yaml
load_options:
  ssp: "ssp2-2.6"
  weather_year: 2013
  prediction_year: 2030
  scale: 1
```

The `weather_year` value corresponds to the weather data which was used to generate the electricity demand profiles for a selected area while `prediction_year` corresponds to the point of a [Shared Socioeconomic Pathways (SSP)](https://en.wikipedia.org/wiki/Shared_Socioeconomic_Pathways) trajectory. PyPSA-Earth uses SSP2-2.6 scenario within the Shared Socioeconomic Pathways framework, which is characterized by medium challenges to mitigation and adaptation efforts resulting in a global warming of approximately 2.6°C by the end of the 21st century.
The available values for `weather_year` and `prediction_year` can be checked by looking into `pypsa-earth/data/ssp2-2.6` folder. Currently, there are pre-calculated demand data for 2011, 2013, 2018 weather years and for 2030, 2040, 2050, and 2100 scenario prediction years.

## Use custom demand data

It is possible to implement custom demand profiles. It can be done by creating a dedicated custom demand sub-folder in a scenario folder `pypsa-earth/data/ssp2-2.6` and placing there a custom demand file. The name of a custom demand sub-folder should correspond to `weather_year` argument which stands in this case for general identification of a demand input. The name of a demand input file should be a continent name to which belongs a country of initerest. Both csv and nc formats can be used for demand files.

For example, to  `pypsa-earth/data/ssp2-2.6/2013_custom/`

For example, to provide custom inputs for Nigeria, you can put the time-series into `Africa.csv` file and place the file into `pypsa-earth/data/ssp2-2.6/2013_custom/` folder. To make it fetched, you'll need to specify `weather_year: 2013_custom` under `load_options`.

A format of the custom csv demand file should correspond to the csv files supplied with the model: there are `region_code`, `time`, `region_name` and `Electricity demand` columns, while a semicolon is used as a separator.

## Configure `atlite` section

To accurately model the temporal and spatial availability of renewable energy sources, PyPSA-Earth processes historical weather data using [atlite](https://atlite.readthedocs.io/en/latest/).

The workflow derives the annual snapshot range and weather cutout automatically from the configured `weather_year`. For example, `weather_year: 2018` automatically selects the `cutout-2018-era5` weather cutout and snapshots from `2018-01-01` to `2019-01-01`.

The technical cutout parameters are configured under `atlite.cutout`:

```yaml
weather_year: 2018

atlite:
  nprocesses: 4
  cutout:
    module: era5
    dx: 0.3 # cutout resolution
    dy: 0.3 # cutout resolution
```

Pre-built cutouts are currently available for weather year 2013 and can be downloaded by enabling `retrieve_cutout`. Building the 2013 cutout locally is also supported, although the workflow will warn that a pre-built version is available.

For other supported weather years, disable `retrieve_cutout` and enable `build_cutout`. The resulting cutout is created using the same weather year as the snapshots and weather-dependent demand profiles.

Renewable technologies configured with `cutout: auto` automatically use the cutout derived from `weather_year`; no manual replacement of cutout names is required.

The supported weather years depend on the selected demand source. GEGIS supports 2011, 2013, and 2018, while DemandCast supports weather years from 2000 to 2024.

To explore additional settings, refer to the [configuration](../configuration.md) page.
