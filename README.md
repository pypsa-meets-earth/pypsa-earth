<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# PyPSA-Earth-OSM: enabling synergies between OSM enhansment and energy modelling.

<p align="left">
by
<a href="https://pypsa-meets-earth.github.io">
    <img src="https://github.com/pypsa-meets-earth/pypsa-meets-earth.github.io/raw/main/assets/img/logo.png" width="150">
<a/>
</p>


## Development Status: **In progress**

[![Test workflows](https://github.com/pypsa-meets-earth/pypsa-earth-osm/actions/workflows/test.yml/badge.svg)](https://github.com/pypsa-meets-earth/pypsa-earth-osm/actions/workflows/test.yml)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth-osm)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa-meets-earth/pypsa-earth-osm)](https://api.reuse.software/info/github.com/pypsa-meets-earth/pypsa-earth-osm)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth-osm/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)

That is an experimental modeling implementation intended to be used for bridging OSM mapping work with enhansement of the energy modelling workflow.

To run the stable modelling workflow for any country of the workd, please refer to [PyPSA-Earth repository](https://github.com/pypsa-meets-earth/pypsa-earth).

## How Earth-OSM Works

PyPSA-Earth-OSM leverages OpenStreetMap (OSM) data to build power grid topology models for energy system analysis. The workflow extracts, processes, and validates power infrastructure data from OSM, enabling researchers to analyze electricity networks using crowd-sourced geographic information.

### Data Sources and Workflow

The Earth-OSM workflow supports three different approaches for obtaining OpenStreetMap data:

#### 1. Latest OSM Data (Default)
Downloads the most current version of OSM data directly from OpenStreetMap servers. This ensures access to the latest updates from the OSM community, including recent infrastructure additions and corrections.

#### 2. Historical OSM Data
Enables analysis of past grid configurations by downloading OSM data from a specific date. This is valuable for:
- Temporal analysis and grid evolution studies
- Validation against historical records
- Reproducibility of previous research
- Comparison of infrastructure development over time

Historical data is available from approximately 2012 onwards, depending on OSM data availability.

#### 3. Custom OSM Data
Allows users to provide their own OpenStreetMap data files (in `.pbf` format). This option is useful for:
- Working with locally cached OSM data
- Using pre-processed or validated OSM datasets
- Offline workflows or environments with limited internet access
- Testing with specific OSM data snapshots

### Configuration Options

The `osm_data` section in the configuration file (`config.default.yaml` or `config.yaml`) controls how Earth-OSM retrieves and processes OpenStreetMap data:

```yaml
osm_data:
  source: "latest"  # Options: "latest", "historical", or "custom"

  # For historical data
  target_date: "2020-01-01"  # Format: YYYY-MM-DD

  # For custom data
  custom_path:
    pbf: "data/custom/osm/pbf"      # Path to custom .pbf files (required)
    power: "data/custom/osm/power"  # Path to processed power files (optional)
```

**Key Configuration Parameters:**

- **`source`**: Determines the OSM data source
  - `"latest"`: Downloads current OSM data (default)
  - `"historical"`: Downloads OSM data from a specific date
  - `"custom"`: Uses user-provided OSM files

- **`target_date`**: Specifies the date for historical data (format: YYYY-MM-DD)
  - Only used when `source: "historical"`
  - Example dates: `"2018-01-01"`, `"2020-06-15"`, `"2022-12-31"`

- **`custom_path`**: Defines paths for custom OSM data files
  - **`pbf`**: Directory containing raw OSM data files (`.osm.pbf` format)
    - Files must follow naming pattern: `<country>-latest.osm.pbf`
    - Example: `bolivia-latest.osm.pbf`, `nigeria-latest.osm.pbf`
  - **`power`**: Directory for pre-processed power infrastructure JSON files (optional)

### Data Processing Pipeline

Once OSM data is obtained, Earth-OSM processes it through several stages:

1. **Download/Load**: Retrieves OSM data based on the configured source
2. **Filter**: Extracts power infrastructure elements (substations, lines, cables, generators)
3. **Clean**: Validates and standardizes the data using configurable thresholds and rules
4. **Build Network**: Constructs the power network topology, connecting buses and lines
5. **Validate**: Performs quality checks and consistency validation


## Collaborators

<!-- https://github.com/marketplace/actions/contribute-list -->

<!-- readme: collaborators,contributors,restyled-commits/- -start -->
<table>
	<tbody>
		<tr>
            <td align="center">
                <a href="https://github.com/FabianHofmann">
                    <img src="https://avatars.githubusercontent.com/u/19226431?v=4" width="100;" alt="FabianHofmann"/>
                    <br />
                    <sub><b>FabianHofmann</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/fneum">
                    <img src="https://avatars.githubusercontent.com/u/29101152?v=4" width="100;" alt="fneum"/>
                    <br />
                    <sub><b>fneum</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/ekatef">
                    <img src="https://avatars.githubusercontent.com/u/30229437?v=4" width="100;" alt="ekatef"/>
                    <br />
                    <sub><b>ekatef</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/euronion">
                    <img src="https://avatars.githubusercontent.com/u/42553970?v=4" width="100;" alt="euronion"/>
                    <br />
                    <sub><b>euronion</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Justus-coded">
                    <img src="https://avatars.githubusercontent.com/u/44394641?v=4" width="100;" alt="Justus-coded"/>
                    <br />
                    <sub><b>Justus-coded</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/mnm-matin">
                    <img src="https://avatars.githubusercontent.com/u/45293386?v=4" width="100;" alt="mnm-matin"/>
                    <br />
                    <sub><b>mnm-matin</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/GbotemiB">
                    <img src="https://avatars.githubusercontent.com/u/48842684?v=4" width="100;" alt="GbotemiB"/>
                    <br />
                    <sub><b>GbotemiB</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/martacki">
                    <img src="https://avatars.githubusercontent.com/u/53824825?v=4" width="100;" alt="martacki"/>
                    <br />
                    <sub><b>martacki</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/LukasFrankenQ">
                    <img src="https://avatars.githubusercontent.com/u/55196140?v=4" width="100;" alt="LukasFrankenQ"/>
                    <br />
                    <sub><b>LukasFrankenQ</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/pz-max">
                    <img src="https://avatars.githubusercontent.com/u/61968949?v=4" width="100;" alt="pz-max"/>
                    <br />
                    <sub><b>pz-max</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/davide-f">
                    <img src="https://avatars.githubusercontent.com/u/67809479?v=4" width="100;" alt="davide-f"/>
                    <br />
                    <sub><b>davide-f</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/koen-vg">
                    <img src="https://avatars.githubusercontent.com/u/74298901?v=4" width="100;" alt="koen-vg"/>
                    <br />
                    <sub><b>koen-vg</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/Eddy-JV">
                    <img src="https://avatars.githubusercontent.com/u/75539255?v=4" width="100;" alt="Eddy-JV"/>
                    <br />
                    <sub><b>Eddy-JV</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/hazemakhalek">
                    <img src="https://avatars.githubusercontent.com/u/87850910?v=4" width="100;" alt="hazemakhalek"/>
                    <br />
                    <sub><b>hazemakhalek</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/energyLS">
                    <img src="https://avatars.githubusercontent.com/u/89515385?v=4" width="100;" alt="energyLS"/>
                    <br />
                    <sub><b>energyLS</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/AnasAlgarei">
                    <img src="https://avatars.githubusercontent.com/u/101210563?v=4" width="100;" alt="AnasAlgarei"/>
                    <br />
                    <sub><b>AnasAlgarei</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/yerbol-akhmetov">
                    <img src="https://avatars.githubusercontent.com/u/113768325?v=4" width="100;" alt="yerbol-akhmetov"/>
                    <br />
                    <sub><b>yerbol-akhmetov</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/doneachh">
                    <img src="https://avatars.githubusercontent.com/u/132910766?v=4" width="100;" alt="doneachh"/>
                    <br />
                    <sub><b>doneachh</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/danielelerede-oet">
                    <img src="https://avatars.githubusercontent.com/u/175011591?v=4" width="100;" alt="danielelerede-oet"/>
                    <br />
                    <sub><b>danielelerede-oet</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/DeniseGiub">
                    <img src="https://avatars.githubusercontent.com/u/113139589?v=4" width="100;" alt="DeniseGiub"/>
                    <br />
                    <sub><b>DeniseGiub</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/finozzifa">
                    <img src="https://avatars.githubusercontent.com/u/167071962?v=4" width="100;" alt="finozzifa"/>
                    <br />
                    <sub><b>finozzifa</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/virio-andreyana">
                    <img src="https://avatars.githubusercontent.com/u/114650479?v=4" width="100;" alt="virio-andreyana"/>
                    <br />
                    <sub><b>virio-andreyana</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Tomkourou">
                    <img src="https://avatars.githubusercontent.com/u/5240283?v=4" width="100;" alt="Tomkourou"/>
                    <br />
                    <sub><b>Tomkourou</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Eric-Nitschke">
                    <img src="https://avatars.githubusercontent.com/u/152230633?v=4" width="100;" alt="Eric-Nitschke"/>
                    <br />
                    <sub><b>Eric-Nitschke</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/GridGrapher">
                    <img src="https://avatars.githubusercontent.com/u/127969728?v=4" width="100;" alt="GridGrapher"/>
                    <br />
                    <sub><b>GridGrapher</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/drifter089">
                    <img src="https://avatars.githubusercontent.com/u/93286254?v=4" width="100;" alt="drifter089"/>
                    <br />
                    <sub><b>drifter089</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/glenkiely-ieg">
                    <img src="https://avatars.githubusercontent.com/u/99269783?v=4" width="100;" alt="glenkiely-ieg"/>
                    <br />
                    <sub><b>glenkiely-ieg</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/ljansen-iee">
                    <img src="https://avatars.githubusercontent.com/u/47030274?v=4" width="100;" alt="ljansen-iee"/>
                    <br />
                    <sub><b>ljansen-iee</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Emre-Yorat89">
                    <img src="https://avatars.githubusercontent.com/u/62134151?v=4" width="100;" alt="Emre-Yorat89"/>
                    <br />
                    <sub><b>Emre-Yorat89</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/giacfalk">
                    <img src="https://avatars.githubusercontent.com/u/36954873?v=4" width="100;" alt="giacfalk"/>
                    <br />
                    <sub><b>giacfalk</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/Ekaterina-Vo">
                    <img src="https://avatars.githubusercontent.com/u/99509555?v=4" width="100;" alt="Ekaterina-Vo"/>
                    <br />
                    <sub><b>Ekaterina-Vo</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/lkstrp">
                    <img src="https://avatars.githubusercontent.com/u/62255395?v=4" width="100;" alt="lkstrp"/>
                    <br />
                    <sub><b>lkstrp</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/TosinGeorge">
                    <img src="https://avatars.githubusercontent.com/u/78568233?v=4" width="100;" alt="TosinGeorge"/>
                    <br />
                    <sub><b>TosinGeorge</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Ly0n">
                    <img src="https://avatars.githubusercontent.com/u/6413976?v=4" width="100;" alt="Ly0n"/>
                    <br />
                    <sub><b>Ly0n</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Tooblippe">
                    <img src="https://avatars.githubusercontent.com/u/805313?v=4" width="100;" alt="Tooblippe"/>
                    <br />
                    <sub><b>Tooblippe</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/arizeosalac">
                    <img src="https://avatars.githubusercontent.com/u/177637669?v=4" width="100;" alt="arizeosalac"/>
                    <br />
                    <sub><b>arizeosalac</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/SermishaNarayana">
                    <img src="https://avatars.githubusercontent.com/u/156903227?v=4" width="100;" alt="SermishaNarayana"/>
                    <br />
                    <sub><b>SermishaNarayana</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/cpschau">
                    <img src="https://avatars.githubusercontent.com/u/124347782?v=4" width="100;" alt="cpschau"/>
                    <br />
                    <sub><b>cpschau</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/gianvicolux">
                    <img src="https://avatars.githubusercontent.com/u/123154558?v=4" width="100;" alt="gianvicolux"/>
                    <br />
                    <sub><b>gianvicolux</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/choiHenry">
                    <img src="https://avatars.githubusercontent.com/u/51810088?v=4" width="100;" alt="choiHenry"/>
                    <br />
                    <sub><b>choiHenry</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/LucieRC">
                    <img src="https://avatars.githubusercontent.com/u/104382956?v=4" width="100;" alt="LucieRC"/>
                    <br />
                    <sub><b>LucieRC</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/carlosfv92">
                    <img src="https://avatars.githubusercontent.com/u/103258059?v=4" width="100;" alt="carlosfv92"/>
                    <br />
                    <sub><b>carlosfv92</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/rajesh-ieg">
                    <img src="https://avatars.githubusercontent.com/u/120284682?v=4" width="100;" alt="rajesh-ieg"/>
                    <br />
                    <sub><b>rajesh-ieg</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/asolavi">
                    <img src="https://avatars.githubusercontent.com/u/131155817?v=4" width="100;" alt="asolavi"/>
                    <br />
                    <sub><b>asolavi</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/rsparks3">
                    <img src="https://avatars.githubusercontent.com/u/30065966?v=4" width="100;" alt="rsparks3"/>
                    <br />
                    <sub><b>rsparks3</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/ollie-bell">
                    <img src="https://avatars.githubusercontent.com/u/56110893?v=4" width="100;" alt="ollie-bell"/>
                    <br />
                    <sub><b>ollie-bell</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/juli-a-ko">
                    <img src="https://avatars.githubusercontent.com/u/126512394?v=4" width="100;" alt="juli-a-ko"/>
                    <br />
                    <sub><b>juli-a-ko</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/squoilin">
                    <img src="https://avatars.githubusercontent.com/u/4547840?v=4" width="100;" alt="squoilin"/>
                    <br />
                    <sub><b>squoilin</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/siddharth-krishna">
                    <img src="https://avatars.githubusercontent.com/u/10712637?v=4" width="100;" alt="siddharth-krishna"/>
                    <br />
                    <sub><b>siddharth-krishna</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/pitmonticone">
                    <img src="https://avatars.githubusercontent.com/u/38562595?v=4" width="100;" alt="pitmonticone"/>
                    <br />
                    <sub><b>pitmonticone</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/PierreKara1">
                    <img src="https://avatars.githubusercontent.com/u/160237120?v=4" width="100;" alt="PierreKara1"/>
                    <br />
                    <sub><b>PierreKara1</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Netotse">
                    <img src="https://avatars.githubusercontent.com/u/89367243?v=4" width="100;" alt="Netotse"/>
                    <br />
                    <sub><b>Netotse</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/milyas009">
                    <img src="https://avatars.githubusercontent.com/u/144870279?v=4" width="100;" alt="milyas009"/>
                    <br />
                    <sub><b>milyas009</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/kma33">
                    <img src="https://avatars.githubusercontent.com/u/25573938?v=4" width="100;" alt="kma33"/>
                    <br />
                    <sub><b>kma33</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/jome1">
                    <img src="https://avatars.githubusercontent.com/u/49280197?v=4" width="100;" alt="jome1"/>
                    <br />
                    <sub><b>jome1</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/jessLryan">
                    <img src="https://avatars.githubusercontent.com/u/122939887?v=4" width="100;" alt="jessLryan"/>
                    <br />
                    <sub><b>jessLryan</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/jarry7">
                    <img src="https://avatars.githubusercontent.com/u/27745389?v=4" width="100;" alt="jarry7"/>
                    <br />
                    <sub><b>jarry7</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/HanaElattar">
                    <img src="https://avatars.githubusercontent.com/u/87770004?v=4" width="100;" alt="HanaElattar"/>
                    <br />
                    <sub><b>HanaElattar</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/EmreYorat">
                    <img src="https://avatars.githubusercontent.com/u/93644024?v=4" width="100;" alt="EmreYorat"/>
                    <br />
                    <sub><b>EmreYorat</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/AndreCNF">
                    <img src="https://avatars.githubusercontent.com/u/19359510?v=4" width="100;" alt="AndreCNF"/>
                    <br />
                    <sub><b>AndreCNF</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/AlexanderMeisinger">
                    <img src="https://avatars.githubusercontent.com/u/91368938?v=4" width="100;" alt="AlexanderMeisinger"/>
                    <br />
                    <sub><b>AlexanderMeisinger</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/AlessandroPampado99">
                    <img src="https://avatars.githubusercontent.com/u/156424082?v=4" width="100;" alt="AlessandroPampado99"/>
                    <br />
                    <sub><b>AlessandroPampado99</b></sub>
                </a>
            </td>
		</tr>
	<tbody>
</table>
<!-- readme: collaborators,contributors,restyled-commits/- -end -->
