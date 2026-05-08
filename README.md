<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# PyPSA-Zambia. A Flexible Python-based Open Power System Model.

<p align="left">
by
<a href="https://pypsa-meets-earth.github.io">
    <img src="https://github.com/pypsa-meets-earth/pypsa-meets-earth.github.io/raw/main/assets/img/logo.png" width="150">
<a/>
</p>

## Development Status: **Under Development**

[![Test workflows](https://github.com/open-energy-transition/pypsa-zambia/actions/workflows/test_zambia.yml/badge.svg)](https://github.com/open-energy-transition/pypsa-zambia/actions/workflows/test_zambia.yml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-zambia/badge/?version=latest)](https://pypsa-zambia.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/open-energy-transition/pypsa-zambia)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/open-energy-transition/pypsa-zambia)](https://api.reuse.software/info/github.com/open-energy-transition/pypsa-zambia)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/open-energy-transition/pypsa-zambia/main.svg)](https://results.pre-commit.ci/latest/github/open-energy-transition/pypsa-zambia/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.apenergy.2023.121096-blue)](https://doi.org/10.1016/j.apenergy.2023.121096)

## Installation

1. Open your terminal at a location where you want to install pypsa-zambia. Type the following in your terminal to download the package from GitHub:

   ```bash
   .../some/path/without/spaces % git clone https://github.com/open-energy-transition/pypsa-zambia.git
   ```
2. The python package requirements are curated in the `envs/{your operating system}64.lock.yaml` file.
   - On linux, the environment can be installed using:
     ```bash
     .../pypsa-zambia % conda env create -f envs/linux-64.lock.yaml
     ```
   - On newest macOS (arm-based), the environment can be installed using:
     ```bash
     .../pypsa-zambia % conda env create -f envs/osx-arm64.lock.yaml
     ```
     On non-arm macOS, the environment can be installed using:
     ```bash
     .../pypsa-zambia % conda env create -f envs/osx-64.lock.yaml
     ```
    - On Windows, the environment can be installed using:
      ```bash
      .../pypsa-zambia % conda env create -f envs/win-64.lock.yaml
     ```

   If the above takes longer than 30 min, you might want to try mamba for faster installation:

   ```bash
   (base) conda install -c conda-forge mamba

   .../pypsa-zambia % mamba env create -f envs/{{your operating system}}64.lock.yaml
   ```
3. (optional) In step 2, three solvers are installed: HiGHs, glpk and gurobi. HiGHs is the recommended open-source solver. Gurobi is generally faster, but requires a license for full functionality, which is [freely available to academics](https://www.gurobi.com/features/academic-named-user-license/) (see instructions website).

4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

   ```bash
   .../pypsa-zambia % ipython kernel install --user --name=pypsa-earth
   .../pypsa-zambia % jupyter lab
   ```
5. Verify or install a java redistribution from the [official website](https://www.oracle.com/java/technologies/downloads/) or equivalent.
   To verify the successful installation the following code can be tested from bash:

   ```bash
   .../pypsa-zambia % java -version
   ```

   The expected output should resemble the following:

   ```bash
   java version "1.8.0_341"
   Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
   Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)
   ```

## Customisation

To build a regional cutout, you can use pre-compiled configuration files:
- "configs/build_cutout_tutorial_zambia_config.yaml" for a smaller tutorial-styled cutout (good for testing)
- "configs/build_cutout_zambia_config.yaml" for a full-scale cutout

To use them, you need to go through the following steps:
1. Comment-out a line in Snakemake which fetches a default configuration
`configfile: "configs/validation_dispatch_zambia.yaml"``
2. Add a name of a suitable cutout-building configuration file to Smakemake under `#configfile: "configs/validation_dispatch_zambia.yaml"`, e.g. `configfile: "configs/build_cutout_zambia_config.yaml"`
3. Set-up Copernicus API
4. run `build_cutout` rule

## Validation

A set of notebooks which provides ingights on outputs of the model and ways to analyse them are available in [zambia-dev-notebooks](https://github.com/open-energy-transition/zambia-dev-notebooks) repo.

## Test run on tutorial

- In the folder open a terminal/command window to be located at this path `~/pypsa-earth/`
- Activate the environment `conda activate pypsa-earth`
- Rename config.tutorial.yaml to config.yaml. For instance in Linux:
  ```bash
  mv config.tutorial.yaml config.yaml
  ```
- Run a dryrun of the Snakemake workflow by typing simply in the terminal:
  ```bash
  snakemake -j 1 solve_all_networks -n
  ```

  Remove the -n to do a real run. Follow the tutorial of PyPSA-Eur 1 and 2 on [YouTube](https://www.youtube.com/watch?v=ty47YU1_eeQ) to continue with an analysis.

## Questions and Issues

- We are happy to answer questions and help with issues **if they are public**. Through being public the wider community can benefit from the raised points. Some tips. **Bugs** and **feature requests** should be raised in the [**GitHub Issues**](https://github.com/open-energy-transition/pypsa-zambia/issues/new/choose). **General workflow** or **user questions** as well as discussion points should be posted at the [**GitHub Discussions**](https://github.com/open-energy-transition/pypsa-zambia/discussions/categories/q-a) tab. Happy coding.

## Documentation

Specific documentation for PyPSA-Zambia is under development [here](https://pypsa-zambia.readthedocs.io/en/latest/index.html). General documentation for PyPSA-Earth is available [here](https://pypsa-earth.readthedocs.io/en/latest/index.html) and can be used to get information on installation and customisation procedures, design of the model and philosophy of the workflow, as well as to access tutorial and support materials.

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
                <a href="https://github.com/Femkemilene">
                    <img src="https://avatars.githubusercontent.com/u/26096675?v=4" width="100;" alt="Femkemilene"/>
                    <br />
                    <sub><b>Femkemilene</b></sub>
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
                <a href="https://github.com/carlosfv92">
                    <img src="https://avatars.githubusercontent.com/u/103258059?v=4" width="100;" alt="carlosfv92"/>
                    <br />
                    <sub><b>carlosfv92</b></sub>
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
                <a href="https://github.com/arizeosalac">
                    <img src="https://avatars.githubusercontent.com/u/177637669?v=4" width="100;" alt="arizeosalac"/>
                    <br />
                    <sub><b>arizeosalac</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/AlexanderMeisinger">
                    <img src="https://avatars.githubusercontent.com/u/91368938?v=4" width="100;" alt="AlexanderMeisinger"/>
                    <br />
                    <sub><b>AlexanderMeisinger</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/asolavi">
                    <img src="https://avatars.githubusercontent.com/u/131155817?v=4" width="100;" alt="asolavi"/>
                    <br />
                    <sub><b>asolavi</b></sub>
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
                <a href="https://github.com/LucieRC">
                    <img src="https://avatars.githubusercontent.com/u/104382956?v=4" width="100;" alt="LucieRC"/>
                    <br />
                    <sub><b>LucieRC</b></sub>
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
                <a href="https://github.com/gianvicolux">
                    <img src="https://avatars.githubusercontent.com/u/123154558?v=4" width="100;" alt="gianvicolux"/>
                    <br />
                    <sub><b>gianvicolux</b></sub>
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
                <a href="https://github.com/milyas009">
                    <img src="https://avatars.githubusercontent.com/u/144870279?v=4" width="100;" alt="milyas009"/>
                    <br />
                    <sub><b>milyas009</b></sub>
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
                <a href="https://github.com/PierreKara1">
                    <img src="https://avatars.githubusercontent.com/u/160237120?v=4" width="100;" alt="PierreKara1"/>
                    <br />
                    <sub><b>PierreKara1</b></sub>
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
                <a href="https://github.com/siddharth-krishna">
                    <img src="https://avatars.githubusercontent.com/u/10712637?v=4" width="100;" alt="siddharth-krishna"/>
                    <br />
                    <sub><b>siddharth-krishna</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/squoilin">
                    <img src="https://avatars.githubusercontent.com/u/4547840?v=4" width="100;" alt="squoilin"/>
                    <br />
                    <sub><b>squoilin</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/juli-a-ko">
                    <img src="https://avatars.githubusercontent.com/u/126512394?v=4" width="100;" alt="juli-a-ko"/>
                    <br />
                    <sub><b>juli-a-ko</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/ollie-bell">
                    <img src="https://avatars.githubusercontent.com/u/56110893?v=4" width="100;" alt="ollie-bell"/>
                    <br />
                    <sub><b>ollie-bell</b></sub>
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
                <a href="https://github.com/saikumarvasa100-hash">
                    <img src="https://avatars.githubusercontent.com/u/228767710?v=4" width="100;" alt="saikumarvasa100-hash"/>
                    <br />
                    <sub><b>saikumarvasa100-hash</b></sub>
                </a>
            </td>
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
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/huyhoang-mike">
                    <img src="https://avatars.githubusercontent.com/u/109945762?v=4" width="100;" alt="huyhoang-mike"/>
                    <br />
                    <sub><b>huyhoang-mike</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/huyhoang-mike">
                    <img src="https://avatars.githubusercontent.com/u/109945762?v=4" width="100;" alt="huyhoang-mike"/>
                    <br />
                    <sub><b>huyhoang-mike</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/HanaElattar">
                    <img src="https://avatars.githubusercontent.com/u/87770004?v=4" width="100;" alt="HanaElattar"/>
                    <br />
                    <sub><b>HanaElattar</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Vamsipriya22">
                    <img src="https://avatars.githubusercontent.com/u/188459113?v=4" width="100;" alt="Vamsipriya22"/>
                    <br />
                    <sub><b>Vamsipriya22</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/flacombe">
                    <img src="https://avatars.githubusercontent.com/u/5690599?v=4" width="100;" alt="flacombe"/>
                    <br />
                    <sub><b>flacombe</b></sub>
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
