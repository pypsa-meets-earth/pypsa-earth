<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# PyPSA-Earth. A Flexible Python-based Open Optimisation Model to Study Energy System Futures around the World.

<p align="left">
by
<a href="https://pypsa-meets-earth.github.io">
    <img src="https://github.com/pypsa-meets-earth/pypsa-meets-earth.github.io/raw/main/assets/img/logo.png" width="150">
<a/>
</p>

## Development Status: **Stable and Active**

[![Test workflows](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/test.yml/badge.svg)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/test.yml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-earth/badge/?version=latest)](https://pypsa-earth.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa-meets-earth/pypsa-earth)](https://api.reuse.software/info/github.com/pypsa-meets-earth/pypsa-earth)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.apenergy.2023.121096-blue)](https://doi.org/10.1016/j.apenergy.2023.121096)


**PyPSA-Earth: A Global Sector-Coupled Open-Source Multi-Energy System Model**

PyPSA-Earth is the first open-source global cross-sectoral energy system model with high spatial and temporal resolution. The workflow provide capabilities for modelling the energy systems of any country in the world, enabling large-scale collaboration and transparent analysis for an inclusive and sustainable energy future. PyPSA-Earth is suitable for both operational studies and capacity expansion studies. Its sector-coupled modeling capabilities enable features for the detailed optimization of multi-energy systems, covering electricity, heating, transport, industry, hydrogen and more.

All the data needed for a simulation are automatically and flexibly retrieved from open sources. This includes, in particular, energy demand across sectors, generation capacities, medium- to high-voltage networks, and renewable energy potentials. Custom datasets can also be integrated as needed, and kept private if required. At the same time, [PyPSA-Earth-Status](https://github.com/pypsa-meets-earth/pypsa-earth-status) provides functionality to share regional insights. If you are willing to contribute your regional expertise, feel free to open an issue there.

PyPSA-Earth is capable of providing the modelling evidence needed to translate the implications of energy scenarios into actionable regional strategies. By making this tool openly available, we aim to foster collaboration, innovation, and informed decision-making to support sustainable and efficient energy solutions worldwide.

Details on the model are available in the following academic publications:

- *power model* M. Parzen et all. "PyPSA-Earth: A new global open energy system optimization model demonstrated in Africa", Applied Energy, 341, 2023. https://doi.org/10.1016/j.apenergy.2023.121096
- *sector-coupled model* H. Abdel-Khalek et al. "PyPSA-Earth sector-coupled: A global open-source multi-energy system model showcased for hydrogen applications in countries of the Global South", Applied Energy, 383, 2025. https://doi.org/10.1016/j.apenergy.2025.125316

**PyPSA meets Earth is an independent research initiative developing a powerful energy system model for Earth.** We work on open data, open source modelling, open source solver support and open communities. Stay tuned and join our mission - We look for users, co-developers and leaders!

<p align="center">
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/ddf041d1b98ca8f8c310f1c6393ec426ab5594cf.png" width=30%>
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/940b2673cfc31c4a6f01b7908f546d39d67df27e.png" width=23.6%>
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/6af089c376b19b72ad148e4e4326c162b94db68f.png" width=34.5%>
</p>

<p align="center"><b> Figure:</b> Example power systems build with PyPSA-Earth.<br>See images of ~193 more countries at <a href="https://zenodo.org/records/10080766">Zenodo</a></p>


The diagram below depicts one representative clustered node for the sector-coupled model with its generation, storage and conversion technologies.

<p align="center">
  <img src="https://ars.els-cdn.com/content/image/1-s2.0-S0306261925000467-gr5_lrg.jpg" width=75%>
</p>

## Livetracker. Most popular global models:

<p align="center">
<a href="https://star-history.com/#pypsa-meets-earth/pypsa-earth&OSeMOSYS/osemosys_global&niclasmattsson/Supergrid&SGIModel/MUSE_OS&etsap-TIMES/TIMES_model&Date">
    <img src="https://api.star-history.com/svg?repos=pypsa-meets-earth/pypsa-earth,OSeMOSYS/osemosys_global,niclasmattsson/Supergrid,SGIModel/MUSE_OS,etsap-TIMES/TIMES_model&type=Date" width="75%">
<a/>

## How to get involved

There are multiple ways to get involved and learn more about our work:
1. **Join our forum** and communication platform on [**PyPSA-meets-Earth**](https://discord.gg/AnuJBk23FU) Discord Server
2. **Chat on Discord with us** in the following open meetings:
    - **General initiative meeting** for project news and [high-level code updates](https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing). Held every [fourth Thursday 16-17:00 (UK time)](https://drive.google.com/file/d/1naH4WwW9drkOkOJ3PLO4fyWdkZQi5-_w/view?usp=share_link) and is a perfect place to meet the community and get a high-level update on PyPSA ecosystem relevant for PyPSA-Earth developments.
    - **Weekly developers meetings**
        - Eastern-Hemisphere friendly *Morning meeting* every [Thursday at 09:00 (UK time)](https://drive.google.com/file/d/1PDdmjsKhzyGRo0_YrP4wPQkn2XTNh6jA/view?usp=share_link).
        - Western-Hemisphere friendly *Evening meeting* every [Thursday 16:00 (UK time)](https://drive.google.com/file/d/1gaLmyV4qGPXsogkeRcAPWjC0ESebUxU-/view?usp=share_link). Every forth Thursday is replaced by the General initiative meeting which has a more high-level perspective, but you can also join to discuss more particular questions.
3. **Look at public materials** at [**google Drive**](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing) to share to minutes, presentations, lists and documents. Feel gree to get a look!
4. **Notify your interest** to on-demand meetings:
    - On-demand meetings
        - Demand creation and prediction meeting
        - AI asset detection meeting
        - Outreach meeting for planning, discussing events, workshops, communication, community activities
5. Join us and **propose your stream**.

## Installation

1. Open your terminal at a location where you want to install pypsa-earth. Type the following in your terminal to download the package from GitHub:

   ```bash
      .../some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth.git
   ```
2. The python package requirements are curated in the `envs/environment.yaml` file.
   The environment can be installed using:

```bash
    .../pypsa-earth % conda env create -f envs/environment.yaml
```

   If the above takes longer than 30min, you might want to try mamba for faster installation:

```bash
    (base) conda install -c conda-forge mamba

    .../pypsa-earth % mamba env create -f envs/environment.yaml
```

3. For running the optimization one has to install the solver. We can recommend the open source HiGHs solver which installation manual is given [here](https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632).
4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

   ```bash
      .../pypsa-earth % ipython kernel install --user --name=pypsa-earth
      .../pypsa-earth % jupyter lab
   ```
5. Verify or install a java redistribution from the [official website](https://www.oracle.com/java/technologies/downloads/) or equivalent.
   To verify the successful installation the following code can be tested from bash:

   ```bash
      .../pypsa-earth % java -version
   ```

   The expected output should resemble the following:

   ```bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)
   ```

## Running the model in previous versions

The model can be run in previous versions by checking out the respective tag. For instance, to run the model in version 0.6.0, which is the last version before the recent PyPSA update, the following command can be used:

```bash
git checkout v0.6.0
```
After checking out the tag, the model can be run as usual. Please make sure to use the environment built for the respective version.



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





## Training

- We recently updated some [hackathon material](https://github.com/pypsa-meets-earth/documentation) for PyPSA-Earth. The hackathon contains jupyter notebooks with exercises. After going through the 1 day theoretical and practical material you should have a suitable coding setup and feel confident about contributing.
- The get a general feeling about the PyPSA functionality, we further recommend going through the [PyPSA](https://github.com/PyPSA/PyPSA/tree/master/examples) and [Atlite](https://github.com/PyPSA/atlite/tree/master/examples) examples.

## Questions and Issues

- We are happy to answer questions and help with issues **if they are public**. Through being public the wider community can benefit from the raised points. Some tips. **Bugs** and **feature requests** should be raised in the [**GitHub Issues**](https://github.com/pypsa-meets-earth/pypsa-earth/issues/new/choose). **General workflow** or **user questions** as well as discussion points should be posted at the [**GitHub Discussions**](https://github.com/pypsa-meets-earth/pypsa-earth/discussions/categories/q-a) tab. Happy coding.

## Documentation

The documentation is available here: [documentation](https://pypsa-earth.readthedocs.io/en/latest/index.html).

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
                <a href="https://github.com/DeniseGiub">
                    <img src="https://avatars.githubusercontent.com/u/113139589?v=4" width="100;" alt="DeniseGiub"/>
                    <br />
                    <sub><b>DeniseGiub</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/danielelerede-oet">
                    <img src="https://avatars.githubusercontent.com/u/175011591?v=4" width="100;" alt="danielelerede-oet"/>
                    <br />
                    <sub><b>danielelerede-oet</b></sub>
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
                <a href="https://github.com/GridGrapher">
                    <img src="https://avatars.githubusercontent.com/u/127969728?v=4" width="100;" alt="GridGrapher"/>
                    <br />
                    <sub><b>GridGrapher</b></sub>
                </a>
            </td>
		</tr>
		<tr>
            <td align="center">
                <a href="https://github.com/drifter089">
                    <img src="https://avatars.githubusercontent.com/u/93286254?v=4" width="100;" alt="drifter089"/>
                    <br />
                    <sub><b>drifter089</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/Eric-Nitschke">
                    <img src="https://avatars.githubusercontent.com/u/152230633?v=4" width="100;" alt="Eric-Nitschke"/>
                    <br />
                    <sub><b>Eric-Nitschke</b></sub>
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
                <a href="https://github.com/Emre-Yorat89">
                    <img src="https://avatars.githubusercontent.com/u/62134151?v=4" width="100;" alt="Emre-Yorat89"/>
                    <br />
                    <sub><b>Emre-Yorat89</b></sub>
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
                <a href="https://github.com/asolavi">
                    <img src="https://avatars.githubusercontent.com/u/131155817?v=4" width="100;" alt="asolavi"/>
                    <br />
                    <sub><b>asolavi</b></sub>
                </a>
            </td>
            <td align="center">
                <a href="https://github.com/rajesh-ieg">
                    <img src="https://avatars.githubusercontent.com/u/120284682?v=4" width="100;" alt="rajesh-ieg"/>
                    <br />
                    <sub><b>rajesh-ieg</b></sub>
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
                <a href="https://github.com/LucieRC">
                    <img src="https://avatars.githubusercontent.com/u/104382956?v=4" width="100;" alt="LucieRC"/>
                    <br />
                    <sub><b>LucieRC</b></sub>
                </a>
            </td>
		</tr>
		<tr>
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
            <td align="center">
                <a href="https://github.com/AlexanderMeisinger">
                    <img src="https://avatars.githubusercontent.com/u/91368938?v=4" width="100;" alt="AlexanderMeisinger"/>
                    <br />
                    <sub><b>AlexanderMeisinger</b></sub>
                </a>
            </td>
		</tr>
	<tbody>
</table>
<!-- readme: collaborators,contributors,restyled-commits/- -end -->
