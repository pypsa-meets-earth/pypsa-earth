# PyPSA meets Africa

## Development Status: **Prototype**

[![GitHub Super-Linter](https://github.com/pz-max/pypsa-meets-africa/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Build Status](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci.yaml/badge.svg)](https://github.com/pypsa-meets-africa/pypsa-africa/actions)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-africa/pypsa-africa)

PyPSA meets Africa is a free and open source software initiative aiming to develop a powerful energy system model for Africa. The tool which is currently under development and will be heavily based on [PyPSA](https://pypsa.readthedocs.io/en/latest/) and [PyPSA-Eur](https://pypsa-eur.readthedocs.io/en/latest/). In 2022 we will focus on Earth wide expansion. Stay tuned and join our mission!

A short presentation about our project and its aims is given on our [website](https://pypsa-meets-africa.github.io/). There you can also sign-up to our Newsletter. Watch our latest discussion with African leaders about [Open Energy System Modelling in Africa: State of the Art and Future Opportunities](https://www.youtube.com/watch?v=E0V0T4U9nmQ). Let's work together for a better future.

<p align="center">
  <img src="doc/img/africa_osm_map.jpeg" width="600">
</p>

## Get involved

There are multiple ways to get involved and learn more about our work. That's how we organise ourselves:

- [**Discord NEW! (Open)**](https://discord.gg/AnuJBk23FU)
  - Chat with the community, team up on features, exchange with developers, code in voice channels
      <p align="left">
        <a href="https://discord.gg/AnuJBk23FU">
          <img src="https://discord.com/assets/cb48d2a8d4991281d7a6a95d2f58195e.svg" width="150">
        <a/>
      </p>
- **General code meeting (Open)**
  - every second Thursday 16-17:00 (UK time) <a href="https://drive.google.com//uc?id=1Xre_N0SioLsehFoMuBS10J4xEWRc-lSW&export=download" >download .ics</a>
  - join for project news and high-level code updates
  - meeting on [Zoom](https://ed-ac-uk.zoom.us/j/89486720170), password: energy101
  - [open agenda](https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing). See what we will discuss. Invited members have edit rights.
- **Buddy talk (Open)**
  - every Friday between 17-18:00 (UK time)
  - book a 20min meeting with Max to discuss anything you like
  - booking link: [app.autobook.me/max-parzen/pypsa-meets-africa](https://app.autobook.me/max-parzen/pypsa-meets-africa) (developed by @mnm-matin)
- **Specific code meeting (by invitation)**
  - current meetings: 1) AI asset detection (WP6) and 2) Data workflow updates (covering WP1-WP5)
  - every week Friday 1) 14-14:45 (UK time), 2) 15-16:00 (UK time)
  - updates, demos, task distribution, weekly targets, Q&A
- **Outreach meeting (by invitation)**
  - every second week
  - planning, discussing events, workshops, communication, community activities
- [**Google Drive (by invitation)**](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)
  - access to presentations, lists, documents


## Installation

1) Open your terminal at a location where you want to install pypsa-africa. Type the following in your terminal to download the package from GitHub:
``` bash
    .../some/path/without/spaces % git clone https://github.com/pypsa-meets-africa/pypsa-africa.git
```
2) The python package requirements are curated in the `envs/environment.yaml` file.
The environment can be installed using:
``` bash
    .../pypsa-africa % conda env create -f envs/environment.yaml
```  
3) For running the optimization one has to install the solver. We can recommend the open source HiGHs solver which installation manual is given [here](https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632).

4) To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:
``` bash
    .../pypsa-africa % ipython kernel install --user --name=pypsa-africa

    .../pypsa-africa % jupyter lab
```
## Test run on tutorial
- In the folder open a terminal/command window to be located at this path `~/pypsa-africa/`
- Activate the environment `conda activate pypsa-africa`
- Rename config.tutorial.yaml to config.yaml. For instance in Linux:
  ``` bash
  mv config.tutorial.yaml config.yaml 
  ```
- Run a dryrun of the Snakemake workflow by typing simply in the terminal:
  ``` bash
  snakemake -j 1 solve_all_networks -n
  ```
  Remove the -n to do a real run. Follow the tutorial of PyPSA-Eur 1 and 2 on [YouTube](https://www.youtube.com/watch?v=ty47YU1_eeQ) to continue with an analysis.
  
## Documentation

The documentation is available here: [documentation](https://pypsa-meets-africa.readthedocs.io/en/latest/index.html).

## Collaborators

<!-- https://github.com/marketplace/actions/contribute-list -->

<!-- readme: collaborators,contributors -start -->
<table>
<tr>
    <td align="center">
        <a href="https://github.com/hazemakhalek">
            <img src="https://avatars.githubusercontent.com/u/26235356?v=4" width="100;" alt="hazemakhalek"/>
            <br />
            <sub><b>Hazemakhalek</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/jarry7">
            <img src="https://avatars.githubusercontent.com/u/27745389?v=4" width="100;" alt="jarry7"/>
            <br />
            <sub><b>Jarrad Wright</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/fneum">
            <img src="https://avatars.githubusercontent.com/u/29101152?v=4" width="100;" alt="fneum"/>
            <br />
            <sub><b>Fabian Neumann</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/euronion">
            <img src="https://avatars.githubusercontent.com/u/42553970?v=4" width="100;" alt="euronion"/>
            <br />
            <sub><b>Euronion</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/mnm-matin">
            <img src="https://avatars.githubusercontent.com/u/45293386?v=4" width="100;" alt="mnm-matin"/>
            <br />
            <sub><b>Mnm-matin</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/desenk">
            <img src="https://avatars.githubusercontent.com/u/48335263?v=4" width="100;" alt="desenk"/>
            <br />
            <sub><b>Desen Kirli</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/pz-max">
            <img src="https://avatars.githubusercontent.com/u/61968949?v=4" width="100;" alt="pz-max"/>
            <br />
            <sub><b>Max Parzen</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Cesare-Caputo">
            <img src="https://avatars.githubusercontent.com/u/62548290?v=4" width="100;" alt="Cesare-Caputo"/>
            <br />
            <sub><b>Cesare-Caputo</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Nayara2020">
            <img src="https://avatars.githubusercontent.com/u/64689686?v=4" width="100;" alt="Nayara2020"/>
            <br />
            <sub><b>Nayara2020</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Ay-Khi">
            <img src="https://avatars.githubusercontent.com/u/65019030?v=4" width="100;" alt="Ay-Khi"/>
            <br />
            <sub><b>Ayman Alkhirbash</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/koen-vg">
            <img src="https://avatars.githubusercontent.com/u/74298901?v=4" width="100;" alt="koen-vg"/>
            <br />
            <sub><b>Koen Van Greevenbroek</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Hazem-IEG">
            <img src="https://avatars.githubusercontent.com/u/87850910?v=4" width="100;" alt="Hazem-IEG"/>
            <br />
            <sub><b>Hazem-IEG</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/energyLS">
            <img src="https://avatars.githubusercontent.com/u/89515385?v=4" width="100;" alt="energyLS"/>
            <br />
            <sub><b>EnergyLS</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/davide-f">
            <img src="https://avatars.githubusercontent.com/u/67809479?v=4" width="100;" alt="davide-f"/>
            <br />
            <sub><b>Davide-f</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/restyled-commits">
            <img src="https://avatars.githubusercontent.com/u/65077583?v=4" width="100;" alt="restyled-commits"/>
            <br />
            <sub><b>Restyled Commits</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Tooblippe">
            <img src="https://avatars.githubusercontent.com/u/805313?v=4" width="100;" alt="Tooblippe"/>
            <br />
            <sub><b>Tobie</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/LukasFrankenQ">
            <img src="https://avatars.githubusercontent.com/u/55196140?v=4" width="100;" alt="LukasFrankenQ"/>
            <br />
            <sub><b>Lukas Franken</b></sub>
        </a>
    </td></tr>
</table>
<!-- readme: collaborators,contributors -end -->
