# PyPSA meets Africa and PyPSA meets Earth

## Development Status: **Stable and Active**

[![Status Linux](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-linux.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-linux.yaml)
[![Status Mac](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-mac.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-mac.yaml)
[![Status Windows](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-windows.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-africa/pypsa-africa/actions/workflows/ci-windows.yaml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-meets-africa/badge/?version=latest)](https://pypsa-meets-africa.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-africa/pypsa-africa)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-africa/pypsa-africa/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-africa/pypsa-africa/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/VHH8TCwn)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)


PyPSA meets Africa is a free and open source software initiative aiming to develop a powerful energy system model for Africa. The tool was first released end of 2022 and is heavily based on [PyPSA](https://pypsa.readthedocs.io/en/latest/) and [PyPSA-Eur](https://pypsa-eur.readthedocs.io/en/latest/). In 2022 we will focus on Earth wide expansion. Stay tuned and join our mission - We look for users, co-developers and leaders!

A short presentation about our project and its aims is given on our [website](https://pypsa-meets-africa.github.io/). There you can also sign-up to our Newsletter. Watch our latest discussion with African leaders about [Open Energy System Modelling in Africa: State of the Art and Future Opportunities](https://www.youtube.com/watch?v=E0V0T4U9nmQ). Let's work together for a better future.

<p align="center">
  <img src="doc/img/africa_osm_map.png" width="600">
</p>

## Get involved

There are multiple ways to get involved and learn more about our work. That's how we organise ourselves:

- [**Discord NEW! (Open)**](https://discord.gg/AnuJBk23FU)
  - chat with the community, team up on features, exchange with developers, code in voice channels
  - registration and usage is for free
      <p align="left">
        <a href="https://discord.gg/AnuJBk23FU">
          <img src="https://discord.com/assets/cb48d2a8d4991281d7a6a95d2f58195e.svg" width="150">
        <a/>
      </p>
- **General code meeting (Open)**
  - every month Thursday 16-17:00 (UTC+0) <a href="https://drive.google.com/file/d/1YY2-oJrBSOGd-mAZV7AmyVhOG3LcM1ry/view?usp=sharing" >download .ics</a>
  - join for project news and high-level code updates
  - meeting hosted on Discord
  - [open agenda](https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing). See what we will discuss. Invited members have edit rights.
- **Buddy talk (Open)**
  - every Friday between 17-18:00 (UTC+0)
  - book a 20min meeting with Max to discuss anything you like
  - booking link: [app.autobook.me/max-parzen/pypsa-meets-africa](https://app.autobook.me/max-parzen/pypsa-meets-africa) (developed by @mnm-matin)
- **Specific code meeting (Open)**
  - meeting hosted on Discord
  - join updates, demos, Q&A's, discussions and the coordination of each work package
    1. Demand creation and prediction meeting, every Wednesday 21:00 UTC+0, <a href="https://drive.google.com/file/d/1OBGBwSWeuYzhssb6tDxyAPIzz3uUhsnv/view?usp=sharing" >download .ics</a>
    2. AI asset detection meeting, every Tuesday 15:30 UTC+0, <a href="https://drive.google.com/file/d/1U9LYvLBezaKn1IuELbKJts0DOGixQHOv/view?usp=sharing" >download .ics</a>
    3. Sector coupling meeting, every Thursday 09:00 UTC+0, <a href="https://drive.google.com/file/d/1TzUcilUdcsnsre7jTEcvyizftET93DQS/view?usp=sharing" >download .ics</a>
    4. Data workflow and architecture meeting, every Thursday 13:30 UTC+0, <a href="https://drive.google.com/file/d/1w9uT6AIC9MIJotFcbglfVFCw1TKwF6bg/view?usp=sharing" >download .ics</a> 
- **Outreach meeting (Open)**
  - every second week, Tuesday 17:00 UTC+0
  - planning, discussing events, workshops, communication, community activities
- [**Google Drive**](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)
  - access to minutes, presentations, lists, documents (access to minutes)

## Installation

1. Open your terminal at a location where you want to install pypsa-africa. Type the following in your terminal to download the package from GitHub:

```bash
    .../some/path/without/spaces % git clone https://github.com/pypsa-meets-africa/pypsa-africa.git
```

2. The python package requirements are curated in the `envs/environment.yaml` file.
   The environment can be installed using:

```bash
    .../pypsa-africa % conda env create -f envs/environment.yaml
```

3. For running the optimization one has to install the solver. We can recommend the open source HiGHs solver which installation manual is given [here](https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632).

4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

```bash
    .../pypsa-africa % ipython kernel install --user --name=pypsa-africa

    .../pypsa-africa % jupyter lab
```

## Test run on tutorial

- In the folder open a terminal/command window to be located at this path `~/pypsa-africa/`
- Activate the environment `conda activate pypsa-africa`
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

- We recently updated some [hackathon material](https://github.com/pypsa-meets-africa/pypsa-africa-hackathon) for PyPSA-Africa. The hackathon contains jupyter notebooks with exercises. After going through the 1 day theoretical and practical material you should have a suitable coding setup and feel confident about contributing.
- The get a general feeling about the PyPSA functionality, we further recommend going through the [PyPSA](https://github.com/PyPSA/PyPSA/tree/master/examples) and [Atlite](https://github.com/PyPSA/atlite/tree/master/examples) examples.

## Questions and Issues

- We are happy to answer questions and help with issues **if they are public**. Through being public the wider community can benefit from the raised points. Some tips. **Bugs** and **feature requests** should be raised in the [**GitHub Issues**](https://github.com/pypsa-meets-africa/pypsa-africa/issues/new/choose). **General workflow** or **user questions** as well as discussion points should be posted at the [**GitHub Discussions**](https://github.com/pypsa-meets-africa/pypsa-africa/discussions/categories/q-a) tab. Happy coding.

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
        <a href="https://github.com/Justus-coded">
            <img src="https://avatars.githubusercontent.com/u/44394641?v=4" width="100;" alt="Justus-coded"/>
            <br />
            <sub><b>Justus Ilemobayo</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/mnm-matin">
            <img src="https://avatars.githubusercontent.com/u/45293386?v=4" width="100;" alt="mnm-matin"/>
            <br />
            <sub><b>Mnm-matin</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/desenk">
            <img src="https://avatars.githubusercontent.com/u/48335263?v=4" width="100;" alt="desenk"/>
            <br />
            <sub><b>Desen Kirli</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/LukasFrankenQ">
            <img src="https://avatars.githubusercontent.com/u/55196140?v=4" width="100;" alt="LukasFrankenQ"/>
            <br />
            <sub><b>Lukas Franken</b></sub>
        </a>
    </td>
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
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/davide-f">
            <img src="https://avatars.githubusercontent.com/u/67809479?v=4" width="100;" alt="davide-f"/>
            <br />
            <sub><b>Davide-f</b></sub>
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
    </td>
    <td align="center">
        <a href="https://github.com/energyLS">
            <img src="https://avatars.githubusercontent.com/u/89515385?v=4" width="100;" alt="energyLS"/>
            <br />
            <sub><b>EnergyLS</b></sub>
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
        <a href="https://github.com/ekatef">
            <img src="https://avatars.githubusercontent.com/u/30229437?v=4" width="100;" alt="ekatef"/>
            <br />
            <sub><b>Ekaterina</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/giacfalk">
            <img src="https://avatars.githubusercontent.com/u/36954873?v=4" width="100;" alt="giacfalk"/>
            <br />
            <sub><b>Giacomo Falchetta</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Tooblippe">
            <img src="https://avatars.githubusercontent.com/u/805313?v=4" width="100;" alt="Tooblippe"/>
            <br />
            <sub><b>Jarrad Wright</b></sub>
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
        <a href="https://github.com/squoilin">
            <img src="https://avatars.githubusercontent.com/u/4547840?v=4" width="100;" alt="squoilin"/>
            <br />
            <sub><b>Sylvain Quoilin</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/stephenjlee">
            <img src="https://avatars.githubusercontent.com/u/11340470?v=4" width="100;" alt="stephenjlee"/>
            <br />
            <sub><b>Stephen J Lee</b></sub>
        </a>
    </td></tr>
</table>
<!-- readme: collaborators,contributors -end -->
