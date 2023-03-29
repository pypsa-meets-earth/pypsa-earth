

# PyPSA-Earth-Sec: A Sector-Coupled Open Optimisation Model of the Global Energy System

## Development Status: Version 0.0.2

[![Status Linux](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-linux.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-linux.yaml)
[![Status Mac](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-mac.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-mac.yaml)
[![Status Windows](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-windows.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth-sec/actions/workflows/ci-windows.yaml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-meets-earth/badge/?version=latest)](https://pypsa-meets-earth.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth-sec)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth-sec/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth-sec/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/VHH8TCwn)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)

**Disclaimer: PyPSA-Earth-Sec is still under development.**

The workflow is adaped to work smoothly for the following countries: Morocco, Namibia, Nigeria and Benin. The spatial and temporal resolution of the model are felxible. It's advisable to use more than 3 nodes per country and a timestep not smaller than 3-hours. 


Currently, no real sectoral demand data is used for the country inspected, instead, we use dummy data. The collection, compilation and processing of real data is underway.


The model now includes the following energy carriers: **electricity**, **hydrogen**, **gas**, *oil** and biomass as well as **carbon** as gas emissions and feedstock for the synthesis of the carriers.

The demand sectors covered are: **residential**, **industry**, **land transport**, **aviation**, **shipping**, **services** and agriculture.

The diagram below depicts one representative clustered node showing the combination of carriers and sectors covered in the model as well as the generation and conversion technologies included. 

![alt text](https://github.com/pypsa-meets-earth/pypsa-earth-sec/blob/main/docs/pes_v0.0.2.png?raw=true)



## Installation

1. Open your terminal at a location where you want to install pypsa-earth-sec. Type the following in your terminal to download the package and the dependency (pypsa-earth) from GitHub:

```bash
    .../some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth.git
    .../some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth-sec.git
```

2. The python package requirements are curated in the `envs/environment.yaml` file.
   The environment can be installed using:

```bash
    .../some/path/without/spaces % cd pypsa-earth
    .../pypsa-earth % conda env create -f envs/environment.yaml
```

3. For running the optimization one has to install the solver. We can recommend the open source HiGHs solver which installation manual is given [here](https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632).

4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

```bash
    .../pypsa-earth % ipython kernel install --user --name=pypsa-earth

    .../pypsa-earth % jupyter lab
```

## Test run

- In the folder *pypsa-earth* open a terminal window to be located at this path `~/pypsa-earth/`
- Rename config.default.yaml to config.yaml. For instance in Linux:
`mv config.default.yaml config.yaml`
- Open the file `config.yaml` and 
  - choose the country you want to model. For example
    `countries: ["MA"]`
  - Make sure `build_cutout` is set to `true`
- In the terminal run the following commands:
  - Activate the environment: `conda activate pypsa-earth`
  - Create the base network for the country by running : `snakemake -j 1 base_network`


- In the folder *pypsa-earth-sec* open a terminal/command window to be located at this path `~/pypsa-earth-sec/`
- Rename config.default.yaml to config.yaml. For instance in Linux:
`mv config.default.yaml config.yaml`
- Open the file `config.yaml` and follow the steps done before in pypsa-earth
  - choose the country you want to model. For example
    `countries: ["MA"]`
 - Open the file `config.pypsa-earth.yaml` and make sure `build_cutout` is set to `true`

- Run a dryrun of the Snakemake workflow by typing simply in the terminal:
  ```bash
  snakemake -j 1 solve_all_networks -n
  ```
  Remove the -n to do a real run. Follow the tutorial of PyPSA-Eur 1 and 2 on [YouTube](https://www.youtube.com/watch?v=ty47YU1_eeQ) to continue with an analysis.
