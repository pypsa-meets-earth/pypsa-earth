..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Clone the Repository
====================

First of all, clone the `PyPSA meets Africa repository <https://github.com/pz-max/pypsa_meets_africa>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces!
If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/pypsa-meets-africa/pypsa-africa.git

.. _deps:

Install Python Dependencies
===============================

PyPSA meets Africa relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the envs/environment.yaml file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa_africa % conda env create -f envs/environment.yaml

    .../pypsa_africa % conda activate pypsa_meets_africa
    
To use jupyter lab (new jupyter notebooks) **continue** with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks>`_ and test if your jupyter lab works:
    
.. code:: bash

    .../pypsa_africa % ipython kernel install --user --name=pypsa-africa

    .../pypsa_africa % jupyter lab 

