..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

Software requirements
=====================

PyPSA meets Africa builds on top of several open-source packages, which are here recalled together with recommended sources to better understand the main dependencies.
In the case the user needs to revise your knowledge of some of the requirements, please, check out our learning material in the tutorial section.

Programs and programming language
---------------------------------

- `Python <https://www.python.org/>`_ **(mandatory)**: Python is used as our main programming language, thus its knowledge is mandatory.
  To refresh the knowledge, there are plenty of online courses free-of-charge, e.g. `CSDojo playlist <https://www.youtube.com/c/CSDojo/playlists>`_.
  Useful content to watch refer to numpy, pandas
- `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_ **(mandatory)**: in order to use packages in python,
  it is highly recommended to use a conda package manager, such as `Anaconda <https://docs.anaconda.com/>`_ or
  `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ (**recommended**)
- `Git <https://git-scm.com/>`__ **(mandatory)**: Git is a free open source system aimed at tracking changes in the code development 
  and enable to coordinate the parallel software development between many developers.
  It is mandatory to `learn the git basics <https://git-scm.com/doc>`_.
- `IDE Python` **(recommendation)**: in order to write python code, you need an Integrated Development Environment (IDE)
  that is a software used to write code. Any program can be used, however, we recommend `Visual Studio Code <https://code.visualstudio.com/>`_,
  which is freely available online.
  Other alternatives are also viable if you are familiar with them, such as `PyCharm <https://www.jetbrains.com/pycharm/>`_,
  however we recommend Visual Studio Code also given its easy to use interface with Git.
  *Note*: if you decide to use Visual Studio Code, check out the tutorial about how to use 
  `Git <https://code.visualstudio.com/docs/editor/versioncontrol#_git-support>`__ and `Github <https://code.visualstudio.com/docs/editor/github>`__ 
  in Visual Studio Code
- `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ **(mandatory)**: snakemake is a tool to create reproducible and scalable workflow procedures.
  Snakemake provides a set of functions and a reference language where a set of execution blocks can be easily defined
  with specific properties which are executed subsequently in an automated procedure.

Python packages
---------------------
- `PyPSA <https://pypsa.readthedocs.io/en/latest/>`_ **(mandatory)**: Python for Power Systems Analysis is a 
  modelling suite with the goal of providing a set of functions to ease the power systems analysis.
  In particular, it is used as a backbone for formulating the mathematical problem to be solved. 
- `PyPSA-Eur <https://pypsa-eur.readthedocs.io/en/latest/>`_ **(mandatory)**: PyPSA-Eur is an open source model of the European Transmission systems.
  This package has been used a the base for PyPSA meets Africa package, hence, its good knowledge is recommended 
  for successfully working with PyPSA meets Africa.
- `Atlite <https://atlite.readthedocs.io/en/latest/>`_ **(optional)**: atlite is used to generate time series for renewable energy sources,
  such as wind, solar and hydro. The package converts weather datasets into time series to be useful for energy studies


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
Make sure that ``conda`` is already installed on your system or install one of the following two distributions:

- `miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ (recommended), which is a mini version of `Anaconda <https://www.anaconda.com/>`__  
- `Anaconda <https://www.anaconda.com/>`__

For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the envs/environment.yaml file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa-africa % conda env create -f envs/environment.yaml

    .../pypsa-africa % conda activate pypsa_meets_africa
    
To use jupyter lab (new jupyter notebooks) **continue** with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks>`_ and test if your jupyter lab works:
    
.. code:: bash

    .../pypsa-africa % ipython kernel install --user --name=pypsa-africa

    .../pypsa_africa % jupyter lab 


Download data
=============

The entire distribution, including the data for the whole Africa, is very heavy (>40Gb) and it involves a large number of files.
To simplify the installation of the github folder, the main source code is available in the Github folder, whereas the data are stored in cloud.
The rule ``retrieve_databundle_light`` has been specifically developed to set up the raw data, and the procedure below guides in setting up the needed data.

1. Duplicate the file ``config.default.yaml`` and rename the copy as ``config.yaml``
2. Open file ``config.yaml`` using any text editor
3. Make sure that the option ``retrieve_databundle`` is set ``true``
   ``retrieve_databundle: true``

4. Execute the following code on the shell to download initial files. Please, note that around **20Gb zipped files will be downloaded**, 
   so make sure you have a stable connection, time and around 50 Gb available in your system. If no errors show up, then you can proceed.

   .. code:: bash

     .../pypsa-africa % conda activate pypsa_meets_africa

     .../pypsa-africa % snakemake -j1 retrieve_databundle_light --force

5. In the file ``config.yaml`` set the option ``retrieve_databundle`` back to ``false`` and save the file:
   ``retrieve_databundle: false``

Once these tasks have been completed, the package is ready to use.