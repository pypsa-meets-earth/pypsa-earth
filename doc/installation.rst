..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.

System Requirements
===================

Building the model with the scripts in this repository runs on a normal computer.
The entire energy model is relatively heavy and it easily requires > 40Gb of available memory in the disk (HHD/SSD);
the exact space requirements depend on the specific models under interest.


Software requirements
=====================

The complete list of software needed before installing PyPSA Africa is listed below.

- `Python 3 <https://www.python.org/>`_ **(mandatory)**: Python is used as our main programming language, thus its knowledge is mandatory.
  To refresh the knowledge, there are plenty of online courses free-of-charge, e.g. `CSDojo playlist <https://www.youtube.com/c/CSDojo/playlists>`_.
  Useful content to watch refer to numpy, pandas
- `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_ **(mandatory)**: in order to use packages in python,
  it is highly recommended to use a conda package manager, such as `Anaconda <https://docs.anaconda.com/>`_ or
  `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`__ (**recommended**). There are many things you can do wrong with conda. `This article <https://towardsdatascience.com/conda-essential-concepts-and-tricks-e478ed53b5b>`_ provides you a crystal clear explanation of conda (**excellent read**).  
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
- `Solver` **(mandatory)**: an optimization solver is needed to solve the mathematical problem that is build with the automated workflow.
  With the goal of supporting completely open source initiative, we focus on relying on Open-Source solvers, such as `CBC <https://projects.coin-or.org/Cbc>`_ ,
  `GLPK <https://www.gnu.org/software/glpk/>`_, `WinGLPK <http://winglpk.sourceforge.net/>`_ or `HiGHS <https://github.com/ERGO-Code/HiGHS>`_;
  to further improve performances, commercial solvers like `Gurobi <http://www.gurobi.com/>`_ or `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_
  (both commercial licenses with free academic options) can also be used. A recommended instruction to install the HiGHS solver is given `here <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_.
 

.. note::
  Be aware that the list of software listed above is only the prerequisite elements needed to successfully install the PyPSA Africa model.
  The complete list of recommended software and prerequisite needed to enjoy the full PyPSA Africa experience is listed in the 
  `Tutorial section <https://pypsa-meets-africa.readthedocs.io/en/latest/tutorial.html#prerequisites-and-learning-material>`_.
  Most of the dependencies needed will be automatically installed using the conda environments listed below

Clone the Repository
====================

First of all, clone the `PyPSA meets Africa repository <https://github.com/pypsa-meets-africa/pypsa-africa>`_ using the version control system ``git``.
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

    .../pypsa-africa % conda activate pypsa-africa

Sometimes the conda pypsa-africa installation can take more than one our
(indicates some conflicts in the packages). In case a quick installation is necessary,
you might also want to try out ``mamba``. It was observed to take around 5-15min to
install all pypsa-africa dependencies with it.

.. code::bash

    ... conda install -c conda-forge mamba

    .../pypsa-africa % mamba env create -f envs/environment.yaml
    
To use jupyter lab (new jupyter notebooks) **continue** with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks>`_ 
and test if your jupyter lab works:
    
.. code:: bash

    .../pypsa-africa % ipython kernel install --user --name=pypsa-africa

    .../pypsa-africa % jupyter lab

.. note::
  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-africa,
  thereby there is no need to install it manually.
  The snakemake included in the conda environment pypsa-africa installed with the above-mentioned procedure can be executed with the following procedure:

  .. code: bash

    .../pypsa-africa % .../pypsa-africa % conda activate pypsa-africa

    .../pypsa-africa % snakemake < any command here >


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

     .../pypsa-africa % conda activate pypsa-africa

     .../pypsa-africa % snakemake -j1 retrieve_databundle_light --force

5. In the file ``config.yaml`` set the option ``retrieve_databundle`` back to ``false`` and save the file:
   ``retrieve_databundle: false``

Once these tasks have been completed, the package is ready to use.
