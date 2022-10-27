..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _installation:

##########################################
Installation
##########################################

The subsequently described installation steps are demonstrated as shell commands, where the path before the ``%`` sign denotes the
directory in which the commands following the ``%`` should be entered.


Clone the Repository
====================

First of all, clone the `PyPSA-Earth repository <https://github.com/pypsa-meets-earth/pypsa-earth/>`_ using the version control system ``git``.
The path to the directory into which the ``git repository`` is cloned, must **not** have any spaces!
If you do not have ``git`` installed, follow installation instructions `here <https://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`_.

.. code:: bash

    /some/other/path % cd /some/path/without/spaces

    /some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth.git


.. _deps:

Install Dependencies
===============================

Python Dependencies
--------------------------------

PyPSA-Earth relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which is a mini version of `Anaconda <https://www.anaconda.com/>`_ that includes only ``conda`` and its dependencies or make sure ``conda`` is already installed on your system with `conda -V`.
For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the `envs/environment.yaml <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/envs/environment.yaml>`_ file.
The environment can be installed and activated using

.. code:: bash

    .../pypsa-earth % conda env create -f envs/environment.yaml

    .../pypsa-earth % conda activate pypsa-earth

Note that activation is local to the currently open shell!
After opening a new terminal window, one needs to reissue the second command!

.. note::
    If you have troubles with a slow ``conda`` installation, we recommend to install
    `mamba <https://github.com/QuantStack/mamba>`_ as a fast drop-in replacement via

    .. code:: bash

      (base) conda install -c conda-forge mamba

    and then install the environment with

    .. code:: bash

      .../pypsa-earth % mamba env create -f envs/environment.yaml

Java Installation 
--------------------------------       

Verify or install a java redistribution from the [official website](https://www.oracle.com/java/technologies/downloads/) or equivalent.
To verify the successfull installation the following code can be tested from bash:

  .. code:: bash

      .../pypsa-earth % java -version

   The expected output should resemble the following:
   
   .. code:: bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)

System requirements
===================

Building the model with the scripts in this repository runs on a normal computer e.g. 8-16GBRAM.
Depending of the region of interest, different amounts of Gb storage (HHD/SSD) are required.
Africa requires about 40Gb, the world 250Gb, a single country between 1-10Gb.
We have also prepared a tutorial which should be below 10Gb.


Software requirements
=====================

The complete list of software needed before installing PyPSA Earth is listed below.

- `Python 3 <https://www.python.org/>`_ **(mandatory)**: Python is used as our main programming language, thus its knowledge is mandatory.
  To refresh the knowledge, there are plenty of online courses free-of-charge, e.g. `CSDojo playlist <https://www.youtube.com/c/CSDojo/playlists>`_.
  Useful content to watch refer to numpy, pandas
- `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html>`_ **(mandatory)**: in order to use packages in python,
  it is highly recommended to use a conda package manager, such as `Anaconda. <https://docs.anaconda.com/>`_ There are many things you can do wrong with conda. `This article <https://towardsdatascience.com/conda-essential-concepts-and-tricks-e478ed53b5b>`_ provides you a crystal clear explanation of conda (**excellent read**).  
- `Git <https://git-scm.com/>`__ **(mandatory)**: Git is a free open source system aimed at tracking changes in the code development 
  and enable to coordinate the parallel software development between many developers.
  It is mandatory to `learn the git basics <https://git-scm.com/doc>`_.
- `Java <https://www.oracle.com/java/technologies/downloads/>`_ **(mandatory)**: A Java distribution is needed for using `powerplantmatching` package.
  To have a better user experience, please install the redistribution from the website according to your operating system.
- `IDE Python` **(recommendation)**: in order to write python code, you need an Integrated Development Environment (IDE)
  that is a software used to write code. We recommend `Visual Studio Code <https://code.visualstudio.com/>`_, which is freely available online and provides an easy to use interface with Git. Obviously, any alternatives like `PyCharm <https://www.jetbrains.com/pycharm/>`_ or `Sublime <https://www.sublimetext.com/>`_ will work as well.

- `Solver` **(mandatory)**: an optimization solver is needed to solve the mathematical problem that is build with the automated workflow.
  With the goal of supporting completely open source initiative, we focus on relying on Open-Source solvers, such as `CBC <https://projects.coin-or.org/Cbc>`_ ,
  `GLPK <https://www.gnu.org/software/glpk/>`_, `WinGLPK <http://winglpk.sourceforge.net/>`_ or `HiGHS <https://github.com/ERGO-Code/HiGHS>`_;
  to further improve performances, commercial solvers like `Gurobi <http://www.gurobi.com/>`_ or `CPLEX <https://www.ibm.com/analytics/cplex-optimizer>`_
  (both commercial licenses with free academic options) can also be used. A recommended instruction to install the HiGHS solver is given `here <https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632>`_.

.. _deps2:

Install python dependencies
===============================

PyPSA Earth relies on a set of other Python packages to function.
We recommend using the package manager and environment management system ``conda`` to install them.
Make sure that ``conda`` is already installed on your system or install one of the following two distributions:
 
- `Anaconda <https://www.anaconda.com/>`__

For instructions for your operating system follow the ``conda`` `installation guide <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

The python package requirements are curated in the envs/environment.yaml file.
We install only `mamba` in the conda base environment to accelerate the installation.
**Please keep the base environment always clean, meaning don't install anything there!**
The environment can be installed in about 5-15min (reported by users) and activated using

.. code:: bash

    .../pypsa-earth (base) % conda install -c conda-forge mamba

    .../pypsa-earth (base) % mamba env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth

In case mamba did not work for you, you might want to try the traditional conda installation

.. code::bash

    .../pypsa-earth (base) % conda env create -f envs/environment.yaml

    .../pypsa-earth (pypsa-earth) % conda activate pypsa-earth

or use miniconda instead.
    
To use jupyter lab (new jupyter notebooks) **continue** with the `ipython kernel installation <http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks>`_ 
and test if your jupyter lab works:
    
.. code:: bash

    .../pypsa-earth (pypsa-earth) % ipython kernel install --user --name=pypsa-earth

    .../pypsa-earth (pypsa-earth) % jupyter lab

.. note::
  Please, make sure to have properly installed java, from the  `official website <https://www.oracle.com/java/technologies/downloads/>`__ or equivalent.

In linux only, that is possible through the following command.

.. code:: bash

    .../pypsa-earth (pypsa-earth) % conda install -c conda-forge openjdk

To verify the successful installation, you can verify that by using the following code.

.. code:: bash
     
    .../pypsa-earth (pypsa-earth) % java -version

The expected output should resemble the following text:

.. code:: bash
     java version "1.8.0_341"
     Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
     Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)

.. note::
   ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth,
   thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth installed with the above-mentioned procedure can be executed with the following procedure:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % .../pypsa-earth % conda activate pypsa-earth

    .../pypsa-earth (pypsa-earth) % snakemake < any command here >


Download data
=============

The entire distribution, including the data for most parts on Earth, is very heavy (>40Gb for Africa) and it involves a large number of files.
To simplify the installation of the github folder, the main source code is available in the Github folder, whereas the data are stored in cloud.
The rule ``retrieve_databundle_light`` has been specifically developed to set up the raw data, and the procedure below guides in setting up the needed data.

1. Duplicate the file ``config.default.yaml`` and rename the copy as ``config.yaml``
2. Open file ``config.yaml`` using any text editor
3. Make sure that the option ``retrieve_databundle`` is set ``true``
   ``retrieve_databundle: true``

4. Execute the following code on the shell to download initial files. Please, note that around **20Gb zipped files will be downloaded**, 
   so make sure you have a stable connection, time and around 50 Gb available in your system. If no errors show up, then you can proceed.

   .. code:: bash

     .../pypsa-earth (base) % conda activate pypsa-earth

     .../pypsa-earth (pypsa-earth) % snakemake -j1 retrieve_databundle_light --force

5. In the file ``config.yaml`` set the option ``retrieve_databundle`` back to ``false`` and save the file:
   ``retrieve_databundle: false``

Once these tasks have been completed, the package is ready to use.
