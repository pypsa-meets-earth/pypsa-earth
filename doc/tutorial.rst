..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _tutorial:


##########################################
Tutorial
##########################################

You may learn how to get started with PyPSA-Earth, which has a similar structure to PyPSA-EUR, by watching this video:

.. raw:: html

    <iframe width="832" height="468" src="https://www.youtube.com/embed/mAwhQnNRIvs" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

How to customise PyPSA-Earth?
======================================================

The model can be adapted to include any country, even multiple countries (e.g. Nigeria and Benin) or continents
(Currently `Africa` work as a whole continent). Several countries have not been tested yet and might not run smoothly at first:

.. code:: yaml

    countries: ["NG", "BJ"]

Likewise, the example's temporal scope can be restricted (e.g. to 7 days):

.. code:: yaml

    snapshots:
        start: "2013-03-1"
        end: "2013-03-7"
        inclusive: "left" # end is not inclusive

It is also possible to allow less or more carbon-dioxide emissions, while defining the current emissions.
It is possible to model a net-zero target by setting the `co2limit` to zero:

.. code:: yaml

    electricity:
        voltages: [220., 300., 380.]
        co2limit: 1.487e+9
        co2base: 1.487e+9

PyPSA-Earth can generate a database of existing conventional powerplants through open data sources.
It is possible to select which types of powerplants to be included:

.. code:: yaml

    extendable_carriers:
        Generator: [solar, onwind, offwind-ac, offwind-dc, OCGT]
        StorageUnit: [] # battery, H2
        Store: [battery, H2]
        Link: []  # H2 pipeline


To accurately model the temporal and spatial availability of renewables such as wind and solar energy, we rely on historical weather data.
It is advisable to adapt the required range of coordinates to the selection of countries.

.. code:: yaml

    atlite:
        nprocesses: 4
        cutouts:
                africa-2013-era5-tutorial:
                    module: era5
                    dx: 0.3  # cutout resolution
                    dy: 0.3  # cutout resolution
                    # The cutout time is automatically set by the snapshot range.

It is also possible to decide which weather data source should be used to calculate potentials and capacity factor time-series for each carrier.
For example, we may want to use the ERA-5 dataset for solar and not the default SARAH-2 dataset.

.. code:: yaml

    africa-2013-era5-tutorial:
        module: era5

.. code:: yaml

    solar:
        cutout: africa-2013-era5-tutorial

Finally, it is possible to pick a solver. For instance, this tutorial uses the open-source solver glpk and does not rely
on the commercial solvers such as Gurobi or CPLEX (for which free academic licenses are available).

.. code:: yaml

    solver:
        name: glpk


Be mindful that we only noted major changes to the provided default configuration that is comprehensibly documented in :ref:`config`.
There are many more configuration options beyond what is adapted for the tutorial!

How to execute different parts of the workflow?
======================================================

Snakemake is a workflow management tool inherited by PyPSA-Earth from PyPSA-Eur.
Snakemake decomposes a large software process into a set of subtasks, or ’rules’, that are automatically chained to obtain the desired output.

.. note::

  ``Snakemake``, which is one of the major dependencies, will be automatically installed in the environment pypsa-earth, thereby there is no need to install it manually.

The snakemake included in the conda environment pypsa-earth can be used to execute any custom rule with the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake < your custom rule >  

Starting with essential usability features, the implemented PyPSA-Earth `Snakemake procedure <https://github.com/pypsa-meets-earth/pypsa-earth/blob/main/Snakefile>`_ that allows to flexibly execute the entire workflow with various options without writing a single line of python code. For instance, you can model the world energy system or any subset of countries only using the required data. Wildcards, which are special generic keys that can assume multiple values depending on the configuration options, help to execute large workflows with parameter sweeps and various options.

You can execute some parts of the workflow in case you are interested in some specific it's parts.
E.g. power grid topology may be extracted and cleaned with the following command which refers to the script name: 

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean_osm_data

Solar profile for the requested area may be calculated using the output name:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 resources/renewable_profiles/profile_solar.nc 


How to use PyPSA-Earth for your energy problem?
======================================================


The solved networks can be analysed just like any other PyPSA network (e.g. in Jupyter Notebooks).

.. code:: python

    import pypsa
    network = pypsa.Network("results/networks/elec_s_6_ec_lcopt_Co2L-4H.nc")    

We prepared an example notebook of how you could explore the above generated tutorial network in `sample-network-analysis.ipynb`.
You can find the notebook in our `notebooks repository <https://github.com/pypsa-meets-earth/documentation/tree/main/notebooks>`_.
For further inspiration how how to explore networks and also revealing the capabilities of PYPSA you may want to see the
`examples section in the PyPSA framework documentation <https://pypsa.readthedocs.io/en/latest/examples-basic.html>`_.


After playing with the tutorial model and before playing with different fucntions,
it's important to clean-up data in your model folder before to proceed further to avoid data conflicts.
You may use the `clean` rule for making so:

.. code:: bash

    .../pypsa-earth (pypsa-earth) % snakemake -j 1 clean

Generally, it's a good idea to repeat the cleaning procedure every time when the underlying data are changed.

.. note::

  This tutorial only covers Nigeria. To make the workflow run on other regions you need to use the ``config.default.yaml`` as ``config.yaml``.
  To use the model in and outside Africa, you should also read
  `How to create a model for you region of interest with PyPSA-Earth? <https://github.com/pypsa-meets-earth/pypsa-earth/discussions/505>`_
