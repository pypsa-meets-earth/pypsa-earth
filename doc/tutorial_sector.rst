.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _tutorial_sector:

##########################################
Tutorial: Sector-Coupled
##########################################

.. note::
    If you have not yet installed PyPSA-Earth, please refer to the :ref:`installation` section.

In this tutorial, we will show you how to run the sector-coupled model. The sector-coupled model
is a model that considers the energy system as a whole, including the electricity, heat, transport,
and industry sectors. The model is based on the PyPSA-Eur model, which is a model of the European
energy system. The sector-coupled model is a global model that can be used to model any region
of the Earth. This section explains how to run and analyze the tutorial model.


The sector-coupling code can be run as an overnight / greenfield scenario or myopic scenario.
The overnight scenario is a long-term scenario that runs for a year, while the myopic scenario
is a short-term scenario that runs for a day.


Overnight Scenarios
=============================================

Configuration
---------------------------------------------

The configuration for sector-coupled is present in the  ``config.default.yaml`` file.
This foresight parameter is set to overnight by default.

.. code:: yaml

    foresight: overnight

Documentation for all options will be added successively to :ref:`config`.

Scenarios can be defined like for electricity-only studies, but with additional wildcard options.

.. code:: yaml

    scenario:
    simpl: [""]
    ll: ["copt"]
    clusters: [10]
    opts: [Co2L-3h]
    planning_horizons: # investment years for myopic and perfect; or costs year for overnight
    - 2030
    sopts:
    - "144h"
    demand:
    - "AB"

For allowed wildcard values, refer to :ref:`wildcards`.

Execution
---------------------------------------------
To run the tutorial for the sector-coupled model, you need to activate the pypsa-earth environment.
You need to have installed PyPSA-Earth using the instructions provided in the :ref:`installation` section.
Make sure to be in the PyPSA-Earth root directory and run the following command:

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ conda activate pypsa-earth
    .../pypsa-earth (pypsa-earth) $ snakemake solve_sector_networks -j8 --configfile config.tutorial.yaml -n

This covers the retrieval of additional raw data from online resources and preprocessing data about
the transport, industry, and heating sectors as well as additional rules about geological storage
and sequestration potentials, gas infrastructure, and biomass potentials.
This triggers a workflow of multiple preceding jobs that depend on each rule's inputs and outputs:

.. graphviz::
    :class: full-width
    :align: center

    digraph snakemake_dag {
        graph[bgcolor=white, margin=0];
        node[shape=box, style=rounded, fontname=sans,                 fontsize=10, penwidth=2];
        edge[penwidth=2, color=grey];
        0[label = "solve_sector_networks", color = "0.50 0.6 0.85", style="rounded"];
        1[label = "solve_sector_network", color = "0.25 0.6 0.85", style="rounded"];
        2[label = "add_export", color = "0.08 0.6 0.85", style="rounded"];
        3[label = "prepare_ports", color = "0.06 0.6 0.85", style="rounded"];
        4[label = "retrieve_cost_data\nyear: 2030", color = "0.05 0.6 0.85", style="rounded"];
        5[label = "build_ship_profile\nh2export: 10", color = "0.34 0.6 0.85", style="rounded"];
        6[label = "prepare_sector_network", color = "0.28 0.6 0.85", style="rounded"];
        7[label = "override_respot\ndiscountrate: 0.071\nsopts: 144h", color = "0.20 0.6 0.85", style="rounded"];
        8[label = "prepare_network\nll: copt\nopts: Co2L-4H", color = "0.53 0.6 0.85", style="rounded"];
        9[label = "add_extra_components", color = "0.24 0.6 0.85", style="rounded"];
        10[label = "cluster_network\nclusters: 6", color = "0.35 0.6 0.85", style="rounded"];
        11[label = "simplify_network\nsimpl: ", color = "0.44 0.6 0.85", style="rounded"];
        12[label = "add_electricity", color = "0.49 0.6 0.85", style="rounded"];
        13[label = "build_renewable_profiles\ntechnology: onwind", color = "0.52 0.6 0.85", style="rounded"];
        14[label = "build_natura_raster", color = "0.41 0.6 0.85", style="rounded"];
        15[label = "retrieve_databundle_light", color = "0.23 0.6 0.85", style="rounded,dashed"];
        16[label = "build_shapes", color = "0.65 0.6 0.85", style="rounded"];
        17[label = "build_powerplants", color = "0.12 0.6 0.85", style="rounded"];
        18[label = "base_network", color = "0.02 0.6 0.85", style="rounded"];
        19[label = "build_osm_network", color = "0.04 0.6 0.85", style="rounded"];
        20[label = "clean_osm_data", color = "0.14 0.6 0.85", style="rounded"];
        21[label = "download_osm_data", color = "0.40 0.6 0.85", style="rounded"];
        22[label = "build_bus_regions", color = "0.66 0.6 0.85", style="rounded"];
        23[label = "build_renewable_profiles\ntechnology: offwind-ac", color = "0.52 0.6 0.85", style="rounded"];
        24[label = "build_renewable_profiles\ntechnology: offwind-dc", color = "0.52 0.6 0.85", style="rounded"];
        25[label = "build_renewable_profiles\ntechnology: solar", color = "0.52 0.6 0.85", style="rounded"];
        26[label = "build_renewable_profiles\ntechnology: hydro", color = "0.52 0.6 0.85", style="rounded"];
        27[label = "build_demand_profiles", color = "0.39 0.6 0.85", style="rounded"];
        28[label = "prepare_energy_totals\ndemand: AB\nplanning_horizons: 2030", color = "0.29 0.6 0.85", style="rounded"];
        29[label = "build_base_energy_totals", color = "0.10 0.6 0.85", style="rounded"];
        30[label = "prepare_heat_data", color = "0.26 0.6 0.85", style="rounded"];
        31[label = "build_clustered_population_layouts", color = "0.17 0.6 0.85", style="rounded"];
        32[label = "build_population_layouts\nplanning_horizons: 2030", color = "0.11 0.6 0.85", style="rounded"];
        33[label = "prepare_urban_percent", color = "0.13 0.6 0.85", style="rounded"];
        34[label = "build_temperature_profiles", color = "0.60 0.6 0.85", style="rounded"];
        35[label = "build_cop_profiles", color = "0.47 0.6 0.85", style="rounded"];
        36[label = "build_solar_thermal_profiles", color = "0.36 0.6 0.85", style="rounded"];
        37[label = "build_heat_demand", color = "0.38 0.6 0.85", style="rounded"];
        38[label = "prepare_transport_data", color = "0.46 0.6 0.85", style="rounded"];
        39[label = "prepare_transport_data_input", color = "0.31 0.6 0.85", style="rounded"];
        40[label = "build_industry_demand", color = "0.62 0.6 0.85", style="rounded"];
        41[label = "build_industrial_distribution_key", color = "0.19 0.6 0.85", style="rounded"];
        42[label = "build_industrial_database", color = "0.01 0.6 0.85", style="rounded,dashed"];
        43[label = "build_base_industry_totals\ndemand: AB\nplanning_horizons: 2030", color = "0.16 0.6 0.85", style="rounded"];
        44[label = "prepare_airports", color = "0.33 0.6 0.85", style="rounded"];
        45[label = "prepare_gas_network", color = "0.00 0.6 0.85", style="rounded"];
        46[label = "copy_config", color = "0.48 0.6 0.85", style="rounded"];
        1 -> 0
        2 -> 1
        4 -> 1
        46 -> 1
        3 -> 2
        4 -> 2
        5 -> 2
        6 -> 2
        10 -> 2
        7 -> 6
        4 -> 6
        30 -> 6
        38 -> 6
        31 -> 6
        40 -> 6
        28 -> 6
        44 -> 6
        3 -> 6
        10 -> 6
        45 -> 6
        8 -> 7
        28 -> 7
        9 -> 8
        4 -> 8
        10 -> 9
        4 -> 9
        11 -> 10
        16 -> 10
        4 -> 10
        12 -> 11
        4 -> 11
        22 -> 11
        16 -> 11
        13 -> 12
        23 -> 12
        24 -> 12
        25 -> 12
        26 -> 12
        18 -> 12
        4 -> 12
        17 -> 12
        16 -> 12
        27 -> 12
        14 -> 13
        15 -> 13
        16 -> 13
        17 -> 13
        22 -> 13
        15 -> 14
        15 -> 16
        18 -> 17
        20 -> 17
        16 -> 17
        19 -> 18
        16 -> 18
        20 -> 19
        16 -> 19
        21 -> 20
        16 -> 20
        16 -> 22
        18 -> 22
        14 -> 23
        15 -> 23
        16 -> 23
        17 -> 23
        22 -> 23
        14 -> 24
        15 -> 24
        16 -> 24
        17 -> 24
        22 -> 24
        14 -> 25
        15 -> 25
        16 -> 25
        17 -> 25
        22 -> 25
        14 -> 26
        15 -> 26
        16 -> 26
        17 -> 26
        22 -> 26
        18 -> 27
        22 -> 27
        15 -> 27
        16 -> 27
        29 -> 28
        10 -> 30
        28 -> 30
        31 -> 30
        34 -> 30
        35 -> 30
        36 -> 30
        37 -> 30
        32 -> 31
        10 -> 31
        15 -> 31
        16 -> 32
        33 -> 32
        15 -> 32
        32 -> 34
        10 -> 34
        15 -> 34
        34 -> 35
        32 -> 36
        10 -> 36
        15 -> 36
        32 -> 37
        10 -> 37
        15 -> 37
        10 -> 38
        28 -> 38
        39 -> 38
        31 -> 38
        34 -> 38
        41 -> 40
        43 -> 40
        42 -> 40
        4 -> 40
        10 -> 41
        31 -> 41
        42 -> 41
        29 -> 43
        10 -> 45
    }


In the terminal, this will show up as a list of jobs to be run:

.. code:: console

    Building DAG of jobs...
    Job stats:
    job                                   count
    ----------------------------------  -------
    add_electricity                           1
    add_export                                1
    add_extra_components                      1
    base_network                              1
    build_base_energy_totals                  1
    build_base_industry_totals                1
    build_bus_regions                         1
    build_clustered_population_layouts        1
    build_cop_profiles                        1
    build_demand_profiles                     1
    build_heat_demand                         1
    build_industrial_distribution_key         1
    build_industry_demand                     1
    build_natura_raster                       1
    build_osm_network                         1
    build_population_layouts                  1
    build_powerplants                         1
    build_renewable_profiles                  5
    build_shapes                              1
    build_ship_profile                        1
    build_solar_thermal_profiles              1
    build_temperature_profiles                1
    clean_osm_data                            1
    cluster_network                           1
    copy_config                               1
    download_osm_data                         1
    override_respot                           1
    prepare_airports                          1
    prepare_energy_totals                     1
    prepare_gas_network                       1
    prepare_heat_data                         1
    prepare_network                           1
    prepare_ports                             1
    prepare_sector_network                    1
    prepare_transport_data                    1
    prepare_transport_data_input              1
    prepare_urban_percent                     1
    retrieve_cost_data                        1
    retrieve_databundle_light                 1
    simplify_network                          1
    solve_sector_network                      1
    solve_sector_networks                     1
    total                                    46




Myopic Foresight Scenarios
=============================================


Configuration
---------------------------------------------

The configuration to run the tutorial for the myopic foresight scenario is present
in the ``test/config.test_myopic.yaml`` file.

.. code:: yaml

    foresight: myopic

Scenarios can be defined like for electricity-only studies, but with additional
wildcard options. For the myopic foresight mode, the ``{planning_horizons}`` wildcard
defines the sequence of investment horizons.

.. code:: yaml

    scenario:
    simpl:
    - ""
    clusters:
    - 4
    planning_horizons: # investment years for myopic and perfect; or costs year for overnight
    - 2030
    ll:
    - "c1"
    opts:
    - "Co2L"
    sopts:
    - "24H"
    demand:
    - "DF"

For allowed wildcard values, refer to :ref:`wildcards`.
Documentation for all options will be added successively to :ref:`config`.

Execution
---------------------------------------------
To run the tutorial for the sector-coupled model with myopic foresight, you need to activate the
pypsa-earth environment. You need to have installed PyPSA-Earth using the instructions provided in the
:ref:`installation` section. Make sure to be in the PyPSA-Earth root directory and run the following command

.. code:: bash

    .../pypsa-earth (pypsa-earth) $ conda activate pypsa-earth
    .../pypsa-earth (pypsa-earth) $ snakemake solve_sector_networks -j8 --configfile test/config.test_myopic.yaml -n

which will result in additional jobs snakemake wants to run, which translates to the following
workflow diagram which nicely outlines how the sequential pathway optimisation with myopic
foresight is implemented in the workflow:
