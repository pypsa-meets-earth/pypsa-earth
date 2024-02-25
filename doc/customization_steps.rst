.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
..
.. SPDX-License-Identifier: CC-BY-4.0

.. _customization_steps:

#######################################
3. General Guidelines for Modeling
#######################################


1. Configure Model Parameters
-----------------------------

Begin by adjusting the model inputs in the ``config.yaml`` configuration file. This file serves as the cornerstone for tailoring PyPSA-Earth to specific regional requirements.

Pay particular attention to parameters such as:

* ``countries``: This parameter allows you to specify the countries to be included in your model. It's crucial for accurately representing the geographical scope of your analysis.

* ``enable``: This section contains parameters such as `retrieve_databundle`, `retrieve_cost_data`, `download_osm_data`, `build_natura_raster`, and `build_cutouts`. The parameter `retrieve_databundle` is used to download all required technological data. Set it to `false` after initial run to prevent redownload. To build custom cutouts and natura exclusion raster, set `build_cutout` and `build_natura_raster` to `true`, respectively. More details on custom cutouts and natura raster are provided below.

* ``cutouts``: These parameters refer to the climate data archive name used for calculating renewable potential. By defining appropriate cutouts, you ensure that your model is equipped with the necessary weather data for accurate simulations. By default, pre-compiled cutouts have standard name of ``cutout-2013-era5``. Adjust it accordingly, if custom cutouts are utilized.

2. Check if custom cutout is necessary
--------------------------------------

By default, we provide pre-built cutouts for every continent for the weather year 2013. If you require a different year or a more complex region, you may opt to create your own custom cutout using ``build_cutout`` rule. This process involves accessing the Copernicus Climate Data Store and installing the necessary packages for data retrieval as described in :ref:`customization_copernicus`.

.. note::

    The concept of a ``cutout`` is central to climate data management in the PyPSA ecosystem, facilitated through the ``build_cutout`` rule.
    A cutout essentially represents a spatio-temporal subset of topology and weather datasets, allowing you to focus on specific regions and timeframes relevant to your analysis.

3. Assess `natura.tiff` raster
--------------------------------

The ``natura.tiff`` raster file delineates areas designated as protected and reserved nature zones. These zones impose land use restrictions, influencing the calculation of renewable potential within a given geographical region.
The default ``natura.tiff`` is now a CC0 data: Sosa Arango, Chrystian Camilo, 2020, "Protected areas (WDPA)", `<https://doi.org/10.7910/DVN/XIV9BL>`__, Harvard Dataverse, V1.
By activating the ``natura`` option and incorporating the ``natura.tiff`` file, users can account for protected and reserved nature areas, where renewable assets cannot be installed:

.. code:: yaml

    renewable:
        onwind:
            natura: true

.. note::

    By default, the option ``natura: true`` is set for all renewable sources.

It's worth noting that a pre-built ``natura.tiff`` file is included in the dataset and downloaded through ``retrieve_databundle_light`` rule, offering a global coverage of protected areas. This default file serves as a foundational resource for model execution.

While the default file suffices for most analyses, users seeking the latest updates and revisions can opt to enable the ``build_natura_raster`` rule:

.. code:: yaml

    enable:
        build_natura_raster: true

However, it's worth noting that generating a customized ``natura.tiff`` file with the latest data may require a significant amount of time due to processing complexities.

4. Model Validation
-------------------

To validate the data obtained with PyPSA-Earth, we recommend to go through the procedure here detailed. An exampled of the validation procedure is available in the `Nigeria validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb>`_ notebook. Public information on the power system of Nigeria are compared to those obtained from the PyPSA-Earth model.

Simulation procedure
^^^^^^^^^^^^^^^^^^^^

It may be recommended to check the following quantities in the validation:

#. Inputs used by the model:

    #. Network characteristics

    #. Substations

    #. Installed generation by type

#. Outputs of the simulation:

    #. Demand

    #. Energy mix

Where to look for reference data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data availability for many parts of the world is still quite limited. Usually the best sources to compare with are regional data hubs. There is also a collection of harmonized datasets curated by the international organisations. A non-exhaustive list of helpful sources:

* `World Bank <https://energydata.info/>`_;

* International Renewable Energy Agency `IRENA <https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/ELECCAP_2022_cycle2.px/>`_;

* International Energy Agency `IEA <https://www.iea.org/data-and-statistics>`_;

* `BP <https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html>`_ Statistical Review of World Energy;

* `Ember <https://ember-climate.org/data/data-explorer/>`_ Data Explorer.


Advanced validation examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following validation notebooks are worth a look when validating your energy model:

1. A detailed `network validation <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb>`_.

2. Analysis of `the installed capacity <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb>`_ for the considered area.

3. Validation of `the power demand <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb>`_ values and profile.

4. Validation of `hydro <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb>`_, `solar and wind <https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb>`_ potentials.
