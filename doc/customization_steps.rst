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

Pay particular attention to regional-dependent parameters such as:

* ``countries``: This parameter allows you to specify the countries to be included in your model. It's crucial for accurately representing the geographical scope of your analysis.

* ``cutouts``: These parameters refer to the climate data archive name used for calculating renewable potential. By defining appropriate cutouts, you ensure that your model is equipped with the necessary weather data for accurate simulations.

2. Generate Custom Cutout
-------------------------

The concept of a `cutout`` is central to climate data management in the PyPSA ecosystem, facilitated through the ``build_cutout`` rule.
A cutout essentially represents a spatio-temporal subset of topology and weather datasets, allowing you to focus on specific regions and timeframes relevant to your analysis.
It's important to understand that the spatial and temporal resolution of your cutout data is determined by the underlying dataset. Typically, a resolution of 30km x 30km grid and hourly data are recommended for optimal results.
While a pre-built cutout for Africa is readily available for the 2013 year, users interested in other regions or updated data can generate their own cutouts by enabling the ``build_cutout`` rule. This process involves accessing the Copernicus Climate Data Store and installing the necessary packages for data retrieval as described in :ref:`customization_copernicus`.

3. Assess `natura.tiff` Raster
--------------------------------

The ``natura.tiff`` raster file delineates areas designated as protected and reserved nature zones. These zones impose land use restrictions, influencing the calculation of renewable potential within a given geographical region.
By activating the ``natura`` option and incorporating the ``natura.tiff`` file, users can account for protected and reserved nature areas, where renewable assets cannot be installed:

.. code:: yaml

    renewable:
        onwind:
            natura: true

.. note::

    By default, the option ``natura: true`` is set for all renewable sources.

It's worth noting that a pre-built ``natura.tiff`` file is included in the dataset and dowloaded through ``retrieve_databundle_light`` rule, offering a global coverage of protected areas. This default file serves as a foundational resource for model execution.

While the default file suffices for most analyses, users seeking the latest updates and revisions can opt to enable the ``build_natura_raster`` rule:

.. code:: yaml

    enable:
        build_natura_raster: true

However, it's worth noting that generating a customized ``natura.tiff`` file with the latest data may require a significant amount of time due to processing complexities.
