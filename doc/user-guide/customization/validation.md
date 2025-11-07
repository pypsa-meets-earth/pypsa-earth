


4. Model Validation


To validate the data obtained with PyPSA-Earth, we recommend to go through the procedure here detailed. An exampled of the validation procedure is available in the [Nigeria validation](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/validation_nigeria.ipynb) notebook. Public information on the power system of Nigeria are compared to those obtained from the PyPSA-Earth model.

### Simulation procedure

It may be recommended to check the following quantities in the validation:

#. Inputs used by the model:

    #. Network characteristics

    #. Substations

    #. Installed generation by type

#. Outputs of the simulation:

    #. Demand

    #. Energy mix

### Where to look for reference data

Data availability for many parts of the world is still quite limited. Usually the best sources to compare with are regional data hubs. There is also a collection of harmonized datasets curated by the international organisations. A non-exhaustive list of helpful sources:

* [World Bank](https://energydata.info/);

* International Renewable Energy Agency [IRENA](https://pxweb.irena.org/pxweb/en/IRENASTAT/);

* International Energy Agency [IEA](https://www.iea.org/data-and-statistics);

* [BP](https://www.bp.com/en/global/corporate/energy-economics/statistical-review-of-world-energy.html) Statistical Review of World Energy;

* [Ember](https://ember-climate.org/data/data-explorer/) Data Explorer.


### Advanced validation examples

The following validation notebooks are worth a look when validating your energy model:

1. A detailed [network validation](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/network_validation.ipynb).

2. Analysis of [the installed capacity](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/capacity_validation.ipynb) for the considered area.

3. Validation of [the power demand](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/demand_validation.ipynb) values and profile.

4. Validation of [hydro](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/hydro_generation_validation.ipynb), [solar and wind](https://github.com/pypsa-meets-earth/documentation/blob/main/notebooks/validation/renewable_potential_validation.ipynb) potentials.
