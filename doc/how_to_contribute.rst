..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Africa authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_contribute:

##########################################
How to contribute
##########################################

Contributions are welcome, and they are greatly appreciated! 
Every little bit helps, and you always earn credits.

You can contribute on the code side in many ways:

- submit feedback,
- add new features,
- report bugs, 
- fix bugs, 
- implement a new cluster/cloud computation backend,
- write documentation

Code related. For linting, formatting and checking your code contributions
against our guidelines (e.g. we use `Black <https://github.com/psf/black>`_ as code style
and aim for `REUSE compliance <https://reuse.software/>`_),
use `pre-commit <https://pre-commit.com/index.html>`_:

1. Installation ``conda install -c conda-forge pre-commit`` or ``pip install pre-commit``
2. Usage:
    * To automatically activate ``pre-commit`` on every ``git commit``: Run ``pre-commit install``
    * To manually run it: ``pre-commit run --all``

Non code related.
If you are interested in non-code related activities then have a look at the `project structure <https://pypsa-meets-africa.readthedocs.io/en/latest/project_structure_and_credits.html>`_ & click on the 'leader' position.
A pdf file will show up, introducing you to other activities in the roam of HR, finance to outreach.


Join Us and Get Involved
========================

Any person/ group is welcome to join us. Be it research leader, researcher, undergraduate, or industry professional.
A simple way to explore this is to join our meeting. All of them are **OPEN**.

- `List of meetings and times <https://github.com/pypsa-meets-africa/pypsa-africa#get-involved>`_

- **Discord**
  
  - Chat with the community, team up on features, exchange with developers, code in voice channels
  - `Discord invitation link <https://discord.gg/AnuJBk23FU>`_


.. _code_wp:

Code work packages
====================


In the code team we are working on 6 different work packages:

- WP1. Demand modelling
- WP2. Conventional generator modelling
- WP3. RES modelling
- WP4. Land coverage constraint modelling
- WP5. Network and substation modelling
- WP6. Data creation and validation

A conference paper describing each of the work packages in more detail introduced [here](https://arxiv.org/abs/2110.10628) 


.. _required_prototype_dev:

Required prototype developments
---------------------------------

The project aims to develop a running prototype model for the African continent until end of this year.
To achieve this we differentiate in the project between essential and optional prototype developents.
You can join the work package development. Some examples of open tasks are:


.. image:: img/africa_osm_map.png
    :width: 60%
    :align: center


.. _essential_prototype:

Essential developments
------------------------

- WP1. Implement GEGIS which applies a machine learning approach based on existing electricity demand time-seriesdata, population densities and spatially resolved income data
- WP2. Support powerplant matching activities to create a rich generator capacity dataset
- WP2. Implement LISFLOOD to add hydro-timerseries
- WP3. Set up Atlite for Africa to create renewable timeseries
- WP3. Add different concentrated solar power (CSP) designs to Atlite
- WP4. Adapt current Atlite land coverage constraint method for African continent. For instance, it may be forbidden to install power plants in protective areas such as national parks or to build wind plants in cities. Atlite can exclude this areas but might need to be adjusted for Africa. 
- WP5. Support the creation of a network topology by applying `various methods <https://github.com/pypsa-meets-africa/pypsa-africa/discussions/15>`_
- WP6. Contribute to the AI satellite image detection for energy asset recognition such `applied for detecting HV lines, substations <https://github.com/pypsa-meets-africa/pypsa-africa/discussions/13>`_ and power plants


.. _optional_prototype:

Optional developments
----------------------


Developments which can

- WP1. Improve and validate GEGIS in different country context
- WP1. Investigate how different demand timeseries could include the state of energy access
- WP2. Improve and validate LISFLOOD in different country context
- WP3. Investigate how and in what quality existing renewable capacities are included
- WP3. Add marine energy to Atlite
- WP4. Validate and extend Atlite capabilities
- WP5. Develop a heuristic to investigate if new east-west or north-west interconnectors within Africa are viable
- WP5. Features that help decision-making on the viability of **'off-grid/mini-grid vs on-grid'**
- WP6. Improve and extend AI satellite image recognition methods
- WP6. Add overall more data and validate datasets

- Connect PyPSA-Africa with PyPSA-Eur-Sec. (Likely after the prototype)


.. _example_studies:

Example case studies 
=====================

Below we list some studies that could be performed after our developments:

- **Long-term capacity expansion planning.** Explore long-term capacity expansion with different renewable energy deployments and different network constraints e.g. business-as-usual, least-cost, RE sub-optimally deployed in other areas/zones to assist just transition
- **Interconnectivity study.** Analysis on improved interconnectivity between African nations or improved interconnectivity between pools.
- **Energy storage study.** Value of short-duration vs long-duration storage in any country that is most appropriate. Could be interesting in any country where high variable renewable energy penetration may already be or is becoming part of the future energy mix.
- **Hydrogen economy.** Potentials of establishing a hydrogen economy in a future energy system. 
- **Energy access.** The impact on changing demand in Africa. Connecting islanded grids to the energy system - a cost and benefit analysis.
- **On-grid vs off-grid study.** Sometimes it could make sense to keep networks isolated or in mini-grid solutions. But when is this the case? Our tool can help to identify regions that are worth keeping isolated. 
- ...

After linking PyPSA-Africa with PyPSA-Eur/PyPSA-Eur-Sec:

- **Intercontinental energy planning study.** The value of collaboration between the EU and the African energy system.
- **Sector coupling.** The benefits of sector coupling (electricity, gas, heat, transport, cooling) in Africa.
- **Electric Vehicles.** Opportunities and pathways to integrate electric vehicles in Africa.
- ...
