Custom Data in PyPSA-Meets-Earth
================================

PyPSA-Meets-Earth allows users to extend the model with **custom data** 
to better reflect local or specialized scenarios. This document guides 
users on how to extend PyPSA-Meets-Earth with custom datasets, ensuring 
smooth integration and reproducibility.

Overview
--------

Custom data can be used to **replace or supplement the default datasets** 
provided by the model. Supported types include:

- Power grids and lines
- Power plants
- Electricity demand
- Heat demand
- Industry demand and industry databases
- Transport demand
- Water costs
- Hydrogen underground storage
- Gas networks
- Export ports
- Airports
- Renewable energy sources (e.g., CSP, rooftop solar, solar PV)

.. note::
   All custom data can remain private if desired. 
   Users are not required to share their data publicly.

Configuration
-------------

The ``config.default.yaml`` file controls which custom data options are enabled. 
Each option can be set to ``true`` or ``false``. When enabled, the model will 
expect corresponding files in the specified directories.

Example Configuration::

   custom_data:
     renewables: []        
     elec_demand: false
     heat_demand: false
     industry_demand: false
     industry_database: false
     transport_demand: false
     water_costs: false
     h2_underground: false
     add_existing: false
     custom_sectors: false
     gas_network: false      
     export_ports: false     
     airports: false    

Required File Locations and Formats
-----------------------------------

+-------------------------------+----------------------------------------+
| Custom Data Type              | Required File Path                     |
+===============================+========================================+
| Gas network                   | ``resources/custom_data/pipelines.csv``|
+-------------------------------+----------------------------------------+
| Export ports                  | ``data/custom/export_ports.csv``       |
+-------------------------------+----------------------------------------+
| Airports                      | ``data/custom/airports.csv``           |
+-------------------------------+----------------------------------------+
| Powerplants                   | ``data/custom_powerplants.csv``        |
+-------------------------------+----------------------------------------+
| Demand/Renewables Time Series | ``.nc`` files compatible with          |
|                               | GEGIS/atlite modules                   |
+-------------------------------+----------------------------------------+

Using Linkers
-------------

Linkers provide a simple interface for integrating custom data:

- Understand and analyze the status of existing data
- Add new data over simple interfaces
- Optionally contribute to global databases (e.g., OpenStreetMap, OpenInfraMap)

Currently, linkers support:

- Power lines
- Custom powerplants
- Demand time series
- Renewable production time series

.. note::
   Additional custom data can be added using **linkers**, which provide 
   simple interfaces for integrating new datasets without modifying core 
   model files.

Reference Data Sources
----------------------

For guidance on sourcing data, refer to the following table:

+------+----------------------------------------+--------+--------+---------+-----+
| Name | Link                                   | Sector | Global | Country | API |
+======+========================================+========+========+=========+=====+
| IEA  | https://www.iea.org/countries/         | All    | ✔      |         | ?   |
+------+----------------------------------------+--------+--------+---------+-----+
| WRI  | https://www.wri.org/                   | All    | ✔      |         | ?   |
+------+----------------------------------------+--------+--------+---------+-----+
| OECD | https://stats.oecd.org/                | All    | ✔      |         | ?   |
+------+----------------------------------------+--------+--------+---------+-----+

.. note::
   This table is continuously updated to include new global and 
   country-level datasets.

Best Practices
--------------

- Keep custom datasets in the recommended directories to avoid conflicts
- Maintain the same format as the default CSV files for seamless integration
- Use descriptive names for new CSVs or NetCDF files
- Document any assumptions or modifications made in custom data for future reproducibility

Additional Notes
----------------

- If using **GADM clustering**, ensure at least one bus per administrative region. 
  Missing buses can be added using a custom CSV created with centroids matching 
  the substation GeoJSON format.
- Private datasets do not need to be shared publicly.
- Users are encouraged to contribute improvements back to the repository 
  following contribution guidelines.

Usage Instructions
------------------

1. Enable the desired options in ``config.default.yaml``.
2. Place required custom CSV/NetCDF files in the specified directories.
3. Use linkers or follow documentation to integrate demand/renewable time series.
4. Run PyPSA-Meets-Earth; the model will automatically use the custom datasets.
