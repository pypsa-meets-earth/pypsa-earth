# Configuration

This chapter compares the default PyPSA-Earth configuration to the PyPSA-ASEAN configuration. 

- [Project Configuration](project-config.md)
- [Scenarios Configuration](scenarios-config.md)

## Methodology

PyPSA-ASEAN is a regional adaptation of the open-source PyPSA-Earth framework [^15]. It builds on PyPSA-Earth’s preprocessing routines and optimisation capabilities, while leveraging ASEAN-specific data on power plants, renewable resources, transmission infrastructure, and electricity demand. Inclusion of customized data, on top of the PyPSA-Earth compatible databases, allowed the model to represent the region’s electricity system and enable the optimization of least-cost investment and dispatch of generation, storage, and transmission assets.

![renewable-pot](../Images/renewable-pot.png)

In this work, input data was compiled from multiple open-source databases to ensure high resolution and consistency.

- Renewable generation potentials were estimated using Atlite, which processes Copernicus ERA5 climate data to generate hourly time series for onshore and offshore wind, solar PV, and hydropower [^16], as shown in the Figure above. 
- Rooftop solar PV potential was estimated using Obane’s methodology [^17], which combines GIS-derived rooftop areas with installation ratios from related literature. 
- Power plant data were gathered using the powerplantmatching Python package [^18], which consolidated and cross-validated datasets, primarily from Global Energy Monitor.
- Renewable energy infrastructure from IRENA [^19] were incorporated and allocated according to regional resource potential. Transmission network data are derived from OpenStreetMap [^20] and national power development plans [^21] [^22] [^23] [^24] [^25] [^26] [^27] [^28], including all 2025 infrastructure (planned or under construction) to ensure internal connectivity within each country, with the exception of Indonesia’s Supergrid, as visualised in Figure 2. 
- Technology cost assumptions combine PyPSA-Earth’s default database and data from AEO8 [^5]. 
- Electricity demand growth followed the AEO8’s Baseline Scenario [^5] with demand projections disaggregated two-fold: first, by country based on UN population projections [^29], and second, nodally within countries based on historical rasterized population and GDP data (adapted from PyPSA-Earth [^15]).

System expansion was simulated from 2025 to 2050 at five-year intervals. With this, investment decisions were determined sequentially using a myopic approach, in which capacity expansion in each period was optimised based on prevailing system conditions without knowledge of future periods. This approach reflected practical planning constraints and captured stage-wise development of generation, storage, and transmission assets.

---
### References

[^5]: ASEAN Centre for Energy (ACE) 2024 8th ASEAN Energy Outlook (AEO8) Vol. 8 ISSN 2963-539X revised November 2024 Jakarta: ASEAN Centre for Energy Available at https://aseanenergy.org/
[^15]: Parzen M, Abdel-Khalek H, Fedotova E, Mahmood M, Frysztacki M M, Hampp J, Franken L, Schumm L, Neumann F, Poli D, Kiprakis A and Fioriti D 2023 PyPSA-Earth. A new global open energy system optimization model demonstrated in Africa Applied Energy 341 121096 
[^16]: Hofmann F, Hampp J, Neumann F, Brown T and Hörsch J 2021 atlite: A Lightweight Python Package for Calculating Renewable Power Potentials and Time Series Journal of Open Source Software 6(62) 3294
[^17]: Obane H 2025 Exploring the technical potential of solar PV and wind energy system in ASEAN  The Institute of Energy Economics Japan
[^18]: Gotzens F Heinrichs H Hörsch J and Hofmann F 2019 Performing energy modelling exercises in a transparent way - The issue of data quality in power plant databases Energy Strategy Reviews 23 1–12
[^19]: IRENA 2025 Renewable energy statistics 2025 International Renewable Energy Agency
[^20]: OpenStreetMap contributors 2025 Planet dump https://www.openstreetmap.org
[^21]: Électricité du Cambodge 2018 Électricité du Cambodge Annual Report. https: //www.edc.com.kh/annually_page/annuallyReport
[^22]: Direktorat Jenderal Ketenagalistrikan 2021 RUPTLN 2021-2030. Direktorat Jenderal Ketenagalistrikan. Sept. 28, 2021. URL: https://gatrik.esdm.go.id/assets/ uploads/download_index/files/38622- ruptl- pln- 2021- 2030.pdf
[^23]: Kimura S and Ueda K 2021 Feasibility Study on the Transmission Highway in ACMECS Economic Research Institute for ASEAN and East Asia 
[^24]: Sarawakenergy. Sarawak’s Power Generation and Network https://www.sarawakenergy.com/assets/images/power-generation/map.png
[^25]: Ministry of Electric Power (MOEP) 2019. Existing Power Grid and Under Construction Projects https://moep.gov.mm/en/ignite/page/641
[^26]: National Grid Corporation of the Philippines 2020 NGCP Transmission Developement Plan 2020-2040 Consultation Draft https://www.ngcp.ph/Attachment- Uploads/TDP%202020- 2040%20Consultation%20Draft%20Volume%201%20Major%20Network%20Development_2020-02-10-17-38-50.pdf
[^27]: Global Energy Network Institue. Singapore Energy Summary https://www.geni.org/globalenergy/library/national_energy_grid/singapore/singaporeannationalelectricitygrid.shtml
[^28]: Energy Comission of Sabah 2023 Sabah Energy Roadmap and Master Plan 2040 https://ecos.gov.my/sites/default/files/uploads/downloads/202309/SABAH%20ENERGY%20ROADMAP%20AND%20MASTER%20PLAN%202040%20%28SERAMP%202040%29.pdf
[^29]: United Nations, Department of Economic and Social Affairs, Population Division 2024 World Population Prospects 2024, Online Edition https://population.un.org/wpp/
