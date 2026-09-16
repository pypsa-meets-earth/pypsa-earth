<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: CC-BY-4.0
-->

# Description of datasets used by workflow

## OpenStreetMap

**Output:** `resources/osm/raw/all_raw_cables.geojson, resources/osm/raw/all_raw_generators.geojson, resources/osm/raw/all_raw_generators.csv, resources/osm/raw/all_raw_lines.geojson, resources/osm/raw/all_raw_substations.geojson`

OpenStreetMap database used as a primary source of power-infrastructure data on substations, lines, cables, generators. OSM data are fetched for requested countries via the earth-osm package (github.com/pypsa-meets-earth/earth-osm) which downloads regional .osm.pbf extracts from Geofabrik and filters them to power-tagged features. For each country, url are resolved dynamically from Geofabrik's index. Cite: © OpenStreetMap contributors, data extracted via Geofabrik (https://download.geofabrik.de/).

## ERA5 hourly reanalysis on single levels from 1940 to present, Copernicus/ECMWF

**Output:** `cutouts/{cutout}.nc`

ERA5 hourly data on single levels from 1940 to present. Atmospheric reanalysis produced by European Centre for Medium-Range Weather Forecasts (ECMWF). Retrieved from Copernicus Climate Change Service (C3S) using atlite package to build weather cutouts for the wind, influx, temperature and runoff features. Cite: Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horanyi, A., Munoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thepaut, J-N. (2023): ERA5 hourly data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: https://doi.org/10.24381/cds.adbb2d47

## Surface Radiation Data Set - Heliosat (SARAH), Edition 3, EUMETSAT CM SAF

**Output:** `cutouts/{cutout}.nc`

SARAH-3 (Surface Radiation Data Set - Heliosat, Edition 3): satellite-based climate data record of surface solar irradiance (SIS/SID) at 0.05 deg x 0.05 deg resolution, with 30-minute/daily/monthly time steps, covering Europe, Africa, parts of South America and adjacent oceans from 1983 to near-present. Produced by EUMETSAT's Satellite Application Facility on Climate Monitoring (CM SAF). The raw files must be obtained manually via the CM SAF Web User Interface or the EUMETSAT Data Store, atlite is used to build a cutout. Cite: Pfeifroth, U., Kothe, S., Drucke, J., Trentmann, J., Schroder, M., Selbach, N., Hollmann, R. (2023): Surface Radiation Data Set - Heliosat (SARAH) - Edition 3. Satellite Application Facility on Climate Monitoring (CM SAF). DOI: https://doi.org/10.5676/EUM_SAF_CM/SARAH/V003

## Global administrative boundaries (GADM), by country

**Output:** `data/gadm/{GADM_filename}/{GADM_filename}.gpkg`

GADM v4.1 database of global administrative boundaries (country, province/state, and county/district levels), used to delineate country and sub-national regions for network clustering; the current, actively maintained GADM release.

## Global administrative boundaries (GADM v3.6)

**Output:** `data/raw/gadm/gadm36_{ISO3}/gadm36_{ISO3}.gpkg`

GADM v3.6 administrative boundaries, fetched from the older biogeo.ucdavis.edu endpoint; still needed by some sector-coupled model outputs even though the function that fetches it directly (_helpers.download_GADM) is not called anywhere in the current codebase.

## Constrained 100m-resolution global population raster (WorldPop) built using Maxar, University of Southampton

**Output:** `data/WorldPop/{WorldPop_filename}`

WorldPop constrained population raster at ~100m resolution, restricting population counts to built-up areas using Ecopia/Maxar building-footprint settlement layers rather than the BSGM growth model; produced by the WorldPop research group, University of Southampton.

## Constrained 100m-resolution global population raster (WorldPop), University of Southampton

**Output:** `data/WorldPop/{WorldPop_filename}`

WorldPop constrained population raster at ~100m resolution, using the Built-Settlement Growth Model (BSGM) to restrict population counts to built-up areas; the standard constrained population product used by default.

## Unconstrained global population raster 100m/3 arc-sec resolution (WorldPop), University of Southampton

**Output:** `data/WorldPop/{WorldPop_filename}`

WorldPop Global Project (wpgp) unconstrained population raster at ~100m/3 arc-second resolution, spreading population across whole administrative units rather than restricting it to built settlement; resolved dynamically per country via the WorldPop REST API, which returns a second, temporary tif download URL.

## Global Dataset of Forecasted Hourly Electricity Demand from 2000 to 2024 Using DemandCast. E. Antonini, P. Goli, K. Steijn

**Output:** `data/demand/forecasts_on_historical_period.parquet`

DemandCast: machine-learning-generated hourly electricity demand forecasts (2000-2024) for 184 countries, combining historical demand data, weather variables and socioeconomic indicators. E. Antonini, P. Goli and K. Steijn.

## GEGIS hourly electricity demand projections by SSP scenario and weather year.

**Output:** `data/ssp2-2.6/2030/era5_2013/{Continent}.nc`

Hourly electricity demand projections by continent, generated by the open-source GlobalEnergyGIS (GEGIS) Julia package, which combines historical demand, population/income data, SSP socioeconomic scenarios, and ERA5 weather data; the pypsa-earth bundle currently ships the SSP2-2.6 scenario, prediction year 2030, weather year 2013.

## World EEZ v11 maritime boundaries worldwide (Flanders Marine Institute / VLIZ)

**Output:** `data/eez/eez_v11.gpkg`

World EEZ v11: global maritime boundary (Exclusive Economic Zone) polygons, published by the Flanders Marine Institute (VLIZ) via marineregions.org.

## GEBCO gridded global bathymetry (2025 grid)

**Output:** `data/gebco/GEBCO_2025_sub_ice.nc`

GEBCO_2025 global gridded bathymetry: seafloor depth data compiled from many trackline sources, largely representative of deep-ocean bathymetry rather than detailed shallow-shelf depths.

## Copernicus Global Land Service 100m land cover (PROBA-V)

**Output:** `data/copernicus/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif`

Copernicus Global Land Service 100m land cover (Collection 3, PROBA-V, epoch 2019): global annual land cover classification.

## World Database on Protected Areas (WDPA), UNEP-WCMC/IUCN

**Output:** `data/landcover/world_protected_areas/*`

World Database on Protected Areas (WDPA): boundaries of protected areas worldwide, maintained by UNEP-WCMC/IUCN and distributed via Protected Planet's monthly-release.

## Rasterized protected-area raster derived from WDPA data generated by PyPSA-Earth. D. Fioriti

**Output:** `data/natura/natura.tiff`

Rasterized protected-area raster derived from WDPA vector data, generated inside PyPSA-Earth by build_natura_raster.py.

## HydroBASINS global watershed boundaries (HydroSHEDS)

**Output:** `data/hydrobasins/hybas_world.shp`

HydroBASINS (HydroSHEDS project): global watershed boundaries and sub-basin delineations at a standard level, fetched per continental-region suffix from data.hydrosheds.org.

## EDGAR gridded fossil CO2 emissions dataset (v6.0), annual gridmaps for 1970-2018

**Output:** `data/co2_emissions/v60_CO2_excl_short-cycle_org_C_1970_2018.xls`

EDGAR (Emissions Database for Global Atmospheric Research) v6.0 gridded inventory of fossil CO2 emissions, excluding short-cycle organic carbon (i.e. biomass burning and LULUCF sources); provided as annual global gridmaps covering 1970-2018, produced by the EU Joint Research Centre.

## IRENA renewable energy statistics

**Output:** `data/IRENA_Statistics_Extract_2025H2.xlsx`

IRENA Renewable Energy Statistics: installed capacity, generation, heat production and related indicators, collected via IRENA's annual member questionnaire and desk research; used in build_renewable_profiles.py to derive hydropower generation potentials.

## Global building footprints, Microsoft

**Output:** `data/global_buildings/{country}_global_buildings_raw.parquet`

Microsoft's Global ML Building Footprints: worldwide building footprint polygons derived from satellite/aerial imagery via machine learning.

## Per-quadkey-partitioned building footprint, Microsoft

**Output:** `data/global_buildings/{country}_global_buildings_raw.parquet`

Per-quadkey-partitioned GeoJSONL building footprint files from Microsoft's Global ML Building Footprints dataset.

## Per-country energy balance dataset, UN

**Output:** `resources/energy_totals_base.csv`

UN Statistics Division per-country energy balance zip files (production, trade, transformation and consumption, in TJ).

## Urban population as a percentage, UNCTAD

**Output:** `resources/urban_percent.csv`

UNCTADstat 'Urban population as percentage of total population' indicator, by country. The dataset is downloaded as a .7z file and contains urban percent for most countries from 1950 and predictions until 2050.

## Location and type of airports worldwide

**Output:** `resources/airports.csv`

Global inventory of airports and airfields from the OurAirports open dataset, including active and closed facilities, heliports, and seaplane bases, with coordinates and basic attributes (name, type, ICAO/IATA codes, elevation). The dataset contains 74844 airports.

## Location and length of airport runways worldwide

**Output:** `resources/airports.csv`

Runway-level detail from the OurAirports open dataset: length, latitude/longitude and heading for each end of every landing surface at the airports listed in the airports dataset.

## Location and technology data of iron and steel plants worldwide, GEM Global Steel Plant Tracker

**Output:** `resources/industrial_database.csv`

Global Energy Monitor's Global Steel Plant Tracker (GSPT): location, operating status (operating/proposed/retired), capacity, and production process (blast furnace, EAF, DRI, etc.) of iron and steel plants worldwide with capacity of 0.5 Mtpa or more.

## Location and technology data of gas pipelines worldwide, GEM Global Gas Infrastructure Tracker

**Output:** `resources/gas_networks/gas_network_elec_s{simpl}_{clusters}.csv`

Global Energy Monitor's Global Gas Infrastructure Tracker (GGIT) gas pipelines dataset: location, length, diameter, capacity, and status of gas transmission pipelines and pipeline projects worldwide. The dataset contains 3144 pipelines.

## SciGRID-gas IGGIELGN dataset: European gas transmission network model

**Output:** `data/gas_network/scigrid-gas/data/IGGIELGN_PipeSegments.geojson`

SciGRID-gas IGGIELGN dataset: a European gas transmission network model combining pipelines, storages, production sites, LNG terminals and interconnection points with capacity, pressure, and diameter attributes. Published on Zenodo by the SciGRID-gas project.

## Location and daily throughput capacity of oil refineries worldwide

**Output:** `resources/industrial_database.csv`

Global Oil Refinery Complex and Daily Capacity: an Esri ArcGIS-hosted point layer with the location and daily throughput capacity of oil refineries worldwide. The dataset contains 536 global Oil refineries.

## Reverse geocoding of refinery coordinates to country codes via OpenStreetMap's Nominatim API

**Output:** `resources/industrial_database.csv`

Reverse geocoding of refinery point coordinates that are missing a country code, resolved to ISO country codes via OpenStreetMap's Nominatim service. Retriven using the default endpoint of the geopy Nominatim client.

## Number of registered motor vehicles, WHO

**Output:** `resources/transport_data.csv`

WHO Global Health Observatory (GHO) 'Registered vehicles' indicator (RS_194, Road Safety theme): number of registered motor vehicles by country. Fetched via the GHO OData API. A few countries are missing from this list (e.g. South Africa, Algeria).

## Number of road motor vehicles per capita

**Output:** `resources/transport_data.csv`

Wikipedia's 'List of countries and territories by motor vehicles per capita': road motor vehicles per 1,000 inhabitants by country/territory, used as a fallback/completion source for countries missing from the WHO vehicles dataset.

## CO2 emissions from transport, World Bank

**Output:** `resources/transport_data.csv`

World Bank indicator EN.CO2.TRAN.ZS ('CO2 emissions from transport, % of total fuel combustion'), used to estimate average land-transport fuel efficiency; the live World Bank API was discontinued in October 2024, so the code fetches an archived Wayback Machine snapshot from 2024-05-21 instead.

## Sea port worldwide, NGA World Port Index

**Output:** `resources/ports.csv`

NGA World Port Index (Publication 150): location, harbor size/type, and facilities of ports, shipping terminals and oil terminals worldwide, published by the US National Geospatial-Intelligence Agency. The dataset is updated monthly and contains 3711 ports.

## Pulp and Paper Mill Database for Latin America, CGFI Spatial Finance Initiative

**Output:** `data/industry/SFI_ALD_Pulp_Paper_Sample_LatAm_Jan_2023.xlsx`

The Spatial Finance Initiative Global Pulp and Paper Mill Database provides information on pulp and paper production facilities around the world. The database contains 3,403 facilities with information about location, operating status, plant type, product type, capacity, fuel, certification status and ownership where available. Dropping the null capacities reduces the dataframe from 3000+  rows to 1672 rows.

## Global Database of Cement Production Assets, CGFI Spatial Finance Initiative

**Output:** `data/industry/SFI-Global-Cement-Database-July-2021.xlsx`

The Spatial Finance Initiative Global Database of Cement Production Assets provides information on cement production facilities worldwide. The database contains 3,117 cement plants with exact geolocation, covering both integrated clinker-producing plants and independent grinding facilities, with ownership, production type, capacity and startup year where available. Cite: McCarten, M., Bayaraa, M., Caldecott, B., Christiaen, C., Foster, P., Hickey, C., Kampmann, D., Layman, C., Rossi, C., Scott, K., Tang, K., Tkachenko, N., and Yoken, D. 2021. Global Database of Cement Production Assets.
