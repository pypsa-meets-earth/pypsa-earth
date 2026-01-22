# Project Configuration

This section details the configuration settings for PyPSA-ASEAN, primarily defined in `configs/config.asean.yaml`. These settings control various aspects of the model. Most of the configurations are available in the [PyPSA-Earth configuration website](https://pypsa-earth.readthedocs.io/en/latest/configuration.html)

## `scenario`

```yaml
foresight: myopic

scenario:
  clusters: [100]
  ll:
  - "v2.0"
  opts:
  - ""
  planning_horizons: # investment years for myopic and perfect; or costs year for overnight
  - 2025
  - 2030
  - 2035
  - 2040
  - 2045
  - 2050
  sopts:
  - "3h"
  demand:
  - "DEC"
```

The PyPSA-ASEAN scenarios are based on the following parameters:

- Foresight: Myopic, with 5-year timesteps (2025–2050).
- Resolution: 3-hourly intervals.
- Demand: DEC (more details in the sector section).

## `countries`

```yaml
countries: [BN, KH, ID, LA, MM, MY, PH, SG, TH, TL, VN]
```

The model consists of 11 countries. Adding more countries requires validation and adjustment that are outside the current scope of the model. It is possible to reduce the number of countries.

## `load_options`

```yaml
load_options:
  substitute: # Load profile not available for LA and TL, use substitute instead
    LA: KH
    TL: ID
  scale:
    LA: 1456.62 # 12.76 TWh / 8760
    TL: 47.35 # 414.76 MWh / 8760
```

The load profile for Laos and Timor Leste is not available in the default version of PyPSA-Earth. Therefore, the load profiles for these countries are substituted with those of Cambodia and Indonesia, respectively. Rough estimates of the total electricity demand for these countries were also made.

## `build_shape_options` and `subregion`

```yaml
build_shape_options:
  simplify_gadm: false

subregion: # remove 'false' if subregion are to be specified
  enable:
    simplify_network: true
    cluster_network: true
  define_by_gadm:
    ID_Java-Bali: [ID.2_1, ID.4_1, ID.7_1, ID.9_1, ID.10_1, ID.11_1, ID.33_1]
    ID_Kalimantan: [ID.12_1, ID.13_1, ID.14_1, ID.34_1, ID.35_1]
    ID_Maluku: [ID.18_1, ID.19_1]
    ID_Nusa-Tenggara: [ID.20_1, ID.21_1]
    ID_Papua: [ID.22_1, ID.23_1]
    ID_Sulawesi: [ID.6_1, ID.25_1, ID.26_1, ID.27_1, ID.28_1, ID.29_1]
    ID_Sumatra: [ID.1_1, ID.3_1, ID.5_1, ID.8_1, ID.16_1, ID.17_1, ID.24_1, ID.30_1, ID.31_1, ID.32_1]
    MY_Peninsular: [MY.1_1, MY.2_1, MY.3_1, MY.4_1, MY.6_1, MY.7_1, MY.8_1, MY.9_1, MY.10_1, MY.11_1, MY.12_1, MY.15_1, MY.16_1]
    MY_Sabah: [MY.5_1, MY.13_1]
    MY_Sarawak: [MY.14_1]
    PH_Luzon: [PH.1_1, PH.5_1, PH.7_1, PH.8_1, PH.10_1, PH.11_1, PH.12_1, PH.13_1, PH.17_1, PH.18_1, PH.19_1, PH.20_1, PH.23_1, PH.24_1, PH.33_1, PH.34_1, PH.35_1, PH.37_1, H.38_1, PH.39_1, PH.40_1, PH.45_1, PH.46_1, PH.47_1, PH.50_1, PH.55_1, PH.56_1, PH.57_1, PH.58_1, PH.59_1, PH.60_1, PH.61_1, PH.62_1, PH.63_1, PH.64_1, PH.65_1, PH.69_1, PH.76_1, PH.78_1]
    PH_Visayas: [PH.4_1, PH.6_1, PH.14_1, PH.15_1, PH.22_1, PH.25_1, PH.31_1, PH.32_1, PH.36_1, PH.43_1, PH.51_1, PH.52_1, PH.54_1, PH.66_1, PH.68_1, PH.71_1]
    PH_Mindanao: [PH.2_1, PH.3_1, PH.9_1, PH.16_1, PH.21_1, PH.26_1, PH.27_1, PH.28_1, PH.29_1, PH.30_1, PH.41_1, PH.42_1, PH.44_1, PH.48_1, PH.49_1, PH.53_1, PH.67_1, PH.70_1, PH.72_1, PH.73_1, PH.74_1, PH.75_1, PH.77_1, PH.79_1, PH.80_1, PH.81_1]
```
The subregion feature in PyPSA-Earth stems from an attempt to model Southeast Asia. The region's archipelago landscape makes clustering and simplifying the model difficult.

In PyPSA-ASEAN, at least three countries are divided into subregions:

- Indonesia: Java-Bali, Kalimantan, Maluku, Nusa Tenggara, Papua, and Sulawesi;
- Malaysia: Peninsular Malaysia, Sabah, and Sarawak
- Philippines: Luzon, Visayas, and Mindanao

## `focus_weights`

```yaml
focus_weights:
  BN: 0.0
  KH: 0.046396
  LA: 0.058770
  MM: 0.169951
  SG: 0.0
  TH: 0.131035
  TL: 0.0
  VN: 0.084013
  ID_Java-Bali: 0.035281
  ID_Kalimantan: 0.135732
  ID_Maluku: 0.0
  ID_Nusa-Tenggara: 0.016782
  ID_Papua: 0.0
  ID_Sulawesi: 0.046963
  ID_Sumatra: 0.120713
  MY_Peninsular: 0.033556
  MY_Sarawak: 0.031638
  MY_Sabah: 0.018709
  PH_Luzon: 0.023487
  PH_Visayas: 0.023487
  PH_Mindanao: 0.023487
```

If left at the default setting, the nodes will be aggregated based on the largest load. The focus weight configuration redistributes the nodes based on the land area of each country and subregion. To focus on a specific country or subregion, increase its focus weight by subtracting from the values of other locations.

## `co2_budget`

```yaml
co2_budget:
  enable: false
  override_co2opt: true
  co2base_value: 1.0e+09 # choose from: [co2limit, co2base, absolute, {float}]
  year:
    2025: 1.0
    2030: 0.82
    2035: 0.64
    2040: 0.46
    2045: 0.28
    2050: 0.1
```

if `enable` is `true`, An ASEAN-wide cap on power-sector emissions was added, tightening from 1000 MtCO₂-eq in 2025 to 100 MtCO₂-eq by 2050, corresponding to a 90% reduction. This pathway reflected a regionally coordinated decarbonisation strategy.

## `electricity` and `transmission_projects`

This configuration is explained in [ASEAN Grid Infrastructure](../Feature/asean-grid-infrastructure.md).

## `cluster_options`

```yaml
cluster_options:
  simplify_network:
    remove_stubs_across_borders: false
    p_threshold_drop_isolated: 0 # [MW] isolated buses are being discarded if bus mean power is below the specified threshold
    p_threshold_merge_isolated: 0 # [MW] isolated buses are being merged into a single isolated bus if a bus mean power is below the specified threshold
    s_threshold_fetch_isolated: 0.1 # [-] a share of the national load for merging an isolated network into a backbone network
```

To maintain as much of the original electricity assumption as possible, buses are not dropped or merged in the simplification process. However, if nodes represent less than 10% of the total country or subregion electricity demand, they are assumed to be connected to the network.

## `costs`

This configuration is explained in [Final Adjustment](../Feature/append-cost-data.md).

## `sector`

```yaml
sector:
  enable:
    heat: true
    biomass: true
    industry: true
    shipping: true
    aviation: true
    land_transport: true
    rail_transport: true
    agriculture: true
    residential: true
    services: true

  biomass_transport: false
  solid_biomass_potential: 360 # TWh/a, technical maximum in Southeast Asia (higher granular data needed). https://doi.org/10.1016/j.energy.2017.06.162
  # biogas_potential: 0.5 # TODO find ASEAN equivalent # TWh/a, Potential of whole modelled area

  hydrogen:
    network: false

  co2_network: false

  land_transport_fuel_cell_share:
    DEC_2020: 0.00
    DEC_2025: 0.00
    DEC_2030: 0.00
    DEC_2035: 0.00
    DEC_2040: 0.00
    DEC_2045: 0.00
    DEC_2050: 0.00
  land_transport_electric_share:
    DEC_2020: 0
    DEC_2025: 0.05
    DEC_2030: 0.2
    DEC_2035: 0.45
    DEC_2040: 0.7
    DEC_2045: 0.85
    DEC_2050: 1
  shipping_hydrogen_share:
    DEC_2020: 0.00
    DEC_2025: 0.00
    DEC_2030: 0.00
    DEC_2035: 0.00
    DEC_2040: 0.00
    DEC_2045: 0.00
    DEC_2050: 0.00

  solar_rooftop: # adds distribution side customer rooftop PV (only work if electricity_distribution_grid: true)
    use_building_size: true
```

In the sector configuration, all sectors are applied, though many components are removed during the final adjustment. Important assumption:

- The solid biomass potential is 360 TWh/a, based on the technical maximum (https://doi.org/10.1016/j.energy.2017.06.162). Currently, the potential is not spatially distributed.
- Hydrogen and CO₂ networks are disabled.
- It is assumed that land transport will be 100% electric by 2050.
- It uses the Global Buildings Dataset to estimate solar rooftop potential.

## `policy_config`, `export` and `solving`

```yaml
policy_config:
  hydrogen:
    temporal_matching: "no_temporal_matching"

export:
  enable: false
  h2export: [0]

solving:
  options:
    load_shedding: false
```

These are minor adjustments to disable:

- Hydrogen-related exports
- Load shedding options

## `final_adjustment`

This configuration is explained in [Final Adjustment](../Feature/final-adjustment.md).




