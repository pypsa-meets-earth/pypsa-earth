# Final Adjustment

This section details the `scripts/final_asean_adjustment.py` script, which performs several critical final modifications to the PyPSA-ASEAN model before optimization, ensuring it accurately reflects specific scenarios related to network scope, demand growth, and transmission infrastructure.

## Readjusting Existing Interconnections

Here's a breakdown of the key configuration options:

```yaml
final_adjustment:
  readjust_existing_interconnections: true
```

If `readjust_existing_interconnections` is set to `true` in the configuration, the script updates the capacity of existing cross-border transmission lines to align with the targets outlined in the ASEAN Interconnection Masterplan Study (AIMS) 2024. This ensures that the model's representation of international electricity trade infrastructure is consistent with regional plans.

> **Note**: If you want all interconnections to exist from the beginning, you can set this to `false`.

## Simplifying to an Electricity-Only Network

```yaml
final_adjustment:
  only_elec_network: true
```

When `only_elec_network` is set to `true` in the configuration, the script prunes the PyPSA network to include only electricity-related buses and carriers. This streamlines the model by removing non-electrical components, focusing the analysis on the electricity system. The specific carriers to keep are defined internally within the script and can be extended by the `electricity` section of the `config.asean.yaml` (e.g., `conventional_carriers`, `renewable_carriers`, `extendable_carriers`).

> **Note**: The current configuration of PyPSA-ASEAN does not have robust assumptions on non-electricity demands, even though the components to model them are available in the default version of PyPSA-Earth.

## Adjusting Total Electricity Demand

```yaml
final_adjustment:
  pop_forecast_path: "data/worldbank_pop_forecast.csv"
  # Electricity consumption per capita from IEA: https://www.iea.org/data-and-statistics/data-tools/energy-statistics-data-browser
  # Units are in MWh/capita
  elec_per_capita:
    BN: 10.67
    KH: 0.93
    ID: 1.45
    LA: 1.67
    MM: 0.37
    MY: 5.08
    PH: 0.93
    SG: 9.75
    TH: 2.97
    TL: 0.29
    VN: 2.59

  total_elec_demand: # Baseline Scenario https://aseanenergy.org/wp-content/uploads/2024/09/8th-ASEAN-Energy-Outlook.pdf
    2025: 1.4014e+09
    2030: 1.6588e+09
    2035: 1.9402e+09
    2040: 2.2646e+09
    2045: 2.6350e+09
    2050: 3.0363e+09
```

This section handles the readjustment of the total electricity demand of the entire model based on the baseline scenarios of the 8th ASEAN Energy Outlook, and distributes it nationally based on a combination of population forecast and current electricity per capita.

*   **`pop_forecast_path`**: This parameter, specified in the configuration, points to the CSV file containing population forecast data. The script uses this data to inform electricity demand growth. If the file is not found locally, it will attempt to download the UN World Population Prospects (WPP 2024) dataset.
*   **`elec_per_capita`**: Defined within the configuration, this dictionary provides country-specific electricity consumption per capita (in MWh/capita). The script uses these values, in conjunction with population forecasts, to scale the electricity demand across various sectors such as residential, industrial, agricultural, and rail transport.
*   **`total_elec_demand`**: An optional configuration, this allows for setting total electricity demand targets for different planning horizons. If provided, the per-capita demand forecasts are adjusted proportionally to meet these overall demand figures, aligning the model with broader energy outlooks like the 8th ASEAN Energy Outlook.

> **Note**: The default PyPSA-Earth has its own electricity demand forecast which can be set in the `load_options`. Readjusting the total electricity demand overwrites the default configuration

## Redistributing Industrial Share of Electricity Demand

```yaml
final_adjustment:
  redistribute_industry: true
  industry_share: # Current IEA industry share https://www.iea.org/data-and-statistics/data-tools/energy-statistics-data-browser
    BN: 0.30
    KH: 0.36
    ID: 0.50
    LA: 0.61
    MM: 0.21
    MY: 0.88
    PH: 0.31
    SG: 0.37
    TH: 0.42
    TL: 0.01  # Data not available
    VN: 0.51
```

If `redistribute_industry` is enabled in the configuration, the script redistributes the total electricity load between industrial and non-industrial sectors. This is based on predefined country-specific industrial shares.

*   **`industry_share`**: This configuration dictionary specifies the proportion of industrial electricity demand relative to the total electricity demand for each country. This data is used by the script to accurately reallocate electricity loads within the model according to sectoral consumption patterns.

> **Note**: This will be useful for studies that focus on the industrial demand of the region.