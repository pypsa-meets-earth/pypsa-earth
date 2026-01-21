# Regionalizing Cost Data

This section details the `append_cost_data.py` script, unique to PyPSA-ASEAN, which adjusts the core technology cost assumptions in your model based on data from the 8th ASEAN Energy Outlook (AEO8). It explains its function and how to use it.

The primary goal of this script is to integrate specific cost assumptions from the AEO8 into PyPSA-ASEAN's technology cost database. This ensures that your model's cost data reflects regional insights and projections.

## Input Data
The script uses several input files, primarily from the `data/AEO8-input/` directory:
*   `AEO8_Table_D15_Cost_Summary.csv`: Summarizes technology costs.
*   `AEO8_Table_D17_Declining_Factor.csv`: Contains declining factors for technology costs over time.
*   `AEO8_Table_D18_Regional_Factor.csv`: Provides a regional adjustment factor for investment costs.
*   `resources/pre_costs_{year}.csv`: This is the base cost data from PyPSA-Earth, which the script modifies.

## Key Operations

To activate the cost data integration from AEO8 and apply the declining factors, set `costs: append_cost_data` to `true` in your configuration file. You can also apply a regional factor by setting `costs: regional_factor` to `true`. An example configuration is shown below:

```yaml
costs:
  append_cost_data: true
  regional_factor: false
```

With these settings, the script performs the following operations:
1.  **Currency Conversion**: Converts all USD-denominated costs from the AEO8 input tables to EUR, as PyPSA-ASEAN (like PyPSA-Earth) uses EUR as its default currency.
2.  **Cost Adjustment**: Applies declining factors from `AEO8_Table_D17_Declining_Factor.csv` to the capital expenditure (Capex) and variable operation and maintenance (VOM) costs of various technologies.
3.  **Regional Factor Application**: If configured, a regional factor from `AEO8_Table_D18_Regional_Factor.csv` (or a user-defined value) is applied to all investment costs.
4.  **Technology Mapping**: Maps "sub-technology" names from the AEO8 data (e.g., "Coal average", "Gas Turbine") to standardized "technology" names used within PyPSA-ASEAN (e.g., "coal", "OCGT").
5.  **Data Integration**: Appends the adjusted investment, fixed operation and maintenance (FOM), and VOM costs to the existing PyPSA-ASEAN cost database. This ensures that the AEO8 assumptions take precedence for overlapping technologies.

## Output
The script generates an updated cost file, typically named `resources/costs_{year}.csv`. This file is then used by subsequent rules in the Snakefile for network building and optimization. The script also logs the changes made to the costs for transparency.
