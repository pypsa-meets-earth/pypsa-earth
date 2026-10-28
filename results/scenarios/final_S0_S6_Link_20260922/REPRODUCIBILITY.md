# Final PyPSA-Earth thesis reproducibility record

## Scope

This directory contains the frozen final national PyPSA-Earth
scenario networks S0-S6 used for the bachelor-thesis results.

The numerical scenario results must not be regenerated or modified
after this freeze unless a documented model error is identified.

## Model framework

- Model: PyPSA-Earth
- PyPSA-Earth basis: v0.8.0 workflow
- Study region: Kazakhstan
- Planning year: 2045
- Meteorological year: 2013
- Temporal resolution: 8760 hourly snapshots
- Spatial resolution: 10 electricity clusters
- Base electricity demand: 187 TWh/a

## Final scenario definitions

- S0: reference electricity system
- S1: zero-direct-CO2 electricity system
- S2: reference + 1 GW flexible BTC mining
- S3: reference + 1 GW PEM electrolysis, 100 kt H2/a
- S4: zero-direct-CO2 + 1 GW flexible BTC mining
- S5: zero-direct-CO2 + 1 GW PEM electrolysis, 100 kt H2/a
- S6: zero-direct-CO2 + both flexible consumers

The zero-direct-CO2 scenarios impose zero represented direct fossil
CO2 emissions in the modeled electricity system. They do not
represent economy-wide or lifecycle carbon neutrality.

## Bitcoin representation

Final BTC scenarios use the validated Link representation:

electricity bus
-> BTC mining Link
-> bookkeeping BTC-service bus
-> accumulator Store

BTC electrical capacity is fixed at 1000 MW.

The legacy negative-Generator representation is retained only in
dedicated regression-validation/history scripts.

## Hydrogen representation

PEM electrolysis is represented by a 1000 MW Link with:

- annual hydrogen target: 100 kt/a
- LHV efficiency: 0.640
- annual electricity requirement: 5.2078125 TWh/a

The H2 accumulator Store is an annual bookkeeping device and does
not represent physical hydrogen storage.

## Software environment

- Python: 3.11.15
- PyPSA: 0.30.3
- pandas: 2.3.3
- NumPy: 1.26.4
- xarray: 2025.1.2
- Linopy: 0.5.8
- NetworkX: 3.6.1
- SciPy: 1.15.2
- Matplotlib: 3.11.1
- Snakemake: 7.32.4
- Gurobi: 13.0.2
- Conda: 26.5.3

Full environment exports are stored in:

- environment_pip_freeze.txt
- environment_conda_explicit.txt
- environment_conda_full.yml

## Solver audit note

At the time of the final reproducibility audit, the current machine
reported a Gurobi HostID mismatch. This prevented verification of a
fresh optimization run on that host.

This does not alter the frozen solved networks. Their integrity is
verified independently using SHA-256 checksums in SHA256SUMS.txt.

## Final validation

The final scenario manifest verifies:

- 8760 snapshots in every scenario
- 187 TWh/a base demand in every scenario
- BTC Link representation only in S2/S4/S6
- 1000 MW BTC capacity in S2/S4/S6
- no BTC component in S0/S1/S3/S5
- 100 kt H2/a in S3/S5/S6
- zero represented direct CO2 in S1/S4/S5/S6

See:

- final_scenario_manifest.csv
- SHA256SUMS.txt

## Result-processing chain

Frozen S0-S6 networks
-> analysis/extract_S0_S6_flexible_kpis.py
-> final S0-S6 KPI / pairwise / interaction / competition CSVs
-> analysis/create_S0_S6_thesis_outputs.py
-> thesis tables and figures

The principal result-processing pipeline uses only the frozen
scenario archive in this directory.
