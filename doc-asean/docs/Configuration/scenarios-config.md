# Scenarios Configuration

This section details the scenarios settings for PyPSA-ASEAN. This is the scenarios excerpt from the paper: 

Three network configurations were examined, ranging from nationally oriented development to a fully meshed regional grid, to explore the impact of regional integration.

**Existing Network**: Included only nationally planned transmission expansions. Subsea interconnections related to Indonesia’s Supergrid were not incorporated, and cross-border exchanges remained limited. This configuration represented a baseline scenario in which countries develop independently.

```yaml
baseline-existing-3H:
  transmission_projects:
    include:
      AIMS: false
      ID_SuperGrid: false
```

**ASEAN Network**: Assumed full internal grid integration within all ASEAN countries, including the implementation of Indonesia’s Supergrid [13]. Cross-border links beyond this level remain at current (2025) capacities. This scenario simulated the benefits of enhanced national-level coordination within countries while maintaining limited regional interconnection.

```yaml
baseline-asean-3H:
  transmission_projects:
    include:
      AIMS: false
      ID_SuperGrid: true # This detail is not necessary. It is set as true by default.
```

**AIMS Network**: Represented a fully interconnected regional system. All cross-border transmission projects proposed in the ASEAN Interconnection Masterplan Study (AIMS) [14], alongside Indonesia’s Supergrid, were included. This scenario captured the potential system-wide benefits of a fully meshed regional network, including optimized renewable integration and reduced system costs.

```yaml
baseline-aims-3H:
  transmission_projects:
    include:
      AIMS: true # This detail is not necessary. It is set as true by default.
      ID_SuperGrid: true # This detail is not necessary. It is set as true by default.
```

3.3 Decarbonisation Policies
Two emissions policy scenarios were applied to assess how policy ambition shapes energy transition pathways:

**Baseline**: No explicit emissions constraint was imposed, representing a continuation of current practices without additional decarbonisation efforts.

**Decarbonised**: An ASEAN-wide cap on power-sector emissions was added, tightening from 1000 MtCO₂-eq in 2025 to 100 MtCO₂-eq by 2050, corresponding to a 90% reduction. This pathway reflected a regionally coordinated decarbonisation strategy.

```yaml
decarbonize-aims-3H:
  co2_budget:
    enable: true
```
