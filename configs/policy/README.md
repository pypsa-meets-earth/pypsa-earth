<!--
SPDX-FileCopyrightText: 2026 Open Energy Transition gGmbH
SPDX-License-Identifier: CC-BY-4.0
-->

# CO2 sector policy

Country-, sector-, and subsector-specific CO2 limits can be configured through
the file referenced by:

```yaml
co2:
  sector_policy:
    policy_file: "configs/policy/co2_sector_limits.csv"
```

The policy is applied when a `Co2L` option is active.

The CSV format is:

```text
year,country,sector,subsector,limit_tco2_per_year
```

The fields have the following semantics:

- `year`: planning horizon to which the limit applies.
- `country`: uppercase country code, for example `US`.
- `sector`: emitting sector.
- `subsector`: optional emitting subsector. Leave blank for an aggregate sector limit.
- `limit_tco2_per_year`: annual CO2 limit in tonnes. Use `inf` to explicitly disable a limit.

A blank `subsector` applies to all subsectors within the selected sector.
A populated `subsector` applies only to that exact subsector. Aggregate sector
and subsector-specific limits may coexist.

Rows for other planning horizons are ignored. If no row exists for a
country/sector combination, no additional sector-specific CO2 constraint is
applied.

Supported sector classifications include:

- `power_and_heat_generation`: `power`, `heat`, `chp`
- `agriculture`
- `residential`
- `services`
- `industry`
- `industry_process_emissions`
- `transport`: `road_transport`, `non_road_rail`, `non_road_aviation`, `non_road_shipping`
- `hydrogen`: `grid_h2`, `green_h2`, `grey_h2`, `blue_h2`

The default `co2_sector_limits.csv` is intentionally header-only.
`co2_sector_limits.example.csv` provides an example policy.

For multi-country models, emitting components must be attributable to a
country. Ambiguous direct-emitting components are not assigned arbitrarily;
the solver raises an error if such a component is covered by an active policy.
