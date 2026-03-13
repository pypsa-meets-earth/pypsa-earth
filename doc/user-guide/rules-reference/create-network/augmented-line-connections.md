<!--
SPDX-FileCopyrightText: PyPSA-Earth and PyPSA-Eur Authors
SPDX-License-Identifier: CC-BY-4.0
-->

# Rule `augmented_line_connections`

This rule ensures that every bus in the network has a minimum number of connections, reinforcing or creating missing line connections to avoid isolated grid components. It runs conditionally — only when `augmented_line_connection.add_to_snakefile` is set to `true` in the configuration.

## When is it used?

The rule is executed after `cluster_network` and before `add_extra_components`. When enabled, `cluster_network` outputs a `_pre_augmentation.nc` file instead of the final network, and this rule augments it with additional lines to produce the final `elec_s{simpl}_{clusters}.nc`.

## Configuration

The relevant configuration section is `augmented_line_connection` in `config.yaml`:

```yaml
augmented_line_connection:
  add_to_snakefile: true  # Set to true to enable
  connectivity_upgrade: 2  # Minimum number of connections per bus
  new_line_type: ["HVAC"]  # Type of new lines to add
```

## Inputs

| Input | Description |
|---|---|
| `tech_costs` | Technology cost data |
| `network` | Pre-augmentation network (`elec_s{simpl}_{clusters}_pre_augmentation.nc`) |
| `regions_onshore` | Onshore bus regions |
| `regions_offshore` | Offshore bus regions |

## Outputs

| Output | Description |
|---|---|
| `network` | Final clustered network (`elec_s{simpl}_{clusters}.nc`) |

## Script Documentation

::: scripts.augmented_line_connections
    options:
        show_root_heading: false
        show_source: false
