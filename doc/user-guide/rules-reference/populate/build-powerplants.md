.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

.. SPDX-License-Identifier: CC-BY-4.0



# Rule `build_powerplants`

!!! note "Workflow Diagram"
    See the complete workflow in the repository.


    digraph snakemake_dag {
        graph [bgcolor=white,
            margin=0,
            size="8,5"
        ];
        node [fontname=sans,
            fontsize=10,
            penwidth=2,
            shape=box,
            style=rounded
        ];
        edge [color=grey,
            penwidth=2
        ];
        4	 [color="0.61 0.6 0.85",
            label=add_electricity];
        6	 [color="0.17 0.6 0.85",
            label=base_network];
        7	 [color="0.58 0.6 0.85",
            fillcolor=gray,
            label=build_powerplants,
            style=filled];
        6 -> 7;
        7 -> 4;
    }

|

!!! info "Source Code Reference"
    For implementation details, see the `build_powerplants` script in the `scripts/` directory.
    :noindex:
