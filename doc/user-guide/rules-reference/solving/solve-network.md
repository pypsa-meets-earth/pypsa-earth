.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

.. SPDX-License-Identifier: CC-BY-4.0



# Rule `solve_network`

!!! note "Workflow Diagram"
    See the complete workflow in the repository.


    digraph snakemake_dag {
        graph [bgcolor=white,
            margin=0,
            size="3,3"
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
        0	 [color="0.64 0.6 0.85",
            fillcolor=gray,
            label=solve_network,
            style=filled];
        1	 [color="0.33 0.6 0.85",
            label=prepare_network];
        1 -> 0;
    }

|

!!! info "Source Code Reference"
    For implementation details, see the `solve_network` script in the `scripts/` directory.
    :noindex:
