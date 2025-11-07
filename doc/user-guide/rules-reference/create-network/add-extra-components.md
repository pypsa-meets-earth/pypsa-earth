.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

.. SPDX-License-Identifier: CC-BY-4.0



# Rule `add_extra_components`

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
        1	 [color="0.56 0.6 0.85",
            label=prepare_network];
        2	 [color="0.47 0.6 0.85",
            fillcolor=gray,
            label=add_extra_components,
            style=filled];
        2 -> 1;
        3	 [color="0.03 0.6 0.85",
            label=cluster_network];
        3 -> 2;
    }

|

!!! info "Source Code Reference"
    For implementation details, see the `add_extra_components` script in the `scripts/` directory.
    :noindex:
