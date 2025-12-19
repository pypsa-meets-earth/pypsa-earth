.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

.. SPDX-License-Identifier: CC-BY-4.0



# Rule `build_natura_raster`

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
        9	 [color="0.22 0.6 0.85",
            label=build_renewable_profiles];
        12	 [color="0.31 0.6 0.85",
            fillcolor=gray,
            label=build_natura_raster,
            style=filled];
        12 -> 9;
    }

|

!!! info "Source Code Reference"
    For implementation details, see the `build_natura_raster` script in the `scripts/` directory.
    :noindex:
