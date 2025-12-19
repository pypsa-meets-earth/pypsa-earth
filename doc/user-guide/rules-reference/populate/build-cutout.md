.. SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

.. SPDX-License-Identifier: CC-BY-4.0



# Rule `build_cutout`

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
        10	 [color="0.44 0.6 0.85",
            label=build_hydro_profile];
        13	 [color="0.17 0.6 0.85",
            fillcolor=gray,
            label=build_cutout,
            style=filled];
        13 -> 9;
        13 -> 10;
    }

|

!!! info "Source Code Reference"
    For implementation details, see the `build_cutout` script in the `scripts/` directory.
    :noindex:
