


Plotting and Summary


# Rule `plot_p_nom_max`


    :align: center

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
        0	 [color="0.42 0.6 0.85",
            fillcolor=gray,
            label=plot_p_nom_max,
            style=filled];
        1	 [color="0.58 0.6 0.85",
            label=cluster_network];
        1 -> 0;
    }

|


# Rule `make_summary`


    :align: center

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
        0	 [color="0.47 0.6 0.85",
            fillcolor=gray,
            label=make_summary,
            style=filled];
        1	 [color="0.11 0.6 0.85",
            label=solve_network];
        1 -> 0;
    }

|


# Rule `plot_summary`


|


# Rule `plot_network`


    :align: center

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
        0	 [color="0.00 0.6 0.85",
            fillcolor=gray,
            label=plot_network,
            style=filled];
        1	 [color="0.50 0.6 0.85",
            label=solve_network];
        1 -> 0;
    }

|



![Image](assets/images/tech-colors.png)
