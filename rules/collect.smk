# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later


rule solve_all_networks:
    input:
        expand(rules.solve_network.output.network, **config["scenario"]),


rule plot_all_p_nom:
    input:
        expand(
            rules.plot_network.output.only_map,
            **config["scenario"],
            attr=["p_nom"],
            ext=["png", "pdf"],
        ),


rule make_all_summaries:
    input:
        expand(
            rules.make_summary.output.summary,
            **config["scenario"],
            country=["all"] + config["countries"],
        ),


rule plot_all_summaries:
    input:
        expand(
            rules.plot_summary.output.plot,
            summary=["energy", "costs"],
            **config["scenario"],
            country=["all"] + config["countries"],
            ext=["png", "pdf"],
        ),


rule prepare_sector_networks:
    input:
        expand(
            rules.prepare_sector_network.output.network,
            **config["scenario"],
            **config["costs"],
        ),


rule solve_sector_networks:
    input:
        expand(sector_postnetwork, **config["scenario"], **config["costs"]),
