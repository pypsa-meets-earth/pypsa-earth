# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later


def input_make_summary(w):
    # It's mildly hacky to include the separate costs input as first entry
    if w.ll.endswith("all"):
        ll = config["scenario"]["ll"]
        if len(w.ll) == 4:
            ll = [l for l in ll if l[0] == w.ll[0]]
    else:
        ll = w.ll
    return [
        rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        )
    ] + expand(
        rules.solve_network.output.network,
        ll=ll,
        **{
            k: config["scenario"][k] if getattr(w, k) == "all" else getattr(w, k)
            for k in ["simpl", "clusters", "opts"]
        },
    )


if config["foresight"] == "overnight":
    sector_postnetwork = rules.solve_sector_network.output.network
elif config["foresight"] == "myopic":
    sector_postnetwork = rules.solve_network_myopic.output.network


rule make_summary:
    params:
        ll=config["scenario"]["ll"],
        scenario=config["scenario"],
    input:
        input_make_summary,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
    output:
        summary=directory(
            "results/"
            + RDIR
            + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}"
        ),
    log:
        "logs/"
        + RDIR
        + "make_summary/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.log",
    script:
        scripts("make_summary.py")


rule plot_summary:
    input:
        rules.make_summary.output.summary,
    output:
        plot="results/"
        + RDIR
        + "plots/summary_{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_summary/{summary}_elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{country}_{ext}.log",
    script:
        scripts("plot_summary.py")


rule plot_network:
    params:
        electricity=config["electricity"],
        plotting=config["plotting"],
    input:
        network=rules.solve_network.output.network,
        extended_country_shape=rules.build_shapes.output.extended_country_shape,
        tech_costs=rules.process_cost_data.output.costs.format(
            year=config["costs"]["year"], scope="elec"
        ),
    output:
        only_map="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}.{ext}",
        ext="results/"
        + RDIR
        + "plots/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_ext.{ext}",
    log:
        "logs/"
        + RDIR
        + "plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{attr}_{ext}.log",
    script:
        scripts("plot_network.py")


rule make_statistics:
    params:
        countries=config["countries"],
        renewable_carriers=config["electricity"]["renewable_carriers"],
        renewable=config["renewable"],
        crs=config["crs"],
        scenario=config["scenario"],
    output:
        stats="results/" + RDIR + "stats.csv",
    threads: 1
    script:
        scripts("make_statistics.py")


rule plot_sector_network:
    input:
        network=sector_postnetwork,
    output:
        map=RESDIR
        + "maps/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RESDIR
            + "benchmarks/plot_network/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
        )
    script:
        scripts("plot_network.py")


rule make_sector_summary:
    params:
        planning_horizons=config["scenario"]["planning_horizons"],
        results_dir=config["results_dir"],
        summary_dir=config["summary_dir"],
        run=run["name"],
        scenario_config=config["scenario"],
        costs_config=config["costs"],
        h2export_qty=config["export"]["h2export"],
        foresight=config["foresight"],
    input:
        networks=expand(sector_postnetwork, **config["scenario"], **config["costs"]),
        costs=rules.process_cost_data.output.costs.format(
            year="{planning_horizons}", scope="sec"
        ),
        plots=expand(
            rules.plot_sector_network.output.map,
            **config["scenario"],
            **config["costs"],
        ),
    output:
        nodal_costs=SDIR + "csvs/nodal_costs.csv",
        nodal_capacities=SDIR + "csvs/nodal_capacities.csv",
        nodal_cfs=SDIR + "csvs/nodal_cfs.csv",
        cfs=SDIR + "csvs/cfs.csv",
        costs=SDIR + "csvs/costs.csv",
        capacities=SDIR + "csvs/capacities.csv",
        curtailment=SDIR + "csvs/curtailment.csv",
        energy=SDIR + "csvs/energy.csv",
        supply=SDIR + "csvs/supply.csv",
        supply_energy=SDIR + "csvs/supply_energy.csv",
        prices=SDIR + "csvs/prices.csv",
        weighted_prices=SDIR + "csvs/weighted_prices.csv",
        market_values=SDIR + "csvs/market_values.csv",
        price_statistics=SDIR + "csvs/price_statistics.csv",
        metrics=SDIR + "csvs/metrics.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        SDIR + "benchmarks/make_summary"
    script:
        scripts("make_summary.py")


rule plot_sector_summary:
    input:
        costs=rules.make_sector_summary.output.costs,
        energy=rules.make_sector_summary.output.energy,
        balances=rules.make_sector_summary.output.supply_energy,
    output:
        costs=SDIR + "graphs/costs.pdf",
        energy=SDIR + "graphs/energy.pdf",
        balances=SDIR + "graphs/balances-energy.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        SDIR + "benchmarks/plot_summary"
    script:
        scripts("plot_summary.py")


rule prepare_db:
    params:
        tech_colors=config["plotting"]["tech_colors"],
    input:
        network=sector_postnetwork,
    output:
        db=RESDIR
        + "summaries/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}-costs-all_{planning_horizons}_{discountrate}.csv",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            RESDIR
            + "benchmarks/prepare_db/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}_{sopts}_{planning_horizons}_{discountrate}"
        )
    script:
        scripts("prepare_db.py")
