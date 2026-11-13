# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Prepare PyPSA network for solving according to :ref:`opts` and :ref:`ll`, such
as.

- adding an annual **limit** of carbon-dioxide emissions,
- adding an exogenous **price** per tonne emissions of carbon-dioxide (or other kinds),
- setting an **N-1 security margin** factor for transmission line capacities,
- specifying an expansion limit on the **cost** of transmission expansion,
- specifying an expansion limit on the **volume** of transmission expansion, and
- reducing the **temporal** resolution by averaging over multiple hours
  or segmenting time series into chunks of varying lengths using ``tsam``.

Relevant Settings
-----------------

.. code:: yaml

    links:

    lines:

    co2:
        limit:
        base:
        emission_price:
        automatic_emission:
        budget:

    electricity:
        max_hours:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`co2_cf`, :ref:`electricity_cf`

Inputs
------

- ``resources/costs.csv.``: The database of cost assumptions for all included technologies for specific years from various sources; e.g. discount rate, lifetime, investment (CAPEX), fixed operation and maintenance (FOM), variable operation and maintenance (VOM), fuel costs, efficiency, carbon-dioxide intensity.
- ``networks/elec_s{simpl}_{clusters}.nc``: confer :ref:`cluster`

Outputs
-------

- ``networks/elec_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc``: Complete PyPSA network that will be handed to the ``solve_network`` rule.

Description
-----------

.. tip::
    The rule :mod:`prepare_all_networks` runs
    for all ``scenario`` s in the configuration file
    the rule :mod:`prepare_network`.
"""

import re

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    configure_logging,
    create_logger,
    sanitize_carriers,
    sanitize_locations,
)
from add_electricity import update_transmission_costs

idx = pd.IndexSlice

logger = create_logger(__name__)


def emission_extractor(
    emission_csv: str, emission_year: int | str, country_names: np.ndarray | list
) -> pd.Series:
    """
    Extracts CO2 emission values for given country codes from the pre-processed
    global emission CSV file.

    Parameters
    ----------
    emission_csv : str
        Path to the CSV file produced by build_co2_emissions. The file has
        columns country_code_a3 (three-letter ISO code), country_code_a2
        (two-letter ISO code), country_name (full country name), and year columns
        named YYYY.
    emission_year : int or str
        Year of CO2 emissions.
    country_names : numpy.ndarray
        Two letter country codes of analysed countries.

    Returns
    -------
    pd.Series
        CO2 emission values (in kt CO2) of studied countries for the given year.
    """
    emission_year = int(emission_year)
    df = pd.read_csv(emission_csv).set_index("country_code_a2")

    # track the available years
    available_years = sorted(
        df.columns.difference(
            ["country_code_a3", "country_code_a2", "country_name"]
        ).map(int)
    )
    # find the closest available year: if present, distance is 0, otherwise the closest year is used
    closest_year = min(available_years, key=lambda y: abs(y - emission_year))
    if emission_year not in available_years:
        logger.warning(
            f"Emission year {emission_year} not found in data. "
            f"Using closest available year {closest_year} instead."
        )

    country_codes = np.atleast_1d(country_names).tolist()
    available_countries = df.index.intersection(country_codes)
    emission_by_country = df.loc[available_countries, str(closest_year)]
    if emission_by_country.index.has_duplicates:
        duplicate_countries = emission_by_country.index[
            emission_by_country.index.duplicated()
        ].unique()
        logger.warning(
            "Multiple emission values found for the following countries: "
            f"{list(duplicate_countries)}. Summing duplicate entries."
        )
        emission_by_country = emission_by_country.groupby(level=0).sum()

    missing_ccs = np.setdiff1d(country_codes, available_countries)
    if missing_ccs.size:
        logger.warning(
            f"The emission value for the following countries has not been found: {missing_ccs}"
        )
    if emission_by_country.empty:
        raise ValueError(
            f"No CO2 emission data could be extracted from '{emission_csv}' for year "
            f"{emission_year} and countries {country_codes}. "
            "The automatic CO2 limit cannot be derived from an empty result. "
            "Please check the emission data file, the requested base year, or the "
            "configured countries."
        )
    return emission_by_country


def add_co2limit(n, annual_emissions, Nyears=1.0):
    n.add(
        "GlobalConstraint",
        "CO2Limit",
        carrier_attribute="co2_emissions",
        sense="<=",
        constant=annual_emissions * Nyears,
    )


def add_gaslimit(n, gaslimit, Nyears=1.0):
    sel = n.carriers.index.intersection(["OCGT", "CCGT", "CHP"])
    n.carriers.loc[sel, "gas_usage"] = 1.0

    n.add(
        "GlobalConstraint",
        "GasLimit",
        carrier_attribute="gas_usage",
        sense="<=",
        constant=gaslimit * Nyears,
    )


def add_emission_prices(n, co2_price=0.0):
    ep = (
        pd.Series({"co2_emissions": co2_price}) * n.carriers.filter(like="_emissions")
    ).sum(axis=1)
    gen_ep = n.generators.carrier.map(ep) / n.generators.efficiency
    n.generators["marginal_cost"] += gen_ep
    su_ep = n.storage_units.carrier.map(ep) / n.storage_units.efficiency_dispatch
    n.storage_units["marginal_cost"] += su_ep


def set_line_s_max_pu(n, s_max_pu):
    n.lines["s_max_pu"] = s_max_pu
    logger.info(f"N-1 security margin of lines set to {s_max_pu}")


def set_transmission_limit(n, ll_type, factor, costs, lines, links):
    links_dc_b = n.links.carrier == "DC" if not n.links.empty else pd.Series()

    _lines_s_nom = (
        np.sqrt(3)
        * n.lines.type.map(n.line_types.i_nom)
        * n.lines.num_parallel
        * n.lines.bus0.map(n.buses.v_nom)
    )
    lines_s_nom = n.lines.s_nom.where(n.lines.type == "", _lines_s_nom)

    col = "capital_cost" if ll_type == "c" else "length"
    ref = (
        lines_s_nom @ n.lines[col]
        + n.links.loc[links_dc_b, "p_nom"] @ n.links.loc[links_dc_b, col]
    )

    update_transmission_costs(n, costs)

    if factor == "opt" or float(factor) > 1.0:
        n.lines["s_nom_min"] = lines_s_nom
        n.lines["s_nom_extendable"] = True

        n.links.loc[links_dc_b, "p_nom_min"] = n.links.loc[links_dc_b, "p_nom"]
        n.links.loc[links_dc_b, "p_nom_extendable"] = True

    if ll_type == "l":
        n.lines["s_nom_max"] = n.lines["s_nom"] * float(factor)
        n.links.loc[links_dc_b, "p_nom_max"] = n.links.loc[links_dc_b, "p_nom"] * float(
            factor
        )
    elif factor != "opt" and float(factor) > 1.0:  # implicitly also ll_type != "l"
        con_type = "expansion_cost" if ll_type == "c" else "volume_expansion"
        rhs = float(factor) * ref
        n.add(
            "GlobalConstraint",
            f"l{ll_type}_limit",
            type=f"transmission_{con_type}_limit",
            sense="<=",
            constant=rhs,
            carrier_attribute="AC, DC",
        )

    set_line_nom_max(n, lines, links)

    return n


def average_every_nhours(n, offset):

    # Please note that a casefold is applied to the offset during the resampling.
    # This is because the usage of offset where the time interval is in capital
    # letters is deprecated in future versions of pandas.
    # For example 24H is deprecated. Instead, 24h is allowed.

    logger.info(f"Resampling the network to {offset}")
    m = n.copy(snapshots=[])

    snapshot_weightings = n.snapshot_weightings.resample(offset.casefold()).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                pnl[k] = df.resample(offset.casefold()).mean()

    return m


def apply_time_segmentation(n, segments, solver_name):
    logger.info(f"Aggregating time series to {segments} segments.")
    try:
        import tsam.timeseriesaggregation as tsam
    except:
        raise ModuleNotFoundError(
            "Optional dependency 'tsam' not found." "Install via 'pip install tsam'"
        )

    p_max_pu_norm = n.generators_t.p_max_pu.abs().max().clip(lower=1e-6)
    p_max_pu = n.generators_t.p_max_pu / p_max_pu_norm

    load_norm = n.loads_t.p_set.abs().max().clip(lower=1e-6)
    load = n.loads_t.p_set / load_norm

    inflow_norm = n.storage_units_t.inflow.abs().max().clip(lower=1e-6)
    inflow = n.storage_units_t.inflow / inflow_norm

    raw = pd.concat([p_max_pu, load, inflow], axis=1, sort=False)

    agg = tsam.TimeSeriesAggregation(
        raw,
        hoursPerPeriod=len(raw),
        noTypicalPeriods=1,
        noSegments=int(segments),
        segmentation=True,
        solver=solver_name,
    )

    segmented = agg.createTypicalPeriods()

    weightings = segmented.index.get_level_values("Segment Duration")
    offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
    snapshots = [n.snapshots[0] + pd.Timedelta(f"{offset}h") for offset in offsets]

    n.set_snapshots(pd.DatetimeIndex(snapshots, name="name"))
    n.snapshot_weightings = pd.Series(
        weightings, index=snapshots, name="weightings", dtype="float64"
    )

    segmented.index = snapshots
    n.generators_t.p_max_pu = segmented[n.generators_t.p_max_pu.columns] * p_max_pu_norm
    n.loads_t.p_set = segmented[n.loads_t.p_set.columns] * load_norm
    n.storage_units_t.inflow = segmented[n.storage_units_t.inflow.columns] * inflow_norm

    return n


def enforce_autarky(n, only_crossborder=False):
    if only_crossborder:
        lines_rm = n.lines.loc[
            n.lines.bus0.map(n.buses.country) != n.lines.bus1.map(n.buses.country)
        ].index
        links_rm = n.links.loc[
            n.links.bus0.map(n.buses.country) != n.links.bus1.map(n.buses.country)
        ].index
    else:
        lines_rm = n.lines.index
        links_rm = n.links.loc[n.links.carrier == "DC"].index
    n.remove("Line", lines_rm)
    n.remove("Link", links_rm)


def set_line_nom_max(n, lines, links):
    s_max, s_min = lines.get("s_nom_max"), lines.get("s_nom_max_min")
    p_max, p_min = links.get("p_nom_max"), links.get("p_nom_max_min")
    n.lines.s_nom_max = n.lines.s_nom_max.clip(lower=s_min, upper=s_max)
    n.links.p_nom_max = n.links.p_nom_max.clip(lower=p_min, upper=p_max)


def add_bitcoin_mining(
    n: pypsa.Network,
    options: dict,
) -> None:
    """
    Add the thesis Bitcoin-mining flexible electricity consumer.

    A fixed-capacity Link consumes electricity from the selected AC bus
    and supplies a dedicated Bitcoin-service bus. A negative-sign
    Generator on that bus represents the gross economic value of the
    mining service.

    Mining can operate continuously between zero and nominal capacity.
    There is no annual energy requirement and curtailed mining is not
    shifted to another hour.
    """

    bus = options["connection_bus"]
    btc = options["bitcoin"]

    p_nom = float(btc["p_nom"])
    marginal_revenue = float(btc["marginal_revenue"])

    if bus not in n.buses.index:
        raise ValueError(f"Bitcoin connection bus '{bus}' does not exist.")

    if p_nom <= 0:
        raise ValueError("Bitcoin p_nom must be greater than zero.")

    if marginal_revenue < 0:
        raise ValueError("Bitcoin marginal_revenue cannot be negative.")

    service_bus = f"{bus} bitcoin"
    link_name = f"{bus} Bitcoin mining"
    service_name = f"{bus} Bitcoin service"

    if any(
        [
            service_bus in n.buses.index,
            link_name in n.links.index,
            service_name in n.generators.index,
        ]
    ):
        raise ValueError("Bitcoin components already exist in the network.")

    logger.info(
        "Adding %.1f MW flexible Bitcoin mining at %s",
        p_nom,
        bus,
    )

    if "bitcoin" not in n.carriers.index:
        n.add("Carrier", "bitcoin")

    if "bitcoin mining" not in n.carriers.index:
        n.add("Carrier", "bitcoin mining")

    n.add(
        "Bus",
        service_bus,
        carrier="bitcoin",
        x=n.buses.at[bus, "x"],
        y=n.buses.at[bus, "y"],
    )

    n.add(
        "Link",
        link_name,
        bus0=bus,
        bus1=service_bus,
        carrier="bitcoin mining",
        p_nom=p_nom,
        p_nom_extendable=False,
        p_min_pu=0.0,
        p_max_pu=1.0,
        efficiency=1.0,
    )

    n.add(
        "Generator",
        service_name,
        bus=service_bus,
        carrier="bitcoin",
        sign=-1,
        p_nom=p_nom,
        p_nom_extendable=False,
        p_min_pu=0.0,
        p_max_pu=1.0,
        marginal_cost=-marginal_revenue,
    )


def add_bitcoin_mining_optimized(
    n: pypsa.Network,
    options: dict,
    costs: pd.DataFrame,
    scenario_opts=None,
) -> None:
    """
    Add endogenous flexible Bitcoin-mining demand.

    A negative-sign Generator directly consumes electricity
    from the existing electricity bus.

    The investment capacity and hourly consumption are
    optimized using standard PyPSA component attributes.
    """
    bus = options["connection_bus"]
    btc = options["bitcoin"]
    technology = "bitcoin mining"

    max_capacity = float(btc.get("p_nom_max", btc["p_nom"]))

    # Optional economic sensitivity, selected by scenario token.
    sensitivities = options.get("btcopt_sensitivities", {})
    tokens = scenario_opts or []

    selected = [token for token in tokens if token in sensitivities]

    sensitivity_tokens = [
        token for token in tokens if token.startswith(("E4R", "E5C", "E6CAP"))
    ]

    if len(selected) > 1 or set(selected) != set(sensitivity_tokens):
        raise ValueError(
            "BTCopt sensitivity token is missing, duplicated, "
            "or not uniquely configured."
        )

    sensitivity = sensitivities[selected[0]] if selected else {}

    max_capacity = float(sensitivity.get("p_nom_max", max_capacity))

    revenue = float(
        sensitivity.get(
            "revenue_eur_per_mwh",
            btc["marginal_revenue"],
        )
    )

    capex_multiplier = float(sensitivity.get("capex_multiplier", 1.0))

    if not np.isfinite(capex_multiplier) or capex_multiplier <= 0:
        raise ValueError("BTCopt CAPEX multiplier must be positive and finite.")

    if bus not in n.buses.index:
        raise ValueError(f"Bitcoin connection bus '{bus}' does not exist.")

    if not np.isfinite(max_capacity) or max_capacity <= 0:
        raise ValueError("BTCopt maximum capacity must be positive and finite.")

    if not np.isfinite(revenue) or revenue < 0:
        raise ValueError("BTCopt revenue must be finite and nonnegative.")

    if technology not in costs.index or "capital_cost" not in costs.columns:
        raise ValueError(
            "BTCopt requires a 'bitcoin mining' entry with "
            "'capital_cost' in the processed PyPSA cost table."
        )

    capital_cost = float(costs.at[technology, "capital_cost"]) * capex_multiplier

    if not np.isfinite(capital_cost) or capital_cost <= 0:
        raise ValueError(
            "BTCopt requires a positive, finite, validated " "annualized capital_cost."
        )

    name = f"{bus} Bitcoin mining optimized"

    if name in n.generators.index:
        raise ValueError(f"Bitcoin component '{name}' already exists.")

    if technology not in n.carriers.index:
        n.add("Carrier", technology)

    n.add(
        "Generator",
        name,
        bus=bus,
        carrier=technology,
        sign=-1,
        p_nom_extendable=True,
        p_nom_min=0.0,
        p_nom_max=max_capacity,
        p_min_pu=0.0,
        p_max_pu=1.0,
        capital_cost=capital_cost,
        marginal_cost=-revenue,
    )

    logger.info(
        "Added BTCopt at %s: capacity limit %.1f MW, "
        "capital cost %.2f EUR/MW, revenue %.2f EUR/MWh",
        bus,
        max_capacity,
        capital_cost,
        revenue,
    )


def add_hydrogen_production(
    n: pypsa.Network,
    options: dict,
    costs: pd.DataFrame | None = None,
    optimize_capacity: bool = False,
) -> None:
    """
    Add the thesis PEM hydrogen-production flexible consumer.

    A fixed-capacity PEM Link converts electricity into H2 product.
    A bookkeeping Store accumulates H2 output and enforces the exact
    annual production target. The Store is an accounting device and
    does not represent physical hydrogen storage.
    """

    bus = options["connection_bus"]
    pem = options["pem"]

    p_nom = float(pem["p_nom"])
    annual_h2_kt = float(pem["annual_h2_kt"])
    efficiency = float(pem["efficiency_lhv"])
    variable_cost = float(pem["variable_cost_eur_per_mwh_el"])
    h2_lhv = float(pem["h2_lhv_kwh_per_kg"])

    capacity_limit = p_nom
    capital_cost = 0.0

    if optimize_capacity:
        capacity_limit = float(pem.get("p_nom_max", p_nom))
        technology = "PEM electrolysis 1 GW"

        if costs is None or technology not in costs.index:
            raise ValueError(
                "PEMopt requires the validated technology "
                "'PEM electrolysis 1 GW' in the processed cost table."
            )

        capital_cost = float(costs.at[technology, "capital_cost"])
        efficiency = float(costs.at[technology, "efficiency"])

        economics = options.get("pemopt_economics", {})
        method = economics.get("stack_replacement_method")

        if method != "full_load_hour_proxy":
            raise ValueError(
                "PEMopt requires the explicit "
                "full_load_hour_proxy replacement method."
            )

        routine_vom = float(economics["routine_vom_eur_per_mwh_el"])
        stack_investment = float(costs.at[technology, "stack_investment"])
        replacement_hours = float(costs.at[technology, "stack_replacement_hours"])

        if not np.isfinite(routine_vom) or routine_vom < 0:
            raise ValueError("Invalid PEMopt routine variable O&M.")

        if not np.isfinite(stack_investment) or stack_investment <= 0:
            raise ValueError("Invalid PEMopt stack investment.")

        if not np.isfinite(replacement_hours) or replacement_hours <= 0:
            raise ValueError("Invalid PEMopt stack replacement interval.")

        # Approximate stack replacement using full-load hours.
        # This is not a discounted replacement-event model.
        variable_cost = routine_vom + stack_investment / replacement_hours

        if not np.isfinite(capital_cost) or capital_cost <= 0:
            raise ValueError(
                "PEMopt requires positive, finite annualized " "capital costs."
            )

        if not np.isfinite(capacity_limit) or capacity_limit <= 0:
            raise ValueError("PEMopt capacity limit must be positive and finite.")

    if bus not in n.buses.index:
        raise ValueError(f"PEM connection bus '{bus}' does not exist.")

    if p_nom <= 0:
        raise ValueError("PEM p_nom must be greater than zero.")

    if annual_h2_kt <= 0:
        raise ValueError("PEM annual_h2_kt must be greater than zero.")

    if not 0 < efficiency <= 1:
        raise ValueError("PEM efficiency_lhv must be in (0, 1].")

    if variable_cost < 0:
        raise ValueError("PEM variable cost cannot be negative.")

    if h2_lhv <= 0:
        raise ValueError("Hydrogen LHV must be greater than zero.")

    h2_kg = annual_h2_kt * 1_000_000.0
    target_mwh_h2 = h2_kg * h2_lhv / 1000.0
    target_mwh_el = target_mwh_h2 / efficiency

    weighted_hours = float(n.snapshot_weightings["generators"].sum())

    maximum_mwh_el = capacity_limit * weighted_hours

    if target_mwh_el > maximum_mwh_el + 1e-6:
        raise ValueError(
            "PEM annual target is infeasible: "
            f"{target_mwh_el:.3f} MWh_el required but "
            f"only {maximum_mwh_el:.3f} MWh_el available "
            "at nominal capacity."
        )

    h2_bus = "THESIS H2 product bus"
    link_name = "THESIS PEM electrolyser"
    store_name = "THESIS H2 annual product accumulator"
    carrier = "thesis H2 product"

    if any(
        [
            h2_bus in n.buses.index,
            link_name in n.links.index,
            store_name in n.stores.index,
        ]
    ):
        raise ValueError("Thesis PEM components already exist in the network.")

    logger.info(
        "Adding %.1f MW PEM at %s with %.3f kt_H2/a target",
        p_nom,
        bus,
        annual_h2_kt,
    )

    if carrier not in n.carriers.index:
        n.add("Carrier", carrier)

    n.add(
        "Bus",
        h2_bus,
        carrier=carrier,
        x=n.buses.at[bus, "x"],
        y=n.buses.at[bus, "y"],
    )

    n.add(
        "Link",
        link_name,
        bus0=bus,
        bus1=h2_bus,
        carrier=carrier,
        p_nom=0.0 if optimize_capacity else p_nom,
        p_nom_extendable=optimize_capacity,
        p_nom_max=capacity_limit if optimize_capacity else float("inf"),
        p_min_pu=0.0,
        p_max_pu=1.0,
        efficiency=efficiency,
        marginal_cost=variable_cost,
        capital_cost=capital_cost,
    )

    e_min_pu = pd.Series(
        0.0,
        index=n.snapshots,
        dtype=float,
    )
    e_max_pu = pd.Series(
        1.0,
        index=n.snapshots,
        dtype=float,
    )

    # Force exact annual H2 accumulation at the final snapshot.
    e_min_pu.iloc[-1] = 1.0
    e_max_pu.iloc[-1] = 1.0

    n.add(
        "Store",
        store_name,
        bus=h2_bus,
        carrier=carrier,
        e_nom=target_mwh_h2,
        e_nom_extendable=False,
        e_initial=0.0,
        e_cyclic=False,
        e_min_pu=e_min_pu,
        e_max_pu=e_max_pu,
        standing_loss=0.0,
        capital_cost=0.0,
        marginal_cost=0.0,
    )

    logger.info(
        "PEM annual requirement: %.3f TWh_H2/a, " "%.3f TWh_el/a",
        target_mwh_h2 / 1e6,
        target_mwh_el / 1e6,
    )


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "prepare_network",
            simpl="",
            clusters="6",
            ll="copt",
            opts="Co2L-4H",
            # configfile="test/config.sector.yaml",
        )

    configure_logging(snakemake)

    opts = snakemake.wildcards.opts.split("-")

    n = pypsa.Network(snakemake.input[0])
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    costs = pd.read_csv(snakemake.input.tech_costs, index_col=0)
    s_max_pu = snakemake.params.lines["s_max_pu"]

    set_line_s_max_pu(n, s_max_pu)

    for o in opts:
        m = re.match(r"^\d+h$", o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break

    for o in opts:
        m = re.match(r"^\d+seg$", o, re.IGNORECASE)
        if m is not None:
            solver_name = snakemake.config["solving"]["solver"]["name"]
            n = apply_time_segmentation(n, m.group(0)[:-3], solver_name)
            break

    for o in opts:
        if "Co2L" in o:
            m = re.findall(r"[0-9]*\.?[0-9]+$", o)
            if snakemake.params.co2["automatic_emission"]["enable"]:
                country_names = n.buses.country.unique()
                emission_year = snakemake.params.co2["automatic_emission"]["base_year"]
                co2limit = emission_extractor(
                    snakemake.input.emissions, emission_year, country_names
                ).sum()
                if len(m) > 0:
                    co2limit = co2limit * float(m[0])
                logger.info("Setting CO2 limit according to emission base year.")
            elif len(m) > 0:
                co2limit = float(m[0]) * float(snakemake.params.co2["base"])
                logger.info("Setting CO2 limit according to wildcard value.")
            else:
                co2limit = float(snakemake.params.co2["limit"])
                logger.info("Setting CO2 limit according to config value.")
            add_co2limit(n, co2limit, Nyears)
            break

    for o in opts:
        if "CH4L" in o:
            m = re.findall(r"[0-9]*\.?[0-9]+$", o)
            if len(m) > 0:
                limit = float(m[0]) * 1e6
                add_gaslimit(n, limit, Nyears)
                logger.info("Setting gas usage limit according to wildcard value.")
            else:
                add_gaslimit(n, snakemake.params.electricity.get("gaslimit"), Nyears)
                logger.info("Setting gas usage limit according to config value.")

    for o in opts:
        oo = o.split("+")
        suptechs = map(lambda c: c.split("-", 2)[0], n.carriers.index)
        if oo[0].startswith(tuple(suptechs)):
            carrier = oo[0]
            # handles only p_nom_max as stores and lines have no potentials
            attr_lookup = {
                "p": "p_nom_max",
                "c": "capital_cost",
                "m": "marginal_cost",
            }
            attr = attr_lookup[oo[1][0]]
            factor = float(oo[1][1:])
            if carrier == "AC":  # lines do not have carrier
                n.lines[attr] *= factor
            else:
                comps = {"Generator", "Link", "StorageUnit", "Store"}
                for c in n.iterate_components(comps):
                    sel = c.df.carrier.str.contains(carrier)
                    c.df.loc[sel, attr] *= factor

    for o in opts:
        if "Ep" in o:
            m = re.findall(r"[0-9]*\.?[0-9]+$", o)
            if len(m) > 0:
                logger.info("Setting emission prices according to wildcard value.")
                add_emission_prices(n, float(m[0]))
            else:
                logger.info("Setting emission prices according to config value.")
                add_emission_prices(n, snakemake.params.co2["emission_price"])
            break

    ll_type, factor = snakemake.wildcards.ll[0], snakemake.wildcards.ll[1:]
    lines = snakemake.params.lines
    links = snakemake.params.links
    set_transmission_limit(n, ll_type, factor, costs, lines, links)

    if "ATK" in opts:
        enforce_autarky(n)
    elif "ATKc" in opts:
        enforce_autarky(n, only_crossborder=True)

    flexible_consumers = snakemake.config["electricity"].get(
        "flexible_consumers",
        {},
    )

    # Thesis flexible electricity consumers are activated through
    # the standard scenario opts wildcard.
    if "BTCopt" in opts:
        add_bitcoin_mining_optimized(n, flexible_consumers, costs, scenario_opts=opts)
    elif "BTC" in opts:
        add_bitcoin_mining(n, flexible_consumers)

    if "PEMopt" in opts:
        add_hydrogen_production(
            n,
            flexible_consumers,
            costs,
            optimize_capacity=True,
        )
    elif "PEM" in opts:
        add_hydrogen_production(n, flexible_consumers)

    sanitize_carriers(n, snakemake.config)
    sanitize_locations(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
