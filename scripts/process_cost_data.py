# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Prepare and process default cost data to fit PyPSA input requirements (capital costs / fixed annualized costs, marginal costs).

The technology cost data is loaded from a CSV file, filtered by scenario and financial case,
and then processed to compute fixed annualized costs and marginal costs for each technology.
Currency conversion is applied to ensure all costs are in a common output currency, using exchange rates for a specified reference year.

The output differs dependind on the wildcard {scope}:
- For "elec", the output contains capital_cost, marginal_cost, and co2_emissions for each technology, and is used in the electricity sector workflow.
- For "sec", the output contains fixed (annualized capital cost), marginal_cost, and co2_emissions for each technology, and is used in the sector-coupling workflow.

Relevant Settings
-----------------
.. code:: yaml

    costs:
        year:
        cost_scenario:
        financial_case:
        output_currency:
        fill_values:
        default_exchange_rate:
        future_exchange_rate_strategy:
        custom_future_exchange_rate:
        investment:  # optional overwrites for specific attributes
        lifetime:
        FOM:
        VOM:
        efficiency:
        fuel:

    electricity:
        max_hours:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at :ref:`costs_cf`,
    :ref:`electricity_cf`
"""
import logging

import pandas as pd
import pypsa
from currency_converter import CurrencyConverter

currency_converter = CurrencyConverter(
    fallback_on_missing_rate=True,
    fallback_on_wrong_date=True,
)

logger = logging.getLogger(__name__)


# PYPSA-EARTH-SEC
def annuity(n: float | pd.Series, r: float | pd.Series) -> float | pd.Series:
    """
    Calculate the annuity factor for an asset with lifetime n years and.

    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6

    parameters
    ----------
    n : int or pd.Series
        Lifetime of the asset in years.
    r : float or pd.Series
        Discount rate (e.g. 0.05 for 5%).

    Returns
    -------
    float or pd.Series
        Annuity factor, i.e. the ratio of fixed annualized cost to initial investment.
    """

    if isinstance(r, pd.Series):
        return pd.Series(1 / n, index=r.index).where(
            r == 0, r / (1.0 - 1.0 / (1.0 + r) ** n)
        )
    elif r > 0:
        return r / (1.0 - 1.0 / (1.0 + r) ** n)
    else:
        return 1 / n


# Single source for the currency reference year (aligned with `technology-data` output files / PyPSA-Earth input cost files).
# Change this value to update the reference year everywhere.
TECH_DATA_REFERENCE_YEAR = 2020

# Simple cache to avoid repeated computations and logging for same (currency, output_currency, year)
_currency_conversion_cache = {}


def get_yearly_currency_exchange_rate(
    initial_currency: str,
    output_currency: str,
    default_exchange_rate: float = None,
    _currency_conversion_cache: dict = None,
    future_exchange_rate_strategy: str = "reference",  # "reference", "latest", "custom"
    custom_future_exchange_rate: float = None,
) -> float:
    """
    Returns the average currency exchange rate for the global reference_year.

    Parameters
    ----------
    initial_currency : str
        Input currency (e.g. "EUR", "USD").
    output_currency : str
        Desired output currency (e.g. "USD").
    default_exchange_rate : float, optional
        Fallback value if no rate data is found.
    _currency_conversion_cache : dict, optional
        Cache for repeated calls.
    future_exchange_rate_strategy : str
        "reference" (use TECH_DATA_REFERENCE_YEAR),
        "latest" (use most recent available year),
        "custom" (use custom_future_exchange_rate).
    custom_future_exchange_rate : float, optional
        Custom exchange rate if strategy is "custom".

    Returns
    -------
    float
        Average exchange rate for the specified year and currency pair.
    """
    if _currency_conversion_cache is None:
        _currency_conversion_cache = {}

    key = (
        initial_currency,
        output_currency,
        TECH_DATA_REFERENCE_YEAR,
        future_exchange_rate_strategy,
    )
    if key in _currency_conversion_cache:
        return _currency_conversion_cache[key]

    # Handle EUR specially (no direct rates, fallback on USD dates)
    if initial_currency == "EUR":
        available_dates = sorted(currency_converter._rates["USD"].keys())
    else:
        if initial_currency not in currency_converter._rates:
            if default_exchange_rate is not None:
                return default_exchange_rate
            raise RuntimeError(f"No data for currency {initial_currency}.")
        available_dates = sorted(currency_converter._rates[initial_currency].keys())

    max_date = available_dates[-1]

    # Decide which year to use
    if future_exchange_rate_strategy == "custom":
        if custom_future_exchange_rate is None:
            raise RuntimeError("Custom strategy selected but no rate provided.")
        avg_rate = custom_future_exchange_rate
    elif future_exchange_rate_strategy == "latest":
        effective_year = max_date.year
        logger.info(
            f"Using latest available year ({effective_year}) for {initial_currency}->{output_currency}."
        )
        dates_to_use = [d for d in available_dates if d.year == effective_year]
        rates = [
            currency_converter.convert(1, initial_currency, output_currency, d)
            for d in dates_to_use
        ]
        avg_rate = sum(rates) / len(rates) if rates else default_exchange_rate
    else:  # "reference": use module-level reference_year
        effective_year = TECH_DATA_REFERENCE_YEAR
        dates_to_use = [d for d in available_dates if d.year == effective_year]
        if not dates_to_use and default_exchange_rate is not None:
            avg_rate = default_exchange_rate
        else:
            rates = [
                currency_converter.convert(1, initial_currency, output_currency, d)
                for d in dates_to_use
            ]
            avg_rate = sum(rates) / len(rates)

    _currency_conversion_cache[key] = avg_rate
    return avg_rate


def build_currency_conversion_cache(
    df: pd.DataFrame,
    output_currency: str,
    default_exchange_rate: float = None,
    future_exchange_rate_strategy: str = "reference",
    custom_future_exchange_rate: float = None,
) -> dict:
    """
    Builds a cache of exchange rates for all unique (currency, output_currency) pairs,
    always using the module-level reference_year.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a 'unit' column with currency information (e.g. "EUR/MW").
    output_currency : str
        Desired output currency (e.g. "USD").
    default_exchange_rate : float, optional
        Fallback value if no rate data is found for a currency pair.
    future_exchange_rate_strategy : str
        "reference" (use TECH_DATA_REFERENCE_YEAR),
        "latest" (use most recent available year),
        "custom" (use custom_future_exchange_rate).
    custom_future_exchange_rate : float, optional
        Custom exchange rate if strategy is "custom".

    Returns
    -------
    dict
        Cache mapping (initial_currency, output_currency, year) to exchange rate.
    """
    currency_list = currency_converter.currencies
    unique_currencies = {
        x["unit"][0:3]
        for _, x in df.iterrows()
        if isinstance(x["unit"], str) and x["unit"][0:3] in currency_list
    }

    _currency_conversion_cache = {}
    for initial_currency in unique_currencies:
        try:
            rate = get_yearly_currency_exchange_rate(
                initial_currency,
                output_currency,
                default_exchange_rate=default_exchange_rate,
                _currency_conversion_cache=_currency_conversion_cache,
                future_exchange_rate_strategy=future_exchange_rate_strategy,
                custom_future_exchange_rate=custom_future_exchange_rate,
            )
            _currency_conversion_cache[
                (initial_currency, output_currency, TECH_DATA_REFERENCE_YEAR)
            ] = rate
        except Exception as e:
            logger.warning(
                f"Failed to get rate for {initial_currency}->{output_currency}: {e}"
            )
            continue

    return _currency_conversion_cache


def apply_currency_conversion(
    cost_dataframe: pd.DataFrame, output_currency: str, cache: dict
) -> pd.DataFrame:
    """
    Applies exchange rates from the cache to convert all cost values and units.

    All rows are assumed to be in `*_reference_year` already (e.g. EUR_2020).

    Parameters
    ----------
    cost_dataframe : pd.DataFrame
        DataFrame containing a 'unit' column with currency information (e.g. "EUR/MW") and a 'value' column with cost values.
    output_currency : str
        Desired output currency (e.g. "USD").
    cache : dict
        Cache mapping (initial_currency, output_currency, year) to exchange rate, as built by `build_currency_conversion_cache`.

    Returns
    -------
    pd.DataFrame
        DataFrame with converted 'value' and 'unit' columns.
    """
    currency_list = currency_converter.currencies

    def convert_row(x):
        unit = x["unit"]
        value = x["value"]

        if not isinstance(unit, str) or "/" not in unit:
            return pd.Series([value, unit])

        currency = unit[:3]

        if currency not in currency_list:
            return pd.Series([value, unit])

        key = (currency, output_currency, TECH_DATA_REFERENCE_YEAR)
        rate = cache.get(key)
        if rate is not None:
            new_value = value * rate
            new_unit = unit.replace(currency, output_currency, 1)
            return pd.Series([new_value, new_unit])
        else:
            logger.warning(f"Missing exchange rate for {key}. Skipping conversion.")
            return pd.Series([value, unit])

    cost_dataframe[["value", "unit"]] = cost_dataframe.apply(convert_row, axis=1)
    return cost_dataframe


def calculate_cost_for_storage_units(
    costs: pd.DataFrame,
    max_hours: dict,
    storage_techs: dict,
    costs_name: str = "capital_cost",
):
    """
    Calculate the capital/fixed costs for storage as defined by storage units.

    Storage units consolidate their charging, discharging, and storage costs into a single combined value.

    Parameters
    ----------
    costs : pd.DataFrame
        DataFrame containing the cost data.
    max_hours : dict
        Dictionary of maximum hours for storage units.
    storage_techs : dict, optional
        A dictionary mapping storage carriers to their respective technology parameters.
    cost_name : str
        Name of the column in the costs DataFrame that represents the capital cost.
    """

    def costs_for_storage(
        store, link1=None, link2=None, max_hours=1.0, costs_name="capital_cost"
    ):
        capital_cost = max_hours * store[costs_name]
        if link1 is not None:
            capital_cost += link1[costs_name]
        if link2 is not None:
            capital_cost += link2[costs_name]
        return pd.Series(
            {
                costs_name: capital_cost,
                "marginal_cost": 0.0,
                "CO2 intensity": 0.0,
            }
        )

    mod_costs_i = costs.index
    missing_store = []
    for k, v in max_hours.items():
        tech = storage_techs.get(k)
        store = tech.get("store") if tech.get("store") in mod_costs_i else None
        bicharger = (
            tech.get("bicharger") if tech.get("bicharger") in mod_costs_i else None
        )
        charger = tech.get("charger") if tech.get("charger") in mod_costs_i else None
        discharger = (
            tech.get("discharger") if tech.get("discharger") in mod_costs_i else None
        )
        if bicharger:
            costs.loc[k] = costs_for_storage(
                costs.loc[store],
                costs.loc[bicharger],
                max_hours=v,
                costs_name=costs_name,
            )
        elif store:
            costs.loc[k] = costs_for_storage(
                costs.loc[store],
                costs.loc[charger] if charger else None,
                costs.loc[discharger] if discharger else None,
                max_hours=v,
                costs_name=costs_name,
            )
        else:
            missing_store += [k]

    if missing_store:
        logger.warning(
            f"No cost data on:\n - "
            + "\n - ".join(missing_store)
            + "\nPlease enable retrieve_cost_data if this storage technology is used"
        )

    return costs


def load_costs(
    tech_costs: str,
    config: dict,
    max_hours: int,
    Nyears: float | int = 1,
    storage_techs: dict = {},
) -> pd.DataFrame:
    """
    Set all asset costs and other parameters. Used only in the electricity sector workflow.

    Parameters
    ----------
    tech_costs : str
        Path to the CSV file containing technology cost data.
    config : dict
        Configuration dictionary containing cost scenario, financial case, output currency, fill values, and any overwrites for specific attributes.
    max_hours : dict
        Dictionary specifying the maximum hours of storage for storage technologies (e.g. {"battery": 4, "H2": 168}).
    Nyears : float or int, optional
        Share of years represented by the model time steps (default is 1, representing one full year).
    storage_techs : dict, optional
        A dictionary mapping storage carriers to their respective technology parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by technology with columns for capital_cost, marginal_cost, co2_emissions, and any other relevant parameters.
    """
    costs = pd.read_csv(tech_costs, index_col=["technology", "parameter"]).sort_index()

    # correct units to MW and output_currency
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.unit = costs.unit.str.replace("/kW", "/MW")
    _currency_conversion_cache = build_currency_conversion_cache(
        costs,
        output_currency=config["output_currency"],
        default_exchange_rate=config["default_exchange_rate"],
        future_exchange_rate_strategy=config.get("future_exchange_rate_strategy"),
        custom_future_exchange_rate=config.get("custom_future_rate", None),
    )
    costs = apply_currency_conversion(
        costs,
        config["output_currency"],
        _currency_conversion_cache,
    )

    # apply filter on financial_case and scenario, if they are contained in the cost dataframe
    wished_cost_scenario = config["cost_scenario"]
    wished_financial_case = config["financial_case"]
    for col in ["scenario", "financial_case"]:
        if col in costs.columns:
            costs[col] = costs[col].replace("", pd.NA)

    if "scenario" in costs.columns:
        costs = costs[
            (costs["scenario"].str.casefold() == wished_cost_scenario.casefold())
            | (costs["scenario"].isnull())
        ]

    if "financial_case" in costs.columns:
        costs = costs[
            (costs["financial_case"].str.casefold() == wished_financial_case.casefold())
            | (costs["financial_case"].isnull())
        ]

    costs = costs.value.unstack().fillna(config["fill_values"])

    for attr in ("investment", "lifetime", "FOM", "VOM", "efficiency", "fuel"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            logger.info(
                f"Overwriting {attr} of {overwrites.index} to {overwrites.values}"
            )

    costs["capital_cost"] = (
        (annuity(costs["lifetime"], costs["discount rate"]) + costs["FOM"] / 100.0)
        * costs["investment"]
        * Nyears
    )

    costs.at["OCGT", "fuel"] = costs.at["gas", "fuel"]
    costs.at["CCGT", "fuel"] = costs.at["gas", "fuel"]

    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})
    # rename because technology data & pypsa earth costs.csv use different names
    # TODO: rename the technologies in hosted tutorial data to match technology data
    costs = costs.rename(
        {
            "hydrogen storage": "hydrogen storage tank",
            "hydrogen storage tank": "hydrogen storage tank",
            "hydrogen storage tank type 1": "hydrogen storage tank",
            "hydrogen underground storage": "hydrogen storage underground",
        },
    )

    costs.at["OCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]
    costs.at["CCGT", "co2_emissions"] = costs.at["gas", "co2_emissions"]

    costs.at["solar", "capital_cost"] = (
        config["rooftop_share"] * costs.at["solar-rooftop", "capital_cost"]
        + (1 - config["rooftop_share"]) * costs.at["solar-utility", "capital_cost"]
    )
    costs.loc["csp"] = costs.loc["csp-tower"]

    costs = calculate_cost_for_storage_units(
        costs, max_hours, storage_techs, costs_name="capital_cost"
    )

    for attr in ("marginal_cost", "capital_cost"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites
            logger.info(
                f"Overwriting {attr} of {overwrites.index} to {overwrites.values}"
            )

    return costs


def prepare_costs(
    cost_file: str,
    config: dict,
    output_currency: str,
    fill_values: dict,
    Nyears: float | int = 1,
    default_exchange_rate: float = None,
    future_exchange_rate_strategy: str = "latest",
    custom_future_exchange_rate: float = None,
    max_hours: dict = {},
    storage_techs: dict = {},
) -> pd.DataFrame:
    """
    Loads and processes cost data, converting units and currency to a common format.

    Applies currency conversion, fills missing values, and computes fixed annualized costs.
    Always uses the module-level reference_year.

    Used only in the sector-coupling workflow.

    Parameters
    ----------
    cost_file : str
        Path to the CSV file containing cost data with columns for 'technology', 'parameter', 'value', 'unit', and optionally 'scenario' and 'financial_case'.
    config : dict
        Configuration dictionary containing cost scenario, financial case, and any overwrites for specific attributes.
    output_currency : str
        Desired output currency for all cost values (e.g. "USD").
    fill_values : dict
        Dictionary specifying fill values for missing cost parameters (e.g. {"FOM": 0, "VOM": 0}).
    Nyears : int or float, optional
        Share of years represented by the model time steps (default is 1, representing one full year).
    default_exchange_rate : float, optional
        Fallback exchange rate if no data is available for a currency pair.
    future_exchange_rate_strategy : str
        Strategy for handling exchange rates when the reference year is in the future relative to available data:
        "reference" (use TECH_DATA_REFERENCE_YEAR),
        "latest" (use most recent available year),
        "custom" (use custom_future_exchange_rate).
    custom_future_exchange_rate : float, optional
        Custom exchange rate to use if future_exchange_rate_strategy is "custom".
    storage_techs : dict, optional
        A dictionary mapping storage carriers to their respective technology parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by technology with columns for fixed (annualized capital cost), marginal_cost, co2_emissions, and any other relevant parameters.
    """
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3

    # apply filter on financial_case and scenario, if they are contained in the cost dataframe
    wished_cost_scenario = config["cost_scenario"]
    wished_financial_case = config["financial_case"]
    for col in ["scenario", "financial_case"]:
        if col in costs.columns:
            costs[col] = costs[col].replace("", pd.NA)

    if "scenario" in costs.columns:
        costs = costs[
            (costs["scenario"].str.casefold() == wished_cost_scenario.casefold())
            | (costs["scenario"].isnull())
        ]

    if "financial_case" in costs.columns:
        costs = costs[
            (costs["financial_case"].str.casefold() == wished_financial_case.casefold())
            | (costs["financial_case"].isnull())
        ]

    if costs["currency_year"].isnull().any():
        logger.warning(
            "Some rows are missing 'currency_year' and will be skipped in currency conversion."
        )

    # Build a shared cache for exchange rates using the global reference_year
    _currency_conversion_cache = build_currency_conversion_cache(
        costs,
        output_currency,
        default_exchange_rate=default_exchange_rate,
        future_exchange_rate_strategy=future_exchange_rate_strategy,
        custom_future_exchange_rate=custom_future_exchange_rate,
    )

    modified_costs = apply_currency_conversion(
        costs, output_currency, _currency_conversion_cache
    )

    modified_costs = (
        modified_costs.loc[:, "value"]
        .unstack(level=1)
        .groupby("technology")
        .sum(min_count=1)
    )
    modified_costs = modified_costs.fillna(fill_values)

    for attr in ("investment", "lifetime", "FOM", "VOM", "efficiency", "fuel"):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            modified_costs.loc[overwrites.index, attr] = overwrites
            logger.info(
                f"Overwriting {attr} of {overwrites.index} to {overwrites.values}"
            )

    def annuity_factor(v):
        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    modified_costs["fixed"] = [
        annuity_factor(v) * v["investment"] * Nyears
        for _, v in modified_costs.iterrows()
    ]
    modified_costs["marginal_cost"] = (
        modified_costs["VOM"] + modified_costs["fuel"] / modified_costs["efficiency"]
    )

    modified_costs = calculate_cost_for_storage_units(
        modified_costs, max_hours, storage_techs, costs_name="fixed"
    )

    return modified_costs


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("process_cost_data", scope="sec")

    n = pypsa.Network(snakemake.input.network)
    Nyears = n.snapshot_weightings.generators.sum() / 8760.0

    if snakemake.wildcards.scope == "elec":
        costs = load_costs(
            snakemake.input.costs,
            snakemake.params.costs,
            snakemake.params.max_hours,
            Nyears,
            snakemake.params.storage_techs,
        )

    if snakemake.wildcards.scope == "sec":
        costs = prepare_costs(
            snakemake.input.costs,
            snakemake.params.costs,
            snakemake.params.costs["output_currency"],
            snakemake.params.costs["fill_values"],
            Nyears,
            snakemake.params.costs["default_exchange_rate"],
            snakemake.params.costs["future_exchange_rate_strategy"],
            snakemake.params.costs["custom_future_exchange_rate"],
            snakemake.params.max_hours,
            snakemake.params.storage_techs,
        )

    costs.to_csv(snakemake.output[0])
