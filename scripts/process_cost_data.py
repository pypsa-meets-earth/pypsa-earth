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
def annuity(n, r):
    """
    Calculate the annuity factor for an asset with lifetime n years and.

    discount rate of r, e.g. annuity(20, 0.05) * 20 = 1.6
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
):
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
    df,
    output_currency,
    default_exchange_rate=None,
    future_exchange_rate_strategy: str = "reference",
    custom_future_exchange_rate: float = None,
):
    """
    Builds a cache of exchange rates for all unique (currency, output_currency) pairs,
    always using the module-level reference_year.
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


def apply_currency_conversion(cost_dataframe, output_currency, cache):
    """
    Applies exchange rates from the cache to convert all cost values and units.

    All rows are assumed to be in `*_reference_year` already (e.g. EUR_2020).
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


# ==============
# Overnight Runs
# ==============


def load_costs(tech_costs, config, max_hours, Nyears=1):
    """
    Set all asset costs and other parameters.
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

    def costs_for_storage(store, link1, link2=None, max_hours=1.0):
        capital_cost = link1["capital_cost"] + max_hours * store["capital_cost"]
        if link2 is not None:
            capital_cost += link2["capital_cost"]
        return pd.Series(
            dict(capital_cost=capital_cost, marginal_cost=0.0, co2_emissions=0.0)
        )

    costs.loc["battery"] = costs_for_storage(
        costs.loc["battery storage"],
        costs.loc["battery inverter"],
        max_hours=max_hours["battery"],
    )
    costs.loc["H2"] = costs_for_storage(
        costs.loc["hydrogen storage tank"],
        costs.loc["fuel cell"],
        costs.loc["electrolysis"],
        max_hours=max_hours["H2"],
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


# ===================
# Myopic/Perfect Runs
# ===================


def prepare_costs(
    cost_file: str,
    config: dict,
    output_currency: str,
    fill_values: dict,
    Nyears: float | int = 1,
    default_exchange_rate: float = None,
    future_exchange_rate_strategy: str = "latest",
    custom_future_exchange_rate: float = None,
):
    """
    Loads and processes cost data, converting units and currency to a common format.

    Applies currency conversion, fills missing values, and computes fixed annualized costs.
    Always uses the module-level reference_year.
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
        )

    costs.to_csv(snakemake.output[0])
