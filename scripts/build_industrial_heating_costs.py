# -*- coding: utf-8 -*-
"""Builds a file including technology cost and efficiency data for industrial heating.
The data is separate from other technology costs because they are specifically compiled
to provide technological competition for Enhanced Geothermal Systems (EGS) in the north
American context.
"""

import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
from _helpers import configure_logging

if __name__ == "__main__":

    configure_logging(snakemake)

    techdata = pd.read_csv(snakemake.input["costs"], index_col=[0, 1])

    cost_year = int(snakemake.params["cost_year"])

    temperature_bands = {
        "low": "<= 80C",
        "medium": "80C - 150C",
        "high": "150C - 250C",
    }

    idx = pd.IndexSlice

    manual_costs = pd.DataFrame(
        columns=techdata.columns,
        index=pd.MultiIndex.from_tuples((), names=techdata.index.names),
    )

    # solar thermal (parabolic trough delivering heat at 150C-250C):
    manual_costs.loc[idx["solar thermal parabolic trough", "investment"], "value"] = (
        800_000
    )
    manual_costs.loc[idx["solar thermal parabolic trough", "investment"], "unit"] = (
        "$/MWth"
    )
    manual_costs.loc[idx["solar thermal parabolic trough", "investment"], "source"] = (
        "https://proceedings.ises.org/conference/eurosun2024/papers/eurosun2024-0042-Akar.pdf"
    )
    manual_costs.loc[
        idx["solar thermal parabolic trough", "investment"], "further description"
    ] = "Varies between 300-500$/kWth depending on assumptions on building materials"

    manual_costs.loc[idx["solar thermal parabolic trough", "lifetime"], "value"] = 25
    manual_costs.loc[idx["solar thermal parabolic trough", "lifetime"], "unit"] = (
        "years"
    )
    manual_costs.loc[idx["solar thermal parabolic trough", "lifetime"], "source"] = (
        "https://proceedings.ises.org/conference/eurosun2024/papers/eurosun2024-0042-Akar.pdf"
    )

    manual_costs.loc[idx["solar thermal parabolic trough", "VOM"], "value"] = 0.012
    manual_costs.loc[idx["solar thermal parabolic trough", "VOM"], "unit"] = "$/MWth"
    manual_costs.loc[idx["solar thermal parabolic trough", "VOM"], "source"] = (
        "https://proceedings.ises.org/conference/eurosun2024/papers/eurosun2024-0042-Akar.pdf"
    )

    # oil-based storage unit coupled to solar thermal
    manual_costs.loc[idx["oil-based storage unit", "investment"], "value"] = (
        150_000 / 6
    )  # 6 is the E-P ratio. investment is expressed as power-capacity, because it is implemented as a storage unit in PyPSA
    manual_costs.loc[idx["oil-based storage unit", "investment"], "unit"] = "$/MWth"
    manual_costs.loc[idx["oil-based storage unit", "investment"], "source"] = (
        "https://docs.nrel.gov/docs/fy10osti/47605.pdf"
    )
    manual_costs.loc[
        idx["oil-based storage unit", "investment"], "further description"
    ] = "This assumes the heat-transfer fluid of the solar field is also used for storage"

    manual_costs.loc[idx["oil-based storage unit", "efficiency_store"], "value"] = np.sqrt(
        0.9
    )
    manual_costs.loc[idx["oil-based storage unit", "efficiency_store"], "unit"] = "per unit"
    manual_costs.loc[idx["oil-based storage unit", "efficiency_store"], "source"] = (
        "https://docs.nrel.gov/docs/fy10osti/47605.pdf"
    )

    manual_costs.loc[idx["oil-based storage unit", "efficiency_dispatch"], "value"] = np.sqrt(
        0.9
    )
    manual_costs.loc[idx["oil-based storage unit", "efficiency_dispatch"], "unit"] = "per unit"
    manual_costs.loc[idx["oil-based storage unit", "efficiency_dispatch"], "source"] = (
        "https://docs.nrel.gov/docs/fy10osti/47605.pdf"
    )

    manual_costs.loc[idx["oil-based storage unit", "lifetime"], "value"] = 25
    manual_costs.loc[idx["oil-based storage unit", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["oil-based storage unit", "lifetime"], "source"] = (
        "https://docs.nrel.gov/docs/fy10osti/47605.pdf"
    )

    manual_costs.loc[idx["oil-based storage unit", "max_hours"], "value"] = (
        6  # 6 is the E-P ratio. investment is expressed as power-capacity, because it is implemented as a storage unit in PyPSA
    )
    manual_costs.loc[idx["oil-based storage unit", "max_hours"], "unit"] = "hours"
    manual_costs.loc[idx["oil-based storage unit", "max_hours"], "source"] = (
        "https://docs.nrel.gov/docs/fy10osti/47605.pdf"
    )

    # add heat pump costs
    manual_costs = pd.concat(
        [
            manual_costs,
            techdata.loc[
                idx[
                    [
                        "industrial heat pump high temperature",
                        "industrial heat pump medium temperature",
                    ],
                    :,
                ]
            ],
        ]
    )

    # adjusting COP based on NREL data
    manual_costs.loc[idx["industrial heat pump high temperature", "efficiency"], "value"] = 2.5
    manual_costs.loc[idx["industrial heat pump high temperature", "efficiency"], "unit"] = "per unit"
    manual_costs.loc[idx["industrial heat pump high temperature", "efficiency"], "source"] = (
        "https://docs.nrel.gov/docs/fy23osti/84560.pdf"
    )

    # solar thermal for 80-150C temperature band: Linear Fresnel Reflectors
    manual_costs.loc[idx["linear fresnel reflector", "investment"], "value"] = 700_000
    manual_costs.loc[idx["linear fresnel reflector", "investment"], "unit"] = "$/MWth"
    manual_costs.loc[idx["linear fresnel reflector", "investment"], "source"] = (
        "https://docs.nrel.gov/docs/fy22osti/81147.pdf#page=1.00&gsr=0"
    )

    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "value"] = 25
    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "source"] = (
        "https://task64.iea-shc.org/Data/Sites/1/publications/IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf"
    )

    manual_costs.loc[idx["linear fresnel reflector", "VOM"], "value"] = 0.015
    manual_costs.loc[idx["linear fresnel reflector", "VOM"], "unit"] = "$/MWthh"
    manual_costs.loc[idx["linear fresnel reflector", "VOM"], "source"] = (
        "https://docs.nrel.gov/docs/fy25osti/93381.pdf"
    )

    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "value"] = 25
    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["linear fresnel reflector", "lifetime"], "source"] = (
        "https://docs.nrel.gov/docs/fy25osti/93381.pdf"
    )

    # for storage the same thermal-oil based storage unit is used

    # Solar for 50-80C temperature band: glazed flat plate collectors
    manual_costs.loc[idx["glazed flat plate collector", "investment"], "value"] = (
        550_000
    )
    manual_costs.loc[idx["glazed flat plate collector", "investment"], "unit"] = (
        "$/MWth"
    )
    manual_costs.loc[idx["glazed flat plate collector", "investment"], "source"] = (
        "https://solarthermalworld.org/news/cost-comparison-of-industrial-heat-from-solar-thermal-and-pv/"
    )

    manual_costs.loc[idx["glazed flat plate collector", "VOM"], "value"] = 1.8
    manual_costs.loc[idx["glazed flat plate collector", "VOM"], "unit"] = "$/MWthh"
    manual_costs.loc[idx["glazed flat plate collector", "VOM"], "source"] = (
        "https://solarthermalworld.org/news/cost-comparison-of-industrial-heat-from-solar-thermal-and-pv/"
    )

    manual_costs.loc[idx["glazed flat plate collector", "lifetime"], "value"] = 25
    manual_costs.loc[idx["glazed flat plate collector", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["glazed flat plate collector", "lifetime"], "source"] = (
        "https://task64.iea-shc.org/Data/Sites/1/publications/IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf"
    )

    # Steel or concrete water tank for 6h 50-80C storage
    manual_costs.loc[idx["steel or concrete water tank", "investment"], "value"] = (
        5_000 / 6
    )  # 6 is the E-P ratio. investment is expressed as power-capacity, because it is implemented as a storage unit in PyPSA
    manual_costs.loc[idx["steel or concrete water tank", "investment"], "unit"] = (
        "$/MWth"
    )
    manual_costs.loc[idx["steel or concrete water tank", "investment"], "source"] = (
        "https://iea-es.org/wp-content/uploads/public/FactSheet_Thermal_Sensible_Water_2024-07-10.pdf"
    )
    manual_costs.loc[
        idx["steel or concrete water tank", "investment"], "further description"
    ] = "Investment is expressed as power-capacity, because it is implemented as a StorageUnit in PyPSA"

    manual_costs.loc[idx["steel or concrete water tank", "efficiency_store"], "value"] = (
        np.sqrt(0.9)
    )
    manual_costs.loc[idx["steel or concrete water tank", "efficiency_store"], "unit"] = (
        "per unit"
    )
    manual_costs.loc[idx["steel or concrete water tank", "efficiency_store"], "source"] = (
        "https://iea-es.org/wp-content/uploads/public/FactSheet_Thermal_Sensible_Water_2024-07-10.pdf"
    )

    manual_costs.loc[idx["steel or concrete water tank", "efficiency_dispatch"], "value"] = (
        np.sqrt(0.9)
    )
    manual_costs.loc[idx["steel or concrete water tank", "efficiency_dispatch"], "unit"] = (
        "per unit"
    )
    manual_costs.loc[idx["steel or concrete water tank", "efficiency_dispatch"], "source"] = (
        "https://iea-es.org/wp-content/uploads/public/FactSheet_Thermal_Sensible_Water_2024-07-10.pdf"
    )

    manual_costs.loc[idx["steel or concrete water tank", "max_hours"], "value"] = 6
    manual_costs.loc[idx["steel or concrete water tank", "max_hours"], "unit"] = "hours"
    manual_costs.loc[idx["steel or concrete water tank", "max_hours"], "source"] = (
        "https://iea-es.org/wp-content/uploads/public/FactSheet_Thermal_Sensible_Water_2024-07-10.pdf"
    )

    # add biogas costs
    manual_costs.loc[idx["biogas", "fuel cost"], "value"] = 79.0
    manual_costs.loc[idx["biogas", "fuel cost"], "unit"] = "$/MWh"
    manual_costs.loc[idx["biogas", "fuel cost"], "source"] = (
        "IEA: https://www.iea.org/reports/outlook-for-biogas-and-biomethane-prospects-for-organic-growth/sustainable-supply-potential-and-costs"
    )
    manual_costs.loc[idx["biogas", "fuel cost"], "further description"] = (
        "From regionalised supply curve for North America, assuming the higher end will dictate the price."
    )

    # add respective steam boiler costs
    investment_costs = {
        2020: 60_000,  # EUR/MW
        2030: 50_000,  # EUR/MW
        2035: 50_000,  # EUR/MW
        2040: 50_000,  # EUR/MW
        2050: 50_000,  # EUR/MW
    }
    manual_costs.loc[idx["steam boiler gas cond", "investment"], "value"] = (
        investment_costs[cost_year]
    )
    manual_costs.loc[idx["steam boiler gas cond", "investment"], "unit"] = "EUR/MWth"
    manual_costs.loc[idx["steam boiler gas cond", "investment"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["steam boiler gas cond", "investment"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW"
    manual_costs.loc[idx["steam boiler gas cond", "investment"], "currency_year"] = (
        2019.0
    )

    manual_costs.loc[idx["steam boiler gas cond", "VOM"], "value"] = 1.0
    manual_costs.loc[idx["steam boiler gas cond", "VOM"], "unit"] = "EUR/MWhth"
    manual_costs.loc[idx["steam boiler gas cond", "VOM"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[idx["steam boiler gas cond", "VOM"], "further description"] = (
        "Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW"
    )
    manual_costs.loc[idx["steam boiler gas cond", "VOM"], "currency_year"] = 2019.0

    manual_costs.loc[idx["steam boiler gas cond", "efficiency"], "value"] = 103.0
    manual_costs.loc[idx["steam boiler gas cond", "efficiency"], "unit"] = "%"
    manual_costs.loc[idx["steam boiler gas cond", "efficiency"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["steam boiler gas cond", "efficiency"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW"
    manual_costs.loc[idx["steam boiler gas cond", "efficiency"], "currency_year"] = (
        2019.0
    )

    manual_costs.loc[idx["steam boiler gas cond", "lifetime"], "value"] = 25
    manual_costs.loc[idx["steam boiler gas cond", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["steam boiler gas cond", "lifetime"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["steam boiler gas cond", "lifetime"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW"
    manual_costs.loc[idx["steam boiler gas cond", "lifetime"], "currency_year"] = 2019.0

    # add respective steam boiler costs
    investment_costs = {
        2020: 50_000,  # EUR/MW
        2030: 40_000,  # EUR/MW
        2035: 40_000,  # EUR/MW
        2040: 40_000,  # EUR/MW
        2050: 40_000,  # EUR/MW
    }
    manual_costs.loc[idx["hot water boiler gas cond", "investment"], "value"] = (
        investment_costs[cost_year]
    )
    manual_costs.loc[idx["hot water boiler gas cond", "investment"], "unit"] = (
        "EUR/MWth"
    )
    manual_costs.loc[idx["hot water boiler gas cond", "investment"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["hot water boiler gas cond", "investment"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW"
    manual_costs.loc[
        idx["hot water boiler gas cond", "investment"], "currency_year"
    ] = 2019.0

    manual_costs.loc[idx["hot water boiler gas cond", "VOM"], "value"] = 1.0
    manual_costs.loc[idx["hot water boiler gas cond", "VOM"], "unit"] = "EUR/MWhth"
    manual_costs.loc[idx["hot water boiler gas cond", "VOM"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[idx["hot water boiler gas cond", "VOM"], "further description"] = (
        "Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW"
    )
    manual_costs.loc[idx["hot water boiler gas cond", "VOM"], "currency_year"] = 2019.0

    manual_costs.loc[idx["hot water boiler gas cond", "efficiency"], "value"] = 103.0
    manual_costs.loc[idx["hot water boiler gas cond", "efficiency"], "unit"] = "%"
    manual_costs.loc[idx["hot water boiler gas cond", "efficiency"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["hot water boiler gas cond", "efficiency"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW"
    manual_costs.loc[
        idx["hot water boiler gas cond", "efficiency"], "currency_year"
    ] = 2019.0

    manual_costs.loc[idx["hot water boiler gas cond", "lifetime"], "value"] = 25.0
    manual_costs.loc[idx["hot water boiler gas cond", "lifetime"], "unit"] = "years"
    manual_costs.loc[idx["hot water boiler gas cond", "lifetime"], "source"] = (
        "Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx"
    )
    manual_costs.loc[
        idx["hot water boiler gas cond", "lifetime"], "further description"
    ] = "Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW"
    manual_costs.loc[idx["hot water boiler gas cond", "lifetime"], "currency_year"] = (
        2019.0
    )

    # lower temperature thermal storage
    # manual_costs = pd.concat(
    #     [manual_costs, techdata.loc[idx[["central water tank storage"], :]]]
    # )

    manual_costs.to_csv(snakemake.output["industrial_heating_costs"])
