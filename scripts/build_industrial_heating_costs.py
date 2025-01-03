# -*- coding: utf-8 -*-
"""Builds a file including technology cost and efficiency data for industrial heating.
The data is separate from other technology costs because they are specifically compiled
to provide technological competition for Enhanced Geothermal Systems (EGS) in the north
American context.
"""

import logging

logger = logging.getLogger(__name__)

import pandas as pd

from _helpers import configure_logging


if __name__ == '__main__':

    configure_logging(snakemake)

    techdata = pd.read_csv(snakemake.input['costs'], index_col=[0,1])

    cost_year = int(snakemake.params['cost_year'])

    temperature_bands = {
        'low': '<= 150C',
        'medium': '150C - 300C',
    }

    idx = pd.IndexSlice

    manual_costs = pd.DataFrame(
        columns=techdata.columns,
        index=pd.MultiIndex.from_tuples((), names=techdata.index.names)
    )

    # add molten salt storage costs
    manual_costs.loc[idx['low-temp molten salt store', 'investment'], 'value'] = 41_600
    manual_costs.loc[idx['low-temp molten salt store', 'investment'], 'unit'] = '$/MWhth'
    manual_costs.loc[idx['low-temp molten salt store', 'investment'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt store', 'investment'], 'further description'] = 'Salt Media + Storage Tank, for a system with 1000MWh capacity'

    manual_costs.loc[idx['low-temp molten salt store', 'FOM'], 'value'] = 1.5
    manual_costs.loc[idx['low-temp molten salt store', 'FOM'], 'unit'] = '%'
    manual_costs.loc[idx['low-temp molten salt store', 'FOM'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt store', 'FOM'], 'further description'] = 'Salt Media + Storage Tank, for a system with 1000MWh capacity'

    manual_costs.loc[idx['low-temp molten salt store', 'lifetime'], 'value'] = 35
    manual_costs.loc[idx['low-temp molten salt store', 'lifetime'], 'unit'] = 'years'
    manual_costs.loc[idx['low-temp molten salt store', 'lifetime'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt store', 'lifetime'], 'further description'] = 'Taken from a typical range of 30-50 years'

    manual_costs.loc[idx['low-temp molten salt store', 'efficiency'], 'value'] = 99.95
    manual_costs.loc[idx['low-temp molten salt store', 'efficiency'], 'unit'] = '%/hour'
    manual_costs.loc[idx['low-temp molten salt store', 'efficiency'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt store', 'efficiency'], 'further description'] = 'Typical value assumed for molten salt storage with a daily thermal loss of 1-2% per day'


    manual_costs.loc[idx['low-temp molten salt discharger', 'investment'], 'value'] = 30_000
    manual_costs.loc[idx['low-temp molten salt discharger', 'investment'], 'unit'] = '$/MWth'
    manual_costs.loc[idx['low-temp molten salt discharger', 'investment'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt discharger', 'investment'], 'further description'] = 'Assumes output is heat only, and does not need conversion to AC'
    manual_costs.loc[idx['low-temp molten salt discharger', 'investment'], 'currency_year'] = 2022.

    manual_costs.loc[idx['low-temp molten salt discharger', 'FOM'], 'value'] = 1.
    manual_costs.loc[idx['low-temp molten salt discharger', 'FOM'], 'unit'] = '%'
    manual_costs.loc[idx['low-temp molten salt discharger', 'FOM'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt discharger', 'FOM'], 'further description'] = 'Assumes output is heat only, and does not need conversion to AC'
    manual_costs.loc[idx['low-temp molten salt discharger', 'FOM'], 'currency_year'] = 2022.

    manual_costs.loc[idx['low-temp molten salt discharger', 'lifetime'], 'value'] = 35
    manual_costs.loc[idx['low-temp molten salt discharger', 'lifetime'], 'unit'] = 'years'
    manual_costs.loc[idx['low-temp molten salt discharger', 'lifetime'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt discharger', 'lifetime'], 'further description'] = 'Taken from a typical range of 30-50 years'
    manual_costs.loc[idx['low-temp molten salt discharger', 'lifetime'], 'currency_year'] = 2022.

    manual_costs.loc[idx['low-temp molten salt discharger', 'efficiency'], 'value'] = 98.
    manual_costs.loc[idx['low-temp molten salt discharger', 'efficiency'], 'unit'] = '%'
    manual_costs.loc[idx['low-temp molten salt discharger', 'efficiency'], 'source'] = 'Viswanathan_2022 - Energy Storage Grand Challenge Cost and Performance Assessment 2022'
    manual_costs.loc[idx['low-temp molten salt discharger', 'efficiency'], 'further description'] = 'Typical value for heat exchangers'
    manual_costs.loc[idx['low-temp molten salt discharger', 'efficiency'], 'currency_year'] = 2022.

    # add SHIP costs
    manual_costs.loc[idx['solar heat for industrial processes', 'investment'], 'value'] = 350_000
    manual_costs.loc[idx['solar heat for industrial processes', 'investment'], 'unit'] = '$/MWth'
    manual_costs.loc[idx['solar heat for industrial processes', 'investment'], 'source'] = 'IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf'
    manual_costs.loc[idx['solar heat for industrial processes', 'investment'], 'further description'] = 'For temperatures between 100-400C; Assumes Solar Collectors Coupled with Thermal Storage; Value reverse-engineered from LCOH; kWth is nameplate capacity'
    manual_costs.loc[idx['solar heat for industrial processes', 'investment'], 'currency_year'] = 2022.

    manual_costs.loc[idx['solar heat for industrial processes', 'FOM'], 'value'] = 2
    manual_costs.loc[idx['solar heat for industrial processes', 'FOM'], 'unit'] = '%'
    manual_costs.loc[idx['solar heat for industrial processes', 'FOM'], 'source'] = 'IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf'
    manual_costs.loc[idx['solar heat for industrial processes', 'FOM'], 'further description'] = 'Typical value for solar thermal'
    manual_costs.loc[idx['solar heat for industrial processes', 'FOM'], 'currency_year'] = 2022.

    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-mexico'], 'value'] = 30
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-mexico'], 'unit'] = '%'
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-mexico'], 'source'] = 'IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf'
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-mexico'], 'further description'] = 'Data point number one for capacity factor in Mexico. Should be used to interpolate capacity factor for other locations' 
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-mexico'], 'currency_year'] = 2022.

    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-germany'], 'value'] = 20
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-germany'], 'unit'] = '%'
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-germany'], 'source'] = 'IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf'
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-germany'], 'further description'] = 'Data point number two for capacity factor in Germany/Austria. Should be used to interpolate capacity factor for other locations' 
    manual_costs.loc[idx['solar heat for industrial processes', 'efficiency-germany'], 'currency_year'] = 2022.

    manual_costs.loc[idx['solar heat for industrial processes', 'lifetime'], 'value'] = 25
    manual_costs.loc[idx['solar heat for industrial processes', 'lifetime'], 'unit'] = 'years'
    manual_costs.loc[idx['solar heat for industrial processes', 'lifetime'], 'source'] = 'IEA-SHC-Task64-SubtaskE-D.E2-D.E3.pdf'
    manual_costs.loc[idx['solar heat for industrial processes', 'lifetime'], 'further description'] = 'Typical value for solar thermal'
    manual_costs.loc[idx['solar heat for industrial processes', 'lifetime'], 'currency_year'] = 2022.

    # add heat pump costs
    manual_costs = pd.concat([
        manual_costs,
        techdata.loc[
            idx[['industrial heat pump high temperature', 'industrial heat pump medium temperature'],:]
        ]
        ])

    # add solar thermal costs
    manual_costs = pd.concat([
        manual_costs,
        techdata.loc[
            idx[['csp-tower'], :]
        ]
        ])

    # add biogas costs
    manual_costs.loc[idx['biogas', 'fuel cost'], 'value'] = 79.
    manual_costs.loc[idx['biogas', 'fuel cost'], 'unit'] = '$/MWh'
    manual_costs.loc[idx['biogas', 'fuel cost'], 'source'] = 'IEA: https://www.iea.org/reports/outlook-for-biogas-and-biomethane-prospects-for-organic-growth/sustainable-supply-potential-and-costs'
    manual_costs.loc[idx['biogas', 'fuel cost'], 'further description'] = 'From regionalised supply curve for North America, assuming the higher end will dictate the price.'

    # add respective steam boiler costs
    investemt_costs = {
        2020: 60_000, # EUR/MW
        2030: 50_000, # EUR/MW
        2040: 50_000, # EUR/MW
        2050: 50_000, # EUR/MW
    }
    manual_costs.loc[idx['steam boiler gas cond', 'investment'], 'value'] = investemt_costs[cost_year]
    manual_costs.loc[idx['steam boiler gas cond', 'investment'], 'unit'] = 'EUR/MWth'
    manual_costs.loc[idx['steam boiler gas cond', 'investment'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['steam boiler gas cond', 'investment'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW'
    manual_costs.loc[idx['steam boiler gas cond', 'investment'], 'currency_year'] = 2019.

    manual_costs.loc[idx['steam boiler gas cond', 'VOM'], 'value'] = 1.
    manual_costs.loc[idx['steam boiler gas cond', 'VOM'], 'unit'] = 'EUR/MWhth'
    manual_costs.loc[idx['steam boiler gas cond', 'VOM'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['steam boiler gas cond', 'VOM'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW'
    manual_costs.loc[idx['steam boiler gas cond', 'VOM'], 'currency_year'] = 2019.

    manual_costs.loc[idx['steam boiler gas cond', 'efficiency'], 'value'] = 103.
    manual_costs.loc[idx['steam boiler gas cond', 'efficiency'], 'unit'] = '%'
    manual_costs.loc[idx['steam boiler gas cond', 'efficiency'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['steam boiler gas cond', 'efficiency'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW'
    manual_costs.loc[idx['steam boiler gas cond', 'efficiency'], 'currency_year'] = 2019.

    manual_costs.loc[idx['steam boiler gas cond', 'lifetime'], 'value'] = 25
    manual_costs.loc[idx['steam boiler gas cond', 'lifetime'], 'unit'] = 'years'
    manual_costs.loc[idx['steam boiler gas cond', 'lifetime'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['steam boiler gas cond', 'lifetime'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular steam boilers; typical size 20MW'
    manual_costs.loc[idx['steam boiler gas cond', 'lifetime'], 'currency_year'] = 2019.

    # add respective steam boiler costs
    investemt_costs = {
        2020: 50_000, # EUR/MW
        2030: 40_000, # EUR/MW
        2040: 40_000, # EUR/MW
        2050: 40_000, # EUR/MW
    }
    manual_costs.loc[idx['hot water boiler gas cond', 'investment'], 'value'] = investemt_costs[cost_year]
    manual_costs.loc[idx['hot water boiler gas cond', 'investment'], 'unit'] = 'EUR/MWth'
    manual_costs.loc[idx['hot water boiler gas cond', 'investment'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['hot water boiler gas cond', 'investment'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW'
    manual_costs.loc[idx['hot water boiler gas cond', 'investment'], 'currency_year'] = 2019.

    manual_costs.loc[idx['hot water boiler gas cond', 'VOM'], 'value'] = 1.
    manual_costs.loc[idx['hot water boiler gas cond', 'VOM'], 'unit'] = 'EUR/MWhth'
    manual_costs.loc[idx['hot water boiler gas cond', 'VOM'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['hot water boiler gas cond', 'VOM'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW'
    manual_costs.loc[idx['hot water boiler gas cond', 'VOM'], 'currency_year'] = 2019.

    manual_costs.loc[idx['hot water boiler gas cond', 'efficiency'], 'value'] = 103.
    manual_costs.loc[idx['hot water boiler gas cond', 'efficiency'], 'unit'] = '%'
    manual_costs.loc[idx['hot water boiler gas cond', 'efficiency'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['hot water boiler gas cond', 'efficiency'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW'
    manual_costs.loc[idx['hot water boiler gas cond', 'efficiency'], 'currency_year'] = 2019.

    manual_costs.loc[idx['hot water boiler gas cond', 'lifetime'], 'value'] = 25.
    manual_costs.loc[idx['hot water boiler gas cond', 'lifetime'], 'unit'] = 'years'
    manual_costs.loc[idx['hot water boiler gas cond', 'lifetime'], 'source'] = 'Technology Data for Industrial Process Heat; https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_industrial_process_heat_0.xlsx'
    manual_costs.loc[idx['hot water boiler gas cond', 'lifetime'], 'further description'] = 'Condensed version chosen for simplicity, marginal difference to regular hot water boilers; typical size 20MW'
    manual_costs.loc[idx['hot water boiler gas cond', 'lifetime'], 'currency_year'] = 2019.

    # lower temperature thermal storage
    manual_costs = pd.concat([
        manual_costs,
        techdata.loc[
            idx[['central water tank storage'], :]
            ]
        ])

    manual_costs.to_csv(snakemake.output['industrial_heating_costs'])
