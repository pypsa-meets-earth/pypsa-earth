# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Build COP time series for air- or ground-sourced heat pumps.
"""

import xarray as xr


def coefficient_of_performance(delta_T, source="air"):
    """
    COP is function of temp difference source to sink.

    The quadratic regression is based on Staffell et al. (2012)
    https://doi.org/10.1039/C2EE22653G.
    """
    if source == "air":
        return 6.81 - 0.121 * delta_T + 0.000630 * delta_T**2
    elif source == "soil":
        return 8.77 - 0.150 * delta_T + 0.000734 * delta_T**2
    else:
        raise NotImplementedError("'source' must be one of  ['air', 'soil']")


# TODO Translate wet-bulb temperature into the relative humidity


def eir_air_conditioner(t_ewb, t_odb, unit_type):
    """
    Source: Tables 17 & 19 of NREL Report Cutler et al (2013)
    "Improved Modeling of Residential Air Conditioners and Heat Pumps for
     Energy Calculations" available via
     https://www.nrel.gov/docs/fy13osti/56354.pdf

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odb is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = -0.30428
        b = 0.11805
        c = -0.00342
        d = -0.00626
        e = 0.00070
        f = -0.00047
    if unit_type == "two stage, low speed":
        a = -0.42738
        b = 0.14191
        c = -0.00412
        d = -0.01406
        e = 0.00083
        f = -0.00043
    if unit_type == "two stage, high speed":
        a = 0.04232
        b = 0.07892
        c = -0.00238
        d = -0.00304
        e = 0.00053
        f = -0.00032

    EIR = (
        a
        + b * t_ewb
        + c * (t_ewb * t_ewb)
        + d * t_odb
        + e * (t_odb * t_odb)
        + f * t_ewb * t_odb
    )

    return EIR


def eir_heat_pump_cooling(t_ewb, t_odb, unit_type):
    """
    Source: Equation (7) and Tables 17 & 19 of Cutler et al (2023) NREL Report
    "Improved Modeling of Residential Air Conditioners and Heat Pumps for
     Energy Calculations" available via
     https://www.nrel.gov/docs/fy13osti/56354.pdf

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odb is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = -0.350448
        b = 0.116810
        c = -0.003400
        d = -0.001226
        e = 0.000601
        f = -0.000467
    if unit_type == "two stage, low speed":
        a = -0.582916
        b = 0.158101
        c = -0.004398
        d = -0.020335
        e = 0.001080
        f = -0.000640
    if unit_type == "two stage, high speed":
        a = -0.488196
        b = 0.099162
        c = -0.002370
        d = 0.019503
        e = 0.000430
        f = -0.001097

    EIR = (
        a
        + b * t_ewb
        + c * (t_ewb * t_ewb)
        + d * t_odb
        + e * (t_odb * t_odb)
        + f * t_ewb * t_odb
    )

    return EIR


def eir_heat_pump_heating(t_ewb, t_odb, unit_type):
    """
    Source: Tables 17 & 19 of Cutler et al (2023) NREL Report
    "Improved Modeling of Residential Air Conditioners and Heat Pumps for
     Energy Calculations" available via
     https://www.nrel.gov/docs/fy13osti/56354.pdf

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odb is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = 0.704658
        b = 0.008767
        c = 0.000625
        d = -0.009037
        e = 0.000738
        f = -0.001025
    if unit_type == "two stage, low speed":
        a = 0.551837
        b = 0.020380
        c = 0.000546
        d = -0.009638
        e = 0.000785
        f = -0.001250
    if unit_type == "two stage, high speed":
        a = 0.815840
        b = -0.006150
        c = 0.001021
        d = -0.001301
        e = 0.001083
        f = -0.001487

    EIR = (
        a
        + b * t_ewb
        + c * (t_ewb * t_ewb)
        + d * t_odb
        + e * (t_odb * t_odb)
        + f * t_ewb * t_odb
    )

    return EIR


def capft_absorption_chiller_air_cool(t_chws, t_odb, unit_type):
    """
    Source: Table 74 of PNNL (2016) ANSI-ASHRAE-IES Standard Performance Rating Method

    Share of the available cooling capacity from the rated value

    t_chws [F] is the chilled water supply temperature
    t_odb [F] is the outdoor air dry-bulb temperature (or condenser water supply temperature)
    """
    if unit_type == "single stage absorption":
        a = 0.723412
        b = 0.079006
        c = -0.000897
        d = -0.025285
        e = -0.000048
        f = 0.000276
    if unit_type == "double stage absorption":
        a = -0.816039
        b = -0.038707
        c = 0.000450
        d = 0.071491
        e = -0.000636
        f = 0.000312
    if unit_type == "engine driven chiller":
        a = 0.573597
        b = 0.0186802
        c = 0.000000
        d = -0.00465325
        e = 0.000000
        f = 0.000000

    # !!!!!!!!!!!!!!!!!!!!!
    # TODO Transform F -> C
    # !!!!!!!!!!!!!!!!!!!!!
    CAP_ft = (
        a
        + b * t_chws
        + c * (t_chws ^ 2)
        + d * t_odb
        + e * (t_odb * t_odb)
        + f * t_chws * t_odb
    )

    return CAP_ft


def cop_air_conditioner(t_ewb, t_odb, unit_type):
    1 / eir_air_conditioner(t_ewb, t_odb, unit_type)


def cop_heat_pump_heating(t_ewb, t_odb, unit_type):
    return 1 / eir_heat_pump_heating(t_ewb, t_odb, unit_type)


def cop_heat_pump_cooling(t_ewb, t_odb, unit_type):
    1 / eir_heat_pump_cooling(t_ewb, t_odb, unit_type)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles", simpl="", clusters=10, planning_horizons=2030
        )

    t_wet_bulb = snakemake.params.t_wet_bulb

    # heating performance
    for area in ["total", "urban", "rural"]:
        # assuming that the heat source for heat pumps is outdoor air or soil
        for source in ["air", "soil"]:
            source_T = xr.open_dataarray(snakemake.input[f"temp_{source}_{area}"])

            delta_T = snakemake.params.heat_pump_sink_T - source_T

            cop = coefficient_of_performance(delta_T, source)

            cop.to_netcdf(snakemake.output[f"cop_{source}_{area}"])

    # currently we don't distinguish between rural and urban areas
    # which is different from how heating is being treated where `area` is in `["total", "urban", "rural"]`
    for area in ["total"]:
        # assuming the bulb point inside at some reasonable level
        # TODO check if K->C transformations is needed
        outdoor_T = xr.open_dataarray(snakemake.input[f"temp_air_{area}"])

        cop_hp_cooling = 1 / eir_heat_pump_cooling(
            t_ewb=t_wet_bulb, t_odb=outdoor_T, unit_type="single stage"
        )
        cop_ac_cooling = 1 / eir_air_conditioner(
            t_ewb=t_wet_bulb, t_odb=outdoor_T, unit_type="single stage"
        )

        capft_abch_cooling = capft_absorption_chiller_air_cool(
            t_chws=t_wet_bulb, t_odb=outdoor_T, unit_type="single stage absorption"
        )

        cop_hp_cooling.to_netcdf(snakemake.output[f"cop_hp_cooling_{area}"])
        cop_ac_cooling.to_netcdf(snakemake.output[f"cop_ac_cooling_{area}"])
        capft_abch_cooling.to_netcdf(snakemake.output[f"capft_abch_cooling_{area}"])
