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
    Source: Tables 17 & 19 of Cutler et al (2023) NREL Report

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odw is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = –0.30428
        b =  0.11805
        c = –0.00342
        d = –0.00626
        e =  0.00070
        f = –0.00047
    if unit_type == "two stage, low speed":
        a = –0.42738
        b =  0.14191
        c = –0.00412
        d = –0.01406
        e =  0.00083
        f = –0.00043
    if unit_type == "two stage, high speed":
        a =  0.04232
        b =  0.07892
        c = –0.00238
        d = –0.00304
        e =  0.00053
        f = –0.00032

    EIR = a + b * t_ewb + c * (t_ebw^2) + d * t_odb + e * (t_odb^2) + f * t_ewb * t_odb

    return EIR

def eir_heat_pump_cooling(t_ewb, t_odw, unit_type):
    """
    Source: Tables 17 & 19 of Cutler et al (2023) NREL Report

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odw is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = –0.350448
        b =  0.116810 
        c = –0.003400
        d = –0.001226
        e =  0.000601 
        f = –0.000467
    if unit_type == "two stage, low speed":
        a = –0.582916
        b = 0.158101 
        c = –0.004398
        d = –0.020335
        e = 0.001080 
        f = –0.000640
    if unit_type == "two stage, high speed":
        a = –0.488196
        b = 0.099162 
        c = –0.002370
        d = 0.019503 
        e = 0.000430 
        f = –0.001097

    EIR = a + b * t_ewb + c * (t_ebw^2) + d * t_odb + e * (t_odb^2) + f * t_ewb * t_odb

    return EIR

def eir_heat_pump_heating(t_ewb, t_odw, unit_type):
    """
    Source: Tables 17 & 19 of Cutler et al (2023) NREL Report

    EIR = 1/COP
    t_ewb is the entering wet-bulb temperature
    t_odw is the outdoor dry-bulb temperature
    """
    if unit_type == "single stage":
        a = 0.704658 
        b = 0.008767 
        c = 0.000625 
        d = –0.009037
        e = 0.000738 
        f = –0.001025
    if unit_type == "two stage, low speed":
        a = 0.551837 
        b = 0.020380 
        c = 0.000546 
        d = –0.009638
        e = 0.000785 
        f = –0.001250
    if unit_type == "two stage, high speed":
        a = 0.815840
        b = –0.006150
        c = 0.001021
        d = –0.001301
        e = 0.001083
        f = –0.001487

    EIR = a + b * t_ewb + c * (t_ebw^2) + d * t_odb + e * (t_odb^2) + f * t_ewb * t_odb

    return EIR    

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_cop_profiles",
            simpl="",
            clusters=4,
        )

    for area in ["total", "urban", "rural"]:
        for source in ["air", "soil"]:
            source_T = xr.open_dataarray(snakemake.input[f"temp_{source}_{area}"])

            delta_T = snakemake.params.heat_pump_sink_T - source_T

            cop = coefficient_of_performance(delta_T, source)

            cop.to_netcdf(snakemake.output[f"cop_{source}_{area}"])
