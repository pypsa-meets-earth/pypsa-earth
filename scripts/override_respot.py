# -*- coding: utf-8 -*-

from itertools import dropwhile
import os
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pypsa
import pytz
import xarray as xr
from helpers import override_component_attrs, sets_path_to_root, mock_snakemake
# from helpers import (
#     create_dummy_data,
#     create_network_topology,
#     cycling_shift,
#     locate_bus,
#     mock_snakemake,
#     override_component_attrs,
#     prepare_costs,
#     three_2_two_digits_country,
#     two_2_three_digits_country,
# )

if __name__ == "__main__":
    if "snakemake" not in globals():
            os.chdir(os.path.dirname(os.path.abspath(__file__)))

            snakemake = mock_snakemake(
                "override_respot",
                simpl="",
                clusters="1002",
                ll="c1.0",
                opts="Co2L",
                planning_horizons="2030",
                sopts="3H",
                ir="low",
            )
            sets_path_to_root("pypsa-earth-sec")

    techs = snakemake.config["custom_data"]["renewables"]
    years = snakemake.config["scenario"]["planning_horizons"]#['2030', '2050']  # , 'offwind', 'onwind', 'csp']
    drs = snakemake.config["costs"]["discountrate"]  # , 'offwind', 'onwind', 'csp']

    suff={'sopv': 'solar', 'csp': 'csp', 'pvr': 'rooftop pv'}

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    m=n.copy()
    
    def override_values(tech, year, dr):
        buses=list(n.buses[n.buses.carrier=='AC'].index)
        
        # enertile_res_pot=pd.read_csv('postprocessed/{0}_{1}_{2}_potential.csv'.format(
        #     tech, year,scenario), index_col=0, parse_dates=True).filter(
        #         buses, axis=1)#.add_suffix(' ' + suff[tech])
        
        enertile_res_pot = pd.read_csv(snakemake.input["custom_res_pot_{0}_{1}_{2}".format(tech, year, dr)]
                           , index_col=0, parse_dates=True).filter(
                buses, axis=1)
    
        enertile_installable=pd.read_csv(snakemake.input["custom_res_ins_{0}_{1}_{2}".format(tech, year, dr)]
                             , index_col=0).filter(buses, axis=0).reset_index()
        
        #enertile_installable.columns=['Generator', 'p_nom_max']
        enertile_installable.rename(columns={'region':'Generator', 'potstepsizeMW':'p_nom_max'}
                                    , inplace=True)
        enertile_installable['Generator']=enertile_installable['Generator'].apply(lambda x: x+' '+suff[tech])#.rename('Generator')
        enertile_installable = enertile_installable.set_index('Generator')
        # def override_solar_techs(n, tech, respot, install):
        
        # tech_cols=n.generators_t.p_max_pu.filter(regex="solar$").columns 
        if suff[tech] in n.generators.carrier.unique():

            n.generators_t.p_max_pu.update(enertile_res_pot)
            n.generators.update(pd.Series(enertile_installable['p_nom_max']))

        else:
            n.madd(
                "Generator",
                buses,
                " " + suff[tech],
                bus=buses,
                carrier=suff[tech],
                p_nom_extendable=True,
                p_nom_max=enertile_installable["p_nom_max"].values,
                #weight=ds["weight"].to_pandas(),
                marginal_cost=enertile_installable['fixedomEuroPKW']*1000,
                capital_cost=enertile_installable['investmentEuroPKW']*1000,
                efficiency=1.0,
                p_max_pu=enertile_res_pot,
            )

    for tech in techs:
        for year in years:
            for dr in drs:
                override_values(tech, year, dr)

    n.export_to_netcdf(snakemake.output[0])    

    # override_values(techs[1], years[0], irs[0])                
                
    # for option in options:

    #     if option in n.generators.carrier:
    #         gen_ind = n.generators[n.generators.carrier=='solar'].index
    #         n.generators[gen_ind]= custom_res['']
    #     else
