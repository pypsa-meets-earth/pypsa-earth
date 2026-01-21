import pypsa

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import country_converter as coco

elec_bus_carrier = ["AC","DC","low voltage"]

battery_pair = {
    "BEV charger": "land transport EV",
    "V2G": "land transport EV",
    "battery charger": "battery",
    "battery discharger": "battery",
    "home battery charger": "home battery",
    "home battery discharger": "home battery",
    "biomass EOP": "biomass"
}

order_carrier = [
    'electricity distribution grid',
    'PHS',
    'DC',
    'AC',
    'nuclear',         # Baseload, always on
    'lignite',         # Baseload, slow to ramp
    'coal',            # Dispatchable, used as baseload or mid-merit
    'oil',             # Peaker, dispatchable, expensive
    'CCGT',            # Flexible, mid-merit
    'OCGT',            # Peaker, fast-ramping
    'biomass',         # Dispatchable, renewable, relatively stable
    'geothermal',      # Very stable, constant output
    'hydro',           # Controllable to some extent, seasonal
    'ror',             # Run-of-river, weather-dependent
    'offwind-dc',      # Offshore wind, variable but often more stable
    'offwind-ac',      # Same as above, different transmission
    'onwind',          # Land-based wind, variable
    'solar',           # Daylight only, weather-sensitive
    'solar rooftop',   # Most variable, decentralized, small-scale
    'battery',
    'home battery',
    'lignite fuel',
    'coal fuel',
    'gas fuel',
    'solid biomass',
]

nice_name_rename = {
    "DC":"DC",
    "AC":"AC",
    "oil":"Oil",
    "gas":"Gas",
    "lignite":"Lignite",
    "coal":"Coal",
    "solar rooftop":"Solar Rooftop",
    'lignite fuel': 'Lignite Fuel',
    'coal fuel': 'Coal Fuel',
    'gas fuel': 'Gas Fuel',
    "biomass":"Biomass",
    "solid biomass":"Solid Biomass",
    "home battery":"Home Battery",
    "electricity distribution grid": "Distribution Grid"
}

country_prefered_order = [
    'MM',
    'TH',
    'MY',
    "Sumatra (ID)",
    "Peninsular (MY)",
    'SG',
    'LA',
    'KH',
    'VN',
    "Java-Bali (ID)",
    'ID',
    "Kalimantan (ID)",
    "Sarawak (MY)",
    'BN',
    "Sabah (MY)",
    "Nusa-Tenggara (ID)",
    'TL'
    "Sulawesi (ID)",
    'PH',
    "Luzon (PH)",
    "Visayas (PH)",
    "Mindanao (PH)",
    "Maluku (ID)",
    "Papua (ID)",
]

subregions = {
    "ID_Java-Bali": "Java-Bali (ID)",
    "ID_Kalimantan": "Kalimantan (ID)",
    "ID_Maluku": "Maluku (ID)",
    "ID_Nusa-Tenggara": "Nusa-Tenggara (ID)",
    "ID_Papua": "Papua (ID)",
    "ID_Sulawesi": "Sulawesi (ID)",
    "ID_Sumatra": "Sumatra (ID)",
    "MY_Peninsular": "Peninsular (MY)",
    "MY_Sabah": "Sabah (MY)",
    "MY_Sarawak": "Sarawak (MY)",
    "PH_Luzon": "Luzon (PH)",
    "PH_Visayas": "Visayas (PH)",
    "PH_Mindanao": "Mindanao (PH)",
    "PH": "Luzon (PH)",
}

def assign_region_or_country(idx):
    for region in subregions.keys():
        if region in idx:
            return subregions[region]
    # If no region matches, use the country code (first two letters before '_')
    return None

def clean_carrier(n):

    append_color = {
        "low voltage": '#110d63',
        "electricity distribution grid": '#97ad8c',
        'rail transport electricity': '#494778', 
        'industry electricity': '#2d2a66',
        'lignite fuel': n.carriers.loc["lignite", "color"],
        'coal fuel': n.carriers.loc["coal", "color"],
        'gas fuel': n.carriers.loc["gas", "color"],
    }

    for key, value in append_color.items():
        n.carriers.loc[key, "color"] = value

    for key, value in nice_name_rename.items():
        n.carriers.loc[key,"nice_name"] = value

    return n

def strip_network(n, scenario=""):

    m = n.copy()
    earth_bus = m.buses[m.buses.location.isin(["Earth"])].index
    m.remove("Bus", earth_bus)
    
    non_dc = m.links[m.links.bus0.isin(earth_bus)].index
    m.remove("Link", non_dc)

    if scenario == "no-AIMS":
        links_project = m.links[(m.links.build_year >= 2025) & ~(m.links.bus0.map(m.buses.country).isin(["ID"]) & m.links.bus1.map(m.buses.country).isin(["ID"]))].index
        lines_project = m.lines[(m.lines.build_year >= 2025) & ~(m.lines.bus0.map(m.buses.country).isin(["ID"]) & m.lines.bus1.map(m.buses.country).isin(["ID"]))].index

        m.remove("Link", links_project)
        m.remove("Line", lines_project)

    if scenario == "no-projects":
        links_project = m.links[m.links.build_year >= 2025].index
        lines_project = m.lines[m.lines.build_year >= 2025].index

        m.remove("Link", links_project)
        m.remove("Line", lines_project)

    return m

def get_energy_balance_bus(n):
    
    eb = n.statistics.energy_balance(
        groupby=["bus", "carrier", "bus_carrier"],
        nice_names = False,
        aggregate_across_components=True,
        comps=["Generator", "Link", "StorageUnit"], # , "Load"
    )
    
    df = pd.DataFrame(eb).reset_index()
    df.loc[df.carrier == "biomass EOP","carrier"] = "biomass"
    mask = df.bus_carrier == "low voltage"
    df.loc[mask, "bus"] = df.loc[mask, "bus"].map(n.buses.location)
    df.loc[mask, "bus_carrier"] = "AC"
    
    df = df[(df.bus_carrier == "AC") & df.carrier.isin(order_carrier)]
    
    eb_new = df.groupby(["bus","carrier"]).sum()[0]

    return eb_new

def get_energy_balance_load(n):
    
    eb = n.statistics.energy_balance(
        groupby=["bus", "bus_carrier"],
        nice_names = False,
        aggregate_across_components=True,
        comps=["Load"],
    )
    
    df = pd.DataFrame(eb).reset_index()
    df["bus"] = df["bus"].map(n.buses.location)
    
    eb_new = df.groupby(["bus"]).sum()[0]

    return - eb_new