# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:15:15 2021

@author: haz43975
"""
import matplotlib.pyplot as plt
import pypsa
import os
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from matplotlib.patches import Circle, Ellipse
from matplotlib.legend_handler import HandlerPatch




base_path = os.path.dirname(os.path.realpath(__file__)) 


network_path = os.path.join(base_path, 'Networks')

network = pypsa.Network(os.path.join(
            "../results/test/postnetworks/elec_s_4_2050.nc"))

n = network.copy()

bus_size_factor = 1e5 #Def 1e5
linewidth_factor = 2e4 #Def 1e4
line_lower_threshold = 1e2 # MW below which not drawn. Def 1e3

DC_lines=n.links[n.links.carrier=='DC']

n.links = n.links[n.links.carrier=='H2 pipeline']
n.links.bus0 = n.links.bus0.str.replace(" H2", "")
n.links.bus1 = n.links.bus1.str.replace(" H2", "")

n.lines.append(DC_lines[['bus0', 'bus1']])

n.madd('Line', names=DC_lines.index, bus0=DC_lines.bus0, bus1=DC_lines.bus1)

fig = plt.figure()
fig.set_size_inches(10.5, 9)

n.plot(branch_components=["Link", "Line"], boundaries=(-20, 0, 25, 40),
       color_geomap={'ocean': 'lightblue', 'land': "oldlace"}, line_colors='darkblue',
       link_colors='turquoise', link_widths=3, bus_sizes=0.03, bus_colors='red', line_widths=1)

#Legend
Elec_Circle = plt.Line2D([0], [0], marker='o', color='darkblue',
label='Clustered node', markerfacecolor='red',
markersize=10)
elec_Line = plt.Line2D([0], [0], marker='_', color='darkblue', 
            label='Electricity Lines', markerfacecolor='w', markersize=16, lw=4)

H2_Line = plt.Line2D([0], [0], marker='_', color='turquoise', 
        label='H2 Pipeline',markerfacecolor='w', markersize=16, lw=4)

plt.legend(handles=[Elec_Circle, elec_Line, H2_Line], loc='upper left')
plt.savefig('networks.png', transparent=True,
            bbox_inches="tight", dpi=300)
plt.savefig('networks.pdf', transparent=True,
            bbox_inches="tight", dpi=300)
# #%%
# h2_emap_filled = pd.read_csv('h2_emap_filled.csv')

# h2_emap_filled_ng = h2_emap_filled[h2_emap_filled.color=='chocolate']
# h2_emap_filled_h = h2_emap_filled[h2_emap_filled.color=='cyan']

# n.mremove('Line', names=n.lines.index.tolist())
# n.mremove('Link', names=n.links.index.tolist())

# n.madd('Link', names=h2_emap_filled_h.index, bus0=h2_emap_filled_h.bus0, bus1=h2_emap_filled_h.bus1)
# n.madd('Line', names=h2_emap_filled_ng.index, bus0=h2_emap_filled_ng.bus0, bus1=h2_emap_filled_ng.bus1)

# n.plot(branch_components=["Link", "Line"], boundaries=(-10, 30, 34, 70),
#        color_geomap={'ocean': 'lightblue', 'land':'oldlace'}, line_colors='chocolate',
#        link_colors='cyan', link_widths=1.5, bus_sizes=0.03, bus_colors='darkblue', line_widths=1.5)

# Elec_Circle = plt.Line2D([0], [0], marker='o', color='darkblue',
#                          label='Clustered node', markerfacecolor='darkblue',
#                          markersize=10)

# elec_Line = plt.Line2D([0], [0], marker='_', color='chocolate', label='Based on NG grid', markerfacecolor='w', markersize=16, lw=4)
# H2_Line = plt.Line2D([0], [0], marker='_', color='cyan', label='Based on heuristics',markerfacecolor='w', markersize=16, lw=4)
# plt.legend(handles=[Elec_Circle, elec_Line, H2_Line], loc='upper left')

#%%
def make_handler_map_to_scale_circles_as_in(ax, dont_resize_actively=False):
    fig = ax.get_figure()

    def axes2pt():
        return np.diff(ax.transData.transform([(0, 0), (1, 1)]), axis=0)[
            0] * (72. / fig.dpi)

    ellipses = []
    if not dont_resize_actively:
        def update_width_height(event):
            dist = axes2pt()
            for e, radius in ellipses:
                e.width, e.height = 2. * radius * dist
        fig.canvas.mpl_connect('resize_event', update_width_height)
        ax.callbacks.connect('xlim_changed', update_width_height)
        ax.callbacks.connect('ylim_changed', update_width_height)

    def legend_circle_handler(legend, orig_handle, xdescent, ydescent,
                              width, height, fontsize):
        w, h = 2. * orig_handle.get_radius() * axes2pt()
        e = Ellipse(xy=(0.5 * width - 0.5 * xdescent, 0.5 *
                        height - 0.5 * ydescent), width=w, height=w)
        ellipses.append((e, orig_handle.get_radius()))
        return e
    return {Circle: HandlerPatch(patch_func=legend_circle_handler)}


def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0, 0), radius=(s / scale)**0.5, **kw) for s in sizes]

def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):

        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)

        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue

            names = ifind.index[ifind == i]

            c.df.loc[names, 'location'] = names.str[:i]

#%%
# def plot_h2_map(network):

n = network.copy()

assign_location(n)

bus_size_factor = 1e5
linewidth_factor = 1e4
# MW below which not drawn
line_lower_threshold = 1e1
bus_color = "m"
link_color = "c"

# Drop non-electric buses so they don't clutter the plot
n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

elec = n.links.index[n.links.carrier == "H2 Electrolysis"]

bus_sizes = n.links.loc[elec,"p_nom_opt"].groupby(n.links.loc[elec,"bus0"]).sum() / bus_size_factor

# make a fake MultiIndex so that area is correct for legend
bus_sizes.index = pd.MultiIndex.from_product(
    [bus_sizes.index, ["electrolysis"]])

n.links.drop(n.links.index[n.links.carrier != "H2 pipeline"], inplace=True)

link_widths = n.links.p_nom_opt / linewidth_factor
link_widths[n.links.p_nom_opt < line_lower_threshold] = 0.

n.links.bus0 = n.links.bus0.str.replace(" H2", "")
n.links.bus1 = n.links.bus1.str.replace(" H2", "")



print(link_widths.sort_values())

print(n.links[["bus0", "bus1"]])

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})

fig.set_size_inches(10.5, 9)

n.plot(bus_sizes=bus_sizes,
        bus_colors={"electrolysis": bus_color},
        link_colors=link_color,
        link_widths=link_widths,
        branch_components=["Link"],
        color_geomap={'ocean': 'lightblue', 'land': "oldlace"},
        ax=ax,  boundaries=(-20, 0, 25, 40))

handles = make_legend_circles_for(
        [50000, 10000], scale=bus_size_factor, facecolor=bus_color)
labels = ["{} GW".format(s) for s in (50, 10)]
l2 = ax.legend(handles, labels,
                   loc="upper left", bbox_to_anchor=(0.01, 1.01),
                   labelspacing=0.8,
                   framealpha=1.,
                   title='Electrolyzer capacity',
                   handler_map=make_handler_map_to_scale_circles_as_in(ax))
ax.add_artist(l2)

handles = []
labels = []

for s in (50, 10):
    handles.append(plt.Line2D([0], [0], color=link_color,
                              linewidth=s * 1e3 / linewidth_factor))
    labels.append("{} GW".format(s))
l1_1 = ax.legend(handles, labels,
                  loc="upper left", bbox_to_anchor=(0.32, 1.01),
                  framealpha=1,
                  labelspacing=0.8, handletextpad=1.5,
                  title='H2 pipeline capacity')
ax.add_artist(l1_1)

fig.savefig('h2.png', transparent=True,
            bbox_inches="tight", dpi=300)
fig.savefig('h2.pdf', transparent=True,
            bbox_inches="tight", dpi=300)
# plot_h2_map(n)


#%%###################

preferred_order = pd.Index(["transmission lines","hydroelectricity","hydro reservoir","run of river","pumped hydro storage","solid biomass","biogas","onshore wind","offshore wind","offshore wind (AC)","offshore wind (DC)","solar PV","solar thermal","solar","building retrofitting","ground heat pump","air heat pump","heat pump","resistive heater","power-to-heat","gas-to-power/heat","CHP","OCGT","gas boiler","gas","natural gas","helmeth","methanation","hydrogen storage","power-to-gas","power-to-liquid","battery storage","hot water storage","CO2 sequestration"])


def rename_techs(label):

    prefix_to_remove = ["residential ","services ","urban ","rural ","central ","decentral "]

    rename_if_contains = ["CHP","gas boiler","biogas","solar thermal","air heat pump","ground heat pump","resistive heater","Fischer-Tropsch"]

    rename_if_contains_dict = {"water tanks" : "hot water storage",
                               "retrofitting" : "building retrofitting",
                               "H2" : "hydrogen storage",
                               "battery" : "battery storage",
                               "CCS" : "CCS"}

    rename = {"solar" : "solar PV",
              "Sabatier" : "methanation",
              "offwind" : "offshore wind",
              "offwind-ac" : "offshore wind (AC)",
              "offwind-dc" : "offshore wind (DC)",
              "onwind" : "onshore wind",
              "ror" : "hydroelectricity",
              "hydro" : "hydroelectricity",
              "PHS" : "hydroelectricity",
              "co2 Store" : "DAC",
              "co2 stored" : "CO2 sequestration",
              "AC" : "transmission lines",
              "DC" : "transmission lines",
              "B2B" : "transmission lines"}

    for ptr in prefix_to_remove:
        if label[:len(ptr)] == ptr:
            label = label[len(ptr):]

    for rif in rename_if_contains:
        if rif in label:
            label = rif

    for old,new in rename_if_contains_dict.items():
        if old in label:
            label = new

    for old,new in rename.items():
        if old == label:
            label = new
    return label


def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech or "resistive heater" in tech:
        return "power-to-heat"
    elif tech in ["methanation", "hydrogen storage", "helmeth"]:
        return "power-to-gas"
    elif tech in ["OCGT", "CHP", "gas boiler"]:
        return "gas-to-power/heat"
    elif "solar" in tech:
        return "solar"
    elif tech == "Fischer-Tropsch":
        return "power-to-liquid"
    elif "offshore wind" in tech:
        return "offshore wind"
    else:
        return tech
#%%
def plot_map(network, components=["links", "generators"],
             bus_size_factor=1.7e10, transmission=False):

    n = network.copy()
    assign_location(n)
    # Drop non-electric buses so they don't clutter the plot
    n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

    costs = pd.DataFrame(index=n.buses.index)

    for comp in components:
        df_c = getattr(n, comp)
        df_c["nice_group"] = df_c.carrier.map(rename_techs_tyndp)

        attr = "e_nom_opt" if comp == "stores" else "p_nom_opt"

        costs_c = ((df_c.capital_cost * df_c[attr])
                   .groupby([df_c.location, df_c.nice_group]).sum()
                   .unstack().fillna(0.))
        costs = pd.concat([costs, costs_c], axis=1)

        print(comp, costs)
    costs = costs.groupby(costs.columns, axis=1).sum()

    costs.drop(list(costs.columns[(costs == 0.).all()]), axis=1, inplace=True)

    new_columns = ((preferred_order & costs.columns)
                   .append(costs.columns.difference(preferred_order)))
    costs = costs[new_columns]

    for item in new_columns:
        if item not in tech_colors:
            print("Warning!",item,"not in config/plotting/tech_colors")

    costs = costs.stack()  # .sort_index()

    # hack because impossible to drop buses...
    #n.buses.loc["EU gas", ["x", "y"]] = n.buses.loc["MA 0", ["x", "y"]]

    n.links.drop(n.links.index[(n.links.carrier != "DC") & (
        n.links.carrier != "B2B")], inplace=True)

    # drop non-bus
    to_drop = costs.index.levels[0] ^ n.buses.index
    if len(to_drop) != 0:
        print("dropping non-buses", to_drop)
        costs.drop(to_drop, level=0, inplace=True, axis=0)

    # make sure they are removed from index
    costs.index = pd.MultiIndex.from_tuples(costs.index.values)

    # PDF has minimum width, so set these to zero
    line_lower_threshold = 500.
    line_upper_threshold = 1e4
    linewidth_factor = 2e3
    ac_color = "gray"
    dc_color = "m"

    #if snakemake.wildcards["lv"] == "1.0":
        # should be zero
    line_widths = n.lines.s_nom_opt - n.lines.s_nom
    link_widths = n.links.p_nom_opt - n.links.p_nom
    title = "Technologies"

    if transmission:
        line_widths = n.lines.s_nom_opt
        link_widths = n.links.p_nom_opt
        linewidth_factor = 2e3
        line_lower_threshold = 0.
        title = "Technologies"
    # else:
    #     line_widths = n.lines.s_nom_opt - n.lines.s_nom_min
    #     link_widths = n.links.p_nom_opt - n.links.p_nom_min
    #     title = "Transmission reinforcement"

    #     if transmission:
    #         line_widths = n.lines.s_nom_opt
    #         link_widths = n.links.p_nom_opt
    #         title = "Total transmission"

    line_widths[line_widths < line_lower_threshold] = 0.
    link_widths[link_widths < line_lower_threshold] = 0.

    line_widths[line_widths > line_upper_threshold] = line_upper_threshold
    link_widths[link_widths > line_upper_threshold] = line_upper_threshold

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()})
    fig.set_size_inches(10.5, 9)

    n.plot(bus_sizes=costs / bus_size_factor,
           bus_colors=tech_colors,
           line_colors=ac_color,
           link_colors=dc_color,
           line_widths=line_widths / linewidth_factor,
           link_widths=link_widths / linewidth_factor,
           ax=ax,  boundaries=(-20, 0, 25, 40),
           color_geomap={'ocean': 'lightblue', 'land': "oldlace"})

    handles = make_legend_circles_for(
        [5e9, 1e9], scale=bus_size_factor, facecolor="gray")
    labels = ["{} b€/a".format(s) for s in (5, 1)]
    l2 = ax.legend(handles, labels,
                   loc="upper left", bbox_to_anchor=(0.33, 1.005),
                   labelspacing=1.0,
                   framealpha=1.,
                   title='System cost',
                   fontsize=12,
                   handler_map=make_handler_map_to_scale_circles_as_in(ax))
    ax.add_artist(l2)

    handles = []
    labels = []

    for s in list(plot_labeles.keys()):
        handles.append(plt.Line2D([0], [0], color=tech_colors[s],
                                  linewidth=5))
        labels.append("{}".format(s))

    l1_1 = ax.legend(handles, labels,
                      loc="upper left", bbox_to_anchor=(0.001, 1.002),
                      framealpha=1,
                      labelspacing=0.4, handletextpad=1.5,
                      fontsize=10)

    ax.add_artist(l1_1)
    import matplotlib.patches as mpatches

    # red_patch = mpatches.Patch(color='red', label='The red data')
    # plt.legend(handles=[red_patch])

    plt.show()
    fig.savefig('plot_map.pdf', transparent=True,
                bbox_inches="tight")
    fig.savefig('plot_map.png', transparent=True,
            bbox_inches="tight", dpi=300)
    
plot_labeles={
    'onshore wind': 'b',
    'offshore wind': 'c',
    'hydroelectricity' :'',
    'solar': 'y',
    'power-to-gas':'#FF1493',
    "gas-to-power/heat" : "orange",
    'power-to-heat':'',
    'power-to-liquid':'',
    'DAC':'',
    'electricity distribution grid':''
    }


tech_colors={
    "CO2 pipeline" : "gray",
    "onwind" : "b",
    "Africa gas" : "b",
    "Africa oil" : "b",
    "onshore wind" : "dodgerblue",                                               #needed
    'offwind' : "c",
    'offshore wind' : "c",                                                    #needed
    'offwind-ac' : "c",
    'offshore wind (AC)' : "b",
    'offwind-dc' : "b",
    'offshore wind (DC)' : "b",
    'wave' : "b",
    "hydro" : "b",
    "hydro reservoir" : "b",
    "ror" : "b",
    "run of river" : "b",
    'hydroelectricity' : 'blue',                                                                #needed
    'solar' : "orange",                                                         #needed
    'solar PV' : "red",
    'solar thermal' : 'b',
    'solar rooftop' : 'b',
    "OCGT" : "wheat",
    "OCGT marginal" : "b",
    "OCGT-heat" : "b",
    "gas boiler" : "b",
    "gas boilers" : "b",
    "gas boiler marginal" : "b",
    "gas-to-power/heat" : "brown",                                               #needed
    "gas" : "b",
    "natural gas" : "b",
    "SMR" : "#4F4F2F",
    "oil" : "red",
    "oil boiler" : "#B5A677",
    "lines" : "k",
    "transmission lines" : "k",
    "H2" : "m",
    "hydrogen storage" : "m",
    "battery" : "slategray",
    "battery storage" : "slategray",
    "home battery" : "b",
    "home battery storage" : "b",
    "Nuclear" : "b",
    "Nuclear marginal" : "b",
    "nuclear" : "b",
    "uranium" : "b",
    "Coal" : "k",
    "coal" : "k",
    "Coal marginal" : "k",
    "Lignite" : "grey",
    "lignite" : "grey",
    "Lignite marginal" : "grey",
    "CCGT" : "b",
    "CCGT marginal" : "orange",
    "heat pumps" : "b", #light green
    "heat pump" : "b", #light green
    "air heat pump" : "b", #light green
    "ground heat pump" : "b", #dark green
    "power-to-heat" : "red",                                                   #needed
    "resistive heater" : "b",
    "Sabatier" : "#b",
    "methanation" : "#b",
    "power-to-gas" : "purple",                                                       #needed
    "power-to-liquid" : "darkgreen",                                            #needed
    "helmeth" : "b",
    "helmeth" : "b",
    "DAC" : "deeppink",                                                           #needed
    "co2 stored" : "b",
    "CO2 sequestration" : "b",
    "CCS" : "b",
    "co2" : "b",
    "co2 vent" : "b",
    "solid biomass for industry co2 from atmosphere" : "b",
    "solid biomass for industry co2 to stored": "b",
    "gas for industry co2 to atmosphere": "b",
    "gas for industry co2 to stored": "b",
    "Fischer-Tropsch" : "b",
    "kerosene for aviation": "b",
    "naphtha for industry" : "b",
    "water tanks" : "b",
    "hot water storage" : "b",
    "hot water charging" : "b",
    "hot water discharging" : "b",
    "CHP" : "b",
    "CHP heat" : "b",
    "CHP electric" : "b",
    "PHS" : "b",
    "Ambient" : "b",
    "Electric load" : "b",
    "Heat load" : "b",
    "Transport load" : "b",
    "heat" : "b",
    "rural heat" : "b",
    "central heat" : "b",
    "decentral heat" : "b",
    "low-temperature heat for industry" : "b",
    "process heat" : "b",
    "heat demand" : "b",
    "electric demand" : "b",
    "Li ion" : "b",
    "district heating" : "b",
    "retrofitting" : "b",
    "building retrofitting" : "b",
    "BEV charger" : "b",
    "V2G" : "b",
    "transport" : "b",
    "electricity" : "b",
    "gas for industry" : "b",
    "solid biomass for industry" : "b",
    "industry electricity" : "b",
    "industry new electricity" : "b",
    "process emissions to stored" : "b",
    "process emissions to atmosphere" : "b",
    "process emissions" : "b",
    "transport fuel cell" : "b",
    "biogas" : "b",
    "solid biomass" : "b",
    "today" : "b",
    "shipping" : 'b',
    "electricity distribution grid" : "y",                                  #needed
    'SMR CC': 'black',
    'gas': 'black',
    'process emissions CC': 'yellow',
    'biomass for industry CC': 'yellow',
    'gas for industry CC':'b',
    'solid biomass for industry CC': 'b'}           

nice_names={
    # OCGT: "Gas",
    # OCGT marginal: "Gas (marginal)",
    'offwind': "offshore wind",
    'onwind': "onshore wind",
    'battery': "Battery storage",
    'lines': "Transmission lines",
    'AC line': "AC lines",
    'AC-AC': "DC lines",
    'ror': "Run of river"}

nice_names_n={
    'offwind': "offshore\nwind",
    'onwind': "onshore\nwind",
    # OCGT: "Gas"
    'H2': "Hydrogen\nstorage",
    # OCGT marginal: "Gas (marginal)"
    'lines': "transmission\nlines",
    'ror': "run of river"}

plot_map(network, transmission=True)
#%%



dfc=pd.DataFrame(data=tech_colors, index=['color']).T
dfc['length'] = 255

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.barh(y=[50]*len(list(dfc.iloc[:30, 1])), width=0.052, color=dfc.color.iloc[:30])
ax1.set_yticklabels(list(dfc.iloc[:30, 1].index))
# ax2=dfc.iloc[30:60, 1].plot.barh(color=dfc.iloc[30:60].color)
# ax3=dfc.iloc[60:, 1].plot.barh(color=dfc.iloc[60:].color)