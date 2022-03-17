#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 15:50:42 2022

@author: user
"""
import matplotlib.pyplot as plt
import pypsa
import os
import cartopy.crs as ccrs
import pandas as pd
import numpy as np
from matplotlib.patches import Circle, Ellipse
from matplotlib.legend_handler import HandlerPatch


#%%
def assign_location(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):

        ifind = pd.Series(c.df.index.str.find(" ", start=4), c.df.index)

        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue

            names = ifind.index[ifind == i]

            c.df.loc[names, 'location'] = names.str[:i]

network = pypsa.Network(os.path.join(
            "../results/test/postnetworks/elec_s_4_2050.nc"))
# def plot_h2_map(network):

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
#%%
n = network.copy()

assign_location(n)

bus_size_factor = 1e5
linewidth_factor = 1e4
# MW below which not drawn
line_lower_threshold = 1e1
bus_color = "m"
link_color = "c"

n.links.loc[:,'p_nom_opt']=n.links.loc[:,'p_nom_opt']*1e7

# n.links.loc[n.links.carrier == "H2 Electrolysis"].p_nom_opt

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