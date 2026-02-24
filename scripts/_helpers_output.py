# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Set of helper functions to deal with human-redable outputs.
"""

import logging
from pathlib import Path

import pandas as pd
import yaml

logger = logging.getLogger(__name__)

# definitions of output blocks ------------------------------------------------
# model_keys = [
#   "version", "tutorial", "foresight", "countries",
#   "scenario", "snapshots", "enable", "run",
#   "clean_osm_data_options", "build_osm_network", "base_network",
#   "build_shape_options",
#   "electricity", "lines", "links", "transformers",
#   "load_options", "demand_data",
#   "augmented_line_connection",
#   "cluster_options",
#   "atlite", "renewable", "solar_thermal",
#   "policy_config", "fossil_reserves",
#   "sector", "export", "industry"
#   "custom_data", "existing_capacities", "costs",
#   "monte_carlo", "solving",
#   "plotting", "crs",
# ]

major_run_keys = [
    "version",
    "tutorial",
    "foresight",
    "countries",
    "enable",
    "run",
]
major_model_keys = [
    "scenario",
    "snapshots",
]
grid_model_keys = [
    "clean_osm_data_options",
    "build_osm_network",
    "base_network",
]
electrycity_model_keys = ["electricity"]

# visual settings -------------------------------------------------------------
style = """<style type='text/css'>
html {
  font-family: Georgia;
}
r {
  color: #ff0000;
  background-color: coral;
}
g {
  color: #00ff00;
  background-color: coral;
}
b {
  color: #0000ff;
  background-color: grey;
}
white_bg {
  color: #36454F;
  background-color: #FFFFFF;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
blue_bg {
  color: #36454F;
  background-color: #A7C8FF;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
mint_bg {
  color: #36454F;
  background-color: #B7E5D5;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
coral_bg {
  color: #36454F;
  background-color: #FF9E9E;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
lavend_bg {
  color: #36454F;
  background-color: #D8ADEF;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
gold_bg {
  color: #36454F;
  background-color: #F9C6A1;
  width: 100%;
  background-size: contain;
  display:inline-block;
}

blue_hd {
  font-weight: bold;
  color: #36454F;
  background-color: #A7C8FF;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
mint_hd {
  font-weight: bold;
  color: #36454F;
  background-color: #B7E5D5;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
coral_hd {
  font-weight: bold;
  color: #36454F;
  background-color: #FF9E9E;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
lavend_hd {
  font-weight: bold;
  color: #36454F;
  background-color: #D8ADEF;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
gold_hd {
  font-weight: bold;
  color: #36454F;
  background-color: #F9C6A1;
  width: 100%;
  background-size: contain;
  display:inline-block;
}
</style>"""


def make_license_text(html=False):
    # Avoid confusing REUSE as in https://reuse.software/faq/#exclude-lines
    # REUSE-IgnoreStart
    license_str_list = [
        "# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors",
        "# SPDX-License-Identifier: AGPL-3.0-or-later ",
    ]
    # REUSE-IgnoreEnd
    if html:
        license_str = "<br />".join(license_str_list)
    else:
        license_str = "\n\r".join(license_str_list)
    return license_str


def write_html(style, fl, type, str_):
    fl.write("<%(type)s>%(str)s</%(type)s>" % {"type": type, "str": str_})


def write_dict_key(
    header, dict, style_text, style_header, style_def=style, fl_name="out.html"
):

    # define styles to be used in htmls generated below
    white_bg = "white_bg"
    coral_bg = "coral_bg"
    mint_bg = "mint_bg"
    blue_bg = "blue_bg"

    content_str = yaml.dump(dict)

    clean_str = (
        content_str.replace("\t", "     ")
        .replace("\n", "<br />")
        .replace("\r", "<br />")
    )

    f = open(fl_name, "a")

    f.write("<html>")
    f.write(style)

    write_html(
        style=style_def,
        fl=f,
        type=style_text,
        str_="------------------------------------------ <br />",
    )
    write_html(style=style_def, fl=f, type=style_header, str_=header + "<br />")
    write_html(
        style=style_def,
        fl=f,
        type=style_text,
        str_="------------------------------------------ <br />",
    )
    write_html(style=style_def, fl=f, type=style_text, str_=clean_str)

    write_html(style=style_def, fl=f, type=mint_bg, str_="<br />")

    f.write("</html>")


def get_vals(test_dict, key_list):
    df = pd.json_normalize(test_dict, sep=" ")
    d_flat = df.to_dict(orient="records")[0]

    res = [(key, value) for key, value in d_flat.items() if key_list in key.lower()]

    return res


def update_cutout(config, fl="config.yaml"):
    """
    Process configuration inputs for a cutout file to be used in a model.
    Updates the Snakemake config variable according to the parameters
    set by a user in the actual configuration file for a simulation accounting
    for the need to provide multiple cutout entries to the Snakemake config
    variable.

    In case the custom configuration file contains only one cutout entry, this
    entry will be used to update the Snakemake config variable for all the
    cutout entries.

    If the custom configuration file contains multiple entries for a cutout
    they will be merged with all the other cutout entries with a normally-used
    Snakemake approach.

    Parameters
    ----------
    config : dictionary
        A dictionary which corresponds to the Snakemake config variable
        to be used in a simulation
    fl : string
        A name of the yaml file which contains custom user-defined configuration
        parameters to be prioritised over default parameters
    """
    snakem_config_cutouts = [
        d_value.get("cutout") for tc, d_value in config.get("renewable").items()
    ] + list(config.get("atlite").get("cutouts"))

    # A user-defined configuration file is loaded directly to avoid merges
    # made by Snakemake
    local_config = yaml.safe_load(Path(fl).read_text())
    local_config_cutouts = [
        d_value.get("cutout")
        # "renewable" and "atlite"" can be missing from the the local config
        for tc, d_value in local_config.get("renewable", {}).items()
    ] + list(local_config.get("atlite", {}).get("cutouts", {}))
    # A local config file can miss some or all keys
    local_config_cutouts = [x for x in local_config_cutouts if x is not None]

    # A user must have an opportunity to over-write all the cutouts
    # with a single value defined in `config.yaml`
    if len(set(local_config_cutouts)) == 1:

        cutout_name = local_config_cutouts[0]

        setup_cutout_dict = {
            "atlite": {"cutouts": cutout_name},
            # TODO Read keys from the config
            "renewable": {
                "onwind": {"cutout": cutout_name},
                "offwind-ac": {"cutout": cutout_name},
                "offwind-dc": {"cutout": cutout_name},
                "solar": {"cutout": cutout_name},
                "hydro": {"cutout": cutout_name},
                "csp": {"cutout": cutout_name},
            },
        }

        config.update(setup_cutout_dict)

    return config


def check_config_keys(config, fl="config.yaml"):
    """
    Check if all the keys of `config.yaml` present in the default configs

    Parameters
    ----------
    config : dictionary
        A dictionary which corresponds to the Snakemake config variable
        to be used in a simulation
    fl : string
        A name of the yaml file which contains custom user-defined configuration
        parameters to be prioritised over default parameters
    """
    config.keys()

    # A user-defined configuration file is loaded directly to avoid merges
    # made by Snakemake
    local_config = yaml.safe_load(Path(fl).read_text())

    config_default_flatten = pd.json_normalize(config, sep=": ")
    config_local_flatten = pd.json_normalize(local_config, sep=": ")

    # elements of index not in other
    # config_local_flatten - config_default_flatten
    config_diff = config_local_flatten.keys().difference(config_default_flatten.keys())

    if not config_diff.empty:
        config_discrepancy_string = "<br />    - ".join(config_diff) + "<br />"
    else:
        config_discrepancy_string = None

    return config_discrepancy_string


def parse_config(config, fl_name=None, style_def=style):
    """
    Outputs a config summary to highlight points which require attention
    from a modeler

    Currently, the following checks are implemented:
    1. List all the year-related parameters used in a simulation
    2. Checks which cutouts are used in a simulation
    3. Check if any keys from `config.yaml` are missed from /config/*.yaml
    """

    year_items = dict(get_vals(config, "year"))

    # 1. Outline all the years used in a simulation
    # The config contains
    year_config_params = [None] * len(year_items)
    k = 0
    for i, j in year_items.items():
        year_config_params[k] = str(i) + ": " + str(j)
        k += 1

    years_html_string = "<br /> - ".join(year_config_params) + "<br />"

    # 2. Extract names of the cutouts used
    # cutouts_used = [
    #     d_value.get("cutout") for tc, d_value in config.get("renewable").items()
    # ] + list(config.get("atlite").get("cutouts"))

    # Different cutout names can be defined in `atlite` section
    # to make possible using different weather datasets
    # Merging of the configs can lead to adding more cutouts
    # without an intention to do so and having a check may be a good idea
    cutouts_atlite = [k for k in config.get("atlite").get("cutouts").keys()]

    cutouts_renewable = [
        d_value.get("cutout") for tc, d_value in config.get("renewable").items()
    ]

    cutout_atlite_string = (
        "<br />" + "- to build new cutouts: " + " ".join(set(cutouts_atlite)) + "<br />"
    )
    cutout_renew_string = (
        "<br />"
        + "- to evaluate the renewable potential: "
        + " ".join(set(cutouts_renewable))
        + "<br />"
    )

    # 3. Check if there is any discrepancy between the default config
    # and a local one
    structure_check_string = check_config_keys(config=config, fl="config.yaml")

    # define styles to be used in htmls generated below
    white_bg = "white_bg"
    coral_bg = "coral_bg"
    mint_bg = "mint_bg"
    blue_bg = "blue_bg"

    f = open("config_check.html", "w")

    f.write("<html>")
    f.write(style)

    write_html(
        style="white_bg",
        fl=f,
        type="white_bg",
        str_="------------------------------------------ <br />",
    )
    f.write(make_license_text(html=True))
    write_html(
        style="white_bg",
        fl=f,
        type="white_bg",
        str_="------------------------------------------ <br />",
    )

    if structure_check_string is not None:
        write_html(
            style=coral_bg,
            fl=f,
            type=coral_bg,
            str_="------------------------------------------ <br />",
        )
        write_html(
            style=coral_bg,
            fl=f,
            type=coral_bg,
            str_="The following configuration parameters do not present in config.yaml:<br />",
        )
        write_html(
            style=coral_bg,
            fl=f,
            type=coral_bg,
            str_=structure_check_string,
        )
        write_html(
            style=coral_bg,
            fl=f,
            type=coral_bg,
            str_="------------------------------------------ <br />",
        )

    write_html(
        style=mint_bg,
        fl=f,
        type=mint_bg,
        str_="------------------------------------------ <br />",
    )
    write_html(
        style=mint_bg,
        fl=f,
        type=mint_bg,
        str_="The current configuration contains the following year definitions:<br />",
    )
    write_html(
        style=mint_bg,
        fl=f,
        type=mint_bg,
        str_=years_html_string,
    )
    write_html(
        style=mint_bg,
        fl=f,
        type=mint_bg,
        str_="------------------------------------------ <br />",
    )
    write_html(
        style=blue_bg,
        fl=f,
        type=blue_bg,
        str_="------------------------------------------ <br />",
    )
    write_html(
        style=blue_bg,
        fl=f,
        type=blue_bg,
        str_="The current configuration is based on using the following cutouts:<br />",
    )
    write_html(
        style=blue_bg,
        fl=f,
        type=blue_bg,
        str_=cutout_atlite_string,
    )
    write_html(
        style=blue_bg,
        fl=f,
        type=blue_bg,
        str_=cutout_renew_string,
    )
    write_html(
        style=blue_bg,
        fl=f,
        type=blue_bg,
        str_="------------------------------------------ <br />",
    )
