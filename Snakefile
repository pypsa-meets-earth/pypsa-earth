# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import os
import sys
import warnings

sys.path.append("./scripts")

from pathlib import Path
from shutil import copyfile, move, unpack_archive

from _helpers import branch  # Remove if Snakemake >= 8.3.0
from _helpers import (
    BASE_DIR,
    check_config_version,
    content_retrieve,
    copy_default_files,
    create_country_list,
    get_last_commit_message,
    migrate_config,
    update_cutout_config,
    script_path_provider,
)
from build_demand_profiles import get_load_paths_gegis

copy_default_files()


configfile: "config.default.yaml"
configfile: "configs/plotting.default.yaml"
configfile: "configs/solving.default.yaml"
configfile: "configs/bundle_config.yaml"
configfile: "configs/powerplantmatching_config.yaml"
configfile: "config.yaml"


check_config_version(config=config)

config = migrate_config(config)

config.update({"git_commit": get_last_commit_message(".")})

# convert country list according to the desired region
config["countries"] = create_country_list(config["countries"])

print(
    "The PyPSA meets Earth initiative also supports dedicated regional models. "
    "See the documentation at "
    "https://pypsa-earth.readthedocs.io/en/latest/user-guide/customization/basic-setup/"
)

# create a list of iteration steps, required to solve the experimental design
# each value is used as wildcard input e.g. solution_{unc}
config["scenario"]["unc"] = [
    f"m{i}" for i in range(config["monte_carlo"]["options"]["samples"])
]

config = update_cutout_config(config)

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""
CDIR = RDIR if not run.get("shared_cutouts") else ""
SECDIR = run["sector_name"] + "/" if run.get("sector_name") else ""
SDIR = config["summary_dir"].strip("/") + f"/{SECDIR}"
RESDIR = config["results_dir"].strip("/") + f"/{SECDIR}"

ATLITE_NPROCESSES = config["atlite"].get("nprocesses", 4)

scripts = script_path_provider(Path(BASE_DIR))


wildcard_constraints:
    simpl="[a-zA-Z0-9]*|all",
    clusters="[0-9]+(m|flex)?|all|min",
    ll=r"(v|c|l)([0-9\.]+|opt|all)|all",
    opts=r"[-+a-zA-Z0-9\.]*",
    unc=r"[-+a-zA-Z0-9\.]*",
    sopts=r"[-+a-zA-Z0-9\.\s]*",
    discountrate=r"[-+a-zA-Z0-9\.\s]*",
    planning_horizons="20[2-9][0-9]|2100",


if config["custom_rules"] is not []:
    for rule in config["custom_rules"]:

        include: rule


include: "rules/common.smk"
include: "rules/retrieve.smk"
include: "rules/build_electricity.smk"
include: "rules/build_sector.smk"
include: "rules/solve.smk"
include: "rules/postprocess.smk"
include: "rules/collect.smk"


rule clean:
    run:
        try:
            shell("snakemake -j 1 solve_all_networks --delete-all-output")
        except:
            shell("snakemake -j 1 solve_all_networks_monte --delete-all-output")
            pass
        shell("snakemake -j 1 run_all_scenarios --delete-all-output")
