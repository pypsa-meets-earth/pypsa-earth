# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import re


def terminate_if_cutout_exists(w):
    """
    Check if any of the requested cutout files exist.
    If that's the case, terminate execution to avoid data loss.

    Tutorial cutouts should be removed once they are no longer needed.
    """
    cutout_fl = "cutouts/" + CDIR + f"{w.cutout}.nc"

    if os.path.exists(cutout_fl) and not config["tutorial"]:
        raise Exception(
            f"An option `build_cutout` or `retrieve_cutout` is enabled, while a cutout file '{cutout_fl}' "
            "still exists and risks to be overwritten. If this is an intended behavior, "
            "please move, rename or delete this file and re-run the rule. Otherwise, "
            "just disable the `build_cutout` and `retrieve_cutout` rule in the config file."
        )

    return []


def memory(w):
    factor = 3.0
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)h$", o, re.IGNORECASE)
        if m is not None:
            factor /= int(m.group(1))
            break
    for o in w.opts.split("-"):
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            factor *= int(m.group(1)) / 8760
            break
    if w.clusters.endswith("m"):
        return int(factor * (18000 + 180 * int(w.clusters[:-1])))
    elif w.clusters.endswith("flex"):
        return int(factor * (18000 + 180 * int(w.clusters[:-4])))
    elif w.clusters == "all":
        return int(factor * (18000 + 180 * 4000))
    elif w.clusters == "min":
        return int(factor * (18000 + 180 * 20))
    else:
        return int(factor * (10000 + 195 * int(w.clusters)))


rule copy_config:
    params:
        summary_dir=config["summary_dir"],
        run=run,
    output:
        folder=directory(SDIR + "configs"),
        config=SDIR + "configs/config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        SDIR + "benchmarks/copy_config"
    script:
        scripts("copy_config.py")
