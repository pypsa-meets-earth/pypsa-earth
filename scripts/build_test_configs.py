# SPDX-FileCopyrightText: : 2022 PyPSA-Earth Authors
#
# coding: utf-8
"""
Write option files (configs) for the Continous Integration tests

The config.tutorial.yaml has all options.
The test/* config files have only key/value strings that are different from the tutorial config.
The below scripts 'updates' the test configs and adds all options of the tutorial config.
"""
import os
import collections.abc
from pathlib import Path
from ruamel.yaml import YAML


def update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

    #Input paths
    fpbase = Path(Path.cwd(), snakemake.tutorial)
    fp0 = Path(Path.cwd(), snakemake.standard)
    fp1 = Path(Path.cwd(), snakemake.custom)
    fp2 = Path(Path.cwd(), snakemake.monte_carlo)
    fp3 = Path(Path.cwd(), snakemake.landlock)

    #Output paths
    fp0 = Path(Path.cwd(), snakemake.standard)
    fp1 = Path(Path.cwd(), snakemake.custom)
    fp2 = Path(Path.cwd(), snakemake.monte_carlo)
    fp3 = Path(Path.cwd(), snakemake.landlock)

    #Load yaml files
    yaml = YAML()
    with open(fptutorial) as fp:
        data_tutorial = yaml.load(fp)   
    with open(fp0) as fp:
        data0 = yaml.load(fp)
    with open(fp1) as fp:
        data1 = yaml.load(fp)
    with open(fp2) as fp:
        data2 = yaml.load(fp)
    with open(fp3) as fp:
        data3 = yaml.load(fp)
    
    #Modify yamls
    a = update(data_tutorial, data0)
    b = update(data_tutorial, data1)
    c = update(a, data2)
    d = update(a, data3)

    #Save files
    yaml.dump(a, fp0)
    yaml.dump(b, fp1)
    yaml.dump(c, fp2)
    yaml.dump(d, fp3)

    #Manual output in terminal
    # import sys
    # yaml.dump(data, sys.stdout)
