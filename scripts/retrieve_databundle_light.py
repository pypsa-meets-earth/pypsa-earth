# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors, 2021 PyPSA-Africa
#
# SPDX-License-Identifier: MIT
"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517935.svg
   :target: https://doi.org/10.5281/zenodo.3517935

The data bundle (1.4 GB) contains common GIS datasets like NUTS3 shapes, EEZ shapes, CORINE Landcover, Natura 2000 and also electricity specific summary statistics like historic per country yearly totals of hydro generation, GDP and POP on NUTS3 levels and per-country load time-series.

This rule downloads the data bundle from `zenodo <https://doi.org/10.5281/zenodo.3517935>`_ and extracts it in the ``data`` sub-directory, such that all files of the bundle are stored in the ``data/bundle`` subdirectory.

The :ref:`tutorial` uses a smaller `data bundle <https://zenodo.org/record/3517921/files/pypsa-eur-tutorial-data-bundle.tar.xz>`_ than required for the full model (19 MB)

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517921.svg
    :target: https://doi.org/10.5281/zenodo.3517921

**Relevant Settings**

.. code:: yaml

    tutorial:

.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``cutouts/bundle``: input data collected from various sources

"""
import logging
import os
import tarfile
from pathlib import Path
from numpy import False_
import yaml
from zipfile import ZipFile
import re

from _helpers import _sets_path_to_root
from _helpers import configure_logging
from _helpers import progress_retrieve
from download_osm_data import create_country_list
from google_drive_downloader import GoogleDriveDownloader as gdd

logger = logging.getLogger(__name__)

def load_databundle_config(path):
    "Load databundle configurations from path file"
    with open(path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # parse the "countries" list specified in the file before processing
    for bundle_name in config:
        config[bundle_name]["countries"] = create_country_list(
                                                config[bundle_name]["countries"],
                                                iso_coding=False
                                            )
    
    return config

def download_and_unzip(host, config, rootpath, dest_path):
    """
    Function to download and unzip data depending on the hosting platform.
    Currently, hosts accepted: zenodo and google
    """
    resource="-".join(config["category"])
    file_path=Path(rootpath, "tempfile.zip")

    if host=="zenodo":
        url=config["urls"]["zenodo"]
        progress_retrieve(url, file_path)
        logger.info(f"Extracting resources")
        with ZipFile(file_path, "r") as zipObj:
            # Extract all the contents of zip file in current directory
            zipObj.extractall(path=dest_path)
        os.remove(file_path)
        logger.info(f"Download resource '{resource}' from cloud '{url}'.")
        return True
    elif host=="google":

        url=config["urls"]["google"]
        # retrieve file_id from path
        partition_view = re.split(r"/view|\\view", str(url), 1)  # cut the part before the ending \view
        if len(partition_view) < 2:
            logger.error(f"Resource {resource} cannot be downloaded: \"\\view\" not found in url {url}")
            return False
        
        code_split = re.split(r"\\|/", partition_view[0])  # split url to get the file_id

        if len(code_split) < 2:
            logger.error(f"Resource {resource} cannot be downloaded: character \"\\\" not found in {partition_view[0]}")
            return False
        
        # get file id
        file_id = code_split[-1]

        if os.path.exists(file_path):
            os.remove(file_path)
        gdd.download_file_from_google_drive(
            file_id=file_id,
            dest_path=file_path,
            showsize=True,
            unzip=False,
        )
        with ZipFile(file_path, "r") as zipObj:
            # Extract all the contents of zip file in current directory
            zipObj.extractall(path=dest_path)
        os.remove(file_path)
        logger.info(f"Download resource '{resource}' from cloud '{url}'.")

        return True
    else:
        logger.error(f"Host {host} not implemented")
        return False

def get_best_bundle(country_list, category, config_bundles, tutorial):
    # dictionary with the number of match by configuration for tutorial/non-tutorial configurations
    dict_n_matched = {bname:config_bundles[bname]["n_matched"] for bname in config_bundles
        if config_bundles[bname]["category"] == category and config_bundles[bname].get("tutorial", False) == tutorial
    }

    returned_countries = []

    # check if non-empty dictionary
    if dict_n_matched:
        # if non-empty, then pick bundles until all countries are selected
        # or no mor bundles are found
        dict_sort = sorted(dict_n_matched.items(), key=lambda d: d[1])

        current_matched_countries = []
        remaining_countries = set(country_list)

        for d_val in dict_sort:

            bname = d_val[0]
            cbundle_list = set(config_bundles[bname]["countries"])

            # list of countries in the bundle that are not yet matched
            intersect = cbundle_list.intersection(remaining_countries)

            if intersect:
                current_matched_countries.extend(intersect)
                remaining_countries = remaining_countries.difference(intersect)

                returned_countries.append(bname)

    return returned_countries

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_databundle_light")
    # TODO Make logging compatible with progressbar (see PR #102, PyPSA-Eur)
    configure_logging(snakemake)

    _sets_path_to_root("pypsa-africa")
    

    rootpath = os.getcwd()
    tutorial = snakemake.config["tutorial"]
    countries = snakemake.config["countries"]
    host = "zenodo"  # hard coded for now. Could be a snakemake rule param./attri
    logger.info(f"Retrieving data from {host}.")

    # load databundle configuration
    config_bundles = load_databundle_config(snakemake.input[0])

    # categories of data to download
    categories = list(set([config_bundles[conf]["category"] for conf in config_bundles]))

    # idenfify matched countries for every bundle
    for bname in config_bundles:
        config_bundles[bname]["matched_countries"] = [
            c for c in config_bundles[bname]["countries"] if c in countries]
        n_matched = len(config_bundles[bname]["matched_countries"])
        config_bundles[bname]["n_matched"] = n_matched

    # bundles to download
    bundle_to_download = []

    for cat in categories:
        if tutorial:
            selection_bundles = get_best_bundle(countries, cat, config_bundles, tutorial)

            # check if non-empty dictionary
            if selection_bundles:
                bundle_to_download.extend(selection_bundles)

                if len(selection_bundles) == 1:
                    logger.warning(f"Multiple bundle data for category {cat}: " * ", ".join(selection_bundles))

                continue
            else:
                logger.info(f"Tutorial data for {cat} not found, fall back to non-tutorial data")

        
        selection_bundles = get_best_bundle(countries, cat, config_bundles, False)

        # check if non-empty dictionary
        if selection_bundles:
            # if non-empty, then 
            bundle_to_download.extend(selection_bundles)
    
    # download the selected bundles
    for b_name in bundle_to_download:
        host_list = config_bundles[b_name]["urls"]
        dest_path = os.path.abspath(config_bundles[b_name]["destination"])
        # loop all hosts until data is successfully downloaded
        for host in host_list:
            if download_and_unzip(host, config_bundles[b_name], rootpath, dest_path):
                break

    logger.info("Bundle successfully loaded and unzipped:\n\t" + "\n\t".join(bundle_to_download))
    print("Bundle successfully loaded and unzipped:\n\t" + "\n\t".join(bundle_to_download))
