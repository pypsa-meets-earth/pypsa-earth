# -*- coding: utf-8 -*-
# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText: : 2021-2022 PyPSA-Africa, 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
"""
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5894972.svg
   :target: https://doi.org/10.5281/zenodo.5894972

The data bundles contains common GIS datasets like EEZ shapes, Copernicus Landcover, Hydrobasins
and also electricity specific summary statistics like historic per country yearly totals of hydro generation,
GDP and POP on NUTS3 levels and per-country load time-series.

This rule downloads the data bundle from `zenodo <https://doi.org/10.5281/zenodo.5894972>`_
or `google drive <https://drive.google.com/drive/u/1/folders/1dkW1wKBWvSY4i-XEuQFFBj242p0VdUlM>`_
and extracts it in the ``data``, ``resources`` and ``cutouts`` sub-directory.
Bundle data are then deleted once downloaded and unzipped.

The :ref:`tutorial` uses a smaller `data bundle <https://zenodo.org/record/3517921/files/pypsa-eur-tutorial-data-bundle.tar.xz>`_
than required for the full model (around 500 MB)

The required bundles are downloaded automatically according to the list names, in agreement to
the data bundles specified in the bundle configuration file, typically located in the ``config`` folder.
Each data bundle entry has the following structure:

.. code:: yaml

  bundle_name:  # name of the bundle
    countries: [country code, region code or country list]  # list of countries represented in the databundle
    [tutorial: true/false]  # (optional, default false) whether the bundle is a tutorial or not
    category: common/resources/data/cutouts  # category of data contained in the bundle:
    destination: "."  # folder where to unzip the files with respect to the repository root (\"\" or \".\")
    urls:  # list of urls by source, e.g. zenodo or google
      zenodo: {zenodo url}  # key to download data from zenodo
      gdrive: {google url}  # key to download data from google drive
      protectedplanet: {url}  # key to download data from protected planet
      direct: {url}  # key to download data directly from a url; if unzip option is enabled data are unzipped
      post:  # key to download data using an url post request; if unzip option is enabled data are unzipped
        url: {url}
        [post arguments]
    [unzip: true/false]  # (optional, default false) used in direct download technique to automatically unzip files
    output: [...]  # list of outputs of the databundle
    [disable_by_opt:]  # option to disable outputs from the bundle; it contains a dictionary of options, each one with
                       # each one with its output. When "all" is specified, the entire bundle is not executed
      [{option}: [outputs,...,/all]]  # list of options and the outputs to remove, or "all" corresponding to ignore everything

Depending on the country list that is asked to perform, all needed databundles are downloaded
according to the following rules:

- The databundle shall adhere to the tutorial configuration: when
  the tutorial configuration is running, only the databundles having tutorial flag true
  shall be downloaded
- For every data category, the most suitable bundles are downloaded by order of
  number of countries matched: for every bundles matching the category,
  the algorithm sorts the bundles by the number of countries that are matched and starts
  downloading them starting from those matching more countries till all countries are matched
  or no more bundles are available
- For every bundle to download, it is given priority to the first bundle source,
  as listed in the ``urls`` option of each bundle configuration; when a source fails,
  the following source is used and so on

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3517921.svg
    :target: https://doi.org/10.5281/zenodo.3517921

**Relevant Settings**

.. code:: yaml

    tutorial:  # configuration stating whether the tutorial is needed


.. seealso::
    Documentation of the configuration file ``config.yaml`` at
    :ref:`toplevel_cf`

**Outputs**

- ``data``: input data unzipped into the data folder
- ``resources``: input data unzipped into the resources folder
- ``cutouts``: input data unzipped into the cutouts folder

"""
import logging
import os
import re
from zipfile import ZipFile

import yaml
from _helpers import configure_logging, progress_retrieve, sets_path_to_root
from download_osm_data import create_country_list
from google_drive_downloader import GoogleDriveDownloader as gdd

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_databundle_config(config):
    "Load databundle configurations from path file or dictionary"

    if type(config) is str:
        with open(config) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)["databundles"]
    elif type(config) is not dict:
        logger.error("Impossible to load the databundle configuration")

    # parse the "countries" list specified in the file before processing
    for bundle_name in config:
        config[bundle_name]["countries"] = create_country_list(
            config[bundle_name]["countries"], iso_coding=False
        )

    return config


def download_and_unzip_zenodo(config, rootpath, hot_run=True, disable_progress=False):
    """
        download_and_unzip_zenodo(config, rootpath, dest_path, hot_run=True, disable_progress=False)

    Function to download and unzip the data from zenodo

    Inputs
    ------
    config : Dict
        Configuration data for the category to download
    rootpath : str
        Absolute path of the repository
    hot_run : Bool (default True)
        When true the data are downloaded
        When false, the workflow is run without downloading and unzipping
    disable_progress : Bool (default False)
        When true the progress bar to download data is disabled

    Outputs
    -------
    True when download is successful, False otherwise

    """
    resource = config["category"]
    file_path = os.path.join(rootpath, "tempfile.zip")

    url = config["urls"]["zenodo"]
    if hot_run:
        try:
            logger.info(f"Downloading resource '{resource}' from cloud '{url}'")
            progress_retrieve(url, file_path, disable_progress=disable_progress)
            logger.info(f"Extracting resources")
            with ZipFile(file_path, "r") as zipObj:
                # Extract all the contents of zip file in current directory
                zipObj.extractall(path=config["destination"])
            os.remove(file_path)
            logger.info(f"Downloaded resource '{resource}' from cloud '{url}'.")
        except:
            logger.warning(f"Failed download resource '{resource}' from cloud '{url}'.")
            return False

    return True


def download_and_unzip_protectedplanet(
    config, rootpath, hot_run=True, disable_progress=False
):
    """
        download_and_unzip_protectedplanet(config, rootpath, dest_path, hot_run=True, disable_progress=False)

    Function to download and unzip the data by category from protectedplanet

    Inputs
    ------
    config : Dict
        Configuration data for the category to download
    rootpath : str
        Absolute path of the repository
    hot_run : Bool (default True)
        When true the data are downloaded
        When false, the workflow is run without downloading and unzipping
    disable_progress : Bool (default False)
        When true the progress bar to download data is disabled

    Outputs
    -------
    True when download is successful, False otherwise

    """
    resource = config["category"]
    file_path = os.path.join(rootpath, "tempfile_wpda.zip")

    url = config["urls"]["protectedplanet"]

    if hot_run:
        if os.path.exists(file_path):
            os.remove(file_path)

        try:
            logger.info(f"Downloading resource '{resource}' from cloud '{url}'.")
            progress_retrieve(url, file_path, disable_progress=disable_progress)

            zip_obj = ZipFile(file_path, "r")

            # list of zip files, which contains the shape files
            zip_files = [
                fname for fname in zip_obj.namelist() if fname.endswith(".zip")
            ]

            # extract the nested zip files
            for fzip in zip_files:
                # final path of the file
                inner_zipname = os.path.join(config["destination"], fzip)

                zip_obj.extract(fzip, path=config["destination"])

                with ZipFile(inner_zipname, "r") as nested_zip:
                    nested_zip.extractall(path=config["destination"])

                # remove inner zip file
                os.remove(inner_zipname)

            # remove outer zip file
            os.remove(file_path)

            logger.info(f"Downloaded resource '{resource}' from cloud '{url}'.")
        except:
            logger.warning(f"Failed download resource '{resource}' from cloud '{url}'.")
            return False

    return True


def download_and_unzip_direct(config, rootpath, hot_run=True, disable_progress=False):
    """
        download_and_unzip_direct(config, rootpath, dest_path, hot_run=True, disable_progress=False)

    Function to download the data by category from a direct url with no processing.
    If in the configuration file the unzip is specified True, then the downloaded data is unzipped.

    Inputs
    ------
    config : Dict
        Configuration data for the category to download
    rootpath : str
        Absolute path of the repository
    hot_run : Bool (default True)
        When true the data are downloaded
        When false, the workflow is run without downloading and unzipping
    disable_progress : Bool (default False)
        When true the progress bar to download data is disabled

    Outputs
    -------
    True when download is successful, False otherwise

    """
    resource = config["category"]
    url = config["urls"]["direct"]

    file_path = os.path.join(config["destination"], os.path.basename(url))

    if hot_run:
        if os.path.exists(file_path):
            os.remove(file_path)

        try:
            logger.info(f"Downloading resource '{resource}' from cloud '{url}'.")
            progress_retrieve(url, file_path, disable_progress=disable_progress)

            # if the file is a zipfile and unzip is enabled
            # then unzip it and remove the original file
            if config.get("unzip", False):
                with ZipFile(file_path, "r") as zipfile:
                    zipfile.extractall(config["destination"])

                os.remove(file_path)
            logger.info(f"Downloaded resource '{resource}' from cloud '{url}'.")
        except:
            logger.warning(f"Failed download resource '{resource}' from cloud '{url}'.")
            return False

    return True


def download_and_unzip_post(config, rootpath, hot_run=True, disable_progress=False):
    """
        download_and_unzip_post(config, rootpath, dest_path, hot_run=True, disable_progress=False)

    Function to download the data by category from a post request.

    Inputs
    ------
    config : Dict
        Configuration data for the category to download
    rootpath : str
        Absolute path of the repository
    hot_run : Bool (default True)
        When true the data are downloaded
        When false, the workflow is run without downloading and unzipping
    disable_progress : Bool (default False)
        When true the progress bar to download data is disabled

    Outputs
    -------
    True when download is successful, False otherwise

    """
    resource = config["category"]

    # load data for post method
    postdata = config["urls"]["post"]
    # remove url feature
    url = postdata.pop("url")

    file_path = os.path.join(config["destination"], os.path.basename(url))

    if hot_run:
        if os.path.exists(file_path):
            os.remove(file_path)

        # try:
        logger.info(f"Downloading resource '{resource}' from cloud '{url}'.")

        progress_retrieve(
            url, file_path, data=postdata, disable_progress=disable_progress
        )

        # if the file is a zipfile and unzip is enabled
        # then unzip it and remove the original file
        if config.get("unzip", False):
            with ZipFile(file_path, "r") as zipfile:
                zipfile.extractall(config["destination"])

            os.remove(file_path)
        logger.info(f"Downloaded resource '{resource}' from cloud '{url}'.")
        # except:
        #     logger.warning(f"Failed download resource '{resource}' from cloud '{url}'.")
        #     return False

    return True


def _check_disabled_by_opt(config_bundle, config_enable):
    """
    Checks if the configbundle has conflicts with the enable configuration

    Returns
    -------
    disabled : Bool
        True when the bundle is completely disabled
    """

    disabled_outs = []

    if "disable_by_opt" in config_bundle:
        disabled_config = config_bundle["disable_by_opt"]
        disabled_objs = [
            disabled_outputs
            for optname, disabled_outputs in disabled_config.items()
            if config_enable.get(optname, False)
        ]

        # merge all the lists unique elements
        all_disabled = []
        for tot_outs in disabled_objs:
            for out in tot_outs:
                if out not in all_disabled:
                    all_disabled.append(out)

        if "all" in all_disabled:
            disabled_outs = ["all"]
        elif "output" in config_enable:
            disabled_outs = list(set(all_disabled))

    return disabled_outs


def get_best_bundles(country_list, category, config_bundles, tutorial, config_enable):
    """
        get_best_bundles(country_list, category, config_bundles, tutorial)

    Function to get the best bundles that download the datafor selected countries,
    given category and tutorial characteristics.

    The selected bundles shall adhere to the following criteria:
    - The bundles' tutorial parameter shall match the tutorial argument
    - The bundles' category shall match the category of data to download
    - When multiple bundles are identified for the same set of users,
    the bundles matching more countries are first selected and more bundles
    are added until all countries are matched or no more bundles are available

    Inputs
    ------
    country_list : list
        List of country codes for the countries to download
    category : str
        Category of the data to download
    config_bundles : Dict
        Dictionary of configurations for all available bundles
    tutorial : Bool
        Whether data for tutorial shall be downloaded
    config_enable : dict
        Dictionary of the enabled/disabled scripts

    Outputs
    -------
    returned_bundles : list
        List of bundles to download

    """
    # dictionary with the number of match by configuration for tutorial/non-tutorial configurations
    dict_n_matched = {
        bname: config_bundles[bname]["n_matched"]
        for bname in config_bundles
        if config_bundles[bname]["category"] == category
        and config_bundles[bname].get("tutorial", False) == tutorial
        and _check_disabled_by_opt(config_bundles[bname], config_enable) != ["all"]
    }

    returned_bundles = []

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

                returned_bundles.append(bname)

    return returned_bundles


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_databundle_light")
    # TODO Make logging compatible with progressbar (see PR #102, PyPSA-Eur)
    configure_logging(snakemake)

    sets_path_to_root("pypsa-africa")

    rootpath = os.getcwd()
    tutorial = snakemake.config["tutorial"]
    countries = snakemake.config["countries"]
    logger.info(f"Retrieving data for {len(countries)} countries.")

    disable_progress = not snakemake.config.get("retrieve_databundle", {}).get(
        "show_progress", True
    )

    # load enable configuration
    config_enable = snakemake.config["enable"]
    # load databundle configuration
    config_bundles = load_databundle_config(snakemake.config["databundles"])

    # categories of data to download
    categories = list(
        set([config_bundles[conf]["category"] for conf in config_bundles])
    )

    # idenfify matched countries for every bundle
    for bname in config_bundles:
        config_bundles[bname]["matched_countries"] = [
            c for c in config_bundles[bname]["countries"] if c in countries
        ]
        n_matched = len(config_bundles[bname]["matched_countries"])
        config_bundles[bname]["n_matched"] = n_matched

    # bundles to download
    bundle_to_download = []

    for cat in categories:
        selection_bundles = get_best_bundles(
            countries, cat, config_bundles, tutorial, config_enable
        )

        # check if non-empty dictionary
        if selection_bundles:
            bundle_to_download.extend(selection_bundles)

            if len(selection_bundles) > 1:
                logger.warning(
                    f"Multiple bundle data for category {cat}: "
                    + ", ".join(selection_bundles)
                )

    logger.warning(
        "DISCLAIMER LICENSES: the use of PyPSA-Africa is conditioned \
        to the acceptance of its multiple licenses.\n \
        The use of the code automatically implies that you accept all the licenses.\n \
        See our documentation for more information. \n \
        Link: https://pypsa-meets-africa.readthedocs.io/en/latest/introduction.html#licence"
    )

    # download the selected bundles
    for b_name in bundle_to_download:
        host_list = config_bundles[b_name]["urls"]
        # loop all hosts until data is successfully downloaded
        for host in host_list:
            # try:
            download_and_unzip = globals()[f"download_and_unzip_{host}"]
            if download_and_unzip(
                config_bundles[b_name], rootpath, disable_progress=disable_progress
            ):
                break
            # except KeyError:
            #     logger.error(f"Function for {host} has not been defined")

    logger.info(
        "Bundle successfully loaded and unzipped:\n\t" + "\n\t".join(bundle_to_download)
    )
    # print("Bundle successfully loaded and unzipped:\n\t" +
    #       "\n\t".join(bundle_to_download))
