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
      zenodo: {zenodo url}
      google: {google url}

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
from _helpers import sets_path_to_root
from _helpers import configure_logging
from _helpers import progress_retrieve
from download_osm_data import create_country_list
from google_drive_downloader import GoogleDriveDownloader as gdd

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_databundle_config(path):
    "Load databundle configurations from path file"
    with open(path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # parse the "countries" list specified in the file before processing
    for bundle_name in config:
        config[bundle_name]["countries"] = create_country_list(
            config[bundle_name]["countries"], iso_coding=False)

    return config


def download_and_unzip(host, config, rootpath, hot_run=True):
    """
        download_and_unzip(host, config, rootpath, dest_path, hot_run=True)

    Function to download and unzip the data by category

    Inputs
    ------
    host : str
        Name of the hosting platform: zenodo or google
    config : Dict
        Configuration data for the category to download
    rootpath : str
        Absolute path of the repository
    hot_run : Bool (default True)
        When true the data are downloaded
        When false, the workflow is run without downloading and unzipping

    Outputs
    -------
    True when download is successful, False otherwise

    """
    resource = "-".join(config["category"])
    file_path = os.path.join(rootpath, "tempfile.zip")

    if host == "zenodo":
        url = config["urls"]["zenodo"]
        if hot_run:
            progress_retrieve(url, file_path)
            logger.info(f"Extracting resources")
            with ZipFile(file_path, "r") as zipObj:
                # Extract all the contents of zip file in current directory
                zipObj.extractall(path=config["destination"])
            os.remove(file_path)
            logger.info(f"Download resource '{resource}' from cloud '{url}'.")
        return True
    elif host == "google":

        url = config["urls"]["google"]
        # retrieve file_id from path
        # cut the part before the ending \view
        partition_view = re.split(r"/view|\\view", str(url), 1)
        if len(partition_view) < 2:
            logger.error(
                f'Resource {resource} cannot be downloaded: "\\view" not found in url {url}'
            )
            return False

        # split url to get the file_id
        code_split = re.split(r"\\|/", partition_view[0])

        if len(code_split) < 2:
            logger.error(
                f'Resource {resource} cannot be downloaded: character "\\" not found in {partition_view[0]}'
            )
            return False

        # get file id
        file_id = code_split[-1]

        if hot_run:
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
                zipObj.extractall(path=config["destination"])
            os.remove(file_path)
        logger.info(f"Download resource '{resource}' from cloud '{url}'.")

        return True
    else:
        logger.error(f"Host {host} not implemented")
        return False


def get_best_bundles(country_list, category, config_bundles, tutorial):
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

    # load databundle configuration
    config_bundles = load_databundle_config(snakemake.input[0])

    # categories of data to download
    categories = list(
        set([config_bundles[conf]["category"] for conf in config_bundles]))

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
        selection_bundles = get_best_bundles(countries, cat, config_bundles,
                                             tutorial)

        # check if non-empty dictionary
        if selection_bundles:
            bundle_to_download.extend(selection_bundles)

            if len(selection_bundles) > 1:
                logger.warning(f"Multiple bundle data for category {cat}: " +
                               ", ".join(selection_bundles))

    # download the selected bundles
    for b_name in bundle_to_download:
        host_list = config_bundles[b_name]["urls"]
        # loop all hosts until data is successfully downloaded
        for host in host_list:
            if download_and_unzip(host, config_bundles[b_name], rootpath):
                break

    logger.info("Bundle successfully loaded and unzipped:\n\t" +
                "\n\t".join(bundle_to_download))
    # print("Bundle successfully loaded and unzipped:\n\t" +
    #       "\n\t".join(bundle_to_download))
