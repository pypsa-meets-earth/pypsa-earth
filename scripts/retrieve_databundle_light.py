# Copyright 2019-2020 Fabian Hofmann (FIAS)
# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
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
from _helpers import progress_retrieve, configure_logging
from google_drive_downloader import GoogleDriveDownloader as gdd
import tarfile
from pathlib import Path
from _helpers import _sets_path_to_root

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

        from _helpers import mock_snakemake
        snakemake = mock_snakemake('retrieve_databundle_light')
        rootpath = '..'
    else:
        rootpath = '.'
    configure_logging(snakemake) # TODO Make logging compatible with progressbar (see PR #102)

    _sets_path_to_root("pypsa-africa")


    # BUNDLE 1
    destination = './resources'
    zip_path = destination +'.zip'
    url = "https://drive.google.com/file/d/1nATTTS4t-EBBaAScNbZksVqsfUiL9TTV/view?usp=sharing"
    gdd.download_file_from_google_drive(file_id='1nATTTS4t-EBBaAScNbZksVqsfUiL9TTV',
                                        dest_path=zip_path,
                                        unzip=True)
    os.remove(zip_path)
    logger.info(f"Download data to '{destination}' from cloud '{url}'.")


    # BUNDLE 2
    destination = './data'
    zip_path = destination +'.zip'
    url = "https://drive.google.com/file/d/1lxfxCrUBdnQdFd6Isx3XlGwP6d9KRW2E/view?usp=sharing"
    gdd.download_file_from_google_drive(file_id='1lxfxCrUBdnQdFd6Isx3XlGwP6d9KRW2E',
                                        dest_path=zip_path,
                                        unzip=True)
    os.remove(zip_path)
    logger.info(f"Download data to '{destination}' from cloud '{url}'.")


    # BUNDLE 3
    destination = './cutouts'
    zip_path = destination +'.zip'
    url = "https://drive.google.com/file/d/1FvjV6CmqoQMMUF0otVnG3fBQm2_GGISW/view?usp=sharing"
    gdd.download_file_from_google_drive(file_id='1FvjV6CmqoQMMUF0otVnG3fBQm2_GGISW',
                                        dest_path=zip_path,
                                        unzip=True)
    os.remove(zip_path)
    logger.info(f"Download data to '{destination}' from cloud '{url}'.")