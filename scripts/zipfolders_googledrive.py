# SPDX-FileCopyrightText: : 2021 PyPSA-Africa Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Module to zip the desired folders to be stored in google drive, or equivalent
"""
import os
import zipfile
from os.path import basename

from _helpers import _sets_path_to_root

# Zip the files from given directory that matches the filter


def zipFilesInDir(dirName, zipFileName, filter):
    # create a ZipFile object
    with zipfile.ZipFile(zipFileName, "w",
                         compression=zipfile.ZIP_DEFLATED) as zipObj:
        # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(dirName):
            for filename in filenames:
                if filter(filename):
                    # create complete filepath of file in directory
                    filePath = os.path.join(folderName, filename)
                    # Add file to zip
                    zipObj.write(filePath, filePath)


if __name__ == "__main__":
    # Set path to this file
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # Required to set path to pypsa-africa
    _sets_path_to_root("pypsa-africa")

# Execute zip function
zipFilesInDir("./resources", "resources.zip", lambda x: True)
zipFilesInDir("./data", "data.zip", lambda x: True)
zipFilesInDir("./cutouts", "cutouts.zip", lambda x: True)
