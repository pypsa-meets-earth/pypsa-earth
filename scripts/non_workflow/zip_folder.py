# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Module to zip the desired folders to be stored in google drive, or equivalent.
"""
import os
import zipfile
from os.path import basename
from xml.etree.ElementInclude import include

# Zip the files from given directory that matches the filter


def zipFilesInDir(dirName, zipFileName, filter, include_parent=True):
    # create a ZipFile object
    with zipfile.ZipFile(zipFileName, "w", compression=zipfile.ZIP_DEFLATED) as zipObj:
        # Iterate over all the files in directory
        for folderName, subfolders, filenames in os.walk(dirName):
            for filename in filenames:
                if filter(filename):
                    # create complete filepath of file in directory
                    filePath = os.path.join(folderName, filename)

                    # path of the zip file
                    if include_parent:
                        filePathZip = filePath
                    else:
                        filePathZip = filePath.replace(
                            dirName, ".", 1
                        )  # remove first occurrence of the dirName

                    # Add file to zip
                    zipObj.write(filePath, filePathZip)


if __name__ == "__main__":
    # Set path to this file

    # Execute zip function
    # zipFilesInDir("./resources", "resources.zip", lambda x: True, include_parent=False)
    zipFilesInDir("./data", "data.zip", lambda x: True, include_parent=False)
    # zipFilesInDir("./cutouts", "cutouts.zip", lambda x: True, include_parent=False)
