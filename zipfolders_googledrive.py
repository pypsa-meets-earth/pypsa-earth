"""
Module to zip the desired folders to be stored in google drive, or equivalent
"""

import zipfile
import os
from os.path import basename

# Zip the files from given directory that matches the filter
def zipFilesInDir(dirName, zipFileName, filter):
   # create a ZipFile object
   with zipfile.ZipFile(zipFileName, 'w', compression=zipfile.ZIP_DEFLATED) as zipObj:
       # Iterate over all the files in directory
       for folderName, subfolders, filenames in os.walk(dirName):
           for filename in filenames:
               if filter(filename):
                   # create complete filepath of file in directory
                   filePath = os.path.join(folderName, filename)
                   # Add file to zip
                   zipObj.write(filePath, filePath)


zipFilesInDir("./resources", "resources.zip", lambda x: True)
zipFilesInDir("./data", "data.zip", lambda x: True)
zipFilesInDir("./cutouts", "cutouts.zip", lambda x: True)