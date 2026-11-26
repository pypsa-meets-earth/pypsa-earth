# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Utility to compress PyPSA-Earth output folders for upload to cloud storage.

This non-workflow helper script zips a directory tree (data, resources, or
cutouts) so that the resulting archive can be shared via Google Drive, Zenodo,
or any equivalent service.

**Typical usage**

Run the script directly from the repository root:

```bash
python scripts/non_workflow/zip_folder.py
```

Uncomment the relevant ``zipFilesInDir`` call at the bottom of the file to
choose which folder to compress.

**Inputs**

- A local directory (e.g. ``./data``, ``./resources``, ``./cutouts``) that
  already contains the files to archive.

**Outputs**

- A ``.zip`` archive in the current working directory
  (e.g. ``data.zip``).
"""

import os
import zipfile
from collections.abc import Callable


def zipFilesInDir(
    dirName: str,
    zipFileName: str,
    filter: Callable[[str], bool],
    include_parent: bool = True,
) -> None:
    """
    Compress all files in *dirName* that satisfy *filter* into a ZIP archive.

    Parameters
    ----------
    dirName : str
        Path to the source directory to walk recursively.
    zipFileName : str
        Path (and name) of the output ``.zip`` file to create.
    filter : Callable[[str], bool]
        A function that receives a filename (without its directory path) and
        returns ``True`` if the file should be included in the archive.
    include_parent : bool, optional
        When ``True`` (default), archive paths include the full path from the
        current working directory.  When ``False``, the *dirName* prefix is
        stripped so the archive root starts with ``./``.

    Returns
    -------
    None

    Examples
    --------
    Zip every file in ``./data``, stripping the top-level directory prefix:

    >>> zipFilesInDir("./data", "data.zip", lambda x: True, include_parent=False)
    """
    with zipfile.ZipFile(zipFileName, "w", compression=zipfile.ZIP_DEFLATED) as zipObj:
        for folderName, subfolders, filenames in os.walk(dirName):
            for filename in filenames:
                if filter(filename):
                    filePath = os.path.join(folderName, filename)
                    filePathZip = (
                        filePath
                        if include_parent
                        else filePath.replace(dirName, ".", 1)
                    )
                    zipObj.write(filePath, filePathZip)


if __name__ == "__main__":
    # Uncomment the folder you want to compress before running.
    # zipFilesInDir("./resources", "resources.zip", lambda x: True, include_parent=False)
    zipFilesInDir("./data", "data.zip", lambda x: True, include_parent=False)
    # zipFilesInDir("./cutouts", "cutouts.zip", lambda x: True, include_parent=False)

