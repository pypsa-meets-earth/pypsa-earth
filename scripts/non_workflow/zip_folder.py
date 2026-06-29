# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

"""
Create ZIP archives from selected workflow directories.

This utility script packages workflow outputs and input datasets into ZIP
archives for storage, distribution, or upload to external services such as
Google Drive, Zenodo, or other file repositories.

The script recursively traverses a target directory and adds all files matching
a user-defined filter function to a compressed ZIP archive. Optionally, the
directory structure can be stored either relative to the parent directory or
including the full source directory path.

Typical use cases include packaging:

- ``data``: input datasets used by the workflow
- ``resources``: processed resources generated during preprocessing
- ``cutouts``: weather cutouts and climate datasets

Outputs
-------

- ``*.zip``: compressed archive containing the selected files and directory
  structure

"""

import os
import zipfile
from collections.abc import Callable


def zip_files_in_dir(
    dir_name: str,
    zip_file_name: str,
    file_filter: Callable[[str], bool],
    include_parent: bool = True,
) -> None:
    """
    Create a ZIP archive containing files from a directory.

    The function recursively traverses the specified directory and adds all
    files matching ``file_filter`` to the output ZIP archive.

    Parameters
    ----------
    dir_name : str
        Path to the directory to archive.
    zip_file_name : str
        Name of the ZIP archive to create.
    file_filter : Callable[[str], bool]
        Function receiving a filename and returning ``True`` if the file
        should be included in the archive.
    include_parent : bool, optional
        Whether to preserve the original directory path inside the ZIP archive.
        If ``False``, paths are stored relative to ``dir_name``.
        Default is ``True``.

    Returns
    -------
    None
        The function writes the ZIP archive to disk and does not return
        anything.
    """
    with zipfile.ZipFile(
        zip_file_name, "w", compression=zipfile.ZIP_DEFLATED
    ) as zip_obj:
        for folder_name, _, filenames in os.walk(dir_name):
            for filename in filenames:
                if file_filter(filename):
                    file_path = os.path.join(folder_name, filename)

                    if include_parent:
                        file_path_zip = file_path
                    else:
                        file_path_zip = file_path.replace(dir_name, ".", 1)

                    zip_obj.write(file_path, file_path_zip)


if __name__ == "__main__":
    # Set path to this file

    # Execute zip function
    # zipFilesInDir("./resources", "resources.zip", lambda x: True, include_parent=False)
    zip_files_in_dir("./data", "data.zip", lambda x: True, include_parent=False)
    # zipFilesInDir("./cutouts", "cutouts.zip", lambda x: True, include_parent=False)
