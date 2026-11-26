# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later
"""
Upload or delete dataset files on Zenodo (or sandbox.Zenodo) via the REST API.

This non-workflow helper script is used by maintainers to push release datasets
(ZIP archives, cutouts, etc.) to a Zenodo deposition so that they can be
downloaded by end-users during ``snakemake`` runs.

Requirements
------------
Install the patched ``zenodopy`` client::

    pip install git+https://github.com/pz-max/zenodopy@patch-4

Obtain a personal API token from your Zenodo account and configure
``zenodopy`` as described at:
https://github.com/pz-max/zenodopy/tree/patch-4#using-the-package

For metadata customisation refer to the Zenodo REST API docs:
https://developers.zenodo.org/#introduction

Relevant Settings
-----------------
The following module-level constants control the script behaviour:

.. code-block:: python

    SANDBOX_BOOL        # True → use sandbox.zenodo.org for safe testing
    NEW_PROJECT         # True → create a new deposition; False → reuse EXISTING_PROJECT_ID
    EXISTING_PROJECT_ID # Deposition ID to update when NEW_PROJECT is False
    TYPE                # "upload" or "delete"
    PUBLISH_BOOL        # True → publish after upload (irreversible!)
    ROOT                # Absolute path to the directory containing the files
    UPLOAD_PATHS        # List of filenames (relative to ROOT) to upload
    DELETE_PATHS        # List of filenames (relative to ROOT) to delete
    METADATA            # Dict with Zenodo deposition metadata

Inputs
------
- Local files or directories listed in ``UPLOAD_PATHS``, resolved relative to
  ``ROOT``.  Directories are automatically compressed to ZIP before upload.

Outputs
-------
- Files uploaded to (or deleted from) the specified Zenodo deposition.
- If ``PUBLISH_BOOL`` is ``True`` the deposition is published and becomes
  publicly accessible with a persistent DOI.
"""

from pathlib import Path

import zenodopy

######################
# INPUTS AND OPTIONS #
######################
SANDBOX_BOOL = True  # sandbox should be used for testing
NEW_PROJECT = False  # if False, use existing project ID
EXISTING_PROJECT_ID = 1183583
TYPE = "upload"  # 'delete' or 'upload' files
PUBLISH_BOOL = False  # publish repos CANNOT be deleted
ROOT = "/home/max/Downloads"  # example "/home/max/Downloads"
UPLOAD_PATHS = [
    "tutorial_data_MA.zip",
    "cutouts_NGBJ.zip",
]  # list of files from root e.g. "README.md" in /home/max/README.md, each component will receive a download link
DELETE_PATHS = [
    "README.md",
]  # list of files from root e.g. "README.md" in /home/max/README.md,
METADATA = {
    "title": "PyPSA-Earth (Dataset)",
    "upload_type": "other",
    "description": "Used for model https://github.com/pypsa-meets-earth/pypsa-earth. Multiple data licenses apply https://pypsa-earth.readthedocs.io/en/latest/introduction.html#license",
    "creators": [
        {"name": "PyPSA-Earth Authors", "affiliation": "PyPSA meets Earth"},
    ],
    "access_right": "open",
    "license": {"id": "cc-by-4.0"},
    "keywords": ["Macro Energy Systems", "Power Systems"],
}  # more options visible at Zenodo REST API https://developers.zenodo.org/#introduction


#############
# EXECUTION #
#############
zeno = zenodopy.Client(sandbox=SANDBOX_BOOL)  # test is API key is set
zeno.list_projects


if NEW_PROJECT == True:
    zeno.create_project(title=METADATA["title"])
    zeno.change_metadata(
        dep_id=zeno.deposition_id,
        metadata=METADATA,
    )
    for path in UPLOAD_PATHS:
        path = Path.joinpath(Path(ROOT), path)
        zeno.upload_zip(source_dir=str(path))


if NEW_PROJECT == False:
    zeno.set_project(dep_id=EXISTING_PROJECT_ID)
    zeno.change_metadata(
        dep_id=zeno.deposition_id,
        metadata=METADATA,
    )
    for path in UPLOAD_PATHS:
        path = Path.joinpath(Path(ROOT), path)
        if path.exists():
            if TYPE == "upload":
                if path.exists():
                    try:
                        if path.is_file():
                            zeno.upload_file(str(path))
                        elif path.is_dir():
                            zeno.upload_zip(str(path))
                    except Exception:
                        zeno.update(path)
                        continue

            if TYPE == "delete":
                try:
                    zeno.delete_file(str(path))
                except Exception:
                    print(f"Cannot delete {path}. Repo needs to be in edit mode.")
        else:
            raise FileNotFoundError(f"{path} does not exist")


if PUBLISH_BOOL == True:
    zeno.publish()


zeno.list_projects
