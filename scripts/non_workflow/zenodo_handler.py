# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

# -*- coding: utf-8 -*-
"""
Uploads, updates, deletes, and publishes files in a Zenodo or Zenodo Sandbox
repository.

This utility script provides a simple interface for managing PyPSA-Earth data
bundles and other workflow artifacts stored on Zenodo. It uses the
``zenodopy`` package to interact with the Zenodo REST API and supports both
the creation of new depositions and the modification of existing ones.

Depending on the selected configuration, the script can:

- Create a new Zenodo deposition and assign metadata.
- Upload files or compressed directories to a deposition.
- Update existing files in a deposition.
- Delete files from an editable deposition.
- Publish a deposition, making it publicly accessible and immutable.
- Interact with either the production Zenodo service or the Zenodo Sandbox
  testing environment.

The script is primarily intended for maintaining PyPSA-Earth data bundles,
including tutorial datasets, cutouts, and other resources distributed through
Zenodo.

Requirements
------------

- Install ``zenodopy``:

  .. code:: bash

      pip install git+https://github.com/pz-max/zenodopy@patch-4

- Configure a Zenodo API access token as described in the ``zenodopy``
  documentation.
- Ensure that the target deposition is in edit mode when updating metadata
  or deleting files.

Workflow
--------

The script performs the following steps:

1. Connect to either Zenodo Sandbox or the production Zenodo instance.
2. Create a new deposition or select an existing deposition.
3. Update deposition metadata.
4. Upload, update, or delete files according to the selected operation.
5. Optionally publish the deposition.
6. Display available depositions associated with the authenticated account.

Relevant Settings
-----------------

    SANDBOX_BOOL
        Use the Zenodo Sandbox environment for testing.

    NEW_PROJECT
        Create a new deposition if True; otherwise use an existing deposition.

    EXISTING_PROJECT_ID
        Identifier of the deposition to modify.

    TYPE
        Operation to perform: ``"upload"`` or ``"delete"``.

    PUBLISH_BOOL
        Publish the deposition after modifications.

    ROOT
        Root directory containing files to upload.

    UPLOAD_PATHS
        List of files or directories to upload.

    METADATA
        Metadata associated with the deposition.

Notes
-----

- Published depositions cannot be modified or deleted. A new version must be
  created through the Zenodo interface if further changes are required.
- Directories are automatically compressed before upload.
- Existing files with matching names are updated when possible.
- Zenodo Sandbox should be used for testing before uploading data to the
  production service.

See Also
--------

- Zenodo REST API: https://developers.zenodo.org/
- zenodopy: https://github.com/pz-max/zenodopy
- PyPSA-Earth: https://github.com/pypsa-meets-earth/pypsa-earth

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
                    except:
                        zeno.update(path)
                        continue

            if TYPE == "delete":
                try:
                    zeno.delete_file(str(path))
                except:
                    print(f"Cannot delete {path}. Repo needs to be in edit mode.")
        else:
            raise FileNotFoundError(f"{path} does not exist")


if PUBLISH_BOOL == True:
    zeno.publish()


zeno.list_projects
