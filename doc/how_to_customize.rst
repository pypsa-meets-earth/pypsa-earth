..
  SPDX-FileCopyrightText: 2021 The PyPSA meets Earth authors

  SPDX-License-Identifier: CC-BY-4.0

.. _how_to_customize:

##########################################
How to customize?
##########################################

**1. Set the region in the config**
Better to set as big as it's feasible when generating a cutout. A considered area can be narrowed anytime when building a specific model.

**2. Load all the common data**
The `retrieve_databundle` rule makes the job.
Caution*: mind the flags in the config

**3. Build the custom cutout**
What is cutout
Atlite approach is used 
Copernicis API is needed

**4. Build a natura.tiff raster**
Is used to account for landuse restrictions on the protected and reserved nature areas
There is a pre-built one which is currently valid for Africa
It can be visualized
