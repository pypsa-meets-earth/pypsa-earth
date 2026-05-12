---
name: Data collection issue
about: Add a dataset to the project
title: ""
labels: "data collection"
assignees: ""

---

## A check-list for adding new data sources

- [ ] Download a raw data source (e.g. pdf file)
- [ ] Share the data source in a local storage folder:
    - create a dedicated sub-folder for the particular data source
    - include `Raw` and `Processed` folders inside the sub-folder
    - place the raw source file into `Raw` folder and add `Reference` doc file with url link to the source and a recommended citation if applicable
- [ ] Extract relevant parts of the source into csv format (e.g. tables for the installed generation)
- [ ] Share csv files in `Processed` folder
- [ ] Create a metadata file for csv files and place into `Processed`
- [ ] Upload csv files on Zenodo using the metadata as a description
- [ ] Contribute a source link to the data collection repositories:
   - [Awesome-Electrical-Grid-Mapping](https://github.com/open-energy-transition/Awesome-Electrical-Grid-Mapping)
   - [PyPSA-Earth-Status](https://github.com/pypsa-meets-earth/pypsa-earth-status)
   - [Awesome-Electricity-Demand](https://github.com/open-energy-transition/Awesome-Electricity-Demand)
