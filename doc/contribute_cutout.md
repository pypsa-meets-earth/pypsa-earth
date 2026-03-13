
1. Build a cutout with build_cutout rule
   - snapshots correspond to the standard time frame
   ```
   # full-scale cutout
   snapshots:
     start: "2013-01-01"
     end: "2014-01-01" 
   ```

   ```
   snapshots:
     start: "2013-03-1"
     end: "2013-03-7"   
   ```

   - population and GPD aren't added to the shapes
   ```
   build_shape_options:
       worldpop_method: false
       gdp_method: false
   ```
2. Check the cutout with a cutout checker

3. Compress the cutout and upload a file:
   - use a standard name for the cutout file (cutout-2013-era5.nc for full-scale cutouts and cutout-2013-era5-tutorial.nc for tutorial and testing ones)
   - use zip format for compression   
   - do not include an enclosing folder into the archive
   - change the archive name to a meaningful one
4. Create Zenodo record
   - add a new upload or a new version of the existing one
   - provide a description
5. Add a field to bundle_config.yaml aligning the format with the existing records. An example:
    ```
      bundle_cutouts_tutorial_MA:
        countries: [MA]
        tutorial: true
        category: cutouts
        destination: "cutouts"
        urls:
          zenodo: https://zenodo.org/records/18851102/files/    bundle_cutouts_tutorial_MA.zip?download=1
          gdrive: https://drive.google.com/file/d/1j5v2f4E756jmDMa707QvdNJq3xM4bYUk/    view?usp=drive_link
        output: [cutouts/cutout-2013-era5-tutorial.nc]
        disable_by_opt:
          build_cutout: [all]
    ```   