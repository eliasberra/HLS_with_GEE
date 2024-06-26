# HLS_with_GEE
Harmonized Landsat and Sentinel-2 data with Google Earth Engine

The JavaScript codes here presented are part of a work that aims at developing an ‘all-in-one’ Google Earth Engine (GEE) web-based workflow to produce harmonized surface reflectance data from Landsat-7 (L7) ETM+, Landsat-8 (L8) OLI and Sentinel-2 (S2) MSI top of atmosphere (TOA) reflectance data. Six major processing steps to generate a new source of near-daily Harmonized Landsat and Sentinel (HLS) reflectance observations at 30 m spatial resolution are proposed: band adjustment, atmospheric correction, cloud and cloud shadow masking, view and illumination angle adjustment, co-registration and repro-jection and resampling. The HLS is applied to six equivalent spectral bands, resulting in a surface nadir BRDF-adjusted reflectance (NBAR) time series gridded to a common pixel resolution, map projection, and spatial extent.
The rational behind the code is described in the submitted paper entitled 'Harmonized Landsat and Sentinel-2 data with Google Earth Engine'. Upon acceptance, a link to the paper will be posted here.


## How to assemble and run the code
The codes are organized in two files. The first file ('HLS_code.md') contains the general workflow with the main  parameter settings (e.g. date range) and a call to the modules (functions).
The second file contains the functions themselves ('HLS_modules.md').
Paste both files into GEE and it should work.
