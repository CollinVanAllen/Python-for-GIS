# Python-for-GIS
Final project done for my Python for GIS course Spring 2023

## Project Goal
The goal of this project was to code a GIS toolbox to process occurrence point data and raster data using the arcpy module in python for use in ecological niche modeling (ENM) software. Data used in this project can be sourced as .csv occurrence point data from [this site][gbif] and raster data of bioclimatic variable from [this site][bioclim].  

## Toolbox Description  
This toolbox follows these steps to process the data needed for the running in an ENM software:
1. Import data as a .csv file
2. Select important features (lat,long,species,uncertainty)
3. Clean data for NA values, duplicates, values with no precision (floating points)
4. Points are spatially balanced so no two points are within a certain distance of eachother
5. Creates a buffered mask from occurrence point data
6. Mask is used to process raster data so all raster files are the same shape and coordinates
7. Rasters are translated from .tif to .asc file types
8. Gaussian Kernel Density is run to create a heatmap of species locations

## Notes on the Project
This project requires a large amount of space for all the raster files used and created. The process also takes a lengthy time to carry out given the amount of data it is fed.

[gbif]: https://www.gbif.org/
[bioclim]: https://www.worldclim.org/data/worldclim21.html
