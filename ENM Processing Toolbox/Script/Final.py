# Occurence Data

"""
Created on Wed Apr 12 12:58:53 2023

@author: vanallenc
"""

import arcpy, os, os.path
import pandas as pd
import numpy as np
from arcpy import sa,env
import matplotlib

#workspace stuff

env.workspace = r"E:\Python\Final Project\ProjectActual\Occurrence"
env.overwriteOutput = True

coord = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

# =============================================================================
# Clean Data
# =============================================================================

dataPath = env.workspace + r"\OccurrenceData.csv"

with open(dataPath) as f:
    file_encoding = f.encoding
    print(file_encoding)

dataset = pd.read_csv(dataPath, encoding = file_encoding)

reducedDataset = dataset[["species", "decimalLatitude", "decimalLongitude",
                         "coordinateUncertaintyInMeters", "year"]]

reducedDataset.species = reducedDataset["species"].values[0]

reducedDataset = reducedDataset.dropna(subset = ['decimalLongitude', 'decimalLatitude'])
reducedDataset = reducedDataset[reducedDataset['decimalLongitude' or 'decimalLatitude'] != 0]

reducedDataset = reducedDataset[reducedDataset['decimalLongitude'].apply(lambda x: x.is_integer()) == False]
reducedDataset = reducedDataset[reducedDataset['decimalLatitude'].apply(lambda x: x.is_integer()) == False]

reducedDataset = reducedDataset.drop_duplicates(subset=["decimalLatitude", "decimalLongitude"])

reducedDataset.to_csv(env.workspace + r"\bias.csv")

reducedDataset = reducedDataset.rename(columns={"species": "Species", "decimalLatitude": "decLat", "decimalLongitude": "decLong"})

# =============================================================================
# Create Shape
# =============================================================================
# Define output path and feature class name
def makeShapeFile(out, newFile, geom, dataset):
    outpath = out
    newfc = newFile
    geoType = geom

    # Create new feature class
    arcpy.CreateFeatureclass_management(outpath, newfc, geoType)
    
    newfc2 = outpath + "\\" + newfc
    
    arcpy.AddField_management(newfc2, "Species", "TEXT", field_length=50)
    arcpy.AddField_management(newfc2, "decLong", "FLOAT")
    arcpy.AddField_management(newfc2, "decLat", "FLOAT")


    # Create point object, then use insert cursor to add points to new feature class

    cursor = arcpy.da.InsertCursor(newfc2, ["SHAPE@", "Species", "decLong", "decLat"])
    for index, row in dataset.iterrows():
        point = arcpy.Point(row[dataset.columns.get_loc('decLong')], row[dataset.columns.get_loc('decLat')])
        cursor.insertRow([point, row['Species'], row[dataset.columns.get_loc('decLong')], row[dataset.columns.get_loc('decLat')]])
    del cursor

    arcpy.DefineProjection_management (newfc2, coord)

    arcpy.Describe(newfc2).spatialReference.Name
    
    

makeShapeFile(env.workspace, "bias.shp", "POINT", reducedDataset)


# =============================================================================
# Spatially balance points
# =============================================================================
if arcpy.Exists("bias.shp"):
    arcpy.CopyFeatures_management("bias.shp", env.workspace + r"/points.shp")

infc = env.workspace + r"\points.shp"
fields = ["SHAPE"]
xy_tol = "1.0 Kilometers"

arcpy.management.DeleteIdentical(infc, fields, xy_tol)

arcpy.conversion.ExportTable(infc, env.workspace + r"\points.csv")

# =============================================================================
# Split to test and train
# =============================================================================
points = pd.read_csv(env.workspace + r"/points.csv")

mask = np.random.rand(len(points)) <= 0.6
train = points[mask]
test = points[~mask]

train.to_csv(env.workspace + r"\Sp_train.csv")
test.to_csv(env.workspace + r"\Sp_test.csv")

makeShapeFile(env.workspace, "Sp_train.shp", "POINT", train)
makeShapeFile(env.workspace, "Sp_test.shp", "POINT", test)


# =============================================================================
# Polygon Mask and buffer
# =============================================================================

infc = env.workspace + r"\points.shp"
outfc = env.workspace + r"\MinBoundGeom.shp"

arcpy.management.MinimumBoundingGeometry(infc, outfc, "CONVEX_HULL", "ALL")

#------------------------------------------------------------------------------

infc = env.workspace + r"\MinBoundGeom.shp"
outfc = env.workspace + r"\MinBoundGeomBuff.shp"
buffer = "250 Kilometers"
line_side = "FULL"
dissolve = "ALL"
method = "PLANAR"

arcpy.analysis.Buffer(infc, outfc, buffer, line_side = line_side, dissolve_option = dissolve, method = method)

# =============================================================================
# Process Enviromental Data
# =============================================================================

# =============================================================================
# Mask Bias
# =============================================================================

env.workspace = r"E:\Python\Final Project\ProjectActual"
env.overwriteOutput = True

inRast = env.workspace + r"\BioClim\wc2.1_30s_bio_1.tif"
mask = env.workspace + r"\Occurrence\MinBoundGeomBuff.shp"

arcpy.env.outputCoordinateSystem = arcpy.Describe(inRast).spatialReference
arcpy.env.extent = mask
arcpy.env.cellSize = inRast
arcpy.env.mask = mask
arcpy.env.snapRaster = inRast


outExtractByMask = arcpy.sa.ExtractByMask(inRast, mask, extraction_area = "INSIDE")

outExtractByMask.save(r"E:\Python\Final Project\ProjectActual\MaskBias\rasterMask.tif")

# =============================================================================
# Process Rest of Rasters
# =============================================================================

mask = r"E:\Python\Final Project\ProjectActual\MaskBias\rasterMask.tif"

arcpy.env.outputCoordinateSystem = mask
arcpy.env.extent = mask
arcpy.env.cellSize = mask
arcpy.env.mask = mask
arcpy.env.snapRaster = mask

env.workspace = r"E:\Python\Final Project\ProjectActual\BioClim"

rasters = arcpy.ListRasters("*", "TIF")

index = 1
for raster in rasters:
    outExtractByMask = arcpy.sa.ExtractByMask(raster, mask, extraction_area = "INSIDE")
    filePath = r"E:\Python\Final Project\ProjectActual\BioclimClip"
    outExtractByMask.save(filePath + "\\" + raster)
    print("Raster",index,"complete out of",str(len(rasters)))
    index += 1
    
# =============================================================================
# Raster to ASCII
# =============================================================================

env.workspace = r"E:\Python\FinalProject_CVA\BioclimClip"

rasters = arcpy.ListRasters("*", "TIF")

index = 1
for raster in rasters:
    filename = raster[:-4] + ".asc"
    filepath = r"E:\Python\FinalProject_CVA\BioclimAsc\\"
    output = filepath + filename
    arcpy.conversion.RasterToASCII(raster,output)
    print("Raster",index,"converted out of",str(len(rasters)))
    index += 1

# =============================================================================
# Kernel Density
# =============================================================================
env.workspace = r"E:\Python\Final_Project\ProjectActual"
env.overwriteOutput = True
infc = env.workspace +  r"\Occurrence\bias.shp"
outputcell = env.workspace + r"\MaskBias\rasterMask.tif"
outval = "DENSITIES"
population = "NONE"
method = "PLANAR"

arcpy.env.mask = outputcell
arcpy.env.snapRaster = outputcell

outKDens = arcpy.sa.KernelDensity(in_features = infc, population_field = population, 
                                  cell_size = outputcell, out_cell_values = outval, method = method)
outfc = env.workspace + r"\MaskBias\gkd.tif"
outKDens.save(outfc)


out_raster = arcpy.sa.KernelDensity(
    in_features=r"E:\Python\Final_Project\ProjectActual\Occurrence\bias.shp",
    population_field="NONE",
    cell_size=r"E:\Python\Final_Project\ProjectActual\MaskBias\rasterMask.tif",
    search_radius=None,
    area_unit_scale_factor="SQUARE_MAP_UNITS",
    out_cell_values="DENSITIES",
    method="PLANAR",
    in_barriers=None
)
out_raster.save(r"E:\Python\Final_Project\ProjectActual\MaskBias\gkd.tif")




