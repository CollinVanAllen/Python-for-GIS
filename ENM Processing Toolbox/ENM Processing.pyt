# Toolbox for ENM Modeling

"""
Created on Wed May  3 12:22:46 2023
Purpose: This tool takes point data in csv, cleans it, and converts it to a shape file. 
Once completed the points are used to create buffer and mask shapes to process BioClim rasters,
and provide a finished set of masked rasters in tif and asc format.

@author: Collin Van Allen
"""

import arcpy, os, os.path
import pandas as pd
import numpy as np
from arcpy import sa



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [ENM_Process]


class ENM_Process(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "ENM Processing"
        self.description = "Tool takes point data and converts it shape, as well as returning clipped rasters in tif and asc format"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        # Input Params=========================================================
        param0 = arcpy.Parameter(
            displayName = "Input Dataset",
            name = "in_features",
            datatype = "DETextfile",
            parameterType = "Required",
            direction = "Input")
        param0.filter.list = ["csv"]
        
        param1 = arcpy.Parameter(
            displayName = "XY Tolerance",
            name = "xy_tolerance",
            datatype = "GPLong",
            parameterType = "Required",
            direction = "Input")
        param1.filter.type = "Range"
        param1.filter.list = [0,float('inf')]
        
        param2 = arcpy.Parameter(
            displayName="Train/Test Split",
            name="split_percentage",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        param2.filter.type = "Range"
        param2.filter.list = [0,1]
        
        param3 = arcpy.Parameter(
            displayName="Buffer Distance",
            name="buffer_distance",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        param3.filter.type = "Range"
        param3.filter.list = [0,float('inf')]
        
        param4 = arcpy.Parameter(
            displayName = "Tolerance and Buffer Units",
            name = "distance_unit",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        param4.filter.type = "ValueList"
        param4.filter.list = ["Miles", "Kilometers"]
        
        param5 = arcpy.Parameter(
            displayName="Environmental Data Folder",
            name="raster_folder",
            datatype="DEFolder",
            parameterType="Required")
        
        # Output Params========================================================
        param6 = arcpy.Parameter(
            displayName="Occurrence Data Folder",
            name="occurrence_folder",
            datatype="DEFolder",
            parameterType="Required")
        
        param7 = arcpy.Parameter(
            displayName="Mask Raster Folder",
            name="mask_folder",
            datatype="DEFolder",
            parameterType="Required")
        
        param8 = arcpy.Parameter(
            displayName="Clipped Raster Data Folder",
            name="clipped_mask_folder",
            datatype="DEFolder",
            parameterType="Required")
        
        param9 = arcpy.Parameter(
            displayName="ASCII Raster Data Folder",
            name="ascii_raster_folder",
            datatype="DEFolder",
            parameterType="Required")
        
        params = [param0,param1,param2,param3,param4,param5,param6,param7,param8,param9]
        return params

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        coord = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"

        # =============================================================================
        # Clean Data
        # =============================================================================

        dataPath = parameters[6].valueAsText

        with open(parameters[0].valueAsText) as f:
            file_encoding = f.encoding

        dataset = pd.read_csv(parameters[0].valueAsText, encoding = file_encoding)

        reducedDataset = dataset[["species", "decimalLatitude", "decimalLongitude",
                                 "coordinateUncertaintyInMeters", "year"]]

        reducedDataset.species = reducedDataset["species"].values[0]

        reducedDataset = reducedDataset.dropna(subset = ['decimalLongitude', 'decimalLatitude'])
        reducedDataset = reducedDataset[reducedDataset['decimalLongitude' or 'decimalLatitude'] != 0]

        reducedDataset = reducedDataset[reducedDataset['decimalLongitude'].apply(lambda x: x.is_integer()) == False]
        reducedDataset = reducedDataset[reducedDataset['decimalLatitude'].apply(lambda x: x.is_integer()) == False]

        reducedDataset = reducedDataset.drop_duplicates(subset=["decimalLatitude", "decimalLongitude"])

        reducedDataset.to_csv(dataPath + r"\bias.csv")
        
        arcpy.AddMessage("Data cleaned")

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
            
            

        makeShapeFile(dataPath, "bias.shp", "POINT", reducedDataset)

        arcpy.AddMessage("Bias shapefile created")
        # =============================================================================
        # Spatially balance points
        # =============================================================================
        if arcpy.Exists(dataPath + r"\bias.shp"):
            arcpy.CopyFeatures_management(dataPath + r"\bias.shp", dataPath + r"\points.shp")

        infc = dataPath + r"\points.shp"
        fields = ["SHAPE"]
        xy_tol = str(int(parameters[1].valueAsText)) + " " + parameters[4].valueAsText

        arcpy.management.DeleteIdentical(infc, fields, xy_tol)

        arcpy.conversion.ExportTable(infc, dataPath + r"\points.csv")

        arcpy.AddMessage("Data rarefied")
        # =============================================================================
        # Split to test and train
        # =============================================================================
        points = pd.read_csv(dataPath + r"\points.csv")

        mask = np.random.rand(len(points)) <= float(parameters[2].valueAsText)
        train = points[mask]
        test = points[~mask]

        train.to_csv(dataPath + r"\Sp_train.csv")
        test.to_csv(dataPath + r"\Sp_test.csv")

        makeShapeFile(dataPath, "Sp_train.shp", "POINT", train)
        makeShapeFile(dataPath, "Sp_test.shp", "POINT", test)

        arcpy.AddMessage("Data split into test and train sets")
        # =============================================================================
        # Polygon Mask and buffer
        # =============================================================================

        infc = dataPath + r"\points.shp"
        outfc = dataPath + r"\MinBoundGeom.shp"

        arcpy.management.MinimumBoundingGeometry(infc, outfc, "CONVEX_HULL", "ALL")

        #------------------------------------------------------------------------------

        infc = dataPath + r"\MinBoundGeom.shp"
        outfc = dataPath + r"\MinBoundGeomBuff.shp"
        buffer = str(int(parameters[3].valueAsText)) + " " + parameters[4].valueAsText
        line_side = "FULL"
        dissolve = "ALL"
        method = "PLANAR"

        arcpy.analysis.Buffer(infc, outfc, buffer, line_side = line_side, dissolve_option = dissolve, method = method)

        arcpy.AddMessage("Points have been buffered")
        # =============================================================================
        # Process Enviromental Data
        # =============================================================================
        # =============================================================================
        # Validate File Names        
        # =============================================================================
        
        arcpy.env.workspace = parameters[5].valueAsText

        rasters = arcpy.ListRasters("*", "TIF")

        for raster in rasters:
            inName = raster
            filename = arcpy.ValidateTableName(raster[:-4])
            outName = filename + ".tif"
            os.rename(os.path.join(arcpy.env.workspace, inName), os.path.join(arcpy.env.workspace, outName))
            
        arcpy.AddMessage("Ecological raster file names validated")
        
        # =============================================================================
        # Mask Bias
        # =============================================================================
        arcpy.env.workspace = parameters[5].valueAsText

        inRast = arcpy.ListRasters("*", "TIF")[0]
        mask = dataPath + r"\MinBoundGeomBuff.shp"

        arcpy.env.outputCoordinateSystem = arcpy.Describe(inRast).spatialReference
        arcpy.env.extent = mask
        arcpy.env.cellSize = inRast
        arcpy.env.mask = mask
        arcpy.env.snapRaster = inRast

        outExtractByMask = arcpy.sa.ExtractByMask(inRast, mask)

        outExtractByMask.save(parameters[7].valueAsText + r"\rasterMask.tif")
        
        arcpy.AddMessage("Buffer mask has been made")

        # =============================================================================
        # Process Rest of Rasters
        # =============================================================================
        
        mask = parameters[7].valueAsText + r"\rasterMask.tif"

        arcpy.env.outputCoordinateSystem = mask
        arcpy.env.extent = mask
        arcpy.env.cellSize = mask
        arcpy.env.mask = mask
        arcpy.env.snapRaster = mask

        arcpy.env.workspace = parameters[5].valueAsText

        rasters = arcpy.ListRasters("*", "TIF")

        index = 1
        for raster in rasters:
            outExtractByMask = arcpy.sa.ExtractByMask(raster, mask, extraction_area = "INSIDE")
            filePath = parameters[8].valueAsText
            outExtractByMask.save(filePath + "\\" + raster)
            arcpy.AddMessage("Raster " + str(index) + " complete out of " + str(len(rasters)))
            index += 1
            
        arcpy.AddMessage("All rasters masked")
        # =============================================================================
        # Raster to ASCII
        # =============================================================================
        
        arcpy.env.workspace = parameters[5].valueAsText

        rasters = arcpy.ListRasters("*", "TIF")

        index = 1
        for raster in rasters:
            filename = raster[:-4] + ".asc"
            filepath = parameters[9].valueAsText
            output = filepath + "\\" + filename
            arcpy.conversion.RasterToASCII(raster,output)
            arcpy.AddMessage("Raster " + str(index) + " converted out of " + str(len(rasters)))
            index += 1
            
        arcpy.AddMessage("All rasters converted to ASCII")
        
        # =============================================================================
        # Kernel Density
        # =============================================================================
        
        infc = parameters[6].valueAsText + r"\bias.shp"
        outputcell = parameters[7].valueAsText + r"\rasterMask.tif"
        outval = "DENSITIES"
        population = ""
        method = "PLANAR"

        arcpy.env.mask = outputcell
        arcpy.env.snapRaster = outputcell

        outKDens = arcpy.sa.KernelDensity(in_features = infc, population_field = population, 
                                          cell_size = outputcell, out_cell_values = outval, method = method)
        outfc = parameters[7].valueAsText + r"\gkd.tif"
        outKDens.save(outfc)
        
        arcpy.AddMessage("Gaussian Kernel Density raster created")

        return
