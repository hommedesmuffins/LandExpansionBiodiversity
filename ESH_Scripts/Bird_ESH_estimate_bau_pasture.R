#!/usr/bin/env Rscript
### ESH Script

#-------------------------------------------------------------------------------
# Name:        ESH and agric
# Purpose:      Calculate ESH for species based on range polygons, altitudes
# and land cover associatsions.
#
# Author:      Otter
#
# Created:     17/11/2017
# Copyright:   (c) Otter 2016
#
#-------------------------------------------------------------------------------

# importing libraries
library(raster)
library(gdalUtils)
library(rgdal)

# Directory containing the ESH rasters
esh_folder = "/Users/maclark/Desktop/ESH Models/Birds/SSA_Birds_TIFs/"
# Directory containing ag rasters
agri_folder = "/Users/maclark/Desktop/SSA_Ag_Rasters/BAU/CropMeanExtent/"

# Adding fields----
# I think these creat a table in python
# Execute CreateTable

# arcpy.CreateTable_management(agri_folder, "agri_areas5.dbf")


# Add fields
# arcpy.AddField_management("/Users/maclark/Desktop/SSA_Ag_Rasters/BAU/agri_areas5.dbf", "species", "TEXT")
# arcpy.AddField_management("/Users/maclark/Desktop/SSA_Ag_Rasters/BAU/agri_areas5.dbf", "esh_area", "FLOAT")

# Creating table
esh.table <- data.frame(species = NA,
                        year = NA,
                        esh_area_orig = NA,
                        esh_loss = NA)

# Setting snap raster
Snap.Raster <- raster("/Users/maclark/Desktop/SSA_Ag_Rasters/Pasture_2010_SSA.tif")/100


# Getting list of crop rasters
raster.list <- list.files("/Users/maclark/Desktop/SSA_Ag_Rasters/BAU/",
                          pattern = "mean.tif",
                          full.names = TRUE)
raster.list <- raster.list[grep("Pasture",raster.list)]
raster.list <- raster.list[-grep("Delta", raster.list)]
# Getting 2010 raster file name
crop_2010 <- "/Users/maclark/Desktop/SSA_Ag_Rasters/Pasture_2010_SSA.tif"
# And appending to beginning of list
raster.list <- c(crop_2010,
                 raster.list)
# Getting list of eyars
yrs <- seq(2010, 2060, 5)

# Importing data frame with tifs
species.frame <- read.csv("/Users/maclark/Desktop/ESH Models/Birds/SSA_Birds/SSA_Birds_95percentinSSA.csv",
                          stringsAsFactors = FALSE)

# Looping through years
for(i in 1:length(yrs)) {
  # Checking progress
  print(c("Processing",yrs[i]))
  # Importing raster
  molagraster <- raster(raster.list[i])
  if(i == 1) {
    # Divide by 100 for 2010 data
    # 2010 data is in percent, not proportion of cell
    molagraster <- molagraster/100
  }
  
  # Getting list of .tif rasters
  esh.rasters = species.frame$Raster[grep(".tif",species.frame$Raster)]
  # And looping through species
  for(k in 1:length(esh.rasters)) {
    # Printing update
    if(k/100 == round(k/100)) {
      print(k)
    }
    # Importing esh raster
    species.raster <- raster(esh.rasters[k])
    # Cropping ag raster for a template to ressample
    tmp.crop.raster <- crop(molagraster,
                            species.raster)
    # Resampling species raster
    species.raster <- resample(species.raster, tmp.crop.raster)
    
    # Getting values of rasters
    esh.values <- getValues(species.raster)
    esh.values[is.na(esh.values)] <- 0
    esh.values[esh.values>.95] <- 1
    esh.values[esh.values <= .95] <- 0    
    crop.values <- getValues(tmp.crop.raster)
    # Getting sum of cells with value == 1 in species raster
    # This is value of total cells of suitable habitat without ag
    esh.zoneal <- sum(esh.values, na.rm = TRUE)
    
    # Now getting number of cells of suitable habitat
    # e.g. esh.zoneal - sum(ag) in cells where species.raster == 1
    esh.loss <- sum(crop.values[esh.values == 1], na.rm = TRUE)
    
    # Appending values to species.frame
    # First creating a data table
    tmp.esh.table <- data.frame(species = species.frame$Species[k],
                                year = yrs[i],
                                esh_area_orig = esh.zoneal,
                                esh_loss = esh.loss)
    
    # and updating original table
    esh.table <- rbind(esh.table,
                       tmp.esh.table)
    
    
    
    
  }
  # # Writing csv for 2060
  # if(yrs[i] == 2060) {
  #   esh.table.tmp <- esh.table
  #   esh.table.tmp$Scenario <- "BAU"
  #   write.csv(esh.table.tmp,
  #             "/Users/maclark/Desktop/SSA_Ag_Rasters/Outputs/ESH_Birds_BAU_Pasture_2060.csv",
  #             row.names = FALSE)
  # }
  # # Writing csv for 2010
  # if(yrs[i] == 2010) {
  #   esh.table.tmp <- esh.table
  #   esh.table.tmp$Scenario <- "All"
  #   write.csv(esh.table.tmp,
  #             "/Users/maclark/Desktop/SSA_Ag_Rasters/Outputs/ESH_Birds_BAU_Pasture_2010.csv",
  #             row.names = FALSE)
  # }
  # Need to add rows to the data table that contain the year, species id, and eoo
}


# Adding scenario identifier
esh.table$Scenario <- "BAU"
# And writing csv
write.csv(esh.table,
          "/Users/maclark/Desktop/SSA_Ag_Rasters/Outputs/ESH_Birds_BAU_Pasture_201020.csv",
          row.names = FALSE)

quit(save = "no", status = 0, runLast = TRUE)