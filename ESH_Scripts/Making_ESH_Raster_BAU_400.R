#!/usr/bin/env Rscript
### ESH Script

#-------------------------------------------------------------------------------
# Name:        ESH for SSA
# Purpose:      Calculate average ESH lost by grid cell in SSA
#
# Author:      Mike Clark
#
# Created:     4/4/2018
# Copyright:   (c) Mike Clark
#
#-------------------------------------------------------------------------------


# importing libraries
library(raster)
library(gdalUtils)
library(rgdal)
library(dplyr)

# Data prep ----
# Getting template raster for SSA with row indices
# Making blank raster to add values to
# Importing template raster
ssa.raster <- raster("/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine/Sub-Saharan Africa.tif")
index.raster <- raster("/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine/Sub-Saharan AfricaRow.Index.tif")
# Making raster of row indices
index.vector <- getValues(index.raster)
index.vector2 <- 1:(index.raster@nrows * index.raster@ncols)
index.raster2 <- raster(matrix(index.vector2, nrow = index.raster@nrows, byrow = FALSE))
crs(index.raster2) <- crs(index.raster)
extent(index.raster2) <- extent(index.raster)
# Data prep ----
# Setting up data frame to loop over
# Importing data frame with tifs
# contains
# a) species and sub species identification
# b) full name and path of the raster
species.frame <- read.csv("/Users/maclark/Desktop/ESH Models/Birds/SSA_Birds/SSA_Birds_95percentinSSA.csv",
                          stringsAsFactors = FALSE)
# Importing file that has proportion of remaining ESH lost by species and subspecies
esh_frame <- read.csv("/Users/maclark/Desktop/SSA_Ag_Rasters/Outputs/AOH_Estimates/Plot.Frame_10April.csv",
                      stringsAsFactors = FALSE)
# Updating names for merging
names(esh_frame)[names(esh_frame) == "species_subspecies"] <- "Species"
# Merging in raster list
esh_frame <- left_join(esh_frame,
                       species.frame)
# Getting only 2060 data
yrs.list <- sort(unique(esh_frame$year))
# Getting list of scenarios
scenarios.list <- sort(unique(esh_frame$Scenario))

# Looping through years and scenarios
for(y in 11) {
  # Getting data frame for the specific year
  yrs_frame <- esh_frame[esh_frame$year == yrs.list[y],]
  # Looping through scenarios
  for(s in 3:5) {
    print(c(yrs.list[y], scenarios.list[s]))
    # Getting data frame for th especific scenario
    scenario_frame <- yrs_frame[yrs_frame$Scenario == scenarios.list[s],]
    scenario_frame <- scenario_frame[!is.na(scenario_frame$Percent_AOO_lost),]
    # Creating temporary vectors to hold
    # a) number of species that occur in the cell
    # b) total esh lost in the cell (summed across all species)
    species.vector <- rep(0, length(index.vector))
    esh_vector <- rep(0, length(index.vector))
    
    # Looping through species
    for(esh in 1:400) {
      # getting progress
      print(esh)
      # Getting species raster
      species.raster <- raster(scenario_frame$Raster[esh])
      
      # Cropping ag raster for a template to ressample
      tmp.index.raster <- crop(index.raster2,
                              species.raster)
      # Resampling species raster
      species.raster <- resample(species.raster, tmp.index.raster)
      # updating values in species raster
      # This is surprisingly fast
      species.raster[is.na(species.raster)] <- 0
      # And extending borders
      species.raster <- extend(species.raster, index.raster2)
      
      # Updating values
      species.vector[which(getValues(species.raster) > .95)] <- species.vector[which(getValues(species.raster) > .95)] + 1
      esh_vector[which(getValues(species.raster) > .95)] <- esh_vector[which(getValues(species.raster) > .95)] + scenario_frame$Percent_AOO_lost[esh]
      
      # removing temporary rasters
      rm(species.raster, tmp.index.raster)
      # removing tmp files periodically
      if(esh/25 == round(esh/25)) {
        rasterTmpFile(prefix = 'r_tmp_')
        # showTmpFiles()
        removeTmpFiles(.05)
      }
      
    }
    
    # Making rasters
    # Species first
    species.count.raster <- raster(matrix(species.vector, nrow = index.raster2@nrows, byrow = TRUE))
    crs(species.count.raster) <- crs(index.raster2)
    extent(species.count.raster) <- extent(index.raster2)
    
    # esh second
    # First getting average esh lost by species
    esh_vector[species.vector != 0] <- esh_vector[species.vector != 0] / species.vector[species.vector != 0]
    # Second converting into a raster
    esh_raster <- raster(matrix(esh_vector, nrow = index.raster2@nrows, byrow = TRUE))
    crs(esh_raster) <- crs(index.raster2)
    extent(esh_raster) <- extent(index.raster2)
    
    writeRaster(esh_raster,
                paste("/Users/maclark/Desktop/ESH_Scripts/Making_ESH_Raster/Rasters/ESH_Raster_birds_",
                      scenarios.list[s],
                      "_400_",
                      yrs.list[y],
                      ".tif",
                      sep = ""))
  }
}

quit(save = "no",
     status = 0,
     runLast = FALSE)
