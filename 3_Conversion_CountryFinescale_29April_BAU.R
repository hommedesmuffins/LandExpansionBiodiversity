#!/usr/bin/env Rscript

# Read me -----
# Using MC's model predictions to convert cells in each country to ag land
#
# This version has some major changes from previous versions. 

# Now we get the probability of expansion, contraction or no change in the first step 
#     (previously it was simply "change" vs. "no change")
# We then model expanding and contracting cells separately. 
#
#
# Rough procedure:
# 1) Set constants
# 2) Add urban change into the model outputs
# 3) Add the area of each cell, so that targets can be met 
# 4) Update the modelling data so that it is the latest year for which we have data 
#        (the model will be trained on the five-year period before)
# 5) Run the conversion loop. This version:
#     a) Tries to use targets from the baseline scenarios (email from MC 10th July)
#     b) Converts a number/proportion of cells going the 'wrong' way in each country at 
#        the start of each decade and adds the amount converted to the target
#     c) Then picks only cells going the right way (i.e. for expanding countries pick cells purely
#        based on their expansion probablility)
#     d) Loops multiple times for each country within each iteration to ensure that the 
#        targets can be reached
#     e) Currently does NOT include pasture

### MC Addendum
# update script to forecast land use by region, rather than the world at once
# did this by
# and then putting the whole thing in a loop, where each run of the loop forecasts expansion for a different region

# Packages ----
library(geosphere)
library(raster)
library(plyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(scales)
library(rgdal)
library(maptools)
library(caret)
library(pROC)
library(gbm)
library(e1071)
library(glmnet)
# library(doMC) # For parallel processing - comment out if not on a Mac

library(EBImage)
library(countrycode)
library(gdalUtils)
library(rgdal)
library(countrycode)
library(data.table)
library(matrixStats)

library(doMC)
registerDoMC(cores=5)

# Clean up ----
rm(list=ls())

# Creating counter for number of model iterations to run. 
# Doing this here so that it's easier to find
N.Iterations <- 5

# Who is running the code? ----
name <- "Mike"
# name <- "Dave"
# Load models, functions and data -----
if(name == "Dave"){
  source("Scripts/AnalysisScripts/0_Functions_CountryFinescale_16Feb.R")
  # List for GBM outputs
  file.list.mods_crop <- list.files(pattern = "Crop", path = "Data/FineScaleTests/")
  file.list.mods_pasture <- list.files(pattern = "Pasture", path = "Data/FineScaleTests/")
  # List for Regional data
  file.list.data <- list.files(pattern = "GBM", path = "Data/FineScaleTests/")
  # List for regional raster names
  reg.raster.names <- list.files(path = "Data/FineScaleTests/", pattern = "Row")
  # Urban data
  load("Data/FineScaleTests/UrbanScenarios_HiRes_working.rda")
  # Cells to contract (DRW note: not currently doing this)
  contraction_cells <- read.csv("Data/ModData/FromMC/Prop_Cells_Contract.csv",
                                header = TRUE,
                                stringsAsFactors = FALSE)
  expansion_cells <- read.csv("Data/ModData/FromMC/Prop_Cells_Expand.csv",
                              header = TRUE,
                              stringsAsFactors = FALSE)
} else {
  source("/Users/maclark/Desktop/SSA_Forecast_Runs/0_Functions_CountryFinescale_19March.R")
  # List for GBM outputs
  file.list.mods_crop <- list.files(pattern = "Crop_Coef", path = "/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log")
  file.list.mods_pasture <- list.files(pattern = "Pasture_Coef", path = "/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log")
  # Dropping out models with NOGDD
  # file.list.mods_crop <- file.list.mods_crop[grep("NOGDD", file.list.mods_crop)]
  # file.list.mods_pasture <- file.list.mods_pasture[grep("NOGDD", file.list.mods_pasture)]
  # List for Regional data
  file.list.data <- list.files(pattern = "GBM", path = "/Users/maclark/Desktop/Mike's Files/Forecast Data Files")
  # List for regional raster names
  reg.raster.names <- list.files(path = "/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine_Country/", pattern = "Row")
  # Urban data 
  # load("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/UrbanScenarios_HiRes_working.rda")
  # Cells to contract (DRW note: not currently doing this)
  # contraction_cells <- read.csv("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/PropCellsContract/Prop_Cells_Contract.csv",
  #                               header = TRUE,
  #                               stringsAsFactors = FALSE)
  # expansion_cells <- read.csv("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/PropCellsContract/Prop_Cells_Expand.csv",
  #                             header = TRUE,
  #                             stringsAsFactors = FALSE)
}

# Getting region list to assign names to the model object lists -----
region.names <- gsub("_Gamma_Link_Log_Multinom_Crop_Coefs.rda","",file.list.mods_crop)

# Starting loop -----
t1 <- Sys.time()
for(k in c(1:5,9:12,8,7,6)) {
  # Load targets -----
  # At the moment these are created in "1_DataPrep.R" 
  if(name == "Dave"){
    targets <- read.csv("Data/ModData/Targets_v4.csv",
                        header = TRUE,
                        stringsAsFactors = FALSE)  
    # Importing regional mod, regional data, calculating change in crop prop by cell to getemulate data set format of data_w_preds_crop, then joining with row.index.frame
    # Loading specific model output
    # First for crops
    load(paste("Data/FineScaleTests/",
               file.list.mods_crop[k],
               sep = ""))
  } else {
    targets <- read.csv("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/Land Clearing Targets/Targets_v4.csv",
                        header = TRUE,
                        stringsAsFactors = FALSE)  
    # Importing regional mod, regional data, calculating change in crop prop by cell to get emulate data set format of data_w_preds_crop, then joining with row.index.frame
    # Loading specific model output
    # First for crops
    load(paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log/",
               file.list.mods_crop[k],
               sep = ""))
  }
  # Put objects in list ----
  # Changing names of coefficients
  coefs.crop.increase <- coef.increase
  coefs.crop.decrease <- coef.decrease
  coefs.crop.prob <- coef.multinom
  
  # Removing models
  rm(coef.increase, coef.decrease, coef.multinom)
  
  # And repeating for pastures
  if(name == "Dave"){
    load(paste("Data/FineScaleTests/",
               file.list.mods_pasture[k],
               sep = ""))
  } else {
    load(paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log/",
               file.list.mods_pasture[k],
               sep = ""))  
  }
  
  # Changing names of models to turn into a list
  coefs.pasture.increase <- coef.increase
  coefs.pasture.prob <- coef.multinom
  
  # Removing models
  rm(coef.increase, coef.multinom)
  
  # Loading specifc regional data set -----
  if(name == "Dave"){
    region.data <- read.csv(paste("Data/FineScaleTests/",
                                  file.list.data[k],
                                  sep = ""),
                            stringsAsFactors = FALSE)
  } else {
    region.data <- read.csv(paste("/Users/maclark/Desktop/Mike's Files/Forecast Data Files/",
                                  file.list.data[k],
                                  sep = ""),
                            stringsAsFactors = FALSE)  
  }
  
  
  # ISO3 for Sudan (North part) is missing in our data frame
  # This is because countrycode() doesn't recognize the new iso numeric value (729) as valid
  # Manually updating the iso3n so that it matches with countrycode()
  if(k == 11) {
    region.data$ISO.Numeric[region.data$ISO.Numeric == 729] <- 736
  }
  
  
  # #eliminating row.name column by name, in case order changes
  region.data <- dplyr::select(region.data, -X)
  # Also dropping continuous GDD data and 1st percentile GDD cutoff
  region.data <- dplyr::select(region.data, -GDD.10c)
  region.data <- dplyr::select(region.data, -GDD.1st)
  
  # Emulating layout of data_w_preds_crop
  region.data$prop_to_predict_from_cell_delta_crop <- region.data$Prop.Crop.2010 - region.data$Prop.Crop.2005
  region.data$prop_to_predict_from_cell_delta_pasture <- region.data$Prop.Pasture.2010 - region.data$Prop.Pasture.2005
  
  # Updating value of IUCN.Region.Name
  region.data$IUCN.Region.Name <- region.names[k]
  # Dropping cells for which Iucn.Region is NA
  region.data$IUCN.Region.Name[is.na(region.data$IUCN.Region)] <- NA
  
  # Adding in ISO3 codes to emulate layout of data_w_preds_crop
  region.data$ISO3 <- countrycode(region.data$ISO.Numeric, "iso3n","iso3c")
  
 
  # Renaming columns to match col names in script below
  region.data <- dplyr::rename(region.data,
                        ag_suitability = Ag.suitability,
                        dist_50k = Dist.50k,
                        PA_binary = PA.Binary,
                        GDD_binary = GDD.2.5th,
                        prop_to_predict_from_crop = Prop.Crop.2010,
                        prop_to_predict_from_pasture = Prop.Pasture.2010,
                        prop_adj_to_predict_from_crop = Prop.AdjacentCrop.2010,
                        prop_adj_to_predict_from_pasture = Prop.AdjacentPasture.2010,
                        change_demand_crop = CountryIncrease.Crop.2005.2010,
                        IUCN.Region.Number = IUCN.Region)
  
  
  # Eliminating rounding errors in MODIS data
  # Some cells have total land cover > 1
  # For these cells, dividing proportion of each land cover type by sum of land cover
  region.data$tmp <- region.data$prop_to_predict_from_crop +
    region.data$prop_to_predict_from_pasture + 
    region.data$Prop.NonArable.2010 +
    region.data$Prop.Urban.2010
  # Getting index of cells that have land cover > 1
  tmp_index <- which(region.data$tmp > 1 & !is.na(region.data$tmp))
  # Counter for number of digits to round to
  tmp_counter = 8
  # putting in a while loop
  while(length(tmp_index) > 1 & tmp_counter >= 2) {
    # And dividing land cover of each type by sum of all land cover types for these cells
    region.data$prop_to_predict_from_crop[tmp_index] <- round(region.data$prop_to_predict_from_crop[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
    region.data$prop_to_predict_from_pasture[tmp_index] <- round(region.data$prop_to_predict_from_pasture[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
    region.data$Prop.NonArable.2010[tmp_index] <- round(region.data$Prop.NonArable.2010[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
    region.data$Prop.Urban.2010[tmp_index] <- round(region.data$Prop.Urban.2010[tmp_index] / region.data$tmp[tmp_index], digits = tmp_counter)
    # And updating region.data$tmp
    region.data$tmp <- region.data$prop_to_predict_from_crop +
      region.data$prop_to_predict_from_pasture + 
      region.data$Prop.NonArable.2010 +
      region.data$Prop.Urban.2010
    # and updating tmp_index
    tmp_index <- which(region.data$tmp > 1 & !is.na(region.data$tmp))
    # Updating counter
    tmp_counter = tmp_counter - 1
  }
  
  # Removing tmp_counter and tmp_index
  rm(tmp_counter, tmp_index)
  
  # And removing column tmp from region.data
  region.data <- dplyr::select(region.data, -tmp)
  
  # And recreating values for prop.crop.2010 and prop.pasture.2010
  region.data$Prop.Pasture.2010 <- region.data$prop_to_predict_from_pasture
  region.data$Prop.Crop.2010 <- region.data$prop_to_predict_from_crop
  
  
  
  # Creating iso3c list to import specific row index rasters for each nation -----
  iso3c <- data.frame(iso3c = sort(unique(region.data$ISO3)))
  # Limiting iso3c to only countries for which we have forecasts and converting back into a list
  iso3c <- iso3c[iso3c$iso3c %in% targets$ISO3,]
  
  # Limiting country selection to only those that have cropland change in at least one time interval -----
  targets.tmp = targets
  # Getting targets for scenario we're forecasting
  targets.tmp = targets.tmp[!is.na(targets.tmp$ISO3) & targets.tmp$scenario == "Baseline",]
  
  targets.tmp$Count <- 0
  # Getting countries with no change in crop extent in any time period
  # Doing this to not loop through countries that don't have any change in cropland
  # Looping through time period
  
  for(i in 3:12) {
    # Creating vector of target year
    counter.vector = targets.tmp[,i]
    # Getting rows that don't change
    rows.nochange = which(counter.vector == 0)
    # Adding one to counter
    targets.tmp$Count[rows.nochange] <- targets.tmp$Count[rows.nochange] + 1
  }
  
  # Limiting to only countries that change cropland extent in any time period
  targets.tmp <- targets.tmp[targets.tmp$Count != 10,]
  # And limiting list of countries we're looping through
  # Not doing this anymore now that we're forecasting change in urban extent
  #iso3c <- iso3c[iso3c %in% targets.tmp$ISO3]
  
  # Creating a regional urban data set to avoid joining with a 290 million row data table
  # urban_region <- urban_data[urban_data$ISO3.Numeric %in% countrycode(iso3c, "iso3c", "iso3n"),]
 
  
  # Looping with specific countries -----
  for(c in 1:length(iso3c)) {
    print(iso3c[c])
    # Targets 
    # At the moment these are created in "1_DataPrep.R" 
    if(name == "Dave"){
      targets <- read.csv("Data/ModData/Targets_v4.csv",
                          header = TRUE,
                          stringsAsFactors = FALSE)      
      # Getting name of raster to import
      row.index.name <- list.files(path = "Data/FineScaleTests/",
                                   pattern = "Row")
      # Importing country raster with row.index values
      row.index.raster <- raster(paste("Data/FineScaleTests/",
                                       row.index.name,
                                       sep = ""))
    } else {
      targets <- read.csv("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/Land Clearing Targets/Targets_v4.csv",
                          header = TRUE,
                          stringsAsFactors = FALSE)
      # Getting name of raster to import
      row.index.name <- list.files(path = "/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine_Country",
                                   pattern = "Row")
      row.index.name <- row.index.name[grep(iso3c[c],row.index.name)]
      # Importing country raster with row.index values
      row.index.raster <- raster(paste("/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine_Country/",
                                       row.index.name,
                                       sep = ""))
    }
    
    # Creating row.index.frame that contains only values of row.index associated with the countri-specific raster
    row.index.frame <- data.frame(Row.Index = cbind(c(getValues(row.index.raster))))
    
    # Merging with row.index.frame
    # Limits data frame to only country-level cells
    model_data_all <- left_join(row.index.frame,
                                region.data)
    
    # All non PA areas should have value of one. 
    # This apparently isn't the case for all countries, so inserting manually here
    model_data_all$PA_binary[-which(model_data_all$PA_binary == 1)] <- 0
    # Likewise for GDD data
    model_data_all$GDD_binary[-which(model_data_all$GDD_binary == 1)] <- 0
    
    # Dropping all data not from the country
    row.nums <- grep(iso3c[c], model_data_all$ISO3)
    model_data_all[-row.nums,-1] <- NA
    
    # Removing row.index.frame
    rm(row.index.frame, row.nums)
    
    # 1) Put data into one dataframe and rename if necessary ------
    # Currently not combining pasture and crop as I don't have the pasture data
    # model_data_all <- full_join(data_w_preds_crop,
    #                             data_w_preds_pasture)
    # Rename and convert factors to characters
    
    
    # Changes factor variables into character strings
    model_data_all <- model_data_all %>%
      mutate_if(is.factor, as.character)
    
    # Rename some columns so that they work with the script (easier to change
    #     names than the script)
    model_data_all <- model_data_all %>%
      dplyr::rename(IUCN_region = IUCN.Region.Name,
             row_index = Row.Index, 
             ISO_numeric = ISO.Numeric,
             prop_after_change_crop = Prop.Crop.2010,
             prop_after_change_pasture = Prop.Pasture.2010) %>%
      dplyr::select(row_index, ISO_numeric, IUCN.Region.Number:prop_after_change_crop) %>%
      as.data.frame()
    
    # Order by row_index and rename the rows - this is for checking things later
    model_data_all <- model_data_all[order(model_data_all$row_index),]
    row.names(model_data_all) <- 1:nrow(model_data_all)
    
    # Explanation: Names of variables -----
    # (taken from MC's email, NOTE I'VE CHANGED THESE)
    
    # prop_to_predict_from_crop: proportion of cell in cropland in first decade (originally: Prop.1980)
    # prop_change_cell: Difference in proportion between first and following decade (originally: Delta.Prop)
    # dist_50k: Distance from urban areas (originally: Dist.50k)
    # ag_suitability: Ag suitability index (originally: Ag.suitability)
    # PA_prop: proportion of cell in PAs (originally: PA.mean)
    # has_changed: If a cell has experienced change in crop extent between first and following decades 
    #              Contains values of "No" or "Yes" (originally: Has.Changed)
    # pred_no_change: Probability that a cell has not experienced a change in cropland extent (originally: No)
    # pred_to_change:  Probability that a cell has experienced a change in cropland extent (originally: Yes)
    # likely_to_change: NA if "pred_no_change" and "pred_to_change" are NA; 
    #                   Yes if value in "pred_to_change" is greater than value in "pred_no_change"; 
    #                   No if value in "pred_no_change" is greater than value in "pred_to_change"
    #                   (originally: Pred.Expansion)
    # correct_pred: "YES" if value in has_changed and likely_to_change are the same; 
    #               "NO" if not the same; NA if values in columns are NAs. (originally: Correct.Pred)
    # prop_change_cell: Predicted change in cropland proportion in a cell (originally: predictionsREGION)
    
    
    # Data prep: Set constants: -----
    # predictors, targets, number of model runs for the loop, resolution data, 
    #     convolution filter data
    # Predictors
    predictor_vars_crop <- c('dist_50k_log' ,
                             'dist_50k_log_sq' ,
                             'ag_suitability' ,
                             'PA_binary' ,
                             'prop_adj_to_predict_from_crop' ,
                             'prop_adj_crop_sq' ,
                             'prop_to_predict_from_crop' ,
                             'prop_crop_sq' ,
                             'prop_to_predict_from_pasture' ,
                             'prop_pasture_sq' ,
                             'prop_to_predict_from_cell_delta_crop' ,
                             'prop_cell_delta_crop_sq' ,
                             'ISO_numeric')

    predictor_vars_pasture <- c('dist_50k_log' ,
                                'dist_50k_log_sq' ,
                                'ag_suitability' ,
                                'PA_binary' ,
                                'prop_to_predict_from_crop' ,
                                'prop_crop_sq' ,
                                'prop_to_predict_from_pasture' ,
                                'prop_pasture_sq' ,
                                'prop_to_predict_from_cell_delta_pasture' ,
                                'prop_cell_delta_pasture_sq' ,
                                'prop_adj_to_predict_from_pasture' ,
                                'prop_adj_pasture_sq' ,
                                'ISO_numeric')
    
    # Set scenarios to investigate
    unique(targets$scenario)
    targets <- subset(targets, scenario == "Baseline")
    
    # Model runs - DRW keeping this at 1 for the moment
    N <- N.Iterations
    # model_runs <- subset(urban_scenarios,
    #                      run <= N)
    # model_runs <- model_runs[,c("run", "ISO3.Numeric", 
    #                             names(model_runs)[grep("urb_targ",
    #                                                    names(model_runs))])]
    # 
    # # Changing names of model_runs ISO3.Numeric to ISO_numeric to match col names in model_data_current
    # names(model_runs)[grep("ISO3.Numeric", names(model_runs))] <- "ISO_numeric"
    
    rm(N)
    # Get resolution and area data
    # This could probably be moved outside of the loop
    # Leaving here to minimize number of changes to script
    # area_raster <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/MollweideCountryID_1.5km.tif")
    # area_raster <- raster("/Users/Me/Documents/Datasets/WorldMaps/UN_Numeric_Raster_Mollweide_60km.tif")
    
    # Dimensions are in metres, so convert to square km
    cell_area <- res(row.index.raster)[1]/1000 * res(row.index.raster)[2]/1000
    rows <- nrow(row.index.raster)
    cols <- ncol(row.index.raster)
    
    # Convolution filter  
    adj_cell_mat <- matrix(nrow=3, ncol = 3, data = 1)
    
    # Data prep: Add probabilities of urban change -----
    # (nb. This will change when I tidy things up) 
    # model_data_all$urban_prob <- urban_region$urban_prob[match(model_data_all$row_index,
    #                                                            urban_region$Row.Index)]
    
    # rm(urban_dat)
    # Data prep: Convert data to NAs from countries that we are not forecasting -----
    # We are forecasting country by country, so this removes any data associated with any country that is not the target country
    # Getting list of row numbers associated with the specific target country
    row.nums <- grep(iso3c[c], model_data_all$ISO3)
    # Converting all other rows to NAs
    model_data_all[-row.nums,-1] <- NA
    
    # Data prep: Add area of each cell ----
    model_data_all$cell_area <- cell_area
    # Put NAs in, so it throws an error if necessary
    model_data_all[is.na(model_data_all$prop_to_predict_from_crop), "cell_area"] <- NA
    rm(row.nums)
    
    # Data prep: Update the proportion of ag land with new values -----
    # Update "prop_to_predict_from" to current values
    # Moves next time step data into current time step data to get predictions
    model_data_current <- model_data_all
    model_data_current$prop_to_predict_from_crop <- model_data_current$prop_after_change_crop
    model_data_current$prop_to_predict_from_pasture <- model_data_current$prop_after_change_pasture
    # Add in a reference column for proportion of ag at the start of the loop
    model_data_current$prop_start_loop_crop <- model_data_current$prop_to_predict_from_crop
    model_data_current$prop_start_loop_pasture <- model_data_current$prop_to_predict_from_pasture
    
    # Updating columns to reflect code in lines 408 - 412
    model_data_current$prop_change_cell_crop <- model_data_current$prop_to_predict_from_cell_delta_crop
    model_data_current$prop_change_cell_pasture <- model_data_current$prop_to_predict_from_cell_delta_pasture
    
    
    # Converting iso_numeric to an ordered factor variable
    # model_data_current$ISO_numeric = factor(model_data_current$ISO_numeric, ordered = FALSE)
    # 
    # # And adding in new predictor variables
    # # log(dist), log(dist)^2, etc
    # model_data_current$dist_50k_log <- log(model_data_current$dist_50k)
    # model_data_current$dist_50k_log_sq <- model_data_current$dist_50k_log^2
    # 
    # model_data_current$prop_crop_sq <- model_data_current$prop_to_predict_from_crop^2
    # model_data_current$prop_pasture_sq <- model_data_current$prop_to_predict_from_pasture^2
    # 
    # model_data_current$prop_adj_crop_sq <- model_data_current$prop_adj_to_predict_from_crop^2
    # model_data_current$prop_adj_pasture_sq <- model_data_current$prop_adj_to_predict_from_pasture^2
    # 
    # model_data_current$prop_cell_delta_crop_sq <- model_data_current$prop_to_predict_from_cell_delta_crop^2
    # model_data_current$prop_cell_delta_pasture_sq <- model_data_current$prop_to_predict_from_cell_delta_pasture^2
    # 
    # Some variables go to NAs, this messes things up. Sets everything to NA if a lake or something equivalent
    # Cut out lakes and other areas where not all predictors are present 
    model_data_current[is.na(model_data_current$prop_to_predict_from_crop),
                       c('dist_50k', 'ag_suitability', 'PA_binary',
                         'prop_adj_to_predict_from_crop', 'prop_adj_to_predict_from_pasture',
                         'prop_to_predict_from_crop', 'prop_to_predict_from_pasture',
                         'prop_to_predict_from_cell_delta_crop', 'prop_to_predict_from_cell_delta_pasture',
                         'ISO_numeric',
                         "prop_change_cell_crop", #"prop_change_cell_pasture",
                         "prop_after_change_crop", "prop_after_change_pasture",
                         "cell_area")] <- NA
    
    # Need to run convolution again to update area of cropland in each cell, but have not adjusted adjacent area
    # Can change this to manually update 2007 adjacent crop/pasture with 2012 adjacent crop/pasture
    # Run convolution to get new surrounding values for cropland...
    model_data_current <- adj.ag.f(m = model_data_current, 
                                   r = rows, 
                                   c = cols, 
                                   f = adj_cell_mat, 
                                   crop = "crop")
    
    # ...and pasture
    model_data_current <- adj.ag.f(m = model_data_current, 
                                   r = rows, 
                                   c = cols, 
                                   f = adj_cell_mat, 
                                   crop = "pasture")
    
    # Data prep: Calculate targets based on MODIS data ----
    # Calculate the area of cropland in each country
    # First part before #add ISO3 sums area in cropland in 2012
    # Mutate converts from ISO3 numeric to ISO3 characters
    # Left_join joins current crop area with estimates of proportional increase in cropland
    targets <- model_data_current %>%
      dplyr::select(ISO_numeric, ISO3, prop_start_loop_crop, cell_area) %>%
      mutate(area_cropland_km = prop_start_loop_crop * cell_area) %>%
      dplyr::group_by(ISO_numeric) %>%
      plyr::summarise(area_cropland_km_2010 = sum(area_cropland_km[!is.na(area_cropland_km)])) %>%
      # Add ISO3, so we can bind to targets
      mutate(ISO3 = as.character(iso3c[c])) %>%
      # Bind to proportional increases in cropland
      left_join(., targets %>%
                  dplyr::select(-grep("area", names(.)))) %>%
      as.data.frame()
    
    # Calculate area increases
    # Do first area separately and then loop
    # Multiplies proportional increase by current area in cropland
    targets$target_2010_2015 <- targets$area_cropland_km_2010 * targets$prop_change_2010_2015
    targets$area_cropland_km_2015 <- targets$area_cropland_km_2010 + targets$target_2010_2015
    
    # Loop calculates future targets of cropland by 5 year interval
    for(i in 2:length(grep("prop_change", names(targets)))){
      # Get the interval and variable names
      interval <- names(targets)[grep("prop_change", names(targets))][i]
      n_target <- gsub("prop_change", "target", interval)
      n_area <- paste("area_cropland_km", sub("^.*_", "", n_target), sep = "_")
      n_previous <- paste("area_cropland_km", 
                          gsub("(^.+_)(.*)(_.+$)", "\\2", interval), 
                          sep = "_")
      targets[[n_target]] <- targets[,names(targets)[grep("prop_change", names(targets))][i]] * 
        targets[[n_previous]]
      targets[[n_area]] <- targets[[n_target]] + targets[[n_previous]]
    }
    
    # Can be moved outside of loop
    # Leaving here for now
    # Data prep: Add urban land ----
    # urban_2010 <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 50km/Urban2010.Mollweide.1.5km.tif_Coarse_2Nov.tif")
    # urban_2010 <- raster("/Users/Me/Documents/Datasets/MODIS/Urban2010.Mollweide.5km.tif")
    
    # MC provided a raster in the right resolution, so just add it to the data
    # urban.frame <- data.frame(urban_2010 = cbind(c(getValues(urban_2010))),
    #                           Row.index = 1:259200)
    # urban.frame <- urban.frame[urban.frame$Row.index %in% model_data_current$row_index,]
    
    # Put in baseline (2010) urban proportion
    # model_data_current$urban_2010 <- as.vector(urban.frame$urban_2010) / 100
    
    # This column will get changed during the forecasts to be updated with proportion of cells in urban in future time steps
    # Get dynamic column for current urban proportion
    # model_data_current$urban_prop <- model_data_current$urban_2010
    
    # If you don't have 0s at some point, the code will break
    # Substituting non ISO3s with NAs so that code throws an error
    # Replace NAs with 0 for adding/subtracting, but restrict to cells with data
    # model_data_current[is.na(model_data_current$urban_prop),
    #                    "urban_prop"] <- 0
    # model_data_current[is.na(model_data_current$ISO_numeric),
    #                    "urban_prop"] <- NA
    
    
    # tmp <- raster(matrix(model_data_current$PA_binary,
    #                      nrow = rows))
    # plot(tmp)
    # rm(urban_2010)
    
    # Data prep: Add in proportional land demand ----- 
    # Adds in targets to the big data frame
    model_data_current <- left_join(model_data_current, 
                                    targets[,-grep("area_cropland", names(targets))])
    
    # At this point we have a data frame with all of projection variables updated to 2010
    # Also contains amount of urban land
    # Also contains forecasts of future crop demand for 5 year intervals
    
    # Get new probabilities for a) change b) amount of change -----
    # Cut regions without predictions
    na.data <- model_data_current[is.na(model_data_current$ISO_numeric),]
    model_data_current <- model_data_current[!is.na(model_data_current$ISO_numeric),]
    
    # Add in new demand
    # Updates predictor variable for model forecasts
    # Changed the next line because Dave told me to
    model_data_current$change_demand_crop <- model_data_current$prop_change_2010_2015 + 1
    
    # Creating column named urban_prop to mesh with the functions script
    # Might be able to move this further up when renaming columns
    model_data_current$urban_prop <- model_data_current$Prop.Urban.2010

    # a) Update cropland
    # Calls from pred_fun in 0_Functions
    model_data_current <- model_data_current %>%
      do(pred_fun(p = ., crop = "crop")) %>% 
      as.data.frame()
    
    # Add columns for pasture loss - currently not doing this -----
    # What this will do is take model data, and calculate how much space is available for pasture
    model_data_current$prop_for_pasture <- round(1 - rowSums(model_data_current[,c('urban_prop',
                                                                             'prop_to_predict_from_crop',
                                                                             'Prop.NonArable.2010')],
                                                       na.rm = TRUE),
                                                 digits = 7)
    # Not dealing with pasture outside of country borders...
    model_data_current$prop_for_pasture[is.na(model_data_current$ISO_numeric)] <- NA
    # And really making sure it isn't negative
    # Although the rounding above should take care of that
    model_data_current$prop_for_pasture[model_data_current$prop_for_pasture < 0] <- 0
    
    # b) Update pasture
    # model_data_current <- model_data_current %>%
    #   group_by(IUCN_region) %>%
    #   do(pred_fun(., crop = "pasture")) %>%
    #   as.data.frame()
    
    # Check that nothing mad is being predicted -----
    model_data_current$tmp <- 
      model_data_current$prop_to_predict_from_crop + 
      model_data_current$urban_prop + 
      # model_data_current$prop_to_predict_from_pasture + 
      model_data_current$pred_change_cell_prop_crop_increase
    # range(model_data_current$tmp, na.rm = TRUE)
    
    model_data_current$tmp <- 
      model_data_current$prop_to_predict_from_crop + 
      # model_data_current$prop_to_predict_from_pasture + 
      model_data_current$pred_change_cell_prop_crop_decrease
    # range(model_data_current$tmp, na.rm = TRUE)
    
    model_data_current <- dplyr::select(model_data_current, -tmp)
    
    # c) Update areas without data
    # Adding columns to the NA data set to mesh with the model_data_current data
    na.data$pred_no_change_crop <- NA
    na.data$pred_to_change_crop <- NA
    na.data$area_converted_crop <- NA
    na.data$pred_no_change_pasture <- NA
    na.data$pred_to_change_pasture <- NA
    na.data$area_converted_pasture <- NA
    
    # d) Bind data back together
    model_data_current <- bind_rows(model_data_current, na.data)
    model_data_current <- model_data_current[order(model_data_current$row_index),]
    rm(na.data)
    
    # Limiting columns that not needed for the forecast script
    model_data_current <- dplyr::select(model_data_current,
                                        row_index,
                                        ISO3,
                                        ISO_numeric,
                                        prop_to_predict_from_crop,
                                        prop_to_predict_from_pasture,
                                        Prop.NonArable.2010,
                                        urban_prop,
                                        prop_adj_to_predict_from_crop,
                                        prop_adj_to_predict_from_pasture,
                                        change_demand_crop,
                                        dist_50k,
                                        ag_suitability,
                                        PA_binary,
                                        GDD_binary,
                                        prop_to_predict_from_cell_delta_crop,
                                        prop_to_predict_from_cell_delta_pasture,
                                        pred_to_increase_crop,
                                        pred_to_decrease_crop,
                                        pred_change_cell_prop_crop_decrease,
                                        pred_change_cell_prop_crop_increase,
                                        area_converted_crop_decrease,
                                        area_converted_crop_increase,
                                        prop_for_pasture,
                                        prop_change_2010_2015,
                                        prop_change_2015_2020,
                                        prop_change_2020_2025,
                                        prop_change_2025_2030,
                                        prop_change_2030_2035,
                                        prop_change_2035_2040,
                                        prop_change_2040_2045,
                                        prop_change_2045_2050,
                                        prop_change_2050_2055,
                                        prop_change_2055_2060,
                                        target_2010_2015,
                                        target_2015_2020,
                                        target_2020_2025,
                                        target_2025_2030,
                                        target_2030_2035,
                                        target_2035_2040,
                                        target_2040_2045,
                                        target_2045_2050,
                                        target_2050_2055,
                                        target_2055_2060)
                                        
                                        
                                        
                                        
                                        
    
    
    # Which countries will have a pasture deficit? UNFINISHED -----
    # crop_pasture_totals <- model_data_current
    # crop_pasture_totals$available_for_crop <- crop_pasture_totals$urban_prop
    # crop_pasture_totals[crop_pasture_totals$ag_suitability == 0 & 
    #                       !is.na(crop_pasture_totals$ag_suitability), "available_for_crop"] <- 0
    # 
    # crop_pasture_totals <- crop_pasture_totals %>%
    #   group_by(ISO3) %>%
    #   summarise(crop_total_2010 = sum(prop_to_predict_from_crop * cell_area, na.rm = TRUE),
    #             pasture_total = sum(prop_to_predict_from_pasture * cell_area, na.rm = TRUE),
    #             area_suit_crop = sum(available_for_crop * cell_area, na.rm = TRUE),
    #             area_suit_pasture = sum(urban_prop * cell_area, na.rm = TRUE))
    
    # Plotting to verify things didn't get messed up
    # Doing this by recreating the regional iso3 map
    # temp.dat <- dplyr::select(model_data_current, 
    #                           row_index, 
    #                           ISO_numeric)
    
    # Merging with region.frame
    # country.frame <- data.frame(row_index = cbind(c(getValues(row.index.raster))))
    
    # country.frame <- left_join(country.frame,
    #                           temp.dat)
    # Turning into a matrix
    # country.matrix <- matrix(country.frame$ISO_numeric,
    #                         nrow = row.index.raster@nrows,
    #                         byrow = TRUE)
    # 
    # converting to raster
    # country.raster <- raster(country.matrix)
    # crs(country.raster) <- crs(row.index.raster)
    # extent(country.raster) <- extent(row.index.raster)
    # res(reg.raster) <- res(temp.raster)
    
    # plotting raster
    # plot(country.raster)
    
    # Run the future conversion ----
    # Run the conversion
    # model_runs is the urban scenarios
    # t1 <- Sys.time()
    
    # Setting seed for tourbleshooting to get reproducible results
    #set.seed(1e4)
    
    # Testing quick raster for PA_binary and GDD_binary
    # Used for troubleshooting
    # propincrease.raster <- raster(matrix(model_data_current$pred_to_increase_crop, nrow = row.index.raster@nrows))
    # plot(propincrease.raster)
    # 
    # crop.raster <- raster(matrix(model_data_current$prop_to_predict_from_crop, nrow = row.index.raster@nrows))
    # crop.raster[crop.raster>.1] <-1
    # plot(crop.raster)
    # 
    # suit.raster <- raster(matrix(model_data_current$ag_suitability, nrow = row.index.raster@nrows))
    # plot(suit.raster)
    # 
    # gdd.raster <- raster(matrix(model_data_current$GDD_binary, nrow = row.index.raster@nrows))
    # plot(gdd.raster)
    
    # Subsetting model runs to run == 1
    # Model runs is for urban expansion
    # We're no longer including urban expansion
    # So this does not impact forecast of crop/pasture
    # model_runs <- subset(model_runs, run == 1)
    
    # # Making duplicates of model_runs
    # for(n in 2:N.Iterations) {
    #   tmp.runs <- model_runs
    #   tmp.runs$run <- n
    #   model_runs <- rbind(model_runs,
    #                       tmp.runs)
    # }
    # t1 = Sys.time()
    # print("Starting Conversion Loop")
    # test.list <- dlply(model_runs, 
    #                    .(run),
    #                    .parallel = TRUE,
    #                    year_fun)#,
    # .parallel = TRUE)
    
    # Removing objects we don't need
    rm(model_data_all)
    
    # Sys.time() - t1 # time varies by country
    # Run this if using parallel
    test.list <- c()
    # Making exceptions for geographically large countries
    # Need to do this becasue some countries will take all of the RAM available on the computer if doing 5 iterations
    # And this won't work...as in the computer will turn off
    if(iso3c[c] == 'IND' | iso3c[c] == 'BRA') {
      # First set of forecasts
      # Total of four
      test.list.tmp <- foreach(i = 1:4) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Second set of forecasts
      # Total of eight
      test.list.tmp <- foreach(i = 1:4) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Third set of forecasts
      # Total of twelve
      test.list.tmp <- foreach(i = 1:4) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Fourth set of forecasts
      # Total of sixteen
      test.list.tmp <- foreach(i = 1:4) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Fifth rep of forecasts
      # Total of nineteen
      test.list.tmp <- foreach(i = 1:3) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Sixth rep of forecasts
      # Total of twenty-two
      test.list.tmp <- foreach(i = 1:3) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
      # Seventh rep of forecasts
      # Total of twenty-five
      test.list.tmp <- foreach(i = 1:3) %dopar% year_fun(data.frame(input = "model_runs"))
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
    } else if(iso3c[c] == 'AUS' | iso3c[c] == 'CHN' | iso3c[c] == 'CAN') {
      # Looping through first twelve
      # For total of twenty four model runs
      for(n.iter in 1:12) {
        print(n.iter)
        test.list.tmp <-  foreach(i = 1:2) %dopar% year_fun(data.frame(input = "model_runs"))
        # And appending forecasts to a list
        test.list <- c(test.list, test.list.tmp)
        rm(test.list.tmp)
      }
      # ANd doing the twenty-fifth run
      test.list.tmp <-  foreach(i = 1:1) %dopar% year_fun(data.frame(input = "model_runs"))
      # And appending forecasts to a list
      test.list <- c(test.list, test.list.tmp)
      rm(test.list.tmp)
    } else if(iso3c[c] == 'RUS' | iso3c[c] == 'USA')  {
      print(iso3c[c])
      # Looping through all twenty five iterations
      for(n.iter in 1:25) {
        # if(n.iter == 1) {
        #   test.list <- list()
        # }
        print(n.iter)
        test.list.tmp <-  foreach(i = 1:1) %dopar% year_fun(data.frame(input = "model_runs"))
        # And appending forecasts to a list
        test.list <- c(test.list, test.list.tmp)
        rm(test.list.tmp)
      }
    } else {
      for(n.iter in 1:5) {
        print(c(iso3c[c], n.iter))
        test.list.tmp <-  foreach(i = 1:5) %dopar% year_fun(data.frame(input = "model_runs"))
        # And appending forecasts to a list
        test.list <- c(test.list, test.list.tmp)
        rm(test.list.tmp)
      }
    }    
    
    # print("Conversion Loop Done")
    # Sys.time() - t1
    
    # Need to put this into a matrix  
    # Creating matrix that contains row index, mean cell value, and sd cell value
    # Doing this through a function that is in the 0_Functions script
    # First getting names of columns we need to take mean and sd of
    # Have to do independently for crop and pasture
    column.names.crop = names(test.list[[1]]$LandUse)[grep("prop_crop_20", names(test.list[[1]]$LandUse))]
    # column.names.crop = column.names.crop[-grep("pred", column.names.crop)]
    column.names.pasture = names(test.list[[1]]$LandUse)[grep("prop_pasture_20", names(test.list[[1]]$LandUse))]
    # column.names.pasture = column.names.pasture[-grep("pred", column.names.pasture)]
    # Running functions that creates rasters
    raster.list.crop <- mean.sd.fun(column.names = column.names.crop)
    raster.list.pasture <- mean.sd.fun(column.names = column.names.pasture)
    
    
    # bootstrapping for every 5 countries
    # if(c/5 == round(c/5)) {
    #   bootstrap.frame1 <- bootstrap.sd.function(column.names = column.names.pasture)
    #   bootstrap.frame1 <- reshape(bootstrap.frame1,
    #                               varying = names(bootstrap.frame1)[2:11],
    #                               timevar = 'Interval',
    #                               times = names(bootstrap.frame1)[2:11],
    #                               v.names = 'StandardDeviation',
    #                               direction = 'long')
    #   
    #   bootstrap.frame1 <- bootstrap.frame1[complete.cases(bootstrap.frame1),]
    #   
    #   ggplot(bootstrap.frame1, aes(group = N.Iterations, x = N.Iterations, y = StandardDeviation)) +
    #     geom_boxplot() +
    #     facet_wrap(~Interval)
    #   
    #   ggsave(paste("/users/maclark/Desktop/SSA_Forecast_Runs/Yields/Bootstrap/",
    #                iso3c[c],
    #                "_Pasture_Bootstrap.pdf",
    #                sep = ""),
    #          width = 20,
    #          height = 10)
    #   
    #   bootstrap.frame1 <- bootstrap.sd.function(column.names = column.names.crop)
    #   bootstrap.frame1 <- reshape(bootstrap.frame1,
    #                               varying = names(bootstrap.frame1)[2:11],
    #                               timevar = 'Interval',
    #                               times = names(bootstrap.frame1)[2:11],
    #                               v.names = 'StandardDeviation',
    #                               direction = 'long')
    #   
    #   bootstrap.frame1 <- bootstrap.frame1[complete.cases(bootstrap.frame1),]
    #   
    #   ggplot(bootstrap.frame1, aes(group = N.Iterations, x = N.Iterations, y = StandardDeviation)) +
    #     geom_boxplot() +
    #     facet_wrap(~Interval)
    #   
    #   ggsave(paste("/Users/maclark/Desktop/SSA_Forecast_Runs/Yields/Bootstrap/",
    #                iso3c[c],
    #                "_Crop_Bootstrap.pdf",
    #                sep = ""),
    #          width = 20,
    #          height = 10)
    # }
    # 
    
    # plot(raster.list.crop$prop_crop_2055_2060_mean)
    
    # Checking proportion of a cell in pasture and cropland in 2060
    # Only used for troubleshooting
    # range(test.list$`1`$LandUse$prop_pasture_2055_2060, na.rm = TRUE)
    # range(test.list$`1`$LandUse$prop_crop_2055_2060, na.rm = TRUE)
    # 
    # range(test.list$`1`$LandUse$prop_crop_2055_2060 + test.list$`1`$LandUse$prop_pasture_2055_2060, na.rm = TRUE)
    # 
    # Save it before I do something stupid and clear the workspace -----
    # save(test.list, 
    #      file = paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/Feb_16/BAU/RDAs/",
    #                   iso3c[c],
    #                   region.names[k],
    #                   "_16Feb.rda",
    #                   sep = ""))
    
    # Writing target for 2060 crop into a csv file
    # if(name == "Mike") {
    #   tmp_target_frame <- data.frame(ISO3 = iso3c[c],
    #                                  Target_Met_Crop_2015 = test.list[[1]]$messages$`2010_2015`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2015 = test.list[[1]]$messages$`2010_2015`$Messages.Crop[2],
    #                                  Target_Met_Crop_2020 = test.list[[1]]$messages$`2015_2020`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2020 = test.list[[1]]$messages$`2015_2020`$Messages.Crop[2],
    #                                  Target_Met_Crop_2025 = test.list[[1]]$messages$`2020_2025`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2025 = test.list[[1]]$messages$`2020_2025`$Messages.Crop[2],
    #                                  Target_Met_Crop_2030 = test.list[[1]]$messages$`2025_2030`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2030 = test.list[[1]]$messages$`2025_2030`$Messages.Crop[2],
    #                                  Target_Met_Crop_2035 = test.list[[1]]$messages$`2030_2035`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2035 = test.list[[1]]$messages$`2030_2035`$Messages.Crop[2],
    #                                  Target_Met_Crop_2040 = test.list[[1]]$messages$`2035_2040`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2040 = test.list[[1]]$messages$`2035_2040`$Messages.Crop[2],
    #                                  Target_Met_Crop_2045 = test.list[[1]]$messages$`2040_2045`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2045 = test.list[[1]]$messages$`2040_2045`$Messages.Crop[2],
    #                                  Target_Met_Crop_2050 = test.list[[1]]$messages$`2045_2050`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2050 = test.list[[1]]$messages$`2045_2050`$Messages.Crop[2],
    #                                  Target_Met_Crop_2055 = test.list[[1]]$messages$`2050_2055`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2055 = test.list[[1]]$messages$`2050_2055`$Messages.Crop[2],
    #                                  Target_Met_Crop_2060 = test.list[[1]]$messages$`2055_2060`$Messages.Crop[1],
    #                                  Target_Remaining_Crop_2060 = test.list[[1]]$messages$`2055_2060`$Messages.Crop[2],
    #                                  Target_Met_Pasture_2015 = test.list[[1]]$messages$`2010_2015`$Messages.Pasture,
    #                                  Target_Met_Pasture_2020 = test.list[[1]]$messages$`2015_2020`$Messages.Pasture,
    #                                  Target_Met_Pasture_2025 = test.list[[1]]$messages$`2020_2025`$Messages.Pasture,
    #                                  Target_Met_Pasture_2030 = test.list[[1]]$messages$`2025_2030`$Messages.Pasture,
    #                                  Target_Met_Pasture_2035 = test.list[[1]]$messages$`2030_2035`$Messages.Pasture,
    #                                  Target_Met_Pasture_2040 = test.list[[1]]$messages$`2035_2040`$Messages.Pasture,
    #                                  Target_Met_Pasture_2045 = test.list[[1]]$messages$`2040_2045`$Messages.Pasture,
    #                                  Target_Met_Pasture_2050 = test.list[[1]]$messages$`2045_2050`$Messages.Pasture,
    #                                  Target_Met_Pasture_2055 = test.list[[1]]$messages$`2050_2055`$Messages.Pasture,
    #                                  Target_Met_Pasture_2060 = test.list[[1]]$messages$`2055_2060`$Messages.Pasture)
    #   
    #   # Writing targets as a csv file
    #   write.csv(tmp_target_frame,
    #             paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/Feb_16/Trade/Target_CSVs/",
    #                   region.names[k], "_", iso3c[c], "_Target.csv", sep = ""),
    #             row.names = FALSE)
    #   # And removing target frame
    #   rm(tmp_target_frame)
    # }
    
    
    # And writing rasters
    if(name == "Mike"){
      for(i in 1:length(raster.list.crop)) {
        # Crops
        # png(filename = paste("/Users/maclark/Desktop/LandForecastRuns_29April/BAU/PNGs/",
        #                      iso3c[c],
        #                      region.names[k],
        #                      names(raster.list.crop)[i], 
        #                      ".png",
        #                      sep = ""),
        #     height = 20, width = 40, units = "cm", res = 72)
        # plot(raster.list.crop[[i]], zlim = c(-.1,1.1))
        # dev.off()
        # Writing rasters
        writeRaster(raster.list.crop[[i]],
                    paste("/Users/maclark/Desktop/LandForecastRuns_29April/BAU/TIFs/",
                          iso3c[c],
                          region.names[k],
                          names(raster.list.crop)[i], 
                          ".tif",
                          sep = ""))
        # Pasture
        # png(filename = paste("/Users/maclark/Desktop/LandForecastRuns_29April/BAU/PNGs/",
        #                      iso3c[c],
        #                      region.names[k],
        #                      names(raster.list.pasture)[i], 
        #                      ".png",
        #                      sep = ""),
        #     height = 20, width = 40, units = "cm", res = 72)
        # plot(raster.list.pasture[[i]], zlim = c(-.1,1.1))
        # dev.off()
        # Writing rasters
        writeRaster(raster.list.pasture[[i]],
                    paste("/Users/maclark/Desktop/LandForecastRuns_29April/BAU/TIFs/",
                          iso3c[c],
                          region.names[k],
                          names(raster.list.pasture)[i], 
                          ".tif",
                          sep = ""))
      }
    } else {
      for(i in 1:length(raster.list.crop)) {
        # Crops
        png(filename = paste("Outputs/ConversionMaps_Jan/",
                             iso3c[c],
                             region.names[k],
                             names(raster.list.crop)[i], 
                             "_14Dec_Crop.png",
                             sep = ""),
            height = 20, width = 40, units = "cm", res = 72)
        plot(raster.list.crop[[i]], zlim = c(-.1,1.1))
        dev.off()
        # Writing rasters
        # writeRaster(raster.list.crop[[i]],
        #             paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/Dec_13/Trade/TIFs/",
        #                   iso3c[c],
        #                   region.names[k],
        #                   names(raster.list.crop)[i], 
        #                   "_14Dec_Crop.tif",
        #                   sep = ""))
        # Pasture
        png(filename = paste("Outputs/ConversionMaps_Jan/",
                             iso3c[c],
                             region.names[k],
                             names(raster.list.pasture)[i], 
                             "_14Dec_Pasture.png",
                             sep = ""),
            height = 20, width = 40, units = "cm", res = 72)
        plot(raster.list.pasture[[i]], zlim = c(-.1,1.1))
        dev.off()
        # Writing rasters
        # writeRaster(raster.list.pasture[[i]],
        #             paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/Dec_13/Trade/TIFs/",
        #                   iso3c[c],
        #                   region.names[k],
        #                   names(raster.list)[i], 
        #                   "_14Dec_Pasture.tif",
        #                   sep = ""))
      }
    }
    
    # Removing items to clear space
    rm(raster.list.crop, raster.list.pasture, model_data_current, test.list, row.index.raster)
    
    # removing tmp raster files periodically
    
    rasterTmpFile(prefix = 'r_tmp_')
    # showTmpFiles()
    removeTmpFiles(1)
  }
}


Sys.time() - t1 # 22 mins per rep
# Alot more than 22 mins per rep =)

# Quitting r to reset memory
quit(save = "no", status = 0, runLast = TRUE)