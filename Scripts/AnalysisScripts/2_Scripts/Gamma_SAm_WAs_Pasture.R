#!/usr/bin/env Rscript

# Read me ----
# Modified version of 
# "GBM Model IUCN Region prob expansion outcome Two-step 24.August.R" plus
# "GBM 22.Sept Crop_MC.R" and "GBM 22.Sept Pasture_MC.R"

# I have changed some of the variable names to make them slightly more 
#     intelligible when it comes to the land-clearance loop
# Some of this comes from "ProbabilitiesOfLandClearance_Current.R"
#
# This version (16 October) fits the pasture GBM after the cropland one
#     AND fits a trimmed-down version of the GBM:
# - Does not include ISO Numeric but does include overall land demand
# - Fits a global model (not separate models to countries that have had cropland expansion, 
#   compared to contraction)

# Packages -----
library(plyr)
library(tidyverse)
library(reshape2)
library(caret)
library(pROC)
library(gbm)
library(e1071)
library(glmnet)
library(raster)
library(rgdal)
library(gridExtra)
library(doMC)
library(dplyr)
library(plyr)
library(nnet)
library(DHARMa)
library(SuppDists)
registerDoMC(cores = 4)

# Clean up -----
rm(list = ls())
# change cutoff for where we denote if a cell has experienced a change in crop extent
cutoff <- .025

# getting list of GBM file names
setwd("/Users/maclark/Desktop/Mike's Files/GBM Data Files")
file.list <- list.files(path = ".")

### Prep: creating lists used to assign model outputs
# creating list used to assign values to model output if cells experienced a change in cropland extent
region_list <- gsub("GBM.Data.22.Feb.","",file.list)
region_list <- gsub(".csv","",region_list)
region_list <- as.character(region_list)

region_list_model_change <- paste(gsub(" ","",region_list),"model_change",sep = "_") 
# list used to assign values to model output for model that estimates the change in cropland extent by cell
region_list_model_extentchange_increase <- paste(gsub(" ","",region_list),"model_extentchange_increase",sep = "_")
region_list_model_extentchange_decrease <- paste(gsub(" ","",region_list),"model_extentchange_decrease",sep = "_") 

# Data prep: Set predictor names ----
# I am selecting the columns by name not index to account for possible changes
# This is currently only used to test for scaled predictors
predictor_vars_pasture <- c("dist_50k",
                            "ag_suitability",
                            "PA_binary",
                            "prop_adj_to_predict_from_pasture", 
                            "prop_to_predict_from_pasture", 
                            "prop_to_predict_from_cell_delta_pasture",
                            "prop_to_predict_from_crop")
# And plot variables to test residuals
plot_vars_crop <- c('dist_50k_log' ,
                           'dist_50k_log_sq' ,
                           'ag_suitability' ,
                           'PA_binary',
                           'prop_to_predict_from_crop' ,
                           'prop_crop_sq' ,
                           'prop_to_predict_from_pasture' ,
                           'prop_pasture_sq' ,
                           'prop_to_predict_from_cell_delta_pasture' ,
                           'prop_cell_delta_pasture_sq',
                    "prop_adj_to_predict_from_pasture", 
                    "prop_adj_pasture_sq", 
                           'ISO_numeric')

# creating lists to savemodels
change_list_crop <- list() # first model step
proportion_list_crop_increase <- list() # second model step
proportion_list_crop_decrease <- list() # second model step

# Creating list of names used to save files
file.names.save <- paste(region_list,
                         "_7March_Cutoff_025_Crop.rda",
                         sep = "")
mods.prob <- list()
mods.increase <- list()
mods.decrease <- list()

coef.prob <- list()
coef.increase <- list()
coef.decrease <- list()

t1 = Sys.time()
### big loop to cycle through each file in file.list
for (k in c(10,12)) {
  print(region_list[k])
  
  # importing files
  dat_orig <- read.csv(file.list[k],
                       header = TRUE,
                       stringsAsFactors = FALSE)
  # All non PA areas have value of 0
  dat_orig$PA.Binary[-which(dat_orig$PA.Binary == 1)] <- 0
  
  # creating columns for for change in crop between 2002 and 2007 and 2007 and 2012
  dat_orig$pasture_extent_change_2007.2002 <- dat_orig$Prop.Pasture.2007 - dat_orig$Prop.Pasture.2002
  dat_orig$pasture_extent_change_2012.2007 <- dat_orig$Prop.Pasture.2012 - dat_orig$Prop.Pasture.2007
  # creating column for has_changed_crop between 2007 and 2012
  dat_orig$has_changed_pasture_2012.2007 <- NA
  dat_orig$has_changed_pasture_2012.2007[dat_orig$pasture_extent_change_2012.2007 >= cutoff] <- "Increase"
  dat_orig$has_changed_pasture_2012.2007[dat_orig$pasture_extent_change_2012.2007 <= -cutoff] <- "Decrease"
  dat_orig$has_changed_pasture_2012.2007[abs(dat_orig$pasture_extent_change_2012.2007) < cutoff] <- "aNoChange"
  # Data prep: Remove unnecessary columns, rename predictor_vars ----
  # I am renaming a) for consistency b) so that we can more easily understand
  dat_orig <- dplyr::select(dat_orig,
                            row_index = Row.Index,
                            ISO_numeric = ISO.Numeric, 
                            IUCN_region = IUCN.Region,
                            dist_50k = Dist.50k, 
                            ag_suitability = Ag.suitability, 
                            PA_binary = PA.Binary, 
                            change_demand_crop = CountryIncrease.Crop.2002.2007,
                            prop_adj_to_predict_from_pasture = Prop.AdjacentPasture.2007, 
                            prop_to_predict_from_crop = Prop.Crop.2007,
                            prop_to_predict_from_cell_delta_pasture = pasture_extent_change_2007.2002,
                            pasture_extent_change_2012.2007 = pasture_extent_change_2012.2007,
                            prop_to_predict_from_pasture = Prop.Pasture.2007,
                            prop_pasture_2002 = Prop.Pasture.2002,
                            prop_pasture_2012 = Prop.Pasture.2012,
                            prop_change_cell_pasture_2007.2012 = pasture_extent_change_2012.2007,
                            has_changed_pasture = has_changed_pasture_2012.2007,
                            GDD_10c = GDD_10c)
  # Data prep: Getting model data
  # Limiting to specific variables and complete cases
  dat_mod <- dat_orig[,c(predictor_vars_pasture, "has_changed_pasture", "prop_change_cell_pasture_2007.2012","ISO_numeric")]
  dat_mod <- dat_mod[complete.cases(dat_mod),]

  # Removing original dat
  rm(dat_orig)
  
  # Adding columns for square transformations
  dat_mod$prop_crop_sq <- dat_mod$prop_to_predict_from_crop^2
  dat_mod$prop_pasture_sq <- dat_mod$prop_to_predict_from_pasture^2
  dat_mod$prop_cell_delta_pasture_sq <- dat_mod$prop_to_predict_from_cell_delta_pasture^2
  dat_mod$prop_adj_pasture_sq <- dat_mod$prop_adj_to_predict_from_pasture^2
  # Adding columns for log transformations
  dat_mod$dist_50k_log <- log(dat_mod$dist_50k + .001)
  dat_mod$dist_50k_log_sq <- dat_mod$dist_50k_log^2
  # dat_mod$ISO_numeric <- as.character(dat_mod$ISO_numeric)
  
  # Getting multinomial model
  outcomeName <- "has_changed_pasture" # outcome factor variable
  # Attaching data set

  # Getting multinomial regression model
  # Getting multinomial regression model
  objModel <- multinom(has_changed_pasture ~ 
                         dist_50k_log +
                         dist_50k_log_sq +
                         ag_suitability +
                         PA_binary +
                         prop_to_predict_from_crop +
                         prop_crop_sq +
                         prop_to_predict_from_pasture +
                         prop_pasture_sq +
                         prop_to_predict_from_cell_delta_pasture +
                         prop_cell_delta_pasture_sq +
                         prop_adj_to_predict_from_pasture +
                         prop_adj_pasture_sq +
                         factor(ISO_numeric, ordered = FALSE),
                       dat = dat_mod,
                       model = FALSE,
                       maxit = 250)
  
  coef.multinom <- coef(objModel)
  
  # 
  # coef.prob[[region_list[k]]] <- as.data.frame(exp(coef(objModel)))
  # mods.prob[[region_list[k]]] <- objModel
  # 
  # # Changing outcome name for extent models
  # outcomeName <- "prop_change_cell_crop_2007.2012"
  # Processing: Creating data frames that only contains cells which have increased or decreasd in extent
  dat_region_increase <- dat_mod[dat_mod$prop_change_cell_pasture >= cutoff,]
  dat_region_decrease <- dat_mod[dat_mod$prop_change_cell_pasture <= -cutoff,]
  
  # removing data
  rm(dat_mod)
  
  # Dropping rows with values <0 for increaes and > 0 for decrease
  dat_region_increase <- dat_region_increase[dat_region_increase$prop_change_cell_pasture_2007.2012 >= 0,]
  dat_region_decrease <- dat_region_decrease[dat_region_decrease$prop_change_cell_pasture_2007.2012 <= 0,]
  
  # arcsining response variable
  # dat_region_increase$prop_change_cell_crop_2007.2012 <- asin(sqrt(dat_region_increase$prop_change_cell_crop_2007.2012))
  # dat_region_decrease$prop_change_cell_crop_2007.2012 <- asin(sqrt(abs(dat_region_decrease$prop_change_cell_crop_2007.2012)))
  
  
  
  # Getting model for extent increase
  objModel_increase <- glm(prop_change_cell_pasture_2007.2012 ~ 
                             dist_50k_log +
                             dist_50k_log_sq +
                             ag_suitability +
                             PA_binary +
                             prop_to_predict_from_crop +
                             prop_crop_sq +
                             prop_to_predict_from_pasture +
                             prop_pasture_sq +
                             prop_to_predict_from_cell_delta_pasture +
                             prop_cell_delta_pasture_sq +
                             prop_adj_to_predict_from_pasture +
                             prop_adj_pasture_sq +
                             factor(ISO_numeric, ordered = FALSE),
                            family = Gamma(link = 'log'),
                            dat = dat_region_increase,
                          model = FALSE,
                          control = list(maxit = 250))
  
  coef.increase <- as.data.frame(coef(objModel_increase))
  
  # 
  # 
  # # plot(objModel_increase,
  # #      main = paste(region_list[k],
  # #                   "Increase"))
  # # And getting model for extent decrease
  # 
  # # Getting absolute of decrease in crop extent
  # dat_region_decrease$prop_change_cell_pasture_2007.2012 <- abs(dat_region_decrease$prop_change_cell_pasture_2007.2012)
  # objModel_decrease <- glm(prop_change_cell_pasture_2007.2012 ~ 
  #                            dist_50k_log +
  #                            dist_50k_log_sq +
  #                            ag_suitability +
  #                            PA_binary +
  #                            prop_to_predict_from_crop +
  #                            prop_crop_sq +
  #                            prop_to_predict_from_pasture +
  #                            prop_pasture_sq +
  #                            prop_to_predict_from_cell_delta_pasture +
  #                            prop_cell_delta_pasture_sq +
  #                            prop_adj_to_predict_from_pasture +
  #                            prop_adj_pasture_sq +
  #                            factor(ISO_numeric, ordered = FALSE),
  #                          family = Gamma(link = 'log'),
  #                          dat = dat_region_decrease,
  #                          model = FALSE,
  #                          control = list(maxit = 250))
  # 
  # coef.decrease <- as.data.frame(coef(objModel_decrease))
  # plot(objModel_decrease,
  #      main = paste(region_list[k],
  #                   "Decrease"))
  
  # Simulating results and plotting
  # simincrease  <-  simulateResiduals(objModel_increase,  n=250)
  # 
  # # And plotting simulated residuals
  # graphics.off() # removing plots
  # par(mfrow = c(1,2)) # updating graphical parameters
  # # And plotting
  # png(file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Plots/",
  #                  region_list[k],
  #                  "Sim_Increase_Gamma.png",
  #                  sep = ""), 
  #     height = 10, width = 20, units = 'cm', res = 300)
  # plotSimulatedResiduals(simincrease)
  # dev.off()
  # 
  
  # # And plotting simulated residuals for each predictor variable
  # graphics.off() # removing plots
  # # And plotting individual coefficients
  # for(var.names in 1:length(plot_vars_crop)) {
  #   png(file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Plots/",
  #                    region_list[k],
  #                    plot_vars_crop[var.names],
  #                    "Sim_Increase_Gamma.png",
  #                    sep = ""), 
  #       height = 10, width = 20, units = 'cm', res = 300)
  #   plotResiduals(dat_region_increase[[plot_vars_crop[var.names]]], 
  #                 simincrease$scaledResiduals,
  #                 main = plot_vars_crop[var.names])
  #   dev.off() # Saving plots
  # }
  
  
  
  # removing from memory
  # rm(simincrease, dat_region_increase)
  # 
  # # For decreasing files
  # simdecrease  <-  simulateResiduals(objModel_decrease,  n=250)
  # # And removing plots from R stuido
  # graphics.off()
  # par(mfrow = c(1,2)) # updating graphical parameters
  # 
  # 
  # png(file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Plots/",
  #                  region_list[k],
  #                  "Sim_Decrease_Gamma.png",
  #                  sep = ""), 
  #     height = 10, width = 20, units = 'cm', res = 300)
  # plotSimulatedResiduals(simdecrease)
  # dev.off()
  # 
  # # And plotting simulated residuals for each predictor variable
  # graphics.off() # removing plots
  # # for(var.names in 1:length(plot_vars_crop)) {
  # #   png(file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Plots/",
  # #                    region_list[k],
  # #                    plot_vars_crop[var.names],
  # #                    "Sim_Decrease_Gamma.png",
  # #                    sep = ""), 
  # #       height = 10, width = 20, units = 'cm', res = 300)
  # #   plotResiduals(dat_region_decrease[[plot_vars_crop[var.names]]], 
  # #                 simdecrease$scaledResiduals,
  # #                 main = plot_vars_crop[var.names])
  # #   dev.off() # Saving plots
  # # }
  # # And removing plots to (hopefully) increase speed
  # graphics.off()
  
  
    # Removing from memory
  # rm(simdecrease, dat_region_decrease)
  
  save(objModel_increase, objModel,
       file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log/",
                    region_list[k],
                    "_Gamma_Link_Log_Multinom_Pasture_Model.rda",
                    sep = ""))
  
  rm(objModel_increase)
  # save(simincrease, simdecrease,
  #      file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/",
  #                         region_list[k],
  #                         "_Gamma_Link_Identity_Simulations.rda",
  #                   sep = ""))
  save(coef.increase, coef.multinom,
       file = paste("/Users/maclark/Desktop/Mike's Files/GBM Model Outputs/8_March/Gamma_Log/",
                    region_list[k],
                    "_Gamma_Link_Log_Multinom_Pasture_Coefs.rda",
                    sep = ""))

  
}

Sys.time() - t1

quit(save = 'no')