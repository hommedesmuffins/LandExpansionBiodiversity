# Read me -----
# Quick script to prepare data for use in the analyses
#
# Packages -----
library(plyr)
library(tidyverse)
library(reshape2)
library(ggplot2)

# Clean up ----
rm(list=ls())

# 1) Target setting ----
# 1) Load data ----
# 1) Target setting ----
# 1) Load data ----
targets_orig <- read.csv("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Land Targets/Land Forecasts All Years 9 March.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)

targets_orig <- unique(targets_orig)
# 1) Clean up data ----
# I have to manually add some missing harvested areas:
targets_mod <- targets_orig
targets_mod[targets_mod$ISO3 == "COD", "Area_Harvested_sq_km_2010"] <- 60200.02
targets_mod[targets_mod$ISO3 == "LBY", "Area_Harvested_sq_km_2010"] <- 8078.80

# 1.1) Manually calculatig forecasts for Libya and DRC
# Below chunk of code filters them out because they have NA values
targets_LBY_COD <- targets_mod[targets_mod$ISO3 %in% c("COD", "LBY"),]
# Looping through to update land target values
for(i in 1:8) {
  targets_LBY_COD[,5+i] <- targets_LBY_COD$Area_Harvested_sq_km_2010 * targets_LBY_COD[,13+i] * 100
}

# Dropping LBY and COD from original targets and then adding these back in from targets_LBY_COD
targets_orig <- targets_orig[!(targets_orig$ISO3 %in% c("COD","LBY")),]
targets_orig <- rbind(targets_orig, targets_LBY_COD)

targets_mod <- targets_orig %>%
  # Select only the proportional increases and original area
  dplyr::select(ISO3, 
         year = Year,
         area_km_2010 = Area_Harvested_sq_km_2010,
         Land_Baseline_ha:Land_GroupA.B.Trade_Half.Meat.to.Milk_ClosedYields) %>%
  # Melt data so scenarios are in separate rows
  gather(scenario_orig, area_harvested_ha, -ISO3, -year, -area_km_2010) %>%
  # Convert area harvested into km and year into column name
  mutate(area_harvested_km = area_harvested_ha / 100,
         year_mod = paste("area_km", year, sep = "_")) %>%
  # Rename scenarios
  mutate(scenario_tmp = gsub(paste(c("Land_", "_ha"), collapse = "|"),
                             "", 
                             scenario_orig),
         scenario = gsub("\\.", "_", scenario_tmp)) %>%
  dplyr::select(-scenario_orig, -scenario_tmp, -area_harvested_ha, -year) %>%
  # Cast so that years are in separate columns
  spread(year_mod, area_harvested_km) %>%
  as.data.frame()


# Calculate proportional increases required
# Divide each year by the preceding to get proportional increase
interval <- names(targets_mod)[grep("area_km", names(targets_mod))]

for(i in 2:length(interval)){
  n <- paste("prop_change",
             gsub("area_km_", "", interval[i-1]), 
             gsub("area_km_", "", interval[i]), 
             sep = "_")
  targets_mod[[n]] <- (targets_mod[,interval[i]] / targets_mod[,interval[i - 1]]) - 1
}

# 1) Clean up dataframe -----
targets_mod <- targets_mod %>%
  dplyr::select(ISO3, scenario, grep("prop_change", names(.))) %>%
  as.data.frame()

# 1) Save cleaned data ----
write.csv(targets_mod, 
          file = "/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/Land Clearing Targets/Targets_v4.csv",
          row.names = FALSE)

# 2) Urban expansion ----
# 1) For each country, get the number of cells with non-zero probability 
#    of being urban in 2030.
# 2) Adjust this for cells converted to urban between 2000 and 2010 
# 3) Multiply the number of these cells by the mean probability of conversion in a country 
#        to estimate the number of cells to convert
# 4) Divide this number by the number of five-year periods (4: 2010-2030) to get the number converted each decade
# 5) Pick that number based on the probability of conversion
# 6) All that cell is converted
# 7) Do this within each of the 1000 runs

# 2) Clean up, load data and functions ----
rm(list=ls())
library(raster)
library(rgdal)
# source("Scripts/AnalysisScripts/0_Functions_CountryFinescale.R")
source("/Users/maclark/Desktop/Mike's Files/SourceTree/AgExpansion/Scripts/AnalysisScripts/0_Functions_CountryFinescale.R")
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
goode <- "+proj=igh +ellps=WGS84 +units=m +no_defs"

# Load urban expansion scenarios
# urban_exp <- raster("/Users/Me/Documents/Datasets/UrbanAreas/SetoGuneralpHutyra_UF/seto_uf/w001001.adf",
#                     crs = goode) # Put the CRS in here, so it's ready

urban_exp <- raster("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Datasets/SetoGuneralpHutyra_UF 2/seto_uf/w001001.adf",
                    crs = goode) # Put the CRS in here, so it's ready
# plot(urban_exp)
# Urban data in 2010
# urban2010_raster <- raster("/Users/Me/Documents/Datasets/MODIS/Urban2010.Mollweide.1.5km.tif", 
#                            crs = mollweide)
urban2010_raster <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Urban2010.Mollweide.1.5km.tif",
                           crs = mollweide)

# Country ID
# countryID_raster <- raster("/Users/Me/Documents/Datasets/WorldMaps/MollweideCountryID_1.5km.tif",
#                            crs = mollweide)
countryID_raster <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/MollweideCountryID_1.5km.tif",
                           crs = mollweide)
crs(countryID_raster) <- mollweide
# 2.1) Get data to same resolution and projection ----
# Not sure which is the better way of doing this: either converting the country and urban2010 data
#     to Goode Homolosine or the urban projections to Mollweide. I think the latter may be the 
#     way to go...
urban_exp <- projectRaster(urban_exp,
                           crs = mollweide,
                           res = c(5000,5000),
                           over = TRUE)

# To get the rasters into the same resolution we want to resample the coarse 
#     (urban expansion) data to the same scale that we're going to use - this differs
#     from before, where I aggregated the finer-scale rasters up to the coarser resolution.
# I've kept that code at the end, in case it's useful later.
template <- raster(extent(urban_exp), # The projections map has a smaller extent than others
                   crs = crs(urban2010_raster), 
                   resolution = res(urban2010_raster))

t1 <- Sys.time()
urban_exp <- resample(x = urban_exp, 
                      y = template,
                      method = "bilinear")
Sys.time() - t1 # 4.5 mins
plot(urban_exp)

# 2.2) Extend the urban projections to the extent of the other data -----
t1 <- Sys.time()
urban_exp <- extend(x = urban_exp, 
                    y = extent(urban2010_raster), 
                    value = NA)
Sys.time() - t1 
plot(urban_exp)

# 2.3) Put data into dataframes -----
urban_data <- data.frame(ISO3.Numeric = as.vector(countryID_raster),
                         urban_prop_2010 = as.vector(urban2010_raster),
                         urban_prob = as.vector(urban_exp))

# CHECK THAT THIS MAKES SENSE: AM I NUMBERING THE ROWS IN THE RIGHT WAY?
urban_data$Row.Index <- 1:nrow(urban_data)

# Cap probability of conversion at 100 (i.e. remove 101 - current)
urban_data[urban_data$urban_prob > 100 &
             !is.na(urban_data$urban_prob),"urban_prob"] <- 0

# plot(raster(matrix(urban_data$urban_prob,
#                    nrow = nrow(countryID_raster),
#                    ncol = ncol(countryID_raster),
#                    byrow = TRUE)))
# 2) Clean up -----
rm(countryID_raster, t1, urban_exp)

# 2.4) Count cells to convert and multiply by mean probability ----
# For each country, get the number of cells with non-zero probability 
#    of being urban in 2030.
# Multiply the number of these cells by the mean probability of conversion in a country 
#        to estimate the number of cells to convert

# THIS SECTION IS JUST SUBSETTING SO THAT I CAN RUN THINGS ON MY PUNY COMPUTER.
# COMMENT THEM OUT OR DELETE WHEN RUNNING IT ON THE WHOLE THING
# urban_data_orig <- urban_data
# urban_data <- urban_data_orig
library(data.table) 
# NB. Datatable is FAAAAAAST. But I don't understand its syntax at all, so suggest we keep
#     ploughing on as we were.
# Here I am using it to subset, as working on a dataframe was just killing my computer
#     every time.
 
# t1 <- Sys.time()
# urban_data <- as.data.table(urban_data)
# urban_data <- urban_data[ISO3.Numeric == 8 | ISO3.Numeric == 807]
# urban_data <- as.data.frame(urban_data)
# Sys.time() - t1

urban_target <- urban_data %>%
  # Only want countries
  filter(!is.na(ISO3.Numeric)) %>%
  group_by(ISO3.Numeric) %>%
  # Count cells with probability of urban conversion and get mean probability of conversion
  filter(urban_prob != 0) %>%
  summarise(count = n(),
            mean_prob = mean(urban_prob / 100, na.rm = TRUE)) %>%
  # Multiply number of cells by probability to get number of cells to convert
  mutate(no_to_convert = count * mean_prob,
         # Divide this number by the number of time steps (4: 2010-2030) 
         #     to get the number converted each decade
         per_time_conv_exact = no_to_convert / 4,
         # Get the number to definitely convert, and the number to possibly convert
         pt_definite = floor(per_time_conv_exact),
         pt_possible = per_time_conv_exact - pt_definite) %>%
  as.data.frame()

# 2) Create the scenarios -----
time_step <- c("2010_2015", "2015_2020", "2020_2025", "2025_2030")

t1 <- Sys.time()
urban_scenarios <- data.frame(run = 1:1000) %>%
  group_by(run) %>%
  do(urb_conv_fun(.)) 
Sys.time() - t1
# 13s for Albania and Macedonia

# 2) Save scenarios and data-----
save(urban_scenarios, urban_data,
     # WHATEVER YOUR FILE IS!
     file = "/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/ModData/UrbanScenarios_HiRes_working.rda")

