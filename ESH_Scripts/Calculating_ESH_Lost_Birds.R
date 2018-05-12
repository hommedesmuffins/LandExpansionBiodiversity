# Script calculates percent of AOO lost and changes in RLI classisfication for bird species ----
# Does so by...
# 1) Calculating the extent of crop and pasture in a species original ESH
# 1.b) This is done in the scripts "Rasters_'scenario'.R", where scenario = "yields", "diets", "trade", "all", or "bau"
# 2) Merge calculations from (1) with species habitat preferences from DW
# 3) Calculate AOO remaining by
# 3a) Not penalizing area in crop or pasture (e.g. count all of crop/pasture towards AOO) if these are "suitable" for a species
# 3b) 50% penalty (e.g. count half of crop/pasture towards AOO) on area in crop or pasture if these are "marginal" for a species
# 3c) 100% penalty (e.g. count none of crop/pasture towards AOO) on area in crop or pasture if neither 3a or 3b is true
# 4) calculate percent of remaining habitat lost from 2010 to 2060
# 5) Update estimates of IUCN Red List classification based onc riteria B2 where species are...
# "vulnerable" if AOO < 2000kmsq
# "endangered" if AOO < 500kmsq
# "critically endangered" if AOO < 10kmsq 


# Loading packages
library(raster)
library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# Step 0 ----
# Getting cell area
tmp.raster <- raster("/Users/macuser/Desktop/Ag Expansion and Biodiversity/Crop2010.Mollweide.1.5km.tif")
area <- (res(tmp.raster)[1] * res(tmp.raster)[2])/1e6

# Step 1 ----
# Importing results. These correspond with step (1) at the top of the script 

# For Cropland
# Getting list of files
master.csv <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/ESH_Loss_Crop.csv",
                    stringsAsFactors = FALSE)
# 2010 data is wrong
master.csv <- master.csv[master.csv$year != 2010,]

# Getting 2010 data for baseline AOO left
# Creating data frame with only 2010 data
files.list <- list.files("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Crop_2010",
                         full.names = TRUE)
for(i in 1:length(files.list)) {
  tmp.csv <- read.csv(files.list[i],
                      stringsAsFactors = FALSE)
  if(i == 1) {
    dat_2010 <- tmp.csv
  } else {
    dat_2010 <- rbind(dat_2010, tmp.csv)
  }
}

# and replicating 2010 crop data for other scenarios
scenarios.list <- sort(unique(master.csv$Scenario))

for(i in 1:length(scenarios.list)) {
  tmp.csv <- dat_2010
  tmp.csv$Scenario <- scenarios.list[i]
  if(i == 1) {
    scenarios_2010 <- tmp.csv
  } else {
    scenarios_2010 <- rbind(scenarios_2010, 
                            tmp.csv)
  }
}
# changing name to avoid confusion
names(scenarios_2010)[names(scenarios_2010) == 'esh_loss'] <- 'esh_loss_crop_2010'

scenarios_2010 <- scenarios_2010[!is.na(scenarios_2010$species),]
# joining data sets
master_crop <- left_join(master.csv,
                         dplyr::select(scenarios_2010, -year))

# reorganizing data sets
master_crop <- dplyr::select(master_crop,
                             species,
                             subspecies,
                             Scenario,
                             year,
                             esh_area_orig,
                             esh_loss_crop_2010,
                             esh_loss_crop)

scenarios_2010 <- dplyr::select(scenarios_2010,
                             species,
                             subspecies,
                             Scenario,
                             year,
                             esh_area_orig,
                             esh_loss_crop_2010)
# Adding column to rbind
scenarios_2010$esh_loss_crop <- scenarios_2010$esh_loss_crop_2010


# And combining 2010 with other years of data
master_crop <- rbind(master_crop,
                     scenarios_2010)

# For Pastureland
# Getting list of files
master.csv <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/ESH_Loss_Pasture.csv",
                    stringsAsFactors = FALSE)

# Getting 2010 data for baseline AOO left
# Creating data frame with only 2010 data
dat_2010 <- unique(master.csv[master.csv$year == 2010,])
dat_2010$esh_loss_2010 <- dat_2010$esh_loss
# Merging 2010 data back in
master_pasture <- left_join(master.csv,
                            unique(dplyr::select(dat_2010, species, esh_area_orig, esh_loss_2010)))

# Merging in crop and pasture data
# First renaming columns
names(master_crop)[names(master_crop) == "esh_loss_crop"] <- "esh_loss_crop"
names(master_crop)[names(master_crop) == "esh_loss_2010"] <- "esh_loss_crop_2010"

names(master_pasture)[names(master_pasture) == "esh_loss_pasture"] <- "esh_loss_pasture"
names(master_pasture)[names(master_pasture) == "esh_loss_2010"] <- "esh_loss_pasture_2010"

# And now merging
master_frame <- full_join(master_crop,
                          master_pasture)

# And dropping rows without species identifiers
master_frame <- master_frame[!is.na(master_frame$species),]
# And reorganizing because this will drive me nuts
master_frame <- dplyr::select(master_frame,
                              species,
                              subspecies,
                              Scenario,
                              year,
                              esh_area_orig,
                              esh_loss_crop_2010,
                              esh_loss_pasture_2010,
                              esh_loss_crop,
                              esh_loss_pasture)

# converting NAs to 0s
master_frame$esh_loss_crop_2010[is.na(master_frame$esh_loss_crop_2010)] <- 0
master_frame$esh_loss_pasture_2010[is.na(master_frame$esh_loss_pasture_2010)] <- 0
master_frame$esh_loss_crop[is.na(master_frame$esh_loss_crop)] <- 0
master_frame$esh_loss_pasture[is.na(master_frame$esh_loss_pasture)] <- 0

# Step 2----
# Adding in habitat preferences
# Thanks DW!
bird_preferences <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Habitat_Preferences/Habitats2017.csv",
                             stringsAsFactors = FALSE)

# Getting differences for habitat specialist and generalist
habitat_specialist <- bird_preferences

# Creating counter
# setting to 1 when a habitat is "suitable"
# And setting to .5 when a habitat is "marginal"
habitat_specialist$Suitable[habitat_specialist$Suitability == "Suitable"] <- 1
habitat_specialist$Marginal[habitat_specialist$Suitability == "Marginal"] <- 1
# Changing names of species
names(habitat_specialist)[names(habitat_specialist) == "Scientific.name"] <- "species"
# Creating data frame for habitat preferences
habitat_preferences <- dplyr::select(habitat_specialist,
                                     species,
                                     Habitats.Classification.Scheme.Level.1,
                                     Habitats.Classification.Scheme.Level.2)

# Getting number of rows in data frame
# This corresponds to number of habitats a species can live in
# "Marginal" habitats are penalized
habitat_specialist <- ddply(habitat_specialist,.
                            (species),
                            summarise,
                            Suitable = sum(Suitable, na.rm = TRUE),
                            Marginal = sum(Marginal, na.rm = TRUE))
# And updating format of "species" to mesh with ESH raster data
habitat_specialist$species <- gsub(" ",
                                   "_",
                                   habitat_specialist$species)
# Identifying habitat specialists/generalists
habitat_specialist$Specialist_Genearlist[habitat_specialist$Suitable == 1 & habitat_specialist$Marginal == 0] <- "Specialist"
habitat_specialist$Specialist_Genearlist[habitat_specialist$Suitable == 1 & habitat_specialist$Marginal > 0] <- "Quasi-Specialist"
habitat_specialist$Specialist_Genearlist[habitat_specialist$Suitable > 1] <- "Generalist"
# Getting sum habitats
habitat_specialist$sumHabitats <- habitat_specialist$Suitable + habitat_specialist$Marginal * 0.5


# Getting wide data frame of habitat preferences
habitat_preferences2 <- reshape(habitat_preferences,
                               timevar = "Habitats.Classification.Scheme.Level.1",
                               idvar = "species",
                               direction = "wide")


# Getting habitat preferences
# Value of 1 == species can't survive there (e.g. 0% of land counts towards AOO)
# Value of .5 == species has 50% penalty there (e.g. 50% of land counts towards AOO; corresponds with suitability == "marginal")
# Value of 0 == species has no penalty there (e.g. 100% of land counts towards AOO; corresponds with suitability == "suitable")
# First for cropland
bird_preferences$Cropland <- 1
bird_preferences$Cropland[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Arable Land" &
                            bird_preferences$Suitability == "Marginal"] <- 0.5
bird_preferences$Cropland[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Arable Land" &
                            bird_preferences$Suitability == "Suitable"] <- 0
# Next for pastureland
bird_preferences$Pastureland <- 1
bird_preferences$Pastureland[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Pastureland" &
                               bird_preferences$Suitability == "Marginal"] <- 0.5
bird_preferences$Pastureland[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Pastureland" &
                               bird_preferences$Suitability == "Suitable"] <- 0
# Next for grassland
# Getting row indicator for grassland species

grassland.index_subtropical <- which(gsub("Grassland - Subtropical","",bird_preferences$Habitats.Classification.Scheme.Level.2) != bird_preferences$Habitats.Classification.Scheme.Level.2)
grassland.index_temperate <- which(gsub("Grassland - Temperate","",bird_preferences$Habitats.Classification.Scheme.Level.2) != bird_preferences$Habitats.Classification.Scheme.Level.2)
# And combining
grassland.index <- c(grassland.index_subtropical,
                     grassland.index_temperate)
# Getting uniques only
# There hsould be no duplicates as is
grassland.index <- unique(grassland.index)
# Getting index for "marginal" suitability
grassland.index.marginal <- grassland.index[(bird_preferences$Suitability[grassland.index] == "Marginal")]
# And for suitable
grassland.index.suitable <- grassland.index[(bird_preferences$Suitability[grassland.index] == "Suitable")]

bird_preferences$Grassland <- 1
bird_preferences$Grassland[grassland.index.marginal] <- 0.5
bird_preferences$Grassland[grassland.index.suitable] <- 0

# Getting rid of spaces to merge habitat preferences with ESH data frame
bird_preferences$Scientific.name <- gsub(" ",
                                         "_",
                                         bird_preferences$Scientific.name)

# Adding new column for merging
bird_preferences$species <- bird_preferences$Scientific.name

# getting current IUCN classification
current_iucn <- unique(dplyr::select(bird_preferences,species, X2017.IUCN.Red.List.Category))
# changing names for IUCN classification because its weird
names(current_iucn)[2] <- "IUCN_Current"

# And limiting bird_preferences to only rows with pasture or crop identifiers
# We only penalize extent of area in cropland and pastureland
# Don't need other habitat preferences because of this
bird_preferences_crop <- bird_preferences[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Arable Land",]
bird_preferences_pasture <- bird_preferences[bird_preferences$Habitats.Classification.Scheme.Level.2 == "Artificial/Terrestrial - Pastureland",]
bird_preferences_grassland <- bird_preferences[grassland.index,]

# Getting only cropland and pastureland columns
bird_preferences_crop <- dplyr::select(bird_preferences_crop, species, Cropland)
bird_preferences_pasture <- dplyr::select(bird_preferences_pasture, species, Pastureland)
bird_preferences_grassland <- dplyr::select(bird_preferences_grassland, species, Grassland)

# Merging cropland and pastureland habitat preferences
bird_preferences <- full_join(bird_preferences_crop,
                              bird_preferences_pasture)
bird_preferences <- full_join(bird_preferences,
                              bird_preferences_grassland)

# NA values indicate species don't exist in that habitat
# Changing to 1 to reflecct this
bird_preferences$Cropland[is.na(bird_preferences$Cropland)] <- 1
bird_preferences$Pastureland[is.na(bird_preferences$Pastureland)] <- 1
bird_preferences$Grassland[is.na(bird_preferences$Grassland)] <- 1

# Updating species column in master csv data table 
# Doing this to be able to merge file
# The ESH rasters contain subspecies, whereas the habitat preferences file does not
master_frame$species_subspecies <- master_frame$species
master_frame$species <- gsub("_P1_O1_S1_HabElev.tif","", master_frame$species_subspecies)
master_frame$species <- gsub("_P2_O1_S1_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P3_O1_S1_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P4_O1_S1_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P5_O1_S1_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P6_O1_S1_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P1_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P2_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P3_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P4_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P5_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_P6_O1_S3_HabElev.tif","", master_frame$species)
master_frame$species <- gsub("_S1","", master_frame$species)
master_frame$species <- gsub("_S2","", master_frame$species)
master_frame$species <- gsub("_S3","", master_frame$species)

# Adding in current IUCN status
# Merging ESH data frame with habitat preferences
master_frame <- unique(left_join(master_frame,
                                   dplyr::select(bird_preferences, species, Cropland, Pastureland, Grassland)))

# And merging in current IUCN status
master_frame <- unique(left_join(master_frame,
                                   current_iucn))

# And getting in unique rows to avoid duplication
master_frame <- unique(master_frame)

# And updating cropland and pastureland NAs
# NAs correspond to species that aren't identified as surviving in these areas
# Changing NAs to 1
# 1 Indicates species can't exist in habitat
master_frame$Cropland[is.na(master_frame$Cropland)] <- 1
master_frame$Pastureland[is.na(master_frame$Pastureland)] <- 1
master_frame$Grassland[is.na(master_frame$Grassland)] <- 1

# Adding in organism size
species.size <- read.csv("/Users/macuser/Downloads/allspecies.csv",
                         stringsAsFactors = FALSE)
species.size <- dplyr::select(species.size, Scientific_Name, Combined_Genus_Species_Mass_kg)
names(species.size)[names(species.size) == "Scientific_Name"] <- "species"
names(species.size)[names(species.size) == "Combined_Genus_Species_Mass_kg"] <- "species_mass"
# Updating format of species name to mesh with master_plot
species.size$species <- gsub(" ",
                             "_",
                             species.size$species)

master_frame <- unique(left_join(master_frame,
                                species.size))

# Looking at differences between species mass
master_frame$mass_class[master_frame$species_mass >= 2] <- "large"
master_frame$mass_class[master_frame$species_mass < 0.5] <- "small"
master_frame$mass_class[master_frame$species_mass < 2 & master_frame$species_mass >= 0.5] <- "medium"

# # writing csv
# write.csv(master_frame,
#           "/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Bird_ESH_AOO_3April.csv",
#           row.names = FALSE)
# 
# # Importing data frame
# master_frame <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Bird_ESH_AOO_3April.csv",
#                          stringsAsFactors = FALSE)

# Updating habitat preference for pastureland
# Assming grassland and pastureland is interchangeable habitat
# Essentially taking habitat suitability of grassland or pastureland for whichever habitata species is more tolerant of
master_frame$Pastureland[master_frame$Grassland < master_frame$Pastureland] <- master_frame$Grassland[master_frame$Grassland < master_frame$Pastureland]


# Step 3----
# Calculating (a) percent of remaining habitat lost and (b) AOO remaining

# First getting AOO remaining in 2010
# Calculated as ESH area original - (area lost crop + area lost pasture)
# Note that area los crop and area lost pasture are weighted by species specific habitat tolerance (preferences might be a better word?)
# E.g. species where cropland is "suitable" have no penalty for having cropland in their esh
# "marginal" is a 50% penalty
# And anything else is a 100% penalty
master_frame$AOO_2010 <- master_frame$esh_area_orig - (master_frame$esh_loss_crop_2010 * master_frame$Cropland + master_frame$esh_loss_pasture_2010 * master_frame$Pastureland)

# Second getting AOO remaining in each of the time steps
# Calculated as above, but for each year
# Calculated as ESH area original - area in crop * crop suitability penatly
master_frame$AOO_remaining <- master_frame$esh_area_orig - (master_frame$esh_loss_crop * master_frame$Cropland + master_frame$esh_loss_pasture * master_frame$Pastureland)

# Step 4----
# Calculating percent of AOO lost in each time step
# This is relative to 2010
master_frame$Percent_AOO_lost <- (1 - (master_frame$AOO_remaining / master_frame$AOO_2010)) * 100


# # Writing csv
# write.csv(master_frame,
#           "/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Bird_ESH_AOO_7April.csv",
#           row.names = FALSE)
# 
# # Importing data frame
# master_frame <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Bird_ESH_AOO_7April.csv",
#                          stringsAsFactors = FALSE)
# # Step 5----
# # Getting species that increase in IUCN status
# # Getting our estimates of RLI in 2010
# master_summary$Predicted_IUCN_2010[master_summary$AOO_2010 * area > 3000] <- "LC" 
# master_summary$Predicted_IUCN_2010[master_summary$AOO_2010 * area <= 3000] <- "NT"
# master_summary$Predicted_IUCN_2010[master_summary$AOO_2010 * area <= 2000] <- "VU"
# master_summary$Predicted_IUCN_2010[master_summary$AOO_2010 * area <= 500] <- "EN"
# master_summary$Predicted_IUCN_2010[master_summary$AOO_2010 * area <= 10] <- "CR"
# 
# # And getting estimates in future
# master_summary$Predicted_IUCN[master_summary$AOO_remaining * area > 3000] <- "LC"
# master_summary$Predicted_IUCN[master_summary$AOO_remaining * area <= 3000] <- "NT"
# master_summary$Predicted_IUCN[master_summary$AOO_remaining * area < 2000] <- "VU"
# master_summary$Predicted_IUCN[master_summary$AOO_remaining * area < 500] <- "EN"
# master_summary$Predicted_IUCN[master_summary$AOO_remaining * area < 10] <- "CR"
# 
# # And converting our predictions to a linear scale
# # First for 2010 data
# master_summary$Predicted_IUCN_2010_Linear[master_summary$Predicted_IUCN_2010 == "LC"] <- 1
# master_summary$Predicted_IUCN_2010_Linear[master_summary$Predicted_IUCN_2010 == "NT"] <- 2
# master_summary$Predicted_IUCN_2010_Linear[master_summary$Predicted_IUCN_2010 == "VU"] <- 3
# master_summary$Predicted_IUCN_2010_Linear[master_summary$Predicted_IUCN_2010 == "EN"] <- 4
# master_summary$Predicted_IUCN_2010_Linear[master_summary$Predicted_IUCN_2010 == "CR"] <- 5
# 
# # Second for projected data
# master_summary$Predicted_IUCN_Linear[master_summary$Predicted_IUCN == "LC"] <- 1
# master_summary$Predicted_IUCN_Linear[master_summary$Predicted_IUCN == "NT"] <- 2
# master_summary$Predicted_IUCN_Linear[master_summary$Predicted_IUCN == "VU"] <- 3
# master_summary$Predicted_IUCN_Linear[master_summary$Predicted_IUCN == "EN"] <- 4
# master_summary$Predicted_IUCN_Linear[master_summary$Predicted_IUCN == "CR"] <- 5
# 
# # And calculating change in predicted IUCN status
# master_summary$Linear_Change_IUCN_Status <- master_summary$Predicted_IUCN_Linear - master_summary$Predicted_IUCN_2010_Linear
# 
# # And creating columns for species that increased, decreased, or didn't change in our IUCN predictions
# master_summary$Increase_IUCN_Status[master_summary$Linear_Change_IUCN_Status > 0] <- 1
# master_summary$NoChange_IUCN_Status[master_summary$Linear_Change_IUCN_Status == 0] <- 1
# master_summary$Decrease_IUCN_Status[master_summary$Linear_Change_IUCN_Status < 0] <- 1
# 

# Adding in data for generalist vs specialist species
master_summary <- unique(left_join(master_frame,
                            habitat_specialist))

# Getting a data frame for plotting ----
# Limiting to species that aren't completely wonky
master_plot <- master_summary[master_summary$Percent_AOO_lost <= 100,]
# Limiting to LC, NT, VU, EN, or CR species
master_plot <- master_plot[master_plot$IUCN_Current == "LC" |
                             master_plot$IUCN_Current == "NT" |
                             master_plot$IUCN_Current == "VU" |
                             master_plot$IUCN_Current == "EN" |
                             master_plot$IUCN_Current == "CR" |
                             master_plot$IUCN_Current == "DD",]

# Dropping NAs
master_plot <- master_plot[!is.na(master_plot$IUCN_Current),]
master_plot <- master_plot[!is.na(master_plot$species),]
master_plot <- master_plot[!is.na(master_plot$Scenario),]
master_plot <- master_plot[!is.na(master_plot$year),]
master_plot <- master_plot[!is.na(master_plot$esh_area_orig),]
master_plot <- master_plot[!is.na(master_plot$Percent_AOO_lost),]

# setting levels for IUCN status
master_plot <- transform(master_plot, IUCN_Current = factor(IUCN_Current,
                                                            levels = rev(c("DD","LC","NT","VU","EN","CR"))))



# Getting percent avoided by scenario in 2060
master_plot_2060 <- master_plot[master_plot$year == 2060,]
# individual dfs for different scenarios
master_plot_bau <- master_plot_2060[master_plot_2060$Scenario == "BAU",c(1,18)]
master_plot_all <- master_plot_2060[master_plot_2060$Scenario == "All",c(1,18)]
master_plot_yields <- master_plot_2060[master_plot_2060$Scenario == "Yields",c(1,18)]
master_plot_diets <- master_plot_2060[master_plot_2060$Scenario == "Diets",c(1,18)]
master_plot_trade <- master_plot_2060[master_plot_2060$Scenario == "Trade",c(1,18)]

names(master_plot_bau)[2] <- "BAU"
names(master_plot_all)[2] <- "ALL"
names(master_plot_yields)[2] <- "YIELDS"
names(master_plot_diets)[2] <- "DIETS"
names(master_plot_trade)[2] <- "TRADE"

master_plot_2060 <- left_join(master_plot_bau,
                              master_plot_all)
master_plot_2060 <- left_join(master_plot_2060,
                              master_plot_yields)
master_plot_2060 <- left_join(master_plot_2060,
                              master_plot_diets)
master_plot_2060 <- left_join(master_plot_2060,
                              master_plot_trade)

# Getting percent avoided by species
master_plot_2060$Percent_yields <- master_plot_2060$YIELDS/master_plot_2060$BAU
master_plot_2060$Percent_diets <- master_plot_2060$DIETS/master_plot_2060$BAU
master_plot_2060$Percent_trade <- master_plot_2060$TRADE/master_plot_2060$BAU
master_plot_2060$Percent_all <- master_plot_2060$ALL/master_plot_2060$BAU

master_plot_2060 <- master_plot_2060[!(master_plot_2060$BAU == 0),]


master_plot_2060$Count <- "Count"
# Getting summary
master_plot_2060 <- ddply(master_plot_2060,.
                          (Count),
                          summarise,
                          Percent_yields = mean(Percent_yields, na.rm = TRUE),
                          Percent_diets = mean(Percent_diets, na.rm = TRUE),
                          Percent_trade = mean(Percent_trade, na.rm = TRUE),
                          Percent_all = mean(Percent_all, na.rm = TRUE))



write.csv(master_plot,
          "/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Plot.Frame_10April.csv",
          row.names = FALSE)

master_plot <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Plot.Frame_10April.csv",
                        stringsAsFactors = FALSE)


# View(master_plot)

# Getting number of species that lose more than 10, 25, 50, and 75% of their remaining habitat
master_plot$Count = 1
master_plot$Percent_10[master_plot$Percent_AOO_lost >= 10] = 1
master_plot$Percent_25[master_plot$Percent_AOO_lost >= 25] = 1
master_plot$Percent_50[master_plot$Percent_AOO_lost >= 50] = 1
master_plot$Percent_75[master_plot$Percent_AOO_lost >= 75] = 1
master_plot$Percent_90[master_plot$Percent_AOO_lost >= 90] = 1

# summarizing
esh_loss_percent <- ddply(master_plot,.
                           (year, Scenario),
                           summarise,
                           Count = sum(Count, na.rm = TRUE),
                           Percent_10 = sum(Percent_10, na.rm = TRUE),
                           Percent_25 = sum(Percent_25, na.rm = TRUE),
                           Percent_50 = sum(Percent_50, na.rm = TRUE),
                           Percent_75 = sum(Percent_75, na.rm = TRUE),
                           Percent_90 = sum(Percent_90, na.rm = TRUE))

# Step 6 ----
# Analyses

# Fig 3 a----
# Looking at losses in ESH based on current IUCN status
iucn_status_aoo_lost <- ddply(master_plot,.
                              (Scenario, year, IUCN_Current),
                              summarise,
                              Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE),
                              Count = sum(Count, na.rm = TRUE))
iucn_status_aoo_lost$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(iucn_status_aoo_lost)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == iucn_status_aoo_lost$Scenario[i] & master_plot$year == iucn_status_aoo_lost$year[i] & master_plot$IUCN_Current == iucn_status_aoo_lost$IUCN_Current[i]], na.rm = TRUE)
  # Updating sd
  iucn_status_aoo_lost$sdPercent_AOO_lost[i] <- tmp.sd
}


iucn_status_aoo_lost <- iucn_status_aoo_lost[iucn_status_aoo_lost$IUCN_Current == "LC" |
                                               iucn_status_aoo_lost$IUCN_Current == "NT" |
                                               iucn_status_aoo_lost$IUCN_Current == "VU" |
                                               iucn_status_aoo_lost$IUCN_Current == "EN" |
                                               iucn_status_aoo_lost$IUCN_Current == "CR" |
                                               iucn_status_aoo_lost$IUCN_Current == 'DD',]
iucn_status_aoo_lost <- iucn_status_aoo_lost[!is.na(iucn_status_aoo_lost$IUCN_Current),]

# Transforming to get factors in the correct order
iucn_status_aoo_lost <- transform(iucn_status_aoo_lost,
                                  IUCN_Current = factor(IUCN_Current,
                                                        levels = c("DD","LC","NT","VU","EN","CR")))
# Getting upper and lower estimates
iucn_status_aoo_lost$Lower = iucn_status_aoo_lost$Percent_AOO_lost - iucn_status_aoo_lost$sdPercent_AOO_lost/sqrt(iucn_status_aoo_lost$Count)
iucn_status_aoo_lost$Upper = iucn_status_aoo_lost$Percent_AOO_lost + iucn_status_aoo_lost$sdPercent_AOO_lost/sqrt(iucn_status_aoo_lost$Count)

iucn_status_aoo_lost$Percent_AOO_lost[iucn_status_aoo_lost$year == 2010] <- 0
iucn_status_aoo_lost$sdPercent_AOO_lost[iucn_status_aoo_lost$year == 2010] <- 0

# Line plot
Fig3a.2 <- ggplot(iucn_status_aoo_lost[iucn_status_aoo_lost$Scenario == 'BAU',], aes(x = year, y = Percent_AOO_lost, colour = factor(IUCN_Current))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, color = factor(IUCN_Current)), linetype = 2) + 
  geom_line(aes(x = year, y = Lower, color = factor(IUCN_Current)), linetype = 2) +
  theme_classic() +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2010, y = 46, label = "(a)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  theme(legend.justification=c(0,0), legend.position=c(.27,.83), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=3)) +
  theme(legend.title = element_blank())

# # Line plot
# Fig3a.2 <- ggplot(master_plot[master_plot$Scenario == 'BAU',], aes(x = year, y = Percent_AOO_lost, colour = factor(IUCN_Current))) +
#   geom_smooth(method = 'loess', span=01, se=TRUE, aes(color = factor(IUCN_Current)), alpha=0.3) +
#   theme_classic() +
#   theme_classic() +
#   # facet_wrap(~Scenario) +
#   labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   # scale_y_continuous(limits = c(0, 40), expand = c(0,0), breaks = c(0,10,20,30)) +# scale_y_continuous(limits = c(0, 27), expand = c(0,0), breaks = c(0,5,10,15,20,25)) +
#   # scale_x_reverse() +
#   # coord_flip() +
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
#   theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
#   annotate(geom = 'text', x = 2010, y = 26, label = "(a)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
#   theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(legend.justification=c(0,0), legend.position=c(.27,.83), legend.direction = 'horizontal') +
#   guides(color = guide_legend(ncol=3)) +
#   theme(legend.title = element_blank())



# Bar plot
# Fig3a <- ggplot(iucn_status_aoo_lost[iucn_status_aoo_lost$year == 2060 & iucn_status_aoo_lost$Scenario == "BAU",], aes(x = IUCN_Current, y = Percent_AOO_lost)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   # facet_wrap(~Scenario) +
#   labs(x = NULL, y = NULL) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 27), expand = c(0,0), breaks = c(0,5,10,15,20,25)) +
#   # scale_x_reverse() +
#   # coord_flip() +
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
#   theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
#   annotate(geom = 'text', x = .6, y = 26, label = "(a)", family = "Helvetica", fontface = 2, size = 6, colour = "black")


# Fig 3b----
# c) Looking at differences between habitat specialist vs generalist
master_plot$Specialist_Genearlist[master_plot$Specialist_Genearlist == "Quasi-Specialist"] <- "Specialist"
master_plot$Count = 1

generalist_frame_categorical <- ddply(master_plot[!is.na(master_plot$Specialist_Genearlist),],.
                                      (Scenario, year, Specialist_Genearlist),
                                      summarise,
                                      Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE),
                                      Count = sum(Count, na.rm = TRUE))

# Getting sd
generalist_frame_categorical$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(generalist_frame_categorical)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == generalist_frame_categorical$Scenario[i] & master_plot$year == generalist_frame_categorical$year[i] & master_plot$Specialist_Genearlist == generalist_frame_categorical$Specialist_Genearlist[i]], na.rm = TRUE)
  # Updating sd
  generalist_frame_categorical$sdPercent_AOO_lost[i] <- tmp.sd
}

generalist_frame_categorical$Upper <- generalist_frame_categorical$Percent_AOO_lost + generalist_frame_categorical$sdPercent_AOO_lost / sqrt(generalist_frame_categorical$Count)
generalist_frame_categorical$Lower <- generalist_frame_categorical$Percent_AOO_lost - generalist_frame_categorical$sdPercent_AOO_lost / sqrt(generalist_frame_categorical$Count)
# And plotting this
# Line plot through time

# generalist_frame_categorical$Percent_AOO_lost[generalist_frame_categorical$year == 2010] <- 0

Fig3b.2 <- ggplot(generalist_frame_categorical[generalist_frame_categorical$Scenario == 'BAU',], aes(x = year, y = Percent_AOO_lost, colour = factor(Specialist_Genearlist))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, color = factor(Specialist_Genearlist)), linetype = 2) + 
  geom_line(aes(x = year, y = Lower, color = factor(Specialist_Genearlist)), linetype = 2) +
  theme_classic() +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2010, y = 46, label = "(b)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  theme(legend.justification=c(0,0), legend.position=c(.22,.9), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=3)) +
  theme(legend.title = element_blank())

# # Bar plot in 2060
# Fig3b <- ggplot(generalist_frame_categorical[generalist_frame_categorical$year == 2060 & generalist_frame_categorical$Scenario == "BAU",], aes(x = Specialist_Genearlist, y = Percent_AOO_lost)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   # facet_wrap(~Scenario) +
#   labs(x = NULL, y = NULL) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 27), expand = c(0,0), breaks = c(0,5,10,15,20,25)) +
#   # scale_x_reverse() +
#   # coord_flip() +
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
#   theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
#   annotate(geom = 'text', x = .5, y = 26, label = "(b)", family = "Helvetica", fontface = 2, size = 6, colour = "black")




# Fig 3c----
master_plot$Ag_Suitable <- master_plot$Pastureland + master_plot$Cropland

master_plot$Ag_Suitable_Classification <- "Semi Ag Tolerant"
master_plot$Ag_Suitable_Classification[master_plot$Ag_Suitable == 0] <- "Ag Tolerant"
master_plot$Ag_Suitable_Classification[master_plot$Ag_Suitable == 2] <- "Ag Intolerant"

only_nonag_species <- ddply(master_plot,.
                            (year, Scenario, Ag_Suitable_Classification),
                            summarise,
                            Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE),
                            Count = sum(Count, na.rm = TRUE))

# Getting standard deviations
# Getting sd
only_nonag_species$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(only_nonag_species)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == only_nonag_species$Scenario[i] & master_plot$year == only_nonag_species$year[i] & master_plot$Ag_Suitable_Classification == only_nonag_species$Ag_Suitable_Classification[i]], na.rm = TRUE)
  # Updating sd
  only_nonag_species$sdPercent_AOO_lost[i] <- tmp.sd
}

only_nonag_species$Upper <- only_nonag_species$Percent_AOO_lost + only_nonag_species$sdPercent_AOO_lost / sqrt(only_nonag_species$Count)
only_nonag_species$Lower <- only_nonag_species$Percent_AOO_lost - only_nonag_species$sdPercent_AOO_lost / sqrt(only_nonag_species$Count)

# Only species that can't survive in crop or pasture
# Fig3c <- ggplot(only_nonag_species[only_nonag_species$year == 2060 & only_nonag_species$Scenario == 'BAU' & only_nonag_species$Ag_Suitable_Classification != "Ag Tolerant",], aes(x = Ag_Suitable_Classification, y = Percent_AOO_lost)) +
#   geom_line(size = 1.5) +
#   geom_line(aes(x = year, y = Upper, color = factor(Specialist_Genearlist)), linetype = 2) + 
#   geom_line(aes(x = year, y = Lower, color = factor(Specialist_Genearlist)), linetype = 2) +
#   theme_classic() +
#   theme_classic() +
#   # facet_wrap(~Scenario) +
#   labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
#   # scale_x_reverse() +
#   # coord_flip() +
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
#   theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
#   annotate(geom = 'text', x = 2010, y = 46, label = "(c)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
#   theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(legend.justification=c(0,0), legend.position=c(.27,.83), legend.direction = 'horizontal') +
#   guides(color = guide_legend(ncol=3)) +
  theme(legend.title = element_blank())


Fig3c.2 <- ggplot(only_nonag_species[only_nonag_species$Scenario == 'BAU' & only_nonag_species$Ag_Suitable_Classification != "Ag Tolerant",], aes(x = year, y = Percent_AOO_lost, colour = factor(Ag_Suitable_Classification))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, color = factor(Ag_Suitable_Classification)), linetype = 2) + 
  geom_line(aes(x = year, y = Lower, color = factor(Ag_Suitable_Classification)), linetype = 2) +
  theme_classic() +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2010, y = 46, label = "(c)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  theme(legend.justification=c(0,0), legend.position=c(.3,.83), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=3)) +
  guides(color = guide_legend(ncol=1)) +
  theme(legend.title = element_blank())


# Fig3d----
mass_class_frame <- ddply(master_plot,.
                          (year, Scenario, mass_class),
                          summarise,
                          Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE),
                          Count = sum(Count, na.rm = TRUE))

# Getting standard deviations
# Getting sd
mass_class_frame$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(mass_class_frame)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == mass_class_frame$Scenario[i] & master_plot$year == mass_class_frame$year[i] & master_plot$mass_class == mass_class_frame$mass_class[i]], na.rm = TRUE)
  # Updating sd
  mass_class_frame$sdPercent_AOO_lost[i] <- tmp.sd
}

mass_class_frame$Upper <- mass_class_frame$Percent_AOO_lost + mass_class_frame$sdPercent_AOO_lost / sqrt(mass_class_frame$Count)
mass_class_frame$Lower <- mass_class_frame$Percent_AOO_lost - mass_class_frame$sdPercent_AOO_lost / sqrt(mass_class_frame$Count)



# # Plotting AOO lost vs species size
# Fig3d <- ggplot(mass_class_frame_2060_bau, aes(x = mass_class, y = Percent_AOO_lost)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   # facet_wrap(~Scenario) +
#   labs(x = NULL, y = NULL) +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 27), expand = c(0,0), breaks = c(0,5,10,15,20,25)) +
#   # scale_x_reverse() +
#   # coord_flip() +
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
#   theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
#   theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
#   annotate(geom = 'text', x = .55, y = 26, label = "(d)", family = "Helvetica", fontface = 2, size = 6, colour = "black")

mass_class_frame <- mass_class_frame[!is.na(mass_class_frame$mass_class),]
mass_class_frame$Percent_AOO_lost[mass_class_frame$year == 2010] <- 0

mass_class_frame <- mass_class_frame[!is.na(mass_class_frame$mass_class),]

Fig3d.2 <- ggplot(mass_class_frame[mass_class_frame$Scenario == 'BAU',], aes(x = year, y = Percent_AOO_lost, colour = factor(mass_class))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, color = factor(mass_class)), linetype = 2) + 
  geom_line(aes(x = year, y = Lower, color = factor(mass_class)), linetype = 2) +
  theme_classic() +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2010, y = 46, label = "(d)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  theme(legend.justification=c(0,0), legend.position=c(.23,.9), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=3)) +
  theme(legend.title = element_blank())


# Fig3e ----
# ----
# Getting changes in ESH by habitat preferences
# Thanks DW!
# habitat preferences
bird_preferences <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/Habitat_Preferences/Habitats2017.csv",
                             stringsAsFactors = FALSE)
# habitat frame
master_frame <- read.csv("/Users/macuser/Desktop/Ag Expansion and Biodiversity/ESH_Data/",
                         stringsAsFactors = FALSE)
# Limiting habitat preferences to primary habitat type
bird_preferences <- unique(dplyr::select(bird_preferences, Scientific.name, Habitats.Classification.Scheme.Level.1, Suitability))
# Limiting to only "marginal" and "suitable" habitats, then getting unique values
bird_preferences <- bird_preferences[bird_preferences$Suitability == "Suitable",]
bird_preferences <- unique(dplyr::select(bird_preferences,
                                         Scientific.name,
                                         Habitats.Classification.Scheme.Level.1))
# Changing column names
names(bird_preferences)[grep("Habitats", names(bird_preferences))] <- "Habitat.Preference"
names(bird_preferences)[grep("Scientific", names(bird_preferences))] <- "species"

# Changing format of species to merge with esh data frame
bird_preferences$species <- gsub(" ",
                                 "_",
                                 bird_preferences$species)

# Merging
master_plot <- left_join(master_plot,
                          bird_preferences)

# Narrow down classifications
master_plot$habitat_summary <- master_plot$Habitat.Preference
master_plot$habitat_summary[grep("Marine",master_plot$habitat_summary)] <- "Marine"
master_plot$habitat_summary[grep("Artificial",master_plot$habitat_summary)] <- "Artificial"


habitat_frame <- ddply(master_plot,.
                       (habitat_summary, year, Scenario),
                       summarise,
                       Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE),
                       Count = sum(Count, na.rm = TRUE))


# Getting standard deviations
# Getting sd
habitat_frame$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(habitat_frame)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == habitat_frame$Scenario[i] & master_plot$year == habitat_frame$year[i] & master_plot$habitat_summary == habitat_frame$habitat_summary[i]], na.rm = TRUE)
  # Updating sd
  habitat_frame$sdPercent_AOO_lost[i] <- tmp.sd
}

habitat_frame$Upper <- habitat_frame$Percent_AOO_lost + habitat_frame$sdPercent_AOO_lost / sqrt(habitat_frame$Count)
habitat_frame$Lower <- habitat_frame$Percent_AOO_lost - habitat_frame$sdPercent_AOO_lost / sqrt(habitat_frame$Count)


# changing names
habitat_frame$habitat_summary[grep("Caves", habitat_frame$habitat_summary)] <- "Subterranean"
habitat_frame$habitat_summary[grep("Rocky", habitat_frame$habitat_summary)] <- "Rocky areas"
habitat_frame$habitat_summary[grep("Wetlands", habitat_frame$habitat_summary)] <- "Wetlands"
habitat_frame <- habitat_frame[!is.na(habitat_frame$habitat_summary),]



habitat_frame1 <- habitat_frame[habitat_frame$year == 2060 & habitat_frame$Scenario == 'BAU',]
habitat_frame1 <- habitat_frame1[order(habitat_frame1$Percent_AOO_lost),]

habitat_frame1 <- transform(habitat_frame1,
                            habitat_summary = factor(habitat_summary,
                                                     levels = rev(habitat_frame1$habitat_summary)))


Fig3e <- ggplot(habitat_frame1, aes(x = habitat_summary, y = Percent_AOO_lost)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = NULL) +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 20), expand = c(0,0), breaks = c(0,5,10,15)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black", angle = 90)) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = .75, y = 19.259, label = "(e)", family = "Helvetica", fontface = 2, size = 6, colour = "black")


# habitat_frame$Percent_AOO_lost[habitat_frame$year == 2010] <- 0
Fig3e.2 <- ggplot(habitat_frame[habitat_frame$Scenario == 'BAU',], aes(x = year, y = Percent_AOO_lost, colour = factor(habitat_summary))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, color = factor(habitat_summary)), linetype = 2) + 
  geom_line(aes(x = year, y = Lower, color = factor(habitat_summary)), linetype = 2) +
  theme_classic() +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Remaining Area of \nHabitat Lost (Percent)") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 48), expand = c(0,0), breaks = c(0,10,20,30,40)) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2010, y = 46, label = "(e)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  theme(legend.justification=c(0,0), legend.position=c(.27,.72), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=3)) +
  theme(legend.title = element_blank())



# # And making the entirety of plot 3
# Fig3.A <- ggdraw(plot_grid(Fig3a, Fig3b, Fig3c, Fig3d,
#                   nrow = 2,
#                   align = 'hv'))
# Fig3.B <- ggdraw(plot_grid(Fig3.A,
#                  Fig3e,
#                  nrow = 2,
#                  rel_heights = c(3,2)))
# 
# Fig3.B
# 
# ggsave("/Users/macuser/Desktop/Thesis Chapters/In Progress/Land Expansion Biodiversity/Figures/Fig3.pdf",
#        height = 37.5,
#        width = 30,
#        units = 'cm')

# And making the entirety of plot 3
Fig3.A2 <- ggdraw(plot_grid(Fig3a.2, Fig3b.2, Fig3c.2, Fig3d.2,
                           nrow = 2,
                           align = 'hv'))
Fig3.B2 <- ggdraw(plot_grid(Fig3.A2,
                           Fig3e.2,
                           nrow = 2,
                           rel_heights = c(2,1)))

Fig3.B2

ggsave("/Users/macuser/Desktop/Thesis Chapters/In Progress/Land Expansion Biodiversity/Figures/Fig3.2.pdf",
       height = 37.5,
       width = 30,
       units = 'cm')


# Fig 4----

# Assume species don't get habitat back
# master_plot$Percent_AOO_lost[master_plot$Percent_AOO_lost < 0 & !is.na(master_plot$Percent_AOO_lost)] <- 0

# AOH loss avoided under alternative scenarios
mean_esh_loss <- ddply(master_plot,.
                       (Scenario, year),
                       summarise,
                       esh_loss_mean = mean(Percent_AOO_lost, na.rm = TRUE),
                       Count = sum(Count, na.rm = TRUE))

# # Getting upper and lower limts of esh loss for scenarios
# mean_esh_loss$Upper <- mean_esh_loss$esh_loss_mean + mean_esh_loss$esh_loss_sd
# mean_esh_loss$Lower <- mean_esh_loss$esh_loss_mean - mean_esh_loss$esh_loss_sd


# Getting standard deviations
# Getting sd
mean_esh_loss$sdPercent_AOO_lost = NA
# Manually getting standard deviations
for(i in 1:nrow(mean_esh_loss)) {
  # Getting sd
  tmp.sd <- sd(master_plot$Percent_AOO_lost[master_plot$Scenario == mean_esh_loss$Scenario[i] & master_plot$year == mean_esh_loss$year[i]], na.rm = TRUE)
  # Updating sd
  mean_esh_loss$sdPercent_AOO_lost[i] <- tmp.sd
}

mean_esh_loss$Upper <- mean_esh_loss$esh_loss_mean + mean_esh_loss$sdPercent_AOO_lost
mean_esh_loss$Lower <- mean_esh_loss$esh_loss_mean - mean_esh_loss$sdPercent_AOO_lost

# Getting BAU data and merging back in
mean_esh_loss_bau <- mean_esh_loss[mean_esh_loss$Scenario == "BAU",]
names(mean_esh_loss_bau)[3] <- "meanloss_BAU"
names(mean_esh_loss_bau)[6:7] <- c("Upper_BAU", "Lower_BAU")

mean_esh_loss <- left_join(mean_esh_loss,
                           dplyr::select(mean_esh_loss_bau, year, meanloss_BAU, Upper_BAU, Lower_BAU))


# Getting proportion esh loss avoided
# Mean uppder and lower
mean_esh_loss$Avoided_mean <- (1 - mean_esh_loss$esh_loss_mean / mean_esh_loss$meanloss_BAU) * 100
mean_esh_loss$Avoided_upper <- (1 - mean_esh_loss$Upper / mean_esh_loss$Lower_BAU) * 100
mean_esh_loss$Avoided_lower <- (1 - mean_esh_loss$Lower / mean_esh_loss$Upper_BAU) * 100


# Updating values to 0 for 2010
# mean_esh_loss$Avoided_mean[mean_esh_loss$year == 2010] <- 0
# mean_esh_loss$Avoided_upper[mean_esh_loss$year == 2010] <- 0
# mean_esh_loss$Avoided_lower[mean_esh_loss$year == 2010] <- 0
# 
mean_esh_loss <- transform(mean_esh_loss, 
                           Scenario = factor(Scenario,
                                             levels = c("Diets","Trade","Yields","All","BAU")))






# Plotting esh loss avoided in 2060
Fig4a <- ggplot(mean_esh_loss[mean_esh_loss$Scenario != "BAU" & mean_esh_loss$year == 2060,], aes(x = Scenario, y = Avoided_mean)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(ymax = Avoided_upper, ymin = Avoided_lower), colour = 'black') +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Habitat \nDecline Avoided") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0), breaks = c(0,20,40,60,80)) +
  scale_color_discrete(guide = guide_legend(title = "Scenario", title.theme = element_text(size = 14, family = 'Helvetica', face = 'bold', angle = 0))) +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black", angle = 90)) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = .7, y = 80, label = "(a)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))




# Plotting esh loss  through time
Fig4b <- ggplot(mean_esh_loss[mean_esh_loss$Scenario != "BAU",], aes(x = year, y = esh_loss_mean, colour = factor(Scenario))) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper, colour = factor(Scenario)), linetype = 2) +
  geom_line(aes(x = year, y = Lower, colour = factor(Scenario)), linetype = 2) +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Remaining \nArea of Habitat Lost") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(-10, 25), expand = c(0,0), breaks = c(0,10,20)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_discrete(guide = guide_legend(title = "Scenario", title.theme = element_text(size = 14, family = 'Helvetica', face = 'bold', angle = 0))) +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  annotate(geom = 'text', x = 2014, y = 23.5, label = "(b)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.justification=c(0,0), legend.position=c(.27,.83), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=2)) +
  theme(legend.title = element_blank()) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))


Fig4 <- ggdraw(plot_grid(Fig4a, Fig4b,
                         nrow = 1,
                         align = 'hv'))


Fig4

ggsave("/Users/macuser/Desktop/Thesis Chapters/In Progress/Land Expansion Biodiversity/Figures/Fig4.pdf",
       height = 15,
       width = 30,
       units = 'cm')


# a) Getting median ESH loss through time
median_esh_loss <- ddply(master_plot,.
                         (Scenario, year),
                         summarise,
                         esh_loss_median = median(Percent_AOO_lost, na.rm = TRUE),
                         esh_loss_5th = quantile(Percent_AOO_lost, .05, na.rm = TRUE),
                         esh_loss_95th = quantile(Percent_AOO_lost, .95, na.rm = TRUE))

# b) Getting mean and sd ESH loss through time
mean_esh_loss <- ddply(master_plot,.
                         (Scenario, year),
                         summarise,
                         esh_loss_mean = mean(Percent_AOO_lost, na.rm = TRUE),
                         esh_loss_sd = sd(Percent_AOO_lost, na.rm = TRUE))
mean_esh_loss$Upper <- mean_esh_loss$esh_loss_mean + mean_esh_loss$esh_loss_sd
mean_esh_loss$Lower <- mean_esh_loss$esh_loss_mean - mean_esh_loss$esh_loss_sd
                         

# Looks like no difference in organism sizes
# Possibly because we're not picking up on fragmentation?
# Or absolute decreases in species range?
# Looking at change in RLI based on species size
# Plotting AOO lost vs species size


# Boxplot
tmp.boxplot <- master_plot[master_plot$IUCN_Current == "LC" |
                             master_plot$IUCN_Current == "NT" |
                             master_plot$IUCN_Current == "VU" |
                             master_plot$IUCN_Current == "EN" |
                             master_plot$IUCN_Current == "CR",]
tmp.boxplot <- transform(tmp.boxplot,
                         IUCN_Current = factor(IUCN_Current,
                                               levels = c("LC","NT","VU","EN","CR")))

ggplot(tmp.boxplot[tmp.boxplot$year == 2060,], aes(x = IUCN_Current, y = Percent_AOO_lost, fill = factor(IUCN_Current))) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Remaining Habitat Lost") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(-25, 25))

# # Average increase in IUCN status
# increase_iucn_status <- ddply(master_plot[!is.na(master_plot$IUCN_Current_linear) & master_plot$Linear_Change_IUCN_Status >= 0,],.
#                               (IUCN_Current, year, Scenario),
#                               summarise,
#                               meanChange_Status = mean(Linear_Change_IUCN_Status, na.rm = TRUE),
#                               medianChange_Status = median(Linear_Change_IUCN_Status, na.rm = TRUE))
# # Plotting, bar in 2060
# ggplot(increase_iucn_status[increase_iucn_status$year == 2060,], aes(x = IUCN_Current, y = meanChange_Status, fill = factor(IUCN_Current))) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   facet_wrap(~Scenario) +
#   labs(x = NULL, y = "Percent of Remaining Habitat Lost") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 3.5))
# 
# ggplot(increase_iucn_status[increase_iucn_status$year == 2060,], aes(x = IUCN_Current, y = medianChange_Status, fill = factor(IUCN_Current))) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   facet_wrap(~Scenario) +
#   labs(x = NULL, y = "Percent of Remaining Habitat Lost") +
#   geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 3.5))
# 
# # Proportion of species that increase in IUCN status
# prop_increase_status <- ddply(master_plot,.
#                               (IUCN_Current, year, Scenario),
#                               summarise,
#                               count = sum(Count, na.rm = TRUE),
#                               NIncrease = sum(Increase_IUCN_Status, na.rm = TRUE))
# # Getting percent of species increase that increased in status
# prop_increase_status$PercentIncrease <- (prop_increase_status$NIncrease / prop_increase_status$count) * 100
# 
# # Line plot
# ggplot(prop_increase_status, aes(x = year, y = PercentIncrease, colour = factor(IUCN_Current))) +
#   geom_point() +
#   geom_line() +
#   theme_classic() +
#   facet_wrap(~Scenario) +
#   labs(x = NULL, y = "Percent Species that Increased in IUCN Status") +
#   # geom_hline(yintercept = 0, linetype = 2) +
#   # geom_hline(yintercept = 25, linetype = 2) +
#   # geom_hline(yintercept = 50, linetype = 2) +
#   # geom_hline(yintercept = 75, linetype = 2) +
#   # geom_hline(yintercept = -25, linetype = 2) +
#   # geom_hline(yintercept = -50, linetype = 2) +
#   # geom_hline(yintercept = -75, linetype = 2) +
#   scale_y_continuous(limits = c(0, 100))



# Looking at how changes in ESH vary by current status and body mass size
body.mass_species.status <- ddply(master_plot,.
                                  (Scenario, year, IUCN_Current, mass_class),
                                  summarise,
                                  Percent_AOO_lost = mean(Percent_AOO_lost, na.rm = TRUE))

# Plotting
ggplot(body.mass_species.status[body.mass_species.status$year == 2060 & !is.na(body.mass_species.status$mass_class),], aes(x = IUCN_Current, y = Percent_AOO_lost, fill = factor(mass_class))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() +
  facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Remaining Habitat Lost") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(-25, 25))


# Step something ----
# Plotting!
# Percent of ESH lost relative to 2010
# All species
ggplot(master_plot, aes(x = year, y = Percent_AOO_lost, group = year)) +
  geom_boxplot() + 
  theme_classic() +
  facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Remaining Habitat Lost") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(yintercept = 25, linetype = 2) +
  geom_hline(yintercept = 50, linetype = 2) +
  geom_hline(yintercept = 75, linetype = 2) +
  geom_hline(yintercept = -25, linetype = 2) +
  geom_hline(yintercept = -50, linetype = 2) +
  geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(-100, 100))




dat_upper <- dplyr::select(mean_esh_loss[mean_esh_loss$Scenario == 'BAU',],
                           Scenario,
                           year,
                           Upper)
dat_lower <- dplyr::select(mean_esh_loss[mean_esh_loss$Scenario == 'BAU',],
                           Scenario,
                           year,
                           Lower)
names(dat_upper)[3] <- "SD"
names(dat_lower)[3] <- "SD"

dat2 <- rbind(dat_upper,dat_lower)
dat3 <- data.frame(year = c(2010, 2060, 2060),
                   SD = c(0, 57.199004, 0.6685141))

+geom_polygon(data = dat2, aes(y = bgl, 
                               group = group), alpha = 0.3)


ggplot(mean_esh_loss[mean_esh_loss$Scenario == "BAU",], aes(x = year, y = esh_loss_mean)) +
  geom_line(size = 1.5) +
  geom_line(aes(x = year, y = Upper), linetype = 2, colour = 'black') +
  geom_line(aes(x = year, y = Lower), linetype = 2, colour = 'black') + 
  geom_polygon(data = dat2, aes(x = year, y = SD), alpha = 0.3) +
  geom_polygon(data = dat3, aes(x = year, y = SD), alpha = 0.3) +
  theme_classic() +
  # facet_wrap(~Scenario) +
  labs(x = NULL, y = "Percent of Remaining \nArea of Habitat Lost") +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 25, linetype = 2) +
  # geom_hline(yintercept = 50, linetype = 2) +
  # geom_hline(yintercept = 75, linetype = 2) +
  # geom_hline(yintercept = -25, linetype = 2) +
  # geom_hline(yintercept = -50, linetype = 2) +
  # geom_hline(yintercept = -75, linetype = 2) +
  scale_y_continuous(limits = c(-10, 60), expand = c(0,0), breaks = c(0,20,40)) +
  scale_x_continuous(expand = c(0,0)) +
  # scale_color_discrete(guide = guide_legend(title = "Scenario", title.theme = element_text(size = 14, family = 'Helvetica', face = 'bold', angle = 0))) +
  theme(legend.text = element_text(size = 15,family = "Helvetica",face = "bold",colour = "black")) +
  # scale_x_reverse() +
  # coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(axis.title = element_text(size = 18, family = "Helvetica", face = "bold", colour = "black")) +
  theme(axis.text.x = element_text(size = 18,family = "Helvetica",face = "bold",colour = "black")) +
  theme(axis.text.y = element_text(size = 18,family = "Helvetica",face = "bold", colour = "black")) +
  # annotate(geom = 'text', x = 2014, y = 23.5, label = "(b)", family = "Helvetica", fontface = 2, size = 6, colour = "black") +
  theme(legend.justification=c(0,0), legend.position=c(.27,.83), legend.direction = 'horizontal') +
  guides(color = guide_legend(ncol=2)) +
  theme(legend.title = element_blank()) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))


ggsave("/Users/macuser/Desktop/Presentations/Thesis Seminar/ESH_Loss_BAU.pdf",
       height = 20,
       width = 30,
       units = 'cm')
