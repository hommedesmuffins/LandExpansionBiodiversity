#!/usr/bin/env Rscript

# Read me ----
# I think this code should be quicker than the for-loop extravaganza, but needs 
#     to be speed-checked and also actually checked...

# Packages ---
# Keep these to a bare minimum as they take up space
library(raster)
library(data.table)

# Getting 2010 crop
crop_2010 <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Crop2010.Mollweide.1.5km.tif")
crop_2010_vector <- cbind(c(getValues(crop_2010)))

# Getting row index vectors [DON'T CHANGE THESE FOR EACH SCENARIO] -----
# Rasters for the row indices of individual countries
row_index_rasters <- list.files(path = "/Users/maclark/Desktop/Mike's Files/Region Rasters/Fine_Country", 
                                pattern = "Index",
                                full.names = TRUE)

# Get a global raster of row indices
country_id_raster_global <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/MollweideCountryID_1.5km.tif")

# Raster for original crop / pasture extent [CHANGE FOR CROP VS. PASTURE] ----
original_raster <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Crop2010.Mollweide.1.5km.tif")

# Vector of file names for changed country rasters [CHANGE THESE FOR EACH SCENARIO AND FOR PASTURE!]------
# Get ALL the rasters one year / scenario combination
# Just vary thesefiles for pasture and the different scenario / year combinations
crop_rasters_long <- list.files("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/YieldGaps/TIFs", 
                           pattern = "crop",
                           full.names = TRUE)
# Limiting to mean rasters
crop_rasters_long <- crop_rasters_long[grep('sd', crop_rasters_long)]

# Getting years
yrs <- list.files("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/BAU/TIFs", 
                                pattern = "crop",
                                full.names = TRUE)
# Limiting to mean rasters
yrs <- yrs[grep('sd', yrs)]
yrs <- yrs[grep('AUS', yrs)]
yrs <- gsub("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/BAU/TIFs/AUSOceaniaprop_crop_",
            "",
            yrs)

for(y in 1:length(yrs)) {
  # Getting years of the raster we're looping over
  crop_rasters <- crop_rasters_long[grep(yrs[y], crop_rasters_long)]
  
  world_dt <- data.table(
    cell = 1:289285431,
    crop_orig = rep(0, 289285431))
  
  world_dt$crop_new <- world_dt$crop_orig 
  # This last line is not totally necessary and makes the resulting 
  #     datatable much larger, but it allows me to check the results.
  # When you run this in reality, I suggest you replace it with:
  # world_dt <- data.table(
  #   cell = 1:289285431,
  #   crop_new = as.vector(matrix(getValues(original_raster) / 100,
  #                                nrow = nrow(original_raster),
  #                                byrow = TRUE))
  # )
  
  setkey(world_dt, cell) # This is for putting the values in. I'll be honest, I'm not
  # 100% sure how this works, but you want the datatable ordered by cell.
  # It's also not 100% clear if I need to do this...
  # Sys.time() - t1 # 42 secs, but only needs to be done once for all the scenarios and years
  
  # Put the individual country rasters into this table ----
  # t1 <- Sys.time()
  dt_copy <- world_dt # Create a copy of the data.table, so we don't have to re-make it
  
  # Cycle through the crop rasters...
  for(i in 1:length(crop_rasters)){
    print(i) # Print progress
    # Create a data.table with the proportion of cropland and the row indexes.
    # Note the fucking around with the crop_raster name to get the right file...
    tmp_dt <- data.table(crop = getValues(raster(crop_rasters[[i]])),
                         row_index = getValues(raster(row_index_rasters[grep(substring(gsub("^.*/TIFs/", "", 
                                                                                            crop_rasters[[i]]), 
                                                                                       first = 1, last =3),
                                                                             row_index_rasters)])))
    # Hopefully not needed, but better safe than sorry
    tmp_dt <- tmp_dt[order(tmp_dt$row_index)]
    # Then use insane data.table syntax to replace the values in the data.table.
    # NOTE: DO NOT USE "<-" AS IT SLOWS EVERYTHING WAAAAAAY DOWN
    # http://datatable.r-forge.r-project.org/datatable-intro.pdf
    
    dt_copy[tmp_dt[!is.na(tmp_dt$crop), row_index], 
            # Returns the row_index of the cells in the temporary data.table 
            #     where the cropland cell is not NA
            crop_new := tmp_dt[!is.na(tmp_dt$crop), crop]]
    # the ":=" syntax is some madcap binary search thing which is lightning
  }
  # Making raster
  tmp_new <- raster(matrix(dt_copy$crop_new, 
                           nrow = nrow(country_id_raster_global),
                           byrow = FALSE))
  
  # Because I don't know how to do this the correct way...
  # Converting raster back into a matrix
  tmp.crop.vector <- cbind(c(getValues(tmp_new)))
  
  # Updating NAs
  tmp.crop.vector[is.na(crop_2010_vector)] <- NA
  
  # And also calculating delta crop extent
  # delta.ag.vector <- tmp.crop.vector - crop_2010_vector/100
  
  # And converting back into a raster
  # First absolute crop extent
  tmp_new <- raster(matrix(tmp.crop.vector,
                    nrow = nrow(country_id_raster_global),
                    byrow = TRUE))
  crs(tmp_new) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  extent(tmp_new) <- extent(original_raster)
  # plot(tmp_new)
  # Second delta crop extent
  # tmp_new_delta <- raster(matrix(delta.ag.vector,
  #                   nrow = nrow(country_id_raster_global),
  #                   byrow = TRUE))
  # crs(tmp_new_delta) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  # extent(tmp_new_delta) <- extent(original_raster)
  # plot(tmp_new_delta,
       # zlim = c(-1,1))
  
  # And saving raster
  # First absoluteextent
  writeRaster(tmp_new,
              paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/YieldGaps/World/",
                    "Crop_sd_",
                    yrs[y],
                    sep = ""))
  # Second delta extent
  # writeRaster(tmp_new_delta,
  #             paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/BAU/World/",
  #                   "DeltaCrop_",
  #                   yrs[y],
  #                   sep = ""))
  # 
  
  # And saving png
  # First absolute extent
  png(filename = paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/YieldGaps/World/",
                       "Crop_sd_",
                       gsub(".tif","",yrs[y]),
                       ".png",
                       sep = ""),
      height = 20, width = 40, units = "cm", res = 72)
  plot(tmp_new, zlim = c(-.1,1.1))
  dev.off()
  # Second delta extent
  # png(filename = paste("/Users/maclark/Desktop/Mike's Files/Spatial Land Forecasts/Outputs/Mollweide_1.5km/March_19/BAU/World/",
  #                      "DeltaCrop_",
  #                      gsub(".tif","",yrs[y]),
  #                      ".png",
  #                      sep = ""),
  #     height = 20, width = 40, units = "cm", res = 72)
  # plot(tmp_new_delta, zlim = c(-1,1))
  # dev.off()
  
}

# Save file, repeat for different scenario / year combinations ----
# I'm not sure the best way to save these - it might be to convert them to rasters?
#
# For the different scenarios and years, I would perhaps get a vector with all the 
#     scenario / year combinations and then run through it using d_ply or lapply or
#     something - this is likely to be faster than a for() loop.
# However you do it, obviously clear everything EXCEPT the world_dt after the loop
quit(save = "no", status = 0, runLast = TRUE)