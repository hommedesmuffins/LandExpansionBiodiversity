# Read me ----
# Script to hold the land clearing functions
#
# This version differs from the previous (_v1) because the modelling procedure has changed.
#
# Now we get the probability of expansion, contraction or no change in the first step 
#     (previously it was simply "change" vs. "no change")
# We then model expanding and contracting cells separately. For countries that have 
#     expanding agriculture, we have a proportion of their cells contracting, total up 
#     the area lost and then expand agriculture
#
# pred_fun takes model inputs, assigns probabilities to cells, and how much a cell would expand or contract
# conv_fun picks those cells until target is met
# yr_fun goes through every decade and runs conv_fun until targets are met

# Function for sampling -----
# Much faster than sample()
fastSampleReject <- function(all, n, w){
  out <- numeric(0)
  while(length(out) < n)
    out <- unique(c(out, sample(all, size = n, replace = TRUE, prob = w)))
  out[1:n]
}

# Urban scenario creation -----
# Function to do the picking
urban_scen_fun <- function(x, time_step){
  # Set up a counter for the amount of conversion
  x$converted <- 0
  # Cycle through each time step
  for(i in 1:length(time_step)){
    # Get the name of a column to add
    d <- paste("urb_targ",
               time_step[i], 
               sep = "_")
    # Is any conversion needed?
    if(x$converted < x$no_to_convert){
      # Flip a (weighted) coin to decide if I want to convert the "possible" cells
      #     or not
      if(x$converted + x$pt_definite < x$no_to_convert){
        # The target is the number of definite cells, plus a coin flip for the 
        #     possible cell. The flip is weighted by probability of the possible
        #     cell being converted
        x[[d]] <- x$pt_definite +  rbinom(n = nrow(x), 
                                          size = 1, 
                                          prob=x$pt_possible)
        x$converted <- x$converted + x[[d]]
      } else {
        # If the target doesn't require a possible cell, then don't include it
        x[[d]] <- x$pt_definite
        x$converted <- x$converted + x[[d]]
      }
    } else {
      # If no conversion is needed, set the target to zero
      x[[d]] <- 0
    }
  }
  return(x)
}

# Function to actually convert - this function is just to allow me to cycle through 
#     countries, after the wrapping dplyr loop cycles through runs
urb_conv_fun <- function(y){
  output <- urban_target %>%
    group_by(ISO3.Numeric) %>%
    do(urban_scen_fun(., time_step = time_step))
  return(output)
}


# Proportion of surrounding ag or pasture land ----
convolution_fun <- function(m){
  m_rows <- nrow(m)
  m_cols <- ncol(m)
  m_tmp <- rbind(matrix(0, nrow = 1, ncol = m_cols), m)[-(m_rows + 1),]
  m_tmp <- m_tmp + 
    rbind(matrix(0, nrow = 1, ncol = m_cols),
          cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1])[-(m_rows + 1),]
  m_tmp <- m_tmp + 
    cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1]
  m_tmp <- m_tmp + 
    rbind(cbind(m, matrix(0, nrow = m_rows, ncol = 1))[,-1],
          matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + 
    rbind(m, matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + 
    rbind(cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)],
          matrix(0, nrow = 1, ncol = m_cols))[-1,]
  m_tmp <- m_tmp + cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)]
  m_tmp <- m_tmp + 
    rbind(matrix(0, nrow = 1, ncol = m_cols),
          cbind(matrix(0, nrow = m_rows, ncol = 1), m)[,-(m_cols + 1)])[-(m_rows + 1),]
  return(m_tmp)
}

adj.ag.f <- function(m, r, c, crop){
  if(crop == "crop"){
    # Convert the column to a matrix
    tmp_mat <- matrix(m$prop_to_predict_from_crop, 
                      nrow = r, 
                      ncol = c, 
                      byrow = FALSE)
    # Replace NAs with 0s to allow convolution
    tmp_mat[is.na(tmp_mat)] <- 0
    # Run the convolution
    tmp_mat <- convolution_fun(tmp_mat)
    # Put data in and put NAs back where they should be
    m$prop_adj_to_predict_from_crop <- as.vector(tmp_mat)
    m[is.na(m$prop_to_predict_from_crop), "prop_adj_to_predict_from_crop"] <- NA
    return(m)
  } else {
    # Convert the column to a matrix
    tmp_mat <- matrix(m$prop_to_predict_from_pasture, 
                      nrow = r, 
                      ncol = c, 
                      byrow = FALSE)
    # Replace NAs with 0s to allow convolution
    tmp_mat[is.na(tmp_mat)] <- 0
    # Run the convolution
    tmp_mat <- convolution_fun(tmp_mat)
    # Put data in and put NAs back where they should be
    m$prop_adj_to_predict_from_pasture <- as.vector(tmp_mat)
    m[is.na(m$prop_to_predict_from_pasture), "prop_adj_to_predict_from_pasture"] <- NA
    return(m)
  }
}

# Old version ----
# adj.ag.f <- function(m, r, c, f, crop){
#   if(crop == "crop"){
#     # Convert the column to a matrix
#     tmp_mat <- matrix(m$prop_to_predict_from_crop, 
#                       nrow = r, 
#                       ncol = c, 
#                       byrow = FALSE)
#     # Replace NAs with 0s to allow convolution
#     tmp_mat[is.na(tmp_mat)] <- 0
#     # Run the convolution
#     trans_mat <- filter2(tmp_mat, 
#                          adj_cell_mat, 
#                          boundary = 0)
#     # Take off the value of the focal cell
#     trans_mat <- trans_mat - tmp_mat
#     # Need to round numbers as convolution isn't precise
#     trans_mat <- round(trans_mat, digits = 8)
#     # Put data in and put NAs back where they should be
#     m$prop_adj_to_predict_from_crop <- as.vector(trans_mat)
#     m[is.na(m$prop_to_predict_from_crop), "prop_adj_to_predict_from_crop"] <- NA
#     return(m)
#   } else {
#     # Convert the column to a matrix
#     tmp_mat <- matrix(m$prop_to_predict_from_pasture, 
#                       nrow = r, 
#                       ncol = c, 
#                       byrow = FALSE)
#     # Replace NAs with 0s to allow convolution
#     tmp_mat[is.na(tmp_mat)] <- 0
#     # Run the convolution
#     trans_mat <- filter2(tmp_mat, 
#                          adj_cell_mat, 
#                          boundary = 0)
#     # Take off the value of the focal cell
#     trans_mat <- trans_mat - tmp_mat
#     # Need to round numbers as convolution isn't precise
#     trans_mat <- round(trans_mat, digits = 8)
#     # Put data in and put NAs back where they should be
#     m$prop_adj_to_predict_from_pasture <- as.vector(trans_mat)
#     m[is.na(m$prop_to_predict_from_pasture), "prop_adj_to_predict_from_pasture"] <- NA
#     return(m)
#   }
# }


# Predicting probabilites and amount of conversion ----
# for troubleshooting, p = model_data_current
pred_fun <- function(p, crop){
  if(crop == "crop"){
    
    # And adding in new predictor variables
    # log(dist), log(dist)^2, etc
    p$dist_50k_log <- log(p$dist_50k + .001)
    p$dist_50k_log_sq <- p$dist_50k_log^2
    
    p$prop_crop_sq <- p$prop_to_predict_from_crop^2
    p$prop_pasture_sq <- p$prop_to_predict_from_pasture^2
    
    p$prop_adj_crop_sq <- p$prop_adj_to_predict_from_crop^2
    p$prop_adj_pasture_sq <- p$prop_adj_to_predict_from_pasture^2
    
    p$prop_cell_delta_crop_sq <- p$prop_to_predict_from_cell_delta_crop^2
    p$prop_cell_delta_pasture_sq <- p$prop_to_predict_from_cell_delta_pasture^2
    
    p$ISO_numeric = factor(p$ISO_numeric, ordered = FALSE)
    
    # Cropland -----
    # First: probability of change for crop
    # pred_tmp are temporary predictor values. Gives data frame with probability of either no change, increase, or decrease in crop extent
    # Manual model predictions
    pred_tmp <- multinom.predict.2b(dat = p,
                                    coefs = coefs.crop.prob,
                                    vars = predictor_vars_crop[1:12],
                                    factor_vars = predictor_vars_crop[13])
    
    # Put predictions from pred_tmp into the dataframe
    p$a_pred_to_no_change_crop <- pred_tmp$aNoChange
    p$pred_to_decrease_crop <- pred_tmp$Decrease
    p$pred_to_increase_crop <- pred_tmp$Increase 
    
    # Repeat for decreasing models
    # Proportion to decrease
    # And manually predicting
    # Have to take negative values because the undelrying models can only predict values between 0 and 1
    p$pred_change_cell_prop_crop_decrease <- -glm.predict.2b(dat = p,
                                                            coefs = coefs.crop.decrease,
                                                            vars = predictor_vars_crop[1:12],
                                                            factor_vars = predictor_vars_crop[13])
    
    # Avoiding NAs
    p$pred_change_cell_prop_crop_decrease[is.na(p$pred_change_cell_prop_crop_decrease)] <- 0
    
    # Repeat for increasing models
    # Proportion to increase
    p$pred_change_cell_prop_crop_increase <- glm.predict.2b(dat = p,
                                                             coefs = coefs.crop.increase,
                                                             vars = predictor_vars_crop[1:12],
                                                             factor_vars = predictor_vars_crop[13])
    # And avoiding NAs
    p$pred_change_cell_prop_crop_increase[is.na(p$pred_change_cell_prop_crop_increase)] <- 0
    
    # Check that the amount of cropland plus change plus urban land plus nonarable land is not greater than the area of the cell.
    # Creates temporary cropland to make sure potential future area of cell after crop expansion does not exceed 1
    # Creating temp data frame and converting NAs to 0s to avoid issues summing across columns
    p.tmp <- dplyr::select(p, 
                           urban_prop, 
                           Prop.NonArable.2010, 
                           prop_to_predict_from_crop, 
                           pred_change_cell_prop_crop_increase, 
                           pred_change_cell_prop_crop_decrease)
    # Converting NAs to 0s
    p.tmp[is.na(p.tmp)] <- 0
    # Creating new column to adjust values of predictions so that they cannot exceed 1
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'pred_change_cell_prop_crop_increase')],
                         na.rm = TRUE)



    # Adjust predictions if so
    # If area does exceed 1, then adjust so that predicted cropland change does not cause cell value to exceed 1
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_crop_increase"] <-
      1 - (p.tmp[p.tmp$tmp > 1, "urban_prop"] + 
             p.tmp[p.tmp$tmp > 1, "Prop.NonArable.2010"] + 
             p.tmp[p.tmp$tmp > 1, "prop_to_predict_from_crop"])
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <-     p.tmp$tmp <- rowSums(p.tmp[,c('pred_change_cell_prop_crop_increase',
                                                   'prop_to_predict_from_crop')],
                                          na.rm = TRUE)
    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_crop_increase"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_crop"]
    
    # Repeat for decreasing models
    # Check that the amount of cropland plus change plus urban land plus non arable land is not greater than the area of the cell.
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'pred_change_cell_prop_crop_decrease')],
                         na.rm = TRUE)
    
    # Adjust predictions if so
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_crop_decrease"] <-
      1 - rowSums(p.tmp[p.tmp$tmp>1,c('urban_prop',
                           'Prop.NonArable.2010',
                           'prop_to_predict_from_crop')],
                  na.rm = TRUE)
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <- rowSums(p.tmp[,c('pred_change_cell_prop_crop_decrease',
                                             'prop_to_predict_from_crop')],
                         na.rm = TRUE)

    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_crop_decrease"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_crop"]
    
    # Finally, check that fully urbanised areas can't be converted at all
    # If urban proportion is one, then cropland extent cannot change in these cells
    p[p$urban_prop == 1 & !is.na(p$urban_prop), c("pred_to_increase_crop",
                                                  "pred_to_decrease_crop",
                                                  "pred_change_cell_prop_crop_decrease", 
                                                  "pred_change_cell_prop_crop_increase")] <- 0
    
    # Inserting updated estimates for prop increase and prop decrease into original data frame
    p$pred_change_cell_prop_crop_decrease[!is.na(p$pred_change_cell_prop_crop_decrease)] <- 
      p.tmp$pred_change_cell_prop_crop_decrease[!is.na(p$pred_change_cell_prop_crop_decrease)]
    p$pred_change_cell_prop_crop_increase[!is.na(p$pred_change_cell_prop_crop_increase)] <- 
      p.tmp$pred_change_cell_prop_crop_increase[!is.na(p$pred_change_cell_prop_crop_increase)]
    # p$tmp <- p.tmp$tmp
    # Removing p.tmp
    rm(p.tmp)
    
    # p <- subset(p, select = -c(tmp))
    
    # Make sure areas with a suitability of 0 can't have agricultural expansion, but CAN have contraction
    # Don't need these anymore. 
    # JK. We do

    p[which(p$ag_suitability == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_crop"] <- 0
    
    # Preventing expansion in cells where GDD == 0 and there is no cropland
    # Expansion can occur if GDD == 0 and there is cropland in a cell
    # Contraction can occur in these cells
    p[which(p$GDD_binary == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_crop"] <- 0
    
    # Get the area converted
    p$area_converted_crop_increase <- p$pred_change_cell_prop_crop_increase * cell_area
    p$area_converted_crop_decrease <- p$pred_change_cell_prop_crop_decrease * cell_area
    
    # Round all the numbers to avoid rounding errors
    # Will not meet country-level targets if we don't do this
    p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
          "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
          "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
          "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")] <-
      round(p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
                  "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
                  "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
                  # Re-add this line later
                  "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")],
            digits = 8)
    return(p)
  } else {
    # Pastureland -----
    # First: probability of change
    # pred_tmp are temporary predictor values. 
    # Gives data frame with probability of either no change, increase, or decrease in pasture extent
    pred_tmp <- multinom.predict.2b(dat = p,
                                    coefs = coefs.pasture.prob,
                                    vars = predictor_vars_pasture[1:12],
                                    factor_vars = predictor_vars_pasture[13])
    # Put predictions from pred_tmp into the dataframe
    p$a_pred_to_no_change_pasture <- pred_tmp$aNoChange
    p$pred_to_decrease_pasture <- pred_tmp$Decrease
    p$pred_to_increase_pasture <- pred_tmp$Increase 
    
    # Repeat for increasing models
    # Proportion to increase
    p$pred_change_cell_prop_pasture_increase <- glm.predict.2b(dat = p,
                                                               coefs = coefs.pasture.increase,
                                                               vars = predictor_vars_pasture[1:12],
                                                               factor_vars = predictor_vars_pasture[13])
    # This might cause NAs, converting these to 0s
    p$pred_change_cell_prop_pasture_increase[is.na(p$pred_change_cell_prop_pasture_increase)] <- 0
    
    # Check that the amount of pasture, crop, urban land and nonarable land is not greater 
    #     than the area of the cell.
    # Creates temporary cropland to make sure potential future area of cell after crop expansion 
    #     does not exceed 1
    # Creating temp data frame and converting NAs to 0s to avoid issues summing across columns
    p.tmp <- dplyr::select(p, 
                           urban_prop, 
                           Prop.NonArable.2010, 
                           prop_to_predict_from_crop, 
                           prop_to_predict_from_pasture,
                           pred_change_cell_prop_pasture_increase)
    # Converting NAs to 0s
    # Creating new column to adjust values of predictions so that they cannot exceed 1
    p.tmp$tmp <- rowSums(p.tmp[,c('urban_prop',
                                  'Prop.NonArable.2010',
                                  'prop_to_predict_from_crop',
                                  'prop_to_predict_from_pasture',
                                  'pred_change_cell_prop_pasture_increase')],
                         na.rm = TRUE)

    # Adjust predictions if so
    # If area does exceed 1, then adjust so that predicted pasture change does not cause cell value to exceed 1
    p.tmp[p.tmp$tmp > 1, "pred_change_cell_prop_pasture_increase"] <-
      1 - rowSums(p.tmp[p.tmp$tmp > 1,c('urban_prop',
                           'Prop.NonArable.2010',
                           'prop_to_predict_from_crop',
                           'prop_to_predict_from_pasture')],
                  na.rm = TRUE)
    
    # Repeat so that the area can't be negative
    p.tmp$tmp <- rowSums(p.tmp[,
                               c('pred_change_cell_prop_pasture_increase',
                                 'prop_to_predict_from_pasture')],
                         na.rm = TRUE)
    
    p.tmp[p.tmp$tmp < 0, "pred_change_cell_prop_pasture_increase"] <-
      0 - p.tmp[p.tmp$tmp < 0, "prop_to_predict_from_pasture"]
    
    # Finally, check that fully urbanised areas can't be converted at all
    # If urban proportion is one, then cropland extent cannot change in these cells
    p[p$urban_prop == 1 & !is.na(p$urban_prop), c("pred_to_increase_pasture",
                                                  "pred_to_decrease_pasture",
                                                  "pred_change_cell_prop_pasture_increase")] <- 0
    
    # Inserting updated estimates for prop increase into original data frame
    p$pred_change_cell_prop_pasture_increase[!is.na(p$pred_change_cell_prop_pasture_increase)] <- 
      p.tmp$pred_change_cell_prop_pasture_increase[!is.na(p$pred_change_cell_prop_pasture_increase)]
    # Removing p.tmp
    rm(p.tmp)
    
    # Make sure areas with a suitability of 0 can't have agricultural expansion, but CAN have contraction
    # Unless more than .025 portions of a cell currently have crop
    # Don't need these anymore. 
    # JK we do
    p[which(p$ag_suitability == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_pasture"] <- 0
  
    # Preventing expansion in cells where GDD == 0 and there is no cropland
    # Expansion can occur if GDD == 0 and there is more than .025 cropland in the cell
    # Contraction can occur in these cells
    p[which(p$GDD_binary == 0 & p$prop_to_predict_from_crop < .025), 
      "pred_to_increase_pasture"] <- 0
    
    # Get the area converted
    p$area_converted_pasture_increase <- p$pred_change_cell_prop_pasture_increase * cell_area
    
    # Round all the numbers to avoid rounding errors
    # Will not meet country-level targets if we don't do this
    p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
          "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
          # "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
          "pred_to_increase_pasture", "pred_change_cell_prop_pasture_increase")] <- 
    # "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")] <-
    round(p[, c("prop_adj_to_predict_from_crop", "prop_adj_to_predict_from_pasture",
                "prop_to_predict_from_crop", "prop_to_predict_from_pasture",
                #         "a_pred_to_no_change_crop", "pred_to_increase_crop", "pred_to_decrease_crop", 
                # Re-add this line later
                "pred_to_increase_pasture", "pred_change_cell_prop_pasture_increase")],
          #         "pred_change_cell_prop_crop_increase", "pred_change_cell_prop_crop_decrease")],
          digits = 8)
  }
  return(p)
}

# # Getting predictions, avoiding NA cells (unfinished) -----
# pred_no_na <- function(data, pred_crop){
#   na_data <- data[is.na(data$ISO_numeric),]
#   mod_df <- mod_df[!is.na(mod_df$ISO_numeric),]
#   mod_df <- mod_df %>%
#     group_by(IUCN_region) %>%
#     do(pred_fun(., crop = pred_crop))
#   mod_df <- bind_rows(mod_df, na_data)
#   mod_df <- mod_df[order(mod_df$row_index),]
# }
# 
# Function to convert cells to urban ----
urb_fun <- function(urban, dec){
  # Do any cells need converting?
  if(as.numeric(urban[1, paste("urb", "targ", dec, sep = "_")]) > 0){
    # Sample the cells
    s <- as.numeric(urban[1, paste("urb", "targ", dec, sep = "_")])
    # Make sure I'm not sampling more cells than are available
    if(s < nrow(urban[(1 - urban$urban_prop > 0),])){
      smple <- sample(x = nrow(urban),
                      size = as.numeric(urban[1, paste("urb", "targ", dec, sep = "_")]),
                      replace = FALSE,
                      prob = urban$urban_prob)
      # Subset out these cells
      urban_new <- urban[smple,]
      urban <- urban[-smple,]
      
      # Add the lost agricultural area to the target
      new_crop <- sum(urban_new$prop_to_predict_from_crop * cell_area, na.rm = TRUE)
      
      # Convert the cells to urban
      urban_new$urban_prop <- 1
      
      # Both the area currently pasture / crop and any predicted change should be set to zero
      urban_new$prop_to_predict_from_crop <- 0
      urban_new$prop_to_predict_from_pasture <- 0
      
      urban_new$pred_change_cell_prop_crop <- 0
      urban_new$pred_change_cell_prop_pasture_increase <- 0
      
      # Mark these cells as converted
      urban_new[,paste("conv_urb", dec, sep = "_")] <- 1
      
      # Make sure they can't be picked again
      urban_new$urban_prob <- 0
      
      # Bind data back together
      urban <- bind_rows(urban, urban_new)
      urban <- urban[order(urban$row_index),]
      urban[,paste("area_incr", dec, sep = "_")] <- urban[,paste("area_incr", dec, sep = "_")] + new_crop
      # Return the updated data frame
    } else { # If I need to sample more cells than are available
      smple <- sample(x = nrow(urban),
                      size = nrow(urban[(1 - urban$urban_prop > 0),]),
                      replace = FALSE,
                      prob = urban$urban_prob)
      # Subset out these cells
      urban_new <- urban[smple,]
      urban <- urban[-smple,]
      
      # Add the lost agricultural area to the target
      new_crop <- sum(urban_new$prop_to_predict_from_crop * cell_area, na.rm = TRUE)
      
      # Both the area currently pasture / crop and any predicted change should be set to zero
      urban_new$prop_to_predict_from_crop <- 0
      urban_new$prop_to_predict_from_pasture <- 0
      
      urban_new$pred_change_cell_prop_crop <- 0
      urban_new$pred_change_cell_prop_pasture_increase <- 0
      
      # Mark these cells as converted
      urban_new[,paste("conv_urb", dec, sep = "_")] <- 1
      
      # Make sure they can't be picked again
      urban_new$urban_prob <- 0
      
      # Bind data back together
      urban <- bind_rows(urban, urban_new)
      urban <- urban[order(urban$row_index),]
      urban[,paste("area_incr", dec, sep = "_")] <- urban[,paste("area_incr", dec, sep = "_")] + new_crop
      # Return the updated data frame
    }
    return(urban)
  } else { # If no cells need converting
    return(urban)
  }
}
# Function to pick cells to convert -----
# For troubleshooting, targets = 3 column data frame with ISO3, area of land that needs to be converted, and modulus of target. targets is a 1-row data frame
# Target is specific to 5-year decade we're forecasting
# data = mod_df. Is all of the modelling data we have, is created in another function
# target = target_decade[!is.na(target_decade$ISO3),]
crop_conv_fun <- function(target, data){
  # Get just country's cells, and those with non-NA values
  df <- filter(data, 
               ISO3 == target$ISO3[1]) %>%
    as.data.frame()
  # Exclude NAs
  df <- df %>% 
    filter(!is.na(area_converted_crop_increase),
           !is.na(area_converted_crop_decrease))
  
  # Remove cells where ag_suitability == 0 and prop_crop == 0
  # Prevents expansion in cells where ag_suitability == 0 and there is no current crop production
  df.tmp <- df[which(df$ag_suitability != 0 & df$prop_to_predict_from_crop != 0),]
  
  # Set target. Creating vector from data_frame targets
  t <- target$target_temp
  # If troubleshooting, t = target_decade$target[1]
  
  target_orig <- t # set to target_orig to have something to compare to
  
  # Setting logical value. Used in if loop below
  # If true, multiplies values by -1 
  # Does this so that below script works with both positive and negative targets
  switch = FALSE
  
  # Logical statement to only create prob and area vector if target != 0
  if(t > 0) {
    probs <- as.vector(df$pred_to_increase_crop)
    area <- as.vector(df$pred_change_cell_prop_crop_increase * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
  }
  
  # Logical to make sure target and direction of cell_area are in same direction
  if(t < 0) {
    probs <- as.vector(df$pred_to_decrease_crop)
    area <- -1 * as.vector(df$pred_change_cell_prop_crop_decrease * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
    t = t * -1
    switch = TRUE
  }
  
  # Drop cells with ag_suitability == 0 and prop_crop == 0
  
  # Drop zero probability cells and NA cells
  # Don't want them to be picked
  area <- area[probs > 0 & !is.na(probs)]
  Row.Index <- Row.Index[probs > 0 & !is.na(probs)]
  # GDD <- GDD[probs > 0]
  probs <- probs[probs > 0 & !is.na(probs)]
  
  # Also drop zero area cells
  # Creates shorter vectors and (should) make code faster
  Row.Index <- Row.Index[area!= 0 & !is.na(area)]
  probs <- probs[area!= 0 & !is.na(area)]
  # GDD <- GDD[area!= 0]
  area <- area[area!= 0 & !is.na(area)]
  
  
  # Counter to see if we're double counting cells
  n1 = 0
  # Selecting all cells if sum of cells with correct direction of change is less than target
  if(sum(area[which(area > 0)]) <= t) {
    # And adding row.nums to cells_tmp
    row.num.list <- 1:(length(area))
    # Select all cells with positive values
    cells_tmp <- row.num.list[which(area > 0)]
    # Creating list used to 
    # And removing first value of cell
    # Adjust target
    t = t - sum(area[which(area>0)])
    n1 = n1 + length(which(area > 0))
  } else {
    
    # Need abs(target_orig) for countries that have are forecast to have a contraction in land use
    #if(t == abs(target_orig) & sum(area[which(area > 0)]) != 0) {
    # Extract probabilities of picking a cell and the area to convert in each
    # Gets prob of cell being picked, and area of cell converted if it is picked
    
    # If statement to prevent land from decreasing (or increasing) if abs(mean(area)) is small
    # This would come into play if e.g. target is +100, but most cells are likely to decrease in extent
    # Randomly selecting cells in this case would result in a net decrease in cropland when we need an increase
    # Getting estimate of number of cells needed
    n.cells = round(t/mean(area))
    if(n.cells > length(Row.Index)) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
    }
    
    # n.cells can't be negative
    # Removing cells where area < 0
    # Avoids getting stuck in a while loop
    if(n.cells < 0) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
      # Updating estimate of cells needed for conversion
      n.cells = round(t/mean(area) * .95)
    }
    # Creating list of place in the Row.Index list of cells that have been selected
    # Is used so that these cells are not selected again
    row.num.list <- 1:(length(area) + 1)
    # Need to append values onto probs so that vectors is same length as row.num.list
    # Not doing this causes the fastSampleReject function to crash
    probs <- c(probs, 0)
    
    # Creating list of row indices of cells that were selected for a change in cropland extent
    # First value of cells_tmp is equivalent to last value of row.num.list
    # This value will be dropped
    cells_tmp <- row.num.list[length(row.num.list)]
    # Selecting many cells to meet target rapidly
    while(t > 3) {
      # Getting estimate of number of cells needed
      n.cells = round(t/mean(area) * .95)
      # Adding in exception if n.cells is greater than the number of remaining cells that can be picked
      # Just assuming all of the remaining cells are selected
      # This could result in an overestimate of crop conversion...
      # But this will almost always be a small overestimate
      if(n.cells > length(area[-cells_tmp])) {
        n.cells <- length(area[-cells_tmp])
      }
      
      # Selecting rows
      row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      # Cutting off cells if sum of selected cells is greater than target
      # Put in a while loop so that it loops below code until condition is satisfied
      while(sum(area[row.nums]) > t) {
        # If sum area[row.nums] > t, then reduce number of cells selected and resample until this is met
        n.cells = round(n.cells * .9)
        # Reselect cells
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      }
      
      # Appending selected rows 
      cells_tmp <- c(cells_tmp, row.nums)
      # Adjusting target
      t = t - sum(area[row.nums])
      # Adding counter to verify we're not double-selecting cells
      n1 = n1 + length(row.nums)
    }
    
    # Selecting individual cells to approach target slowly
    while(t > 0) {
      # Selecting rows
      # Need to select first object out because everything else is an NA
      # Putting this in an if statement to prevent the code from freaking out
      if(length(probs[-cells_tmp] < 1) | length(row.num.list[-cells_tmp] < 1)) {
        t = 0
      } else {
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = 1, w = probs[-cells_tmp])[1]
        # Appending selected rows 
        cells_tmp <- c(cells_tmp, row.nums)
        # Adjusting target
        t = t - area[row.nums]
        # Adding counter to verify we're not double-selecting cells
        n1 = n1 + 1
      }
    }
  }
  
  # Resetting switch
  # Also resetting values of area and t
  # Because these were multiplied by -1 above
  if(switch == TRUE) {
    t = t * -1
    area = area * -1
    switch = FALSE
  }
  # Dropping first value from cells_tmp because this is a place holder
  # But only need to do this if counter is not equal to length of cells_tmp
  if(length(cells_tmp) != n1) {
    cells_tmp <- cells_tmp[-1]
  }
  
  # Getting vector of row.indices to tell us what cells to convert
  cells_tmp <- Row.Index[cells_tmp]
  # Data frame containing all of the cells converted, how much area has been converted
  # Data frame can be quite large. A single row for every cell picked, but single value for area converted
  # Data frame can also contain no cells. Need an if statement if this is the case
  if(length(cells_tmp != 0)) {
    out <- data.frame(index = cells_tmp,
                      area_converted = target_orig - t,
                      target_remaining = t,
                      row.names = 1:length(cells_tmp))
    
  }
  # And creating a data frame if no cells were picked for conversion
  # Need to do this so that script continues running
  # Picking an arbitrarily large index for cells_tmp
  # R will return the entirevectors if the index is outside of their range
  if(length(cells_tmp) == 0) {
    out <- data.frame(index = 1e12,
                      area_converted = 0,
                      target_remaining = t,
                      row.names = 1e12)
  }
  # I have no idea why I need to rename things, but this works...
  names(out) <- c("index", "area_converted","target_remaining")
  return(out)
}

# Function to pick cells to convert for pasture (could probably combine this with the one above) ----
pasture_conv_fun <- function(target, data){
  # Filter to country's cells ----
  # Get just country's cells, and those with non-NA values
  df <- filter(data, 
               ISO3 == target$ISO3[1]) %>%
    as.data.frame()
  # Exclude NAs
  df <- df %>% 
    filter(!is.na(area_converted_pasture_increase))
  #!is.na(area_converted_pasture_decrease))
  
  # We DON'T remove cells with a suitability of 0 for pasture
  # Set target. ----
  # Creating vector from data_frame targets
  t <- target$target
  # If troubleshooting, t = target_decade$target[1]
  target_orig <- t # set to target_orig to have something to compare to
  
  # Setting logical value. Used in if loop below
  # If true, multiplies values by -1 
  # Does this so that below script works with both positive and negative target
  switch = FALSE
  
  # Logical statement to only create prob and area vector if target != 0
  if(t > 0) {
    probs <- as.vector(df$pred_to_increase_pasture)
    area <- as.vector(df$pred_change_cell_prop_pasture_increase * cell_area)
    Row.Index <- as.vector(df$row_index)
    # GDD <- as.vector(df$GDD_binary)
  }
  
  # Drop zero probability cells
  # Don't want them to be picked
  area <- area[probs > 0 & !is.na(probs)]
  Row.Index <- Row.Index[probs > 0 & !is.na(probs)]
  # GDD <- GDD[probs > 0]
  probs <- probs[probs > 0 & !is.na(probs)]
  
  # Also drop zero area cells
  # Creates shorter vectors and (should) make code faster
  Row.Index <- Row.Index[area!= 0 & !is.na(area)]
  probs <- probs[area!= 0 & !is.na(area)]
  # GDD <- GDD[area!= 0]
  area <- area[area!= 0 & !is.na(area)]
  
  # Counter to see if we're double counting cells
  n1 = 0
  # Selecting all cells if sum of cells with correct direction of change is less than target
  if(sum(area[which(area > 0)]) <= t) {
    # And adding row.nums to cells_tmp
    row.num.list <- 1:(length(area))
    # Select all cells with positive values
    cells_tmp <- row.num.list[which(area > 0)]
    # Creating list used to 
    # And removing first value of cell
    # Adjust target
    t = t - sum(area[which(area>0)])
    n1 = n1 + length(which(area > 0))
  } else {
    
    # Need abs(target_orig) for countries that have are forecast to have a contraction in land use
    #if(t == abs(target_orig) & sum(area[which(area > 0)]) != 0) {
    # Extract probabilities of picking a cell and the area to convert in each
    # Gets prob of cell being picked, and area of cell converted if it is picked
    
    # If statement to prevent land from decreasing (or increasing) if abs(mean(area)) is small
    # This would come into play if e.g. target is +100, but most cells are likely to decrease in extent
    # Randomly selecting cells in this case would result in a net decrease in cropland when we need an increase
    # Getting estimate of number of cells needed
    n.cells = round(t/mean(area))
    if(n.cells > length(Row.Index)) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
    }
    
    
    # n.cells can't be negative
    # Removing cells where area < 0
    # Avoids getting stuck in a while loop
    if(n.cells < 0 | median(area) < 0) {
      Row.Index <- Row.Index[which(area > 0)]
      probs <- probs[which(area > 0)]
      # GDD <- GDD[which(area > 0)]
      area <- area[which(area > 0)]
      # Updating estimate of cells needed for conversion
      n.cells = round(t/mean(area) * .95)
    }
    
    # Creating list of place in the Row.Index list of cells that have been selected
    # Is used so that these cells are not selected again
    row.num.list <- 1:(length(area) + 1)
    # Need to append values onto probs so that vectors is same length as row.num.list
    # Not doing this causes the fastSampleReject function to crash
    probs <- c(probs, 0)
    
    # Creating list of row indices of cells that were selected for a change in cropland extent
    # First value of cells_tmp is equivalent to last value of row.num.list
    # This value will be dropped
    cells_tmp <- row.num.list[length(row.num.list)]
    # Selecting many cells to meet target rapidly
    while(t > 3) {
      # Getting estimate of number of cells needed
      n.cells = round(t/mean(area[-cells_tmp]) * .95)
      
      # Adding in exception if n.cells is greater than the number of remaining cells that can be picked
      # Just assuming all of the remaining cells are selected
      # This could result in an overestimate of crop conversion...
      # But this will almost always be a small overestimate
      if(n.cells > length(area[-cells_tmp])) {
        n.cells <- length(area[-cells_tmp])
      }
      # Selecting rows
      row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      # Cutting off cells if sum of selected cells is greater than target
      # Put in a while loop so that it loops below code until condition is satisfied
      while(sum(area[row.nums]) > t) {
        # If sum area[row.nums] > t, then reduce number of cells selected and resample until this is met
        n.cells = round(n.cells * .95)
        # Reselect cells
        row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = n.cells, w = probs[-cells_tmp])
      }
      
      # Appending selected rows 
      cells_tmp <- c(cells_tmp, row.nums)
      # Adjusting target
      t = t - sum(area[row.nums])
      # Adding counter to verify we're not double-selecting cells
      n1 = n1 + length(row.nums)
    }
    
    # Selecting individual cells to approach target slowly
    # while(t > 0) {
    #   # Selecting rows
    #   # Need to select first object out because everything else is an NA
    #   row.nums <- fastSampleReject(all = row.num.list[-cells_tmp], n = 1, w = probs[-cells_tmp])[1]
    #   # Appending selected rows 
    #   cells_tmp <- c(cells_tmp, row.nums)
    #   # Adjusting target
    #   t = t - area[row.nums]
    #   # Adding counter to verify we're not double-selecting cells
    #   n1 = n1 + 1
    # }
  }
  
  # Resetting switch
  # Also resetting values of area and t
  # Because these were multiplied by -1 above
  if(switch == TRUE) {
    t = t * -1
    area = area * -1
    switch = FALSE
  }
  # Dropping first value from cells_tmp because this is a place holder
  # But only need to do this if counter is not equal to length of cells_tmp
  if(length(cells_tmp) != n1) {
    cells_tmp <- cells_tmp[-1]
  }
  
  # Getting vector of row.indices to tell us what cells to convert
  cells_tmp <- Row.Index[cells_tmp]
  # Data frame containing all of the cells converted, how much area has been converted
  # Data frame can be quite large. A single row for every cell picked, but single value for area converted
  # Data frame can also contain no cells. Need an if statement if this is the case
  if(length(cells_tmp != 0)) {
    out <- data.frame(index = cells_tmp,
                      area_converted = target_orig - t,
                      target_remaining = t,
                      row.names = 1:length(cells_tmp))
    
  }
  # And creating a data frame if no cells were picked for conversion
  # Need to do this so that script continues running
  # Picking an arbitrarily large index for cells_tmp
  # R will return the entirevectors if the index is outside of their range
  if(length(cells_tmp) == 0) {
    out <- data.frame(index = 1e12,
                      area_converted = 0,
                      target_remaining = t,
                      row.names = 1e12)
  }
  # I have no idea why I need to rename things, but this works...
  names(out) <- c("index", "area_converted","target_remaining")
  return(out)
}

# Function to run through each year, convert cells and update predictions ----
# If troubleshooting x is a subset of model_runs, select model_runs where run == 1
# x <- subset(model_runs, run == 1)
# If running multiple times, probably best to clear it up first
# rm(mod_df, x, decs, dec, i, na.data, target_decade, messages.list,
#    n, y, target_remaining, counter, met, cells, cells_contraction,
#    cells_expansion, out.list)
year_fun <- function(x){
  mod_df <- model_data_current
  # mod_df <- left_join(mod_df, x)
  # # Dropping urban targets to save space
  # mod_df <- dplyr::select(mod_df, 
  #                         -c(run, 
  #                            urb_targ_2010_2015,
  #                            urb_targ_2015_2020,
  #                            urb_targ_2020_2025,
  #                            urb_targ_2025_2030))
  
  # Really 5-yr time intervals. Decs because originally forecasting at 10-yr intervals
  decs <- gsub("target_", 
               "",
               names(mod_df)[grep("target_", names(mod_df))])
  
  # Set up list for messages ----
  # Messages because warnings
  messages.list <- list()
  # Loop through decades -----
  # for(i in 1:1){
  for(i in 1:length(decs)){
    # Set up the decade to use
    dec <- decs[i]
    # Add columns ----
    # Add columns for converted to cropland and pasture or not
    # These chunks of code add columns to the data set that indicate how much land has been converted into cropland
    # mod_df[[paste("conv",
    #               dec,
    #               sep = "_")]] <- 0
    # mod_df[is.na(mod_df$ISO_numeric), 
    #        paste("conv",
    #              dec,
    #              sep = "_")] <- NA
    # mod_df[[paste("conv_pasture",
    #               dec,
    #               sep = "_")]] <- 0
    # mod_df[is.na(mod_df$ISO_numeric), 
    #        paste("conv_pasture",
    #              dec,
    #              sep = "_")] <- NA
    
    # Set demand to correct decade
    mod_df$change_demand_crop <- mod_df[, paste("prop_change", dec, sep = "_")] + 1
    
    # Get targets -----
    target_decade <- as.data.frame(mod_df[!duplicated(mod_df$ISO3),
                                          c("ISO3", paste("target", dec, sep = "_"))])
    names(target_decade)[2] <- "target"
    target_decade <- target_decade[!is.na(target_decade$ISO3),]
    
    # And manually updating target_decade$target if needed
    # For some reasons this can result in an NA for some, but not all, countries
    if(is.na(target_decade$target)) {
      # First getting value for appropriate decade
      tmp.target.dec <- targets[targets$ISO3 == target_decade$ISO3,paste('target_',
                                                                         dec,
                                                                         sep = '')]
      # And second dropping NA value
      tmp.target.dec <- tmp.target.dec[!is.na(tmp.target.dec)]
      
      target_decade$target <- tmp.target.dec
    }
    
    
    # Put warnings into a list
    messages.list[[dec]] <- list()
    
    # Adding remaining target from previous decade to this decade
    if(i != 1) {
      target_decade$target <- target_decade$target + target_remaining
    }
    
    # Get the absolute value of the targets, so I can run a while() loop
    # I don't think is needed anymore
    target_decade$abs.target <- abs(target_decade[,"target"])
    
    # Split target into parts -----
    # Looping through the below while loop n amount of times
    # Each loop hits 1/n part of the total target fo rthe time interval
    # This makes it more likely to for cropland expansion to be concentrated
    # Setting n
    n = 4
    # Logical to fix n if it is less than one
    if(n < 1) {
      n <- 1
    }
    # Setting split target (previously "target_temp" - there is a reason for this...
    # Now doing this outside of the for() loop
    target_decade$split_target <- target_decade$target / n
    # Reset remaining target
    target_remaining <- 0
    
    # Re-run predictions for the decade -----
    # NB. We are doing this only once at the start of each decade, rather than in every
    #     loop of the while() loop, or y-loop.
    # This is because a) It reflects the logic of our model: we are predicting how things
    #     change in a five-year period, based on the situation at the start of that period,
    #     not how it changes throughout the period.
    #     b) Because it speeds things up a lot...
    # First get the adjacent agricultural proportions
    mod_df <- adj.ag.f(m = mod_df,
                       r = rows,
                       c = cols,
                       # f = adj_cell_mat,
                       crop = "crop")  %>%
      as.data.frame()
    mod_df <- adj.ag.f(m = mod_df,
                       r = rows,
                       c = cols,
                       # f = adj_cell_mat,
                       crop = "pasture")  %>%
      as.data.frame()
    
    # Then get new predictions
    na.data <- mod_df[is.na(mod_df$ISO3),]
    mod_df <- mod_df[!is.na(mod_df$ISO3),]
    mod_df <- pred_fun(mod_df, crop = "crop")
    mod_df <- pred_fun(mod_df, crop = "pasture")
    mod_df <- bind_rows(mod_df, na.data)
    mod_df <- mod_df[order(mod_df$row_index),]
    
    # Run through split target ----
    for(y in 1:n) {
      # op <- options(digits.secs = 6)
      # print(c(dec, y, format(Sys.time(), "%a %b %d %X %Y"))) # This will probably slow things down a bit...
      # Adding counter to loop through multiple rounds of conversion
      counter <- 1
      # If statement to avoid looping over countries that do not have any cropland expansion
      if(target_decade$abs.target < .1) {
        # Set counter to 19 to avoid while loop below
        counter <- 19
        # And update target_remaining so that it does not catch on above logical statement
        # Not doing this results in an error
        target_remaining <- target_decade$target
        
        # And updating values for extent conversion and absolute crop for countries that have no expansion
        # mod_df[,paste("conv", dec, sep = "_")] <- 0
        mod_df[[paste("prop_crop", dec, sep = "_")]] <- mod_df$prop_to_predict_from_crop
      }  
      # Get a temporary target: the portion of the overall target that we're trying to meet 
      #     (at the moment 1/4) plus any remaining target from the previous loop
      target_decade$target_temp <- target_decade$split_target + target_remaining
      
      # if(i != 1) {
      #   # Taking away the amount of land that has already been converted
      #   target_decade$target <- target_decade$target - target_decade$converted
      # }
      
      # Pick cells for crop expansion ----
      while(abs(target_decade$target_temp) > 1 & counter <= 15){
        # Make sure that crop predictions can't be too high or low -----
        # This is lifted from pred_fun()
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                      'Prop.NonArable.2010',
                                      'prop_to_predict_from_crop',
                                      'pred_change_cell_prop_crop_increase')],
                             na.rm = TRUE)

        # Adjust predictions if so
        # If area does exceed 1, then adjust so that predicted cropland change does not cause cell value to exceed 1
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_crop_increase"] <-
          1 - rowSums(mod_df[mod_df$tmp>1,c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop')],
                      na.rm = TRUE)

        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_crop_increase',
                                        'prop_to_predict_from_crop')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_crop_increase"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_crop"]
        
        # Repeat for decreasing models
        # Check that the amount of cropland plus change plus urban land plus non arable land is not greater than the area of the cell.
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                        'Prop.NonArable.2010',
                                        'prop_to_predict_from_crop',
                                        'pred_change_cell_prop_crop_decrease')],
                              na.rm = TRUE)
        
        # Adjust predictions if so
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_crop_decrease"] <-
          1 - rowSums(mod_df[mod_df$tmp>1, c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop')],
                      na.rm = TRUE)
        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_crop_decrease',
                                        'prop_to_predict_from_crop')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_crop_decrease"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_crop"]
        
        # Pick cells ----
        # This needs to be split by expansion or contraction. 
        # Some regions won't have any countries that expand / contract, so put in 
        #     little if loops
        # cells_expansion is output of crop_conv_fun
        # Then repeats for countries that have cropland contraction
        if(nrow(target_decade[target_decade$target_temp > 0,]) > 0){
          cells_expansion <- crop_conv_fun(target = target_decade, 
                                           data = mod_df)
        } else {
          cells_expansion <- data.frame(index = NULL,
                                        area_converted = NULL)
        }
        if(nrow(target_decade[target_decade$target_temp < 0,]) > 0){
          cells_contraction <- crop_conv_fun(target = target_decade, 
                                             data = mod_df) 
        } else {
          cells_contraction <- data.frame(index = NULL,
                                          area_converted = NULL)
        }
        
        # Convert chosen cells ----
        # Expansion
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_crop"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_crop"] + 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_crop_increase"]
        # Contraction
        mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_crop"] <- 
          mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_crop"] + 
          mod_df[mod_df$row_index %in% cells_contraction$index, "pred_change_cell_prop_crop_decrease"]
        
        # Update the change in the previous time-step
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_cell_delta_crop"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_crop_increase"]
        # Contraction
        mod_df[mod_df$row_index %in% cells_contraction$index, "prop_to_predict_from_cell_delta_crop"] <- 
          mod_df[mod_df$row_index %in% cells_contraction$index, "pred_change_cell_prop_crop_decrease"]
        
        # Note the proportion in each cell in each year
        # mod_df[mod_df$row_index %in% cells_expansion$index,
        #        paste("conv", dec, sep = "_")] <-  
        #   mod_df[mod_df$row_index %in% cells_expansion$index, paste("conv", dec, sep = "_")] + 
        #   mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_crop_increase"]
        # mod_df[mod_df$row_index %in% cells_contraction$index,
        #        paste("conv", dec, sep = "_")] <-  
        #   mod_df[mod_df$row_index %in% cells_contraction$index, paste("conv", dec, sep = "_")] + 
        #   mod_df[mod_df$row_index %in% cells_contraction$index, "pred_change_cell_prop_crop_decrease"]
        
        # Getting remaining target. e.g. land that could not be converted ----
        # This will be added onto target for future decades
        if(target_decade$target > 0) {
          target_remaining = unique(cells_expansion$target_remaining)
        }
        if(target_decade$target < 0) {
          target_remaining = unique(cells_contraction$target_remaining)
        }
        if(target_decade$target == 0) {
          target_remaining = 0
        }
        
        # Adjust targets -----
        cells <- bind_rows(cells_contraction, 
                           cells_expansion)
        target_decade$converted <- cells$area_converted[1]
        # Take the area converted off the target
        # In an if loop because target_decade$converted is undefined if target_decade$target == 0
        if(target_decade$target != 0) {
          target_decade$target_temp <- target_decade$target_temp - target_decade$converted
          # Recalculate absolute target
          # target_decade$abs.target <- target_decade$abs.target - abs(target_decade$converted)
        }
        
        # }
        # Adding in logical to prevent script from looping through while loop if there are no cells left to expand
        # Also prevents from looping if target has been hit
        if(nrow(cells) == 0) {
          counter <- 16
        } else {if(abs(target_decade$split_target - target_decade$converted) < 1) {
          counter <- 17
        } else {
          counter <- counter + 1
        }
        }
        # Below bracket marks end of while loop that meets 1/n proportion of time interval's crop target
      }
      
      # Pasture target setting -----
      # # Work out how much pasture has been lost:
      # # Area available for pasture is the area without urban or cropland
      mod_df$prop_for_pasture <- 0
      # Getting row index of cells that changed in crop extent
      # tmp_index <- which(mod_df$row_index %in% cells$index)
        
      # Calculating amount for pasture
      mod_df$prop_for_pasture <- 1 - rowSums(mod_df[,c('urban_prop',
                                                                  'prop_to_predict_from_crop',
                                                                  'Prop.NonArable.2010')],
                                                        na.rm = TRUE)
      
      # Calculating amount of pasture lost
      mod_df$pasture_loss_tmp <- 0
      mod_df$pasture_loss_tmp <- mod_df$prop_to_predict_from_pasture -
        mod_df$prop_for_pasture
      
      # Dropping value
      mod_df <- dplyr::select(mod_df,
                              -prop_for_pasture)
      
      # Pasture is not necessarily gained in a cell if cropland extent decreases
      # Negative values correspond with pasture area being gained by a cell
      # So Limiting range of pasture loss values from 0 to 1
      # Negative values in pasture_loss_tmp would indicate a net increase in pasture extent in the cell
      mod_df$pasture_loss_tmp <- squish(mod_df$pasture_loss_tmp, 
                                                   range = c(0,1))
      
      # Adjust the amount of pastureland so it is equal to that available
      mod_df$prop_to_predict_from_pasture <-
        mod_df$prop_to_predict_from_pasture - mod_df$pasture_loss_tmp
      # 
      
      # Round this value. NOTE: I am rounding this to one fewer digit (7 vs 8). This is 
      #     because taking two 8-digit decimals from each other and rounding it to 8 d.p. 
      #     can still give a floating point error
      mod_df$prop_to_predict_from_pasture <- round(mod_df$prop_to_predict_from_pasture, 
                                                   digits = 7)
      
      # # Sum up the total area lost for each country
      pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                   target = sum(mod_df$pasture_loss_tmp * cell_area,
                                                na.rm = TRUE)) %>%
        mutate(abs.target = abs(target))
      
      # Manually updating values of pasture_target in the case that this data frame is not created
      if(nrow(pasture_target) == 0) {
        pasture_target <- data.frame(ISO3 = target_decade$ISO3[1],
                                     target = 0,
                                     abs.target = 0)
      }
      
      # If statement to avoid looping over countries that do not have any cropland expansion
      if(pasture_target$abs.target < .1) {
        # Setting counter to 19
        counter <- 19
        # and updating pasture extent in next time period if no cells are changed
        # And updating values for extent conversion and absolute crop for countries that have no expansion
        # mod_df[[paste("conv_pasture", dec, sep = "_")]] <- 0
        mod_df[[paste("prop_pasture", dec, sep = "_")]] <- mod_df$prop_to_predict_from_pasture
      } else {counter <- 1}
      
      # Pick cells for pasture expansion ----
      while(abs(pasture_target$target) > 1 & counter <= 15){
        # Make sure that pasture predictions can't be too high or low ----
        mod_df$tmp <- rowSums(mod_df[,c('urban_prop',
                                        'Prop.NonArable.2010',
                                        'prop_to_predict_from_crop',
                                        'prop_to_predict_from_pasture',
                                        'pred_change_cell_prop_pasture_increase')],
                              na.rm = TRUE)
        
        # Adjust predictions if so
        # If area does exceed 1, then adjust so that predicted pasture change does not cause cell value to exceed 1
        mod_df[mod_df$tmp > 1, "pred_change_cell_prop_pasture_increase"] <-
          1 - rowSums(mod_df[mod_df$tmp > 1,c('urban_prop',
                                'Prop.NonArable.2010',
                                'prop_to_predict_from_crop',
                                'prop_to_predict_from_pasture')],
                      na.rm = TRUE)
        
        # Repeat so that the area can't be negative
        mod_df$tmp <- rowSums(mod_df[,c('pred_change_cell_prop_pasture_increase',
                                        'prop_to_predict_from_pasture')],
                              na.rm = TRUE)
        
        mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "pred_change_cell_prop_pasture_increase"] <-
          0 - mod_df[mod_df$tmp < 0 & !is.na(mod_df$tmp), "prop_to_predict_from_pasture"]
        
        mod_df <- dplyr::select(mod_df, -tmp)
        
        # Pick cells ----
        cells_expansion <- pasture_conv_fun(target = pasture_target, 
                                            data = mod_df)
        # I'm (naively?) assuming that pasture targets will be met. If they aren't then they aren't, and that's
        #     tough
        # Convert chosen cells ----
        # Expansion
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_pasture"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_pasture"] + 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_pasture_increase"]
        
        # Update the amount changed in the previous time step    
        mod_df[mod_df$row_index %in% cells_expansion$index, "prop_to_predict_from_cell_delta_pasture"] <- 
          mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_pasture_increase"]
        
        # Note the proportion in each cell in each year
        # mod_df[mod_df$row_index %in% cells_expansion$index,
        #        paste("conv_pasture", dec, sep = "_")] <-  
        #   mod_df[mod_df$row_index %in% cells_expansion$index, paste("conv_pasture", dec, sep = "_")] + 
        #   mod_df[mod_df$row_index %in% cells_expansion$index, "pred_change_cell_prop_pasture_increase"]
        
        # Adjust targets -----
        pasture_target$converted <- cells_expansion$area_converted[1]
        # Take the area converted off the target
        # In an if loop because pasture_target$converted is undefined if pasture_target$target == 0
        if(pasture_target$target != 0) {
          pasture_target$target <- pasture_target$target - pasture_target$converted
          # Recalculate absolute target
          # pasture_target$abs.target <- pasture_target$abs.target - abs(pasture_target$converted)
        }
        
        # Subset the targets to only those that remain unmet
        met <- pasture_target[pasture_target$abs.target <= 0,]
        # pasture_target <- pasture_target[pasture_target$abs.target > 0,]
        
        
        # Adding in logical to prevent script from looping through while loop if there are no cells left to expand
        if(nrow(cells) == 0) {
          counter <- 16
        } else {if(abs(pasture_target$target) < 1) {
          counter <- 17
        } else {
          counter <- counter + 1
        }
        }
        # Below bracket marks end of while loop for pasture conversion
      }
      
      # Need to reset the pasture target for each of mini-loops within each time interval
      mod_df$pasture_loss_tmp <- 0
      
      # Small logical that creates target_decade$converted if this column does not already exist
      if(!("converted" %in% colnames(target_decade))) {
        target_decade$converted <- 0
      }
      
      # Also need to adjust the absolute target into make sure that the messages reflect
      #     whether the target is met
      target_decade$abs.target <- target_decade$abs.target - abs(target_decade$converted)
      
      # Below bracket marks end for loop that crop and pasture conversion loops are embedded in
      # The for loop divides up the time interval target into n parts of equal size (assuming each part's target is met in full)
    }
    
    # Print the countries where targets are met
    if(target_decade$abs.target <= 0) {
      messages.list[[dec]][["Messages.Crop"]] <- c("Target Met", NA)
    }
    # Print countries where targets are not met 
    # NEED TO IMPROVE THIS - THE COUNTER COULD BE MET, INSTEAD OF RUNNING OUT OF LAND
    if(target_decade$abs.target > 0) {
      messages.list[[dec]][["Messages.Crop"]] <- c("Target Not Met - Land Remaining", target_decade$abs.target)
    }
    
    # Print the countries where targets are met for pasture
    if(abs(pasture_target$target) < 1) {
      messages.list[[dec]][["Messages.Pasture"]] <- "Pasture target Met"
    }
    # Print countries where targets are not met for pasture
    if(abs(pasture_target$target) >= 1) {
      messages.list[[dec]][["Messages.Pasture"]] <- "Pasture Target Not Met - Not Enough Available Land"
    }
    
    # Get the final area of cropland and pastureland in each cell at the end of the 
    #     time period
    # If statement. These columns are created above if target == 0
    if(!(paste("prop_crop", dec, sep = "_") %in% colnames(mod_df))) {
      mod_df[[paste("prop_crop", dec, sep = "_")]] <- mod_df$prop_to_predict_from_crop
      mod_df[[paste("prop_pasture", dec, sep = "_")]] <- mod_df$prop_to_predict_from_pasture
    }
    # The bracket below marks the end of the for() loop going through each time period
    
    
    # Dropping columns in mod_df to clear space
    # Getting index of columns
    col.identifiers <- c(grep(paste("prop_change_",
                                  dec,
                                  sep = ""),
                            names(mod_df)),
                         grep(paste("target_",
                                    dec,
                                    sep = ""),
                              names(mod_df)))

    mod_df <- mod_df[,-col.identifiers]

  }
  
  # If we only want to get the meaningful cells:
  mod_df <- mod_df[!is.na(mod_df$ISO3),]
  out.list <- list()
  out.list[["LandUse"]] <- mod_df
  out.list[["messages"]] <- messages.list
  return(out.list)
}
# out.list is a dataframe with all the messages created
# Above function returns a list of lists. One item for each model run, and each model run has warning messages and data frames


# This is a script for making rasters
# Will be embedded in in a function
# This function returns one list
# First list is a list of crop extent rasters
# List contains both sd and mean by decade
# Second list is a list of pasture extent rasters
# List contains both sd and mean by decade


# For troubleshooting
# column.names = c('prop_crop_2055_2060')
# x = test.list
mean.sd.fun <- function(column.names) {
  
  # Puting row index into first column in matrix
  plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
  
  # Dropping any random columns that may have ended up in the data frame
  column.names <- column.names[grep('_20', column.names)]
  
  for (i in 1:length(column.names))  {
    # Creating empty matrix
    # Used to hold crop or pasture extent values
    matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = length(test.list))
    # Updating values in matrix.holder
    for(rep in 1:length(test.list)){
      matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
    }
      
    # And taking row means and row Sds of matrix.holder
    plot.table[[paste(column.names[i],
                      "mean",
                      sep = "_")]] <-  rowMeans(matrix.holder,
                                               na.rm = TRUE)
    
    # And putting sds into a column
    plot.table[[paste(column.names[i],
                      "sd",
                      sep = "_")]] <- rowSds(matrix.holder,
                                             na.rm = TRUE)
  }
  # Alright, now have a matrix that contains row sds and row means
  # Next need to join this matrix with a matrix that contains row index values for the entire country
  # This includes water area and ocean area
  row.index.table <- data.table(Row.Index = cbind(c(getValues(row.index.raster))))
  # Changing row names. Don't know why I need to do this
  names(row.index.table)[1] <- "Row.Index"
  # And now joining two matrices together
  plot.table <- left_join(row.index.table,
                          plot.table)
  
  # And ordering by row.index
  plot.table <- setorder(plot.table, 
                         Row.Index)
  raster.list <- list()
  # And now making a list of rasters
  for(i in 2:ncol(plot.table)) {
    raster.list[[names(plot.table)[i]]] <- raster(matrix(plot.table[,i],
                                                         nrow = row.index.raster@nrows,
                                                         byrow = FALSE))
    crs(raster.list[[(i-1)]]) <- crs(row.index.raster)
    extent(raster.list[[(i-1)]]) <- extent(row.index.raster)
  }
  
  return(raster.list)
}


# Functions used to predict from coefficients-----
# Second function part two
# Does not require that variable inputs are in the same order as the model inputs
# Can deal with any number of inputs
# Assumes coefs[1] is the coefficient for the intercept
# Works with factors, but only if data frames that contain a single factor are fed into the function
glm.predict.2b <- function(dat,coefs, vars, factor_vars) {
  # Converting coefs into a matrix
  # First saving row names
  # Will also be used later
  row.names.coefs = row.names(coefs)
  # Second converting into a matrix
  coefs = t(as.matrix(coefs))
  # And third adding column names
  colnames(coefs) <- row.names.coefs
  
  # Ordering coefs data frame to make sure it matches order of the predictor vars
  # There's got to be a simple way of doing this, but I don't know this off the top of my head!
  coefs.tmp <- rep(NA, length(vars) + 1)
  coefs.tmp[1] <- coefs[colnames(coefs) == '(Intercept)']
  # And looping through to reorder coefficients
  for(i in 1:length(vars)) {
    coefs.tmp[i+1] <- coefs[colnames(coefs) == vars[i]]
  }
  
  # Adding in factor coefficients
  # First getting any coefficient for a factor variable
  coefs_factor <- matrix(coefs[grep(factor_vars, colnames(coefs))],
                         nrow = 1)
  # And updating col names
  row.names.coefs <- row.names.coefs[grep(factor_vars, row.names.coefs)]
  colnames(coefs_factor) <- row.names.coefs
  
  # Second getting specific coefficient for a factor variable
  # Doing this in two steps
  # First dropping name of the factor variable from the column name
  colnames(coefs_factor) <- gsub(".*)","",colnames(coefs_factor))
  
  # Saving names of factors for later
  # tmp.factors <- names(coefs_factor)
  
  # And saving location of factor in above string for later
  # tmp.location <- which(tmp.factors == unique(dat[,factor_vars]))
  # Getting iso3 of country
  tmp.iso = unique(dat[,factor_vars])
  tmp.iso = tmp.iso[!is.na(tmp.iso)]
  
  # And second only getting the estimate for the specific factor variable
  coefs_factor <- coefs_factor[which(colnames(coefs_factor) == tmp.iso)]
  
  
  # If categorical outcome we're looking at is the baseline, then intercept is not adjusted
  if(length(coefs_factor) == 0) {
    coefs_factor = 0
  }
  
  # And doing maths!
  return(exp(as.matrix(dat[,vars]) %*% matrix(coefs[2:length(coefs.tmp)]) +
               coefs.tmp[1] + coefs_factor))
}


multinom.predict.2b <- function(dat, coefs, vars, factor_vars) {
  # Ordering coefs data frame to make sure it matches order of the predictor vars
  # There's got to be a simple way of doing this, but I don't know this off the top of my head!
  coefs.tmp <- coefs
  coefs.tmp[,2:ncol(coefs.tmp)] <- NA
  # Dropping coefficients for factor variables
  # Will be added in later
  coefs.tmp <- coefs.tmp[,1:(length(vars) + 1)]
  # And looping through to reorder coefficients
  for(i in 2:ncol(coefs.tmp)) {
    coefs.tmp[,i] <- coefs[,colnames(coefs) == vars[i-1]]
  }
  
  
  # Adding in factor coefficients
  # First getting any coefficient for a factor variable
  coefs_factor <- coefs[,grep(factor_vars, colnames(coefs))]
  
  # Changing names of coefs_factor
  colnames(coefs_factor) <- gsub(".*)","",colnames(coefs_factor))
  
  # Getting iso numeric value
  tmp.iso = unique(dat[,factor_vars])
  tmp.iso = tmp.iso[!is.na(tmp.iso)]
  
  # And getting coefficients for the specific factor being fed to the data frame
  coefs_factor <- coefs_factor[,which(colnames(coefs_factor) == tmp.iso)]
  
  # If cateogrical variable is the baseline categorical variable, then need to create a matrix with 0s
  if(length(coefs_factor) == 0) {
    coefs_factor <- matrix(c(0,0), ncol = 1)
    row.names(coefs_factor) <- row.names(coefs)
  }
  
  # creating empty list 
  # Used to store model estimates
  vectors.list <- list()
  
  # Creating names to store in the list
  cat.names <- c('aNoChange',
                 row.names(coefs))
  
  # First category has coefficients of 0
  # Creating empty vector with these values
  vectors.list[[cat.names[1]]] <- rep(exp(0), nrow(dat))
  
  # And looping through other categories to get variable estimates
  for(n.cats in 1:nrow(coefs)) {
    # Doing matrix multiplication
    # And then exponentiating
    # To get vector of relative model probabilities
    vectors.list[[cat.names[n.cats + 1]]] <- exp(as.matrix(dat[,vars]) %*% matrix(coefs[n.cats,2:ncol(coefs.tmp)]) + 
                                                   coefs.tmp[n.cats, 1] + coefs_factor[n.cats])
  }
  
  # Creating matrix with relative probabilities
  prob.matrix <- matrix(c(vectors.list[[1]],
                          vectors.list[[2]],
                          vectors.list[[3]]),
                        ncol = (nrow(coefs) + 1),
                        byrow = FALSE)
  
  # And then getting relative probabilities
  prob.matrix <- prob.matrix / rowSums(prob.matrix)
  # And adding column names
  colnames(prob.matrix) <- cat.names
  # And returing prob.matrix
  return(as.data.frame(prob.matrix))
}



# Function for getting standard deviation and standard deviation by number of reps
reps.sd.se.function <- function(column.names) {
  
  for(n.its in 2:length(test.list)) {
    # Puting row index into first column in matrix
    plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
    
    # Dropping any random columns that may have ended up in the data frame
    column.names <- column.names[grep('_20', column.names)]
    
    for (i in 1:length(column.names))  {
      # Creating empty matrix
      # Used to hold crop or pasture extent values
      matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = n.its)
      # Updating values in matrix.holder
      for(rep in 1:n.its){
        matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
      }
      
      # # And taking row means and row Sds of matrix.holder
      # plot.table[[paste(column.names[i],
      #                   "mean",
      #                   sep = "_")]] <-  rowMeans(matrix.holder,
      #                                             na.rm = TRUE)
      
      # And putting sds into a column
      plot.table[[paste(column.names[i],
                        "sd",
                        sep = "_")]] <- rowSds(matrix.holder,
                                               na.rm = TRUE)
    }
    # Alright, now have a matrix that contains row sds and row means
    # Next need to get median SD of this matrix
    # But only for cells that don't have a rowsums of 0 for the sds across all years
    # Because this indicates these cells didn't change in crop extent
    # Taking rowsums
    rowsums.tmp <- rowSums(plot.table[,2:ncol(plot.table)], na.rm = TRUE)
    # Setting cell values to NA if cell never experiences a change in crop extent
    plot.table[which(rowsums.tmp == 0),2:ncol(plot.table)] <- NA
    # Getting median row value
    plot.table.median <- colMedians(as.matrix(plot.table[,2:ncol(plot.table)]), 
                                    na.rm = TRUE)
    # Turning into a matrix with row names
    plot.table.median <- as.data.frame(t(as.matrix(plot.table.median)))
    names(plot.table.median) <- column.names
    plot.table.median$N.Iterations <- n.its
    
    if(n.its == 2) {
      median.sd.table <- plot.table.median
    }
    if(n.its > 2) {
      median.sd.table <- rbind(plot.table.median,
                               median.sd.table)
    }
    
  }
  
  # and outputting this table
  return(median.sd.table)
  
}


reps.mean.function <- function(column.names) {
  
  for(n.its in 2:length(test.list)) {
    # Puting row index into first column in matrix
    plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
    
    # Dropping any random columns that may have ended up in the data frame
    column.names <- column.names[grep('_20', column.names)]
    
    for (i in 1:length(column.names))  {
      # Creating empty matrix
      # Used to hold crop or pasture extent values
      matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = n.its)
      # Updating values in matrix.holder
      for(rep in 1:n.its){
        matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
      }
      
      # # And taking row means and row Sds of matrix.holder
      # plot.table[[paste(column.names[i],
      #                   "mean",
      #                   sep = "_")]] <-  rowMeans(matrix.holder,
      #                                             na.rm = TRUE)
      
      # And putting sds into a column
      plot.table[[paste(column.names[i],
                        "sd",
                        sep = "_")]] <- rowSds(matrix.holder,
                                               na.rm = TRUE)
    }
    # Alright, now have a matrix that contains row sds and row means
    # Next need to get median SD of this matrix
    # But only for cells that don't have a rowsums of 0 for the sds across all years
    # Because this indicates these cells didn't change in crop extent
    # Taking rowsums
    rowsums.tmp <- rowSums(plot.table[,2:ncol(plot.table)], na.rm = TRUE)
    # Setting cell values to NA if cell never experiences a change in crop extent
    plot.table[which(rowsums.tmp == 0),2:ncol(plot.table)] <- NA
    # Getting median row value
    plot.table.median <- colMeans(as.matrix(plot.table[,2:ncol(plot.table)]), 
                                    na.rm = TRUE)
    # Turning into a matrix with row names
    plot.table.median <- as.data.frame(t(as.matrix(plot.table.median)))
    names(plot.table.median) <- column.names
    plot.table.median$N.Iterations <- n.its
    
    if(n.its == 2) {
      median.sd.table <- plot.table.median
    }
    if(n.its > 2) {
      median.sd.table <- rbind(plot.table.median,
                               median.sd.table)
    }
    
  }
  
  # and outputting this table
  return(median.sd.table)
  
}


bootstrap.sd.function <- function(column.names) {
  # First looking at sd across all iterations
  # Doing this to limit sample to only cells that experience a change in crop extent in any iteration
  # Don't want to artifically deflate the SD that we're sampling
  # Puting row index into first column in matrix
  plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
  # Looping to get sds across all iterations 
  # for (i in 1:length(column.names))  {
  #   # Creating empty matrix
  #   # Used to hold crop or pasture extent values
  #   matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = length(test.list))
  #   # Updating values in matrix.holder
  #   for(rep in 1:length(test.list)){
  #     matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
  #   }
  #   
  #   # # And taking row means and row Sds of matrix.holder
  #   # plot.table[[paste(column.names[i],
  #   #                   "mean",
  #   #                   sep = "_")]] <-  rowMeans(matrix.holder,
  #   #                                             na.rm = TRUE)
  #   
  #   # And putting sds into a column
  #   plot.table[[paste(column.names[i],
  #                     "sd",
  #                     sep = "_")]] <- rowSds(matrix.holder,
  #                                            na.rm = TRUE)
  # }
  # # Getting row means of the sds
  # # If row mean = 0, then cell didn't experience any change in crop extent in any model iteration
  # tmp.sums <- rowSums(plot.table[,2:ncol(plot.table)])
  # # And getting index to remove these rows from the data table
  # tmp.index <- which(tmp.sums == 0)
  
  # Looping through iterations to make a list of tables
  iteration.list <- list()
  for(n.its in 2:length(test.list)) {
    # Puting row index into first column in matrix
    plot.table <- data.table(Row.Index = test.list[[1]]$LandUse$row_index)
    
    # Dropping any random columns that may have ended up in the data frame
    column.names <- column.names[grep('_20', column.names)]
    for (i in 1:length(column.names))  {
      # Creating empty matrix
      # Used to hold crop or pasture extent values
      matrix.holder <- matrix(NA, nrow = nrow(plot.table), ncol = n.its)
      # Updating values in matrix.holder
      for(rep in 1:n.its){
        matrix.holder[,rep] <- test.list[[rep]]$LandUse[,column.names[i]]
      }
      
      # # And taking row means and row Sds of matrix.holder
      # plot.table[[paste(column.names[i],
      #                   "mean",
      #                   sep = "_")]] <-  rowMeans(matrix.holder,
      #                                             na.rm = TRUE)
      
      # And putting sds into a column
      plot.table[[paste(column.names[i],
                        "sd",
                        sep = "_")]] <- rowSds(matrix.holder,
                                               na.rm = TRUE)
    }
    # Adding plot table to list of tables
    iteration.list[[paste("Iteration_",
                          n.its,
                          sep = "")]] <- plot.table
  }
  
  
  
  sample.vector <- 1:nrow(plot.table)
  sample.size = round(length(sample.vector)/100)
  
  # Looping through bootstrap
  for(bootstrap in 1:100) {
    # Getting subset of cells to sample
    tmp.cells <- fastSampleReject(all = sample.vector,
                                  n = sample.size,
                                  w = rep(1, length(sample.vector)))
    for(z in 1:length(iteration.list)) {
      # Taking subset of cells
      bootstrap.frame <- iteration.list[[z]][tmp.cells,]
      # Adding indicators for the bootstrap and number of iterations
      bootstrap.frame$N.Iterations <- z
      bootstrap.frame$Bootstrap <- bootstrap
      if(z == 1) {
        master.frame <- bootstrap.frame
      } else {
        master.frame <- rbind(master.frame,
                              bootstrap.frame)
      }
    }
    
    if(bootstrap == 1) {
      complete.frame <- master.frame
    }
    if(bootstrap > 1) {
      complete.frame <- rbind(complete.frame,
                              master.frame)
    }
    # Below bracket is end of bootstrap loop
  }
  return(complete.frame)
}
# Bits cut from year_fun() ----
# Bits cut: Urban sections ----
# From underneath:
# Set demand to correct decade
# mod_df$change_demand_crop <- mod_df[, paste("prop_change", dec, sep = "_")] + 1

# # Bits cut: Urban expansion - not currently included ----
# # If it's a decade with urban expansion then convert cells to urban
# if(!is.null(mod_df[[paste("urb", "targ", dec, sep = "_")]])){
#   mod_df[[paste("conv", "urb", dec, sep = "_")]] <- NA
#   na.data <- mod_df[!complete.cases(mod_df[,c("urban_prob",
#                                               paste("urb", "targ", dec, sep = "_"))]),]
#   mod_df <- mod_df[complete.cases(mod_df[,c("urban_prob",
#                                             paste("urb", "targ", dec, sep = "_"))]),]
#   mod_df <- mod_df %>%
#     group_by(ISO3) %>%
#     do(urb_fun(., dec))
#   mod_df <- bind_rows(mod_df, na.data)
#   mod_df <- mod_df[order(mod_df$row_index),]
#   rm(na.data)
# }
# Bits cut: Re-calculate adjacent areas after updating urban areas ----
# mod_df <- adj.ag.f(m = mod_df,
#                    r = rows,
#                    c = cols,
#                    f = adj_cell_mat,
#                    crop = "crop")

# # Bits cut: Re-run predictions after updating urban areas ----
# na.data <- mod_df[is.na(mod_df$ISO_numeric),]
# mod_df <- mod_df[!is.na(mod_df$ISO_numeric),]
# 
# mod_df <- mod_df %>%
#   group_by(IUCN_region) %>%
#   do(pred_fun(., crop = "crop")) %>%
#   as.data.frame()
# # mod_df <- mod_df %>%
# #   group_by(IUCN_region) %>%
# #   do(pred_fun(., crop = "pasture"))
# 
# mod_df <- bind_rows(mod_df, na.data)
# mod_df <- mod_df[order(mod_df$row_index),]
# 
# 
# 
# Little functions to go in year function and make it more readable ----
# NONE OF THESE ARE USED AT THE MOMENT -----
# # year_fun() Adding conversion columns -----
# add_conv_col_fun <- function(df, label){
#   df[[paste("conv",
#             label,
#             sep = "_")]] <- 0
#   df[is.na(df$ISO_numeric), 
#      paste("conv",
#            label,
#            sep = "_")] <- NA
#   df[[paste("conv_pasture",
#             label,
#             sep = "_")]] <- 0
#   df[is.na(df$ISO_numeric), 
#      paste("conv_pasture",
#            label,
#            sep = "_")] <- NA
#   return(df)
# }
# # year_fun() Convert decade's worth of urban cells -----
# urban_dec_fun <- function(df, decade){
#   if(!is.null(df[[paste("urb", "targ", decade, sep = "_")]])){
#     df[[paste("conv", "urb", decade, sep = "_")]] <- NA
#     na_tmp <- 
#       df[!complete.cases(df[,c("urban_prob", paste("urb", "targ", decade, sep = "_"))]),]
#     df <- 
#       df[complete.cases(df[,c("urban_prob", paste("urb", "targ", decade, sep = "_"))]),]
#     df <- df %>%
#       group_by(ISO3) %>%
#       do(urb_fun(., decade))
#     df <- bind_rows(df, na_tmp)
#     df <- df[order(df$row_index),]
#   }
#   return(df)
# }
# 
# # year_fun() Set targets -----
# target_fun <- function(data, decade){
#   target_decade <- as.data.frame(data[!duplicated(data$ISO3),
#                                       c("ISO3", 
#                                         paste("target_crop", 
#                                               decade, 
#                                               sep = "_"))])
#   names(target_decade)[2] <- "target"
#   # Cut countries with no targets, or where the target has been met
#   target_decade <- target_decade[target_decade$target != 0 &
#                                    !is.na(target_decade$target),]
#   # Get the absolute value of the targets, so I can run a while() loop
#   target_decade$abs.target <- abs(target_decade[,"target"])
#   return(target_decade)
# }
# 
# 

# 
# # Summing up the area available for cropland ----