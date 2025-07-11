library(terra)
library(dplyr)
library(randomForestSRC)


data <- read.csv('All_sites.csv')
View(data)

data$DataSource <- as.factor(data$DataSource)
data$Nutrient.Regime <- as.factor(data$Nutrient.Regime)
data$Moisture.Regime <- as.factor(data$Moisture.Regime)
data$Ecosite <- as.factor(data$Ecosite)
data$Natural.Subregion <- as.factor(data$NSRNAME)
data$Site <- as.factor(data$Site)

##Remove Open Water observations, plus swamp, alkali, and not treed for exploratory analyses
data1 <- subset(data, Nutrient.Regime !='Open Water'& Nutrient.Regime!= 'Alkali'& Nutrient.Regime!= 'Swamp'
                & Nutrient.Regime!= 'Not Treed'& Moisture.Regime!="Swamp"& Moisture.Regime!="No Value Provided"&
                  Nutrient.Regime!="No Value Provided" & Nutrient.Regime!=" " & Moisture.Regime !=" ")

View(data1)
# Load spatial and raster data only once
NSR <- vect("W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/GIS/Natural_Sub_Regions_Alberta_2005_shapefile.shp") 

# Load raster layers
raster_paths <- list(
  mean_Clay = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/mean_Clay_0_60.tif",
  mean_Sand = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/mean_Sand_0_60.tif",
  mean_SOC = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/mean_SOC_0_60.tif",
  DEM = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/DEM_100m.tif",
  Slope = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/Slope.tif",
  Aspect = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/Aspect.tif",
  GG = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/GreatGroups_10TM.tif",
  cti = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_cti.tif",
  tpi = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_tpi.tif",
  pcurv = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_pcurv.tif",
  tcurv = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_tcurv.tif",
  tr = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/terrain_roughness_AB.tif",
  MAT = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/MAT_AB.tif",
  MAP = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/MAP_AB.tif",
  pH = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/mean_pH_0_60.tif",
  BD = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/mean_BD_0_60.tif"
)

rasters <- lapply(raster_paths, rast)


# Helper function to assign numeric values
assign_regime_values <- function(df) {
  df %>%
    mutate(Moist_num = case_when(
      Moisture.Regime == 'Subxeric' ~ 0.5,
      Moisture.Regime == 'Xeric' ~ 1.0,
      Moisture.Regime == 'Submesic' ~ 1.5,
      Moisture.Regime == 'Mesic' ~ 2.0,
      Moisture.Regime == 'Subhygric' ~ 2.5,
      Moisture.Regime == 'Hygric' ~ 3.0,
      Moisture.Regime == 'Subhydric' ~ 3.5,
      Moisture.Regime == 'Hydric' ~ 4.0
    ),
    Nut_num = case_when(
      Nutrient.Regime == 'Very Poor' ~ 0.5,
      Nutrient.Regime == 'Poor' ~ 1.0,
      Nutrient.Regime == 'Medium' ~ 2.0,
      Nutrient.Regime == 'Rich' ~ 3.0,
      Nutrient.Regime == 'Very Rich' ~ 3.5
    ))
}




# List of NSR names to loop over
nsr_list <- unique(data1$NSRNAME)

# Predictor variable names (must match raster and data)
predictors <- c("Mean_Clay", "Mean_Sand", "Mean_SOC", "Elevation", "Roughness",
                "MAT", "MAP", "Slope", "Aspect", "Soil_GG", "CTI", "TPI", "pcurv", "tcurv","pH","Bulk_Density")

# Build formula safely
formula_str <- paste("Multivar(Moist_num, Nut_num) ~", paste(predictors, collapse = " + "))
rf_formula <- as.formula(formula_str)
rf_formula

# Create output directory if it doesn't exist
dir.create("Ecosite_output", showWarnings = FALSE, recursive = TRUE)

##Split data into chunks for processing memory

NSR1 <- subset(data1, NSRNAME %in% c("Upper Boreal Highlands","Boreal Subarctic","Lower Boreal Highlands"))
NSR2 <- subset(data1, NSRNAME %in% c("Northern Mixedwood","Kazan Uplands", "Peace River Parkland"))
NSR3 <- subset(data1, NSRNAME %in% c("Foothills Parkland","Lower Foothills","Upper Foothills"))
NSR4 <- subset(data1, NSRNAME %in% c("Alpine","Subalpine","Montane"))
##NSR5 <- subset(data1, NSRNAME %in% c("Central Parkland","Dry Mixedgrass","Foothills Fescue","Mixedgrass","Northern Fescue"))

# Combine into a list of subsets
nsr_groups <- list(NSR1, NSR2, NSR3, NSR4)

# Loop through each group
for (nsr_set in nsr_groups) {
  for (nsr_name in unique(nsr_set$NSRNAME)) {
    message("Processing NSR: ", nsr_name)
    
    tryCatch({
      
      # Subset and convert to spatial
      nsr_data <- subset(nsr_set, NSRNAME == nsr_name)
      nsr_points <- vect(nsr_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")
      nsr_shape <- NSR[NSR$NSRNAME == nsr_name]
      
      # Mask and trim rasters
      masked <- lapply(rasters, function(r) trim(mask(r, nsr_shape)))
      
      # Align all rasters to first raster
      template <- masked[[1]]
      masked_aligned <- lapply(masked, function(r) {
        if (!compareGeom(r, template, stopOnError = FALSE)) {
          resample(r, template, method = "bilinear")
        } else {
          r
        }
      })
      
      stack_rasters <- rast(masked_aligned)
      names(stack_rasters) <- predictors
      
      # Extract raster values to points
      values <- terra::extract(stack_rasters, nsr_points)[, -1]
      nsr_combined <- cbind(nsr_data, values)
      
      # Assign numeric regime values
      nsr_combined <- assign_regime_values(nsr_combined)
      
      # Clean data
      nsr_combined <- nsr_combined[, colSums(is.na(nsr_combined)) < nrow(nsr_combined)]
      nsr_combined[] <- lapply(nsr_combined, function(x) if (is.character(x)) as.factor(x) else x)
      
      # Skip if key variables missing
      if (!("Moist_num" %in% names(nsr_combined)) || !("Nut_num" %in% names(nsr_combined))) {
        message("Skipping ", nsr_name, ": Moist_num or Nut_num missing.")
        next
      }
      
      # Fit model
      model <- rfsrc(rf_formula, data = nsr_combined, importance = TRUE)
      message("\nModel Summary for ", nsr_name, ":")
      print(model)
      plot(model)
      
      # Predict
      pred_df <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
      pred_result <- predict(model, newdata = pred_df, na.action = "na.omit")
      pred_df$Moist_pred <- pred_result$regrOutput$Moist_num$predicted
      pred_df$Nut_pred <- pred_result$regrOutput$Nut_num$predicted
      
      # Create rasters
      moist_rast <- rast(pred_df[, c("x", "y", "Moist_pred")], type = "xyz")
      nut_rast   <- rast(pred_df[, c("x", "y", "Nut_pred")], type = "xyz")
      
      writeRaster(moist_rast, paste0("Ecosite_output/Climate/Moist_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      writeRaster(nut_rast, paste0("Ecosite_output/Climate/Nut_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      
      # Classify
      moist_breaks <- quantile(values(moist_rast), probs = seq(0, 1, length.out = 9), na.rm = TRUE)
      moist_cat <- cbind(moist_breaks[-length(moist_breaks)], moist_breaks[-1], 1:8)
      moist_rast_cat <- classify(moist_rast, moist_cat)
      levels(moist_rast_cat) <- data.frame(
        ID = 1:8,
        Category = c("Subxeric", "Xeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
      )
      
      nut_breaks <- quantile(values(nut_rast), probs = seq(0, 1, length.out = 6), na.rm = TRUE)
      nut_cat <- cbind(nut_breaks[-length(nut_breaks)], nut_breaks[-1], 1:5)
      nut_rast_cat <- classify(nut_rast, nut_cat)
      levels(nut_rast_cat) <- data.frame(
        ID = 1:5,
        Category = c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
      )
      
      ecosite_rast <- (moist_rast_cat * 10) + nut_rast_cat
      names(ecosite_rast) <- paste0("Ecosite_", gsub(" ", "_", nsr_name))
      
      writeRaster(ecosite_rast, paste0("Ecosite_output/Climate/Ecosite_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      
      message("Finished: ", nsr_name)
      
    }, error = function(e) {
      message("Error processing ", nsr_name, ": ", e$message)
    })
  }
}

##Improve memory for large subregions

NSR7 <- subset(data1, NSRNAME %in% c("Central Mixedwood"))
NSR8 <- subset(data1, NSRNAME %in% c("Dry Mixedwood"))##Do one at a time or memory crashes

for (nsr_name in unique(NSR8$Natural.Subregion)) {
  message("Processing NSR: ", nsr_name)
  
  tryCatch({
    
    # Subset and convert to spatial
    nsr_data <- subset(NSR8, Natural.Subregion == nsr_name)
    nsr_points <- vect(nsr_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")
    
    nsr_shape <- NSR[NSR$NSRNAME == nsr_name]
    
    # Mask all rasters
    masked <- lapply(rasters, function(r) trim(mask(r, nsr_shape)))
    
    # Align all masked rasters to a common template (first raster)
    template <- masked[[1]]
    masked_aligned <- lapply(masked, function(r) {
      if (!compareGeom(r, template, stopOnError = FALSE)) {
        resample(r, template, method = "bilinear")
      } else {
        r
      }
    })
    
    # Stack
    stack_rasters <- rast(masked_aligned)
    names(stack_rasters) <- predictors
    
    # Extract raster values to points
    values <- terra::extract(stack_rasters, nsr_points)[, -1]
    nsr_combined <- cbind(nsr_data, values)
    
    # Assign numeric values for regimes
    nsr_combined <- assign_regime_values(nsr_combined)
    
    # Remove unnecessary columns
    nsr_combined <- nsr_combined[, !names(nsr_combined) %in% c("Cardinal_Direction", "Zone", "Transect", "Transect_Location")]
    nsr_combined <- nsr_combined[, colSums(is.na(nsr_combined)) < nrow(nsr_combined)]
    nsr_combined[] <- lapply(nsr_combined, function(x) if (is.character(x)) as.factor(x) else x)
    
    # Skip if key variables are missing
    if (!("Moist_num" %in% names(nsr_combined)) || !("Nut_num" %in% names(nsr_combined))) {
      message("Skipping ", nsr_name, ": Moist_num or Nut_num missing.")
      next
    }
    
    
    # Fit model
    model <- rfsrc(rf_formula, data = nsr_combined, importance = TRUE)
    
    # Show model summary in console
    message("\nModel Summary for ", nsr_name, ":")
    print(model)
    plot(model)
    
    # Define block size
    block_size <- 10000
    pred_df_all <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
    n <- nrow(pred_df_all)
    
    # Initialize results
    moist_vals <- numeric(n)
    nut_vals <- numeric(n)
    
    # Predict in chunks
    for (i in seq(1, n, by = block_size)) {
      i_end <- min(i + block_size - 1, n)
      chunk <- pred_df_all[i:i_end, ]
      
      pred_chunk <- predict(model, newdata = chunk, na.action = "na.omit")
      
      moist_vals[i:i_end] <- pred_chunk$regrOutput$Moist_num$predicted
      nut_vals[i:i_end] <- pred_chunk$regrOutput$Nut_num$predicted
    }
    
    # Attach predictions to full df
    pred_df_all$Moist_pred <- moist_vals
    pred_df_all$Nut_pred <- nut_vals
    
    # Convert to rasters
    moist_rast <- rast(pred_df_all[, c("x", "y", "Moist_pred")], type = "xyz")
    nut_rast <- rast(pred_df_all[, c("x", "y", "Nut_pred")], type = "xyz")
    
    writeRaster(moist_rast, paste0("Ecosite_output/Climate/Moist_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    writeRaster(nut_rast, paste0("Ecosite_output/Climate/Nut_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    
    # Classify
    # Get quantile breaks for Moisture
    moist_breaks <- quantile(values(moist_rast), probs = seq(0, 1, length.out = 9), na.rm = TRUE)
    
    # Create reclassification matrix
    moist_cat <- cbind(
      moist_breaks[-length(moist_breaks)],
      moist_breaks[-1],
      1:8
    )
    
    # Classify
    moist_rast_cat <- classify(moist_rast, moist_cat)
    
    # Assign moisture class labels
    levels(moist_rast_cat) <- data.frame(
      ID = 1:8,
      Category = c("Subxeric", "Xeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
    )
    
    # Get quantile breaks for Nutrient
    nut_breaks <- quantile(values(nut_rast), probs = seq(0, 1, length.out = 6), na.rm = TRUE)
    
    # Create reclassification matrix
    nut_cat <- cbind(
      nut_breaks[-length(nut_breaks)],
      nut_breaks[-1],
      1:5
    )
    
    # Classify
    nut_rast_cat <- classify(nut_rast, nut_cat)
    
    # Assign nutrient class labels
    levels(nut_rast_cat) <- data.frame(
      ID = 1:5,
      Category = c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
    )
    
    ecosite_rast <- (moist_rast_cat * 10) + nut_rast_cat
    names(ecosite_rast) <- paste0("Ecosite_", gsub(" ", "_", nsr_name))
    
    # Save
    writeRaster(ecosite_rast, paste0("Ecosite_output/No_climate/Ecosite_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    
    message("Finished: ", nsr_name)
    
  }, error = function(e) {
    message("Error processing ", nsr_name, ": ", e$message)
  })
}

##Combine Athabasca Plain and Peace-Athabasca Delta (only 9 points in Delta, no model predictions made)

NSR9 <- subset(data1, NSRNAME %in% c("Athabasca Plain","Peace-Athabasca Delta"))

# Define combined subregion name and members
nsr_name <- "Athabasca_Plain__Peace_Athabasca_Delta"
combined_subregions <- c("Athabasca Plain", "Peace-Athabasca Delta")

process_combined_nsr <- function(nsr_name, combined_subregions, NSR9, rasters, predictors, rf_formula) {
  tryCatch({
    message("Processing Combined NSR: ", nsr_name)
    
    # Combine attribute data
    combined_data <- subset(NSR9, NSRNAME %in% combined_subregions)
    if (nrow(combined_data) == 0) {
      message("No data found for ", nsr_name)
      return(NULL)
    }
    
    # Convert to spatial points
    combined_points <- vect(combined_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")
    
    # Combine and dissolve spatial shapes
    combined_shape <- NSR[NSR$NSRNAME %in% combined_subregions, ]
    combined_shape <- aggregate(combined_shape)
    
    # Mask and trim rasters
    masked <- lapply(rasters, function(r) trim(mask(r, combined_shape)))
    
    # Align all rasters to first
    template <- masked[[1]]
    masked_aligned <- lapply(masked, function(r) {
      if (!compareGeom(r, template, stopOnError = FALSE)) {
        resample(r, template, method = "bilinear")
      } else {
        r
      }
    })
    
    # Stack
    stack_rasters <- rast(masked_aligned)
    names(stack_rasters) <- predictors
    
    # Extract to points
    values <- terra::extract(stack_rasters, combined_points)[, -1]
    combined_df <- cbind(combined_data, values)
    
    # Assign numeric regime values
    combined_df <- assign_regime_values(combined_df)
    
    # Clean columns
    combined_df <- combined_df[, !names(combined_df) %in% c("Cardinal_Direction", "Zone", "Transect", "Transect_Location")]
    combined_df <- combined_df[, colSums(is.na(combined_df)) < nrow(combined_df)]
    combined_df[] <- lapply(combined_df, function(x) if (is.character(x)) as.factor(x) else x)
    
    # Check required vars
    if (!("Moist_num" %in% names(combined_df)) || !("Nut_num" %in% names(combined_df))) {
      message("Skipping ", nsr_name, ": Moist_num or Nut_num missing.")
      return(NULL)
    }
    
    # Fit model
    model <- rfsrc(rf_formula, data = combined_df, importance = TRUE)
    
    message("\nModel Summary for ", nsr_name, ":")
    print(model)
    plot(model)
    
    # Predict
    pred_df <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
    pred_result <- predict(model, newdata = pred_df, na.action = "na.omit")
    pred_df$Moist_pred <- pred_result$regrOutput$Moist_num$predicted
    pred_df$Nut_pred   <- pred_result$regrOutput$Nut_num$predicted
    
    # To raster
    moist_rast <- rast(pred_df[, c("x", "y", "Moist_pred")], type = "xyz")
    nut_rast   <- rast(pred_df[, c("x", "y", "Nut_pred")], type = "xyz")
    
    writeRaster(moist_rast, paste0("Ecosite_output/Climate/Moist_", nsr_name, ".tif"), overwrite = TRUE)
    writeRaster(nut_rast,   paste0("Ecosite_output/Climate/Nut_",   nsr_name, ".tif"), overwrite = TRUE)
    
    # Classify
    moist_breaks <- quantile(values(moist_rast), probs = seq(0, 1, length.out = 9), na.rm = TRUE)
    moist_cat <- cbind(moist_breaks[-length(moist_breaks)], moist_breaks[-1], 1:8)
    moist_rast_cat <- classify(moist_rast, moist_cat)
    levels(moist_rast_cat) <- data.frame(
      ID = 1:8,
      Category = c("Subxeric", "Xeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
    )
    
    nut_breaks <- quantile(values(nut_rast), probs = seq(0, 1, length.out = 6), na.rm = TRUE)
    nut_cat <- cbind(nut_breaks[-length(nut_breaks)], nut_breaks[-1], 1:5)
    nut_rast_cat <- classify(nut_rast, nut_cat)
    levels(nut_rast_cat) <- data.frame(
      ID = 1:5,
      Category = c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
    )
    
    ecosite_rast <- (moist_rast_cat * 10) + nut_rast_cat
    names(ecosite_rast) <- paste0("Ecosite_", nsr_name)
    
    writeRaster(ecosite_rast, paste0("Ecosite_output/Climate/Ecosite_", nsr_name, ".tif"), overwrite = TRUE)
    
    message("Finished: ", nsr_name)
    
  }, error = function(e) {
    message("Error processing ", nsr_name, ": ", e$message)
  })
}

# Define combined subregion name and members
nsr_name <- "Athabasca_Plain__Peace_Athabasca_Delta"
combined_subregions <- c("Athabasca Plain", "Peace-Athabasca Delta")

process_combined_nsr(
  nsr_name = "Athabasca_Plain__Peace_Athabasca_Delta",
  combined_subregions = c("Athabasca Plain", "Peace-Athabasca Delta"),
  NSR = NSR,
  rasters = rasters,
  predictors = predictors,
  rf_formula = rf_formula
)

