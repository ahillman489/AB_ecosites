##AB Ecosite Modelling
##Ashley Hillman
##Mar 4, 2026

##This script produces modelled ecosites for separate natural subregions in Alberta's forested region.
##Factor data are converted to continuous, modelled using random forests, and reclassified
##based on optimal cutoff thresholds

##-----------------------
##Load libraries
##-----------------------

library(terra)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)
library(cutpointr)
library(caret)
library(Metrics)


##----------------------------------------------------------------
##Load tabular and spatial data, and set relevant data as factors
##----------------------------------------------------------------

data <- read.csv("C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/All_sites_updated.csv")

data$DataSource <- as.factor(data$DataSource)
data$Nutrient.Regime <- as.factor(data$Nutrient.Regime)
data$Moisture.Regime <- as.factor(data$Moisture.Regime)
data$Ecosite <- as.factor(data$Ecosite)
data$NSRNAME <- as.factor(data$NSRNAME)
data$Site <- as.factor(data$Site)

##Polygons of Natural Subregions
NSR <- vect("W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/GIS/Natural_Sub_Regions_Alberta_2005_shapefile.shp") 

##Predictive rasters (pre-processed)
raster_paths <- list(
  Proportion_sand = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_sand_prop.tif",
  Proportion_clay = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_clay_prop.tif",
  Proportion_silt = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_silt_prop.tif",
  Soil_organic_carbon = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_SOC.tif",
  Surficial_geology = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SurfGeo50m.tif",
  Slope = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/Slope.tif",
  Elevation = "W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/DEM_100m.tif",
  CTI = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_cti.tif",
  TPI = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_tpi.tif",
  Profile_curvature = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_pcurv.tif",
  Tangential_curvature = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/AB_tcurv.tif",
  pH = "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_pH.tif"
)

rasters <- lapply(raster_paths, rast)

##-----------------------------------------
##Prepare input data
##-----------------------------------------

# Helper function to assign numeric values to classified Moisture and Nutrient regimes
assign_regime_values <- function(df) {
  df %>%
    mutate(Moisture_regime = case_when(
      Moisture.Regime == 'Xeric' ~ 1,
      Moisture.Regime == 'Subxeric' ~ 2,
      Moisture.Regime == 'Submesic' ~ 3,
      Moisture.Regime == 'Mesic' ~ 4,
      Moisture.Regime == 'Subhygric' ~ 5,
      Moisture.Regime == 'Hygric' ~ 6,
      Moisture.Regime == 'Subhydric' ~ 7,
      Moisture.Regime == 'Hydric' ~ 8
    ),
    Nutrient_regime = case_when(
      Nutrient.Regime == 'Very Poor' ~ 1,
      Nutrient.Regime == 'Poor' ~ 2,
      Nutrient.Regime == 'Medium' ~ 3,
      Nutrient.Regime == 'Rich' ~ 4,
      Nutrient.Regime == 'Very Rich' ~ 5
    ))
}




# List of NSR names to loop over
nsr_list <- unique(data$NSRNAME)

# Predictor variable names 
predictors <- c("Proportion_sand", "Proportion_clay", "Proportion_silt","Soil_organic_carbon", "Surficial_geology",
                "Slope","Elevation","CTI", "TPI", "Profile_curvature", "Tangential_curvature","pH")


# Build Random Forest formula
##Multivariate using moisture and nutrient regime as response variables
formula_str <- paste("Multivar(Moisture_regime, Nutrient_regime) ~", paste(predictors, collapse = " + "))
rf_formula <- as.formula(formula_str)
rf_formula

# Create output directory for tabular and spatial data
dir.create("Ecosite_output/Thresholds", showWarnings = FALSE, recursive = TRUE)

data$NSRNAME_combined <- as.character(data$NSRNAME)  # make a modifiable copy

# Combine Athabasca Plain + Peace-Athabasca Delta into one NSR
data$NSRNAME_combined[data$NSRNAME %in% c("Athabasca Plain", "Peace-Athabasca Delta")] <- "Athabasca-Delta"

##Split data into chunks to improve processing memory
NSR1 <- subset(data, NSRNAME %in% c("Upper Boreal Highlands","Boreal Subarctic","Lower Boreal Highlands"))
NSR2 <- subset(data, NSRNAME %in% c("Northern Mixedwood","Kazan Uplands", "Peace River Parkland"))
NSR3 <- subset(data, NSRNAME %in% c("Foothills Parkland","Lower Foothills","Upper Foothills"))
NSR4 <- subset(data, NSRNAME %in% c("Alpine","Subalpine","Montane"))
NSR5 <- subset(data, NSRNAME %in% c("Central Mixedwood"))
NSR6 <- subset(data, NSRNAME %in% c("Dry Mixedwood"))
NSR7 <- subset(data, NSRNAME_combined == "Athabasca-Delta")

               

# Combine into a list of subsets
nsr_groups <- list(NSR1,NSR2,NSR3,NSR4,NSR5,NSR6,NSR7)

##-------------------------------------------------------------------------------------
## Function to define optimal cutpoints for moisture and nutrient regimes separately
##-------------------------------------------------------------------------------------
get_cutpoints <- function(df, moist_col = "Moist_pred", nut_col = "Nut_pred") {
  
  results <- list(Moisture = NULL, Nutrient = NULL)  # always initialize
  
  regimes <- list(
    Moisture = list(class_var = "Moisture_regime", pred_var = moist_col),
    Nutrient = list(class_var = "Nutrient_regime", pred_var = nut_col)
  )
  
  for (reg in names(regimes)) {
    
    class_var <- regimes[[reg]]$class_var
    pred_var  <- regimes[[reg]]$pred_var
    
    if (!(class_var %in% names(df)) || !(pred_var %in% names(df))) {
      warning(sprintf("Skipping %s: required columns not found (%s, %s)", reg, class_var, pred_var))
      next
    }
    
    classes <- sort(unique(df[[class_var]]))
    cut_list <- list()
    k <- 1
    
    if (length(classes) < 2) {
      warning(sprintf("Skipping %s: not enough classes to compute cutpoints", reg))
      next
    }
    
    for (i in 1:(length(classes) - 1)) {
      c_low  <- classes[i]
      c_high <- classes[i + 1]
      sub <- df[df[[class_var]] %in% c(c_low, c_high), ]
      
      if (length(unique(sub[[class_var]])) < 2) next
      
      sub$class_bin <- ifelse(sub[[class_var]] == c_high, 1, 0)
      
      cutpt <- tryCatch({
        cutpointr(
          data  = sub,
          x     = !!as.name(pred_var),
          class = class_bin,
          method = maximize_metric,
          metric = sum_sens_spec
        )
      }, error = function(e) {
        warning(sprintf("Cutpoint computation failed for %s between %s and %s: %s", reg, c_low, c_high, e$message))
        return(NULL)
      })
      
      if (is.null(cutpt)) next
      
      cut_list[[k]] <- data.frame(
        low_class   = c_low,
        high_class  = c_high,
        optimal_cut = cutpt$optimal_cutpoint,
        sensitivity = cutpt$sensitivity,
        specificity = cutpt$specificity,
        auc         = cutpt$AUC,
        n_low       = sum(sub$class_bin == 0),
        n_high      = sum(sub$class_bin == 1)
      )
      
      k <- k + 1
    }
    
    if (length(cut_list) > 0) {
      results[[reg]] <- do.call(rbind, cut_list)
    } else {
      results[[reg]] <- NULL
    }
  }
  
  return(results)
}

##----------------------------------------------------------
## Helper function: Fix reversed or non-monotonic cutpoints
##----------------------------------------------------------
fix_cutpoints <- function(cut_table) {
  
  ct <- cut_table
  
  # If only 1 cutpoint, nothing to fix
  if (nrow(ct) < 2) return(ct)
  
  # Extract numeric cutpoints
  cuts <- ct$optimal_cut
  
  # If the sequence is already strictly increasing, return as-is
  if (all(diff(cuts) > 0)) {
    return(ct)
  }
  
  message("⚠ Non-monotonic cutpoints detected — merging affected classes.")
  
  i <- 1
  while (i < nrow(ct)) {
    
    # Bad break if cuts[i] >= cuts[i+1]
    if (ct$optimal_cut[i] >= ct$optimal_cut[i + 1]) {
      
      n_low  <- ct$n_low[i]
      n_high <- ct$n_high[i]
      
      if (n_low < n_high) {
        remove_row <- i
        merged_from <- ct$low_class[i]
        merged_to   <- ct$high_class[i]
      } else {
        remove_row <- i + 1
        merged_from <- ct$high_class[i]
        merged_to   <- ct$low_class[i]
      }
      
      message(sprintf("  → Merging %s → %s", merged_from, merged_to))
      
      # Remove the class boundary
      ct <- ct[-remove_row, ]
      
      # Recompute cuts for next iteration
      next
    }
    
    i <- i + 1
  }
  
  return(ct)
}


library(caret)

# Regression metrics
reg_metrics <- function(obs, pred) {
  data.frame(
    RMSE = RMSE(pred, obs),
    R2   = R2(pred, obs),
    MAE  = mean(abs(pred - obs), na.rm = TRUE)
  )
}


##------------------------------------------------------------------------------
## Run random forest models with train/test validation per NSR
## Thresholds learned on training data and applied to testing data
##------------------------------------------------------------------------------

for (nsr_set in nsr_groups) {
  for (nsr_name in unique(nsr_set$NSRNAME)) {
    
    message("Processing NSR: ", nsr_name)
    
    tryCatch({
      
      ##------------------------------------------------
      ## Subset NSR data and create spatial points
      ##------------------------------------------------
      nsr_data  <- subset(nsr_set, NSRNAME == nsr_name)
      nsr_shape <- NSR[NSR$NSRNAME == nsr_name]
      
      nsr_data$ID <- seq_len(nrow(nsr_data))
      nsr_points <- vect(
        nsr_data,
        geom = c("Easting", "Northing"),
        crs = "EPSG:3400"
      )
      
      masked <- lapply(rasters, function(r) mask(r, nsr_shape))
      stack_rasters <- rast(masked)
      names(stack_rasters) <- predictors
      
      raster_values <- terra::extract(stack_rasters, nsr_points)
      nsr_combined  <- left_join(nsr_data, raster_values, by = "ID")
      
      nsr_combined <- assign_regime_values(nsr_combined)
      nsr_combined <- nsr_combined %>%
        mutate(across(where(is.character), as.factor))
      
      if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(nsr_combined))) {
        message("Skipping ", nsr_name, ": missing regimes.")
        next
      }
      
      ##------------------------------------------------
      ## Split into training and testing
      ##------------------------------------------------
      set.seed(123)
      
      model_df <- nsr_combined %>%
        dplyr::select(
          Site,
          Easting,
          Northing,
          NSRNAME,
          all_of(all.vars(rf_formula))
        ) %>%
        na.omit()
      
      
      train_idx <- createDataPartition(
        y = interaction(
          model_df$Moisture_regime,
          model_df$Nutrient_regime,
          drop = TRUE
        ),
        p = 0.7,
        list = FALSE
      )
      
      train_df <- model_df[train_idx, ]
      test_df  <- model_df[-train_idx, ]
      
      message(
        "Training n = ", nrow(train_df),
        " | Testing n = ", nrow(test_df)
      )
      
      
      ##------------------------------------------------
      ## Fit random forest (TRAINING DATA ONLY)
      ##------------------------------------------------
      model <- rfsrc(
        rf_formula,
        data = train_df,
        importance = "permute",
        na.action = "na.impute"
      )
      
      ##------------------------------------------------
      ## Predict on training and testing
      ##------------------------------------------------
      train_pred <- predict(model, newdata = train_df)$regrOutput
      test_pred  <- predict(model, newdata = test_df)$regrOutput
      
      train_df$Moist_pred <- train_pred$Moisture_regime$predicted
      train_df$Nut_pred   <- train_pred$Nutrient_regime$predicted
      
      test_df$Moist_pred  <- test_pred$Moisture_regime$predicted
      test_df$Nut_pred    <- test_pred$Nutrient_regime$predicted
      
      ##------------------------------------------------
      ## Continuous validation metrics (optional)
      ##------------------------------------------------
      validation_summary <- bind_rows(
        cbind(
          Outcome = "Moisture",
          Dataset = "Training",
          reg_metrics(train_df$Moisture_regime, train_df$Moist_pred)
        ),
        cbind(
          Outcome = "Moisture",
          Dataset = "Testing",
          reg_metrics(test_df$Moisture_regime, test_df$Moist_pred)
        ),
        cbind(
          Outcome = "Nutrient",
          Dataset = "Training",
          reg_metrics(train_df$Nutrient_regime, train_df$Nut_pred)
        ),
        cbind(
          Outcome = "Nutrient",
          Dataset = "Testing",
          reg_metrics(test_df$Nutrient_regime, test_df$Nut_pred)
        )
      )
      
      write.csv(
        validation_summary,
        paste0(
          "Ecosite_output/Thresholds/Validation_",
          gsub(" ", "_", nsr_name),
          ".csv"
        ),
        row.names = FALSE
      )
      
      ##------------------------------------------------
      ## Compute optimal cutpoints (TRAINING DATA ONLY)
      ##------------------------------------------------
      cuts <- get_cutpoints(train_df)
      
      moist_cut_table <- fix_cutpoints(cuts$Moisture)
      nut_cut_table   <- fix_cutpoints(cuts$Nutrient)
      
      write.csv(
        moist_cut_table,
        paste0(
          "Ecosite_output/Thresholds/Cutpoints_Moisture_",
          gsub(" ", "_", nsr_name),
          ".csv"
        ),
        row.names = FALSE
      )
      
      write.csv(
        nut_cut_table,
        paste0(
          "Ecosite_output/Thresholds/Cutpoints_Nutrient_",
          gsub(" ", "_", nsr_name),
          ".csv"
        ),
        row.names = FALSE
      )
      
      ##------------------------------------------------
      ## Define class labels dynamically
      ##------------------------------------------------
      moisture_labels <- c(
        "Xeric", "Subxeric", "Submesic", "Mesic",
        "Subhygric", "Hygric", "Subhydric", "Hydric"
      )
      
      nutrient_labels <- c(
        "Very Poor", "Poor", "Medium", "Rich", "Very Rich"
      )
      
      n_moist_intervals <- length(moist_cut_table$optimal_cut) + 1
      n_nut_intervals   <- length(nut_cut_table$optimal_cut) + 1
      
      moisture_labels_subset <- moisture_labels[1:n_moist_intervals]
      nutrient_labels_subset <- nutrient_labels[1:n_nut_intervals]
      
      ##------------------------------------------------
      ## Classify TRAINING data (diagnostic only)
      ##------------------------------------------------
      train_df$Moist_class_pred <- cut(
        train_df$Moist_pred,
        breaks = c(-Inf, moist_cut_table$optimal_cut, Inf),
        labels = moisture_labels_subset,
        include.lowest = TRUE
      )
      
      
      train_df$Nut_class_pred <- cut(
        train_df$Nut_pred,
        breaks = c(-Inf, nut_cut_table$optimal_cut, Inf),
        labels = nutrient_labels_subset,
        include.lowest = TRUE
      )
      
      ##reference classes from original data
      
      train_df$Moist_class <- factor(
        train_df$Moisture_regime,
        levels = 1:8,
        labels = moisture_labels
      )
      
      train_df$Nut_class <- factor(
        train_df$Nutrient_regime,
        levels = 1:5,
        labels = nutrient_labels
      )
      
      
      ##------------------------------------------------
      ## Classify TESTING data (USED FOR ACCURACY)
      ##------------------------------------------------
      test_df$Moist_class_pred <- cut(
        test_df$Moist_pred,
        breaks = c(-Inf, moist_cut_table$optimal_cut, Inf),
        labels = moisture_labels_subset,
        include.lowest = TRUE
      )
      
      
      test_df$Nut_class_pred <- cut(
        test_df$Nut_pred,
        breaks = c(-Inf, nut_cut_table$optimal_cut, Inf),
        labels = nutrient_labels_subset,
        include.lowest = TRUE
      )
      
      ##reference classes from original data
      
      test_df$Moist_class <- factor(
        test_df$Moisture_regime,
        levels = 1:8,
        labels = moisture_labels
      )
      
      test_df$Nut_class <- factor(
        test_df$Nutrient_regime,
        levels = 1:5,
        labels = nutrient_labels
      )
      
      ##------------------------------------------------
      ## Accuracy assessment (TESTING DATA ONLY)
      ##------------------------------------------------
      
      # Moisture confusion matrix
      moist_cm <- confusionMatrix(
        data      = test_df$Moist_class_pred,
        reference = test_df$Moist_class
      )
      
      moist_tbl <- as.data.frame.matrix(moist_cm$table)
      
      write.csv(
        moist_tbl,
        paste0(
        "Ecosite_output/Thresholds/Validation_Confusion_Moisture_",
        gsub(" ", "_", nsr_name),
        ".csv"
      ),
        row.names = TRUE
      )
      
      # Nutrient confusion matrix
      nut_cm <- confusionMatrix(
        data      = test_df$Nut_class_pred,
        reference = test_df$Nut_class
      )
      
      nut_tbl <- as.data.frame.matrix(nut_cm$table)
      
      write.csv(
        nut_tbl,
        paste0(
          "Ecosite_output/Thresholds/Validation_Confusion_Nutrient_",
          gsub(" ", "_", nsr_name),
          ".csv"
        ),
        row.names = TRUE
      )
      
      
      ##------------------------------------------------
      ## Save point predictions
      ##------------------------------------------------
      write.csv(
        bind_rows(
          mutate(train_df, Dataset = "Train"),
          mutate(test_df,  Dataset = "Test")
          ),
        paste0(
          "Ecosite_output/Thresholds/PointPredictions_",
          gsub(" ", "_", nsr_name),
          ".csv"
        ),
        row.names = FALSE
      )
      
      message(
        "Finished processing ", nsr_name,
        " — thresholds trained on training data, accuracy assessed on testing."
      )
      
    }, error = function(e) {
      message("Error processing ", nsr_name, ": ", e$message)
    })
  }
}


##Combine all confusion matrices

library(caret)
library(dplyr)

# Define the full set of labels (must match the labels used in the NSR-specific confusion matrices)
moisture_labels <- c("Xeric", "Subxeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
nutrient_labels <- c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")

# Empty cumulative confusion tables
full_moist_cm <- matrix(0, nrow = length(moisture_labels), ncol = length(moisture_labels),
                        dimnames = list(Reference = moisture_labels, Prediction = moisture_labels))

full_nut_cm <- matrix(0, nrow = length(nutrient_labels), ncol = length(nutrient_labels),
                      dimnames = list(Reference = nutrient_labels, Prediction = nutrient_labels))

for (nsr_set in nsr_groups) {
  for (nsr_name in unique(nsr_set$NSRNAME)) {
    
    # Read the saved CSV confusion matrices (or use the test_df directly if still in memory)
    moist_file <- paste0("Ecosite_output/Thresholds/Validation_Confusion_Moisture_", gsub(" ", "_", nsr_name), ".csv")
    nut_file   <- paste0("Ecosite_output/Thresholds/Validation_Confusion_Nutrient_", gsub(" ", "_", nsr_name), ".csv")
    
    if (file.exists(moist_file) & file.exists(nut_file)) {
      moist_tbl <- read.csv(moist_file, row.names = 1, check.names = FALSE)
      nut_tbl   <- read.csv(nut_file,   row.names = 1, check.names = FALSE)
      
      # Ensure row/col names match full label sets
      moist_tbl <- moist_tbl[moisture_labels, moisture_labels, drop = FALSE]
      nut_tbl   <- nut_tbl[nutrient_labels, nutrient_labels, drop = FALSE]
      
      # Add to cumulative confusion matrices
      full_moist_cm <- full_moist_cm + as.matrix(moist_tbl)
      full_nut_cm   <- full_nut_cm   + as.matrix(nut_tbl)
    }
  }
}

write.csv(full_moist_cm, "Ecosite_output/Thresholds/Full_Confusion_Moisture.csv", row.names = TRUE)
write.csv(full_nut_cm,   "Ecosite_output/Thresholds/Full_Confusion_Nutrient.csv", row.names = TRUE)


