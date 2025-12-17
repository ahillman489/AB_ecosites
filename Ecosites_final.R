##AB Ecosite Modelling
##Ashley Hillman
##Dec 17, 2025

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

##Split data into chunks to improve processing memory
NSR1 <- subset(data, NSRNAME %in% c("Upper Boreal Highlands","Boreal Subarctic","Lower Boreal Highlands"))
NSR2 <- subset(data, NSRNAME %in% c("Northern Mixedwood","Kazan Uplands", "Peace River Parkland"))
NSR3 <- subset(data, NSRNAME %in% c("Foothills Parkland","Lower Foothills","Upper Foothills"))
NSR4 <- subset(data, NSRNAME %in% c("Alpine","Subalpine","Montane"))
##NSR5 <- subset(data1, NSRNAME %in% c("Central Parkland","Dry Mixedgrass","Foothills Fescue","Mixedgrass","Northern Fescue"))

# Combine into a list of subsets
nsr_groups <- list(NSR1,NSR2,NSR3,NSR4)

##-------------------------------------------------------------------------------------
## Function to define optimal cutpoints for moisture and nutrient regimes separately
##-------------------------------------------------------------------------------------
get_cutpoints <- function(df) {
  
  results <- list()
  lookup <-  list()
  
  regimes <- list(
    Moisture = list(
      class_var = "Moisture_regime",
      pred_var  = "Moist_pred_point"
    ),
    
    Nutrient = list(
      class_var = "Nutrient_regime",
      pred_var  = "Nutrient_pred_point"
    )
  )
  
  for (reg in names(regimes)) {
    
    class_var <- regimes[[reg]]$class_var
    pred_var  <- regimes[[reg]]$pred_var
    
    # Keep only classes present & in the correct order
    classes <- sort(unique(df[[class_var]]))
    
    cut_list <- list()
    k <- 1
    
    for (i in 1:(length(classes) - 1)) {
      
      c_low  <- classes[i]
      c_high <- classes[i + 1]
      
      sub <- df[df[[class_var]] %in% c(c_low, c_high), ]
      
      if (length(unique(sub[[class_var]])) < 2)
        next
      
      sub$class_bin <- ifelse(sub[[class_var]] == c_high, 1, 0)
      
      cutpt <- cutpointr(
        data  = sub,
        x     = !!as.name(pred_var),
        class = class_bin,
        method = maximize_metric,
        metric = sum_sens_spec
      )
      
      cut_list[[k]] <- data.frame(
        low_class       = c_low,
        high_class      = c_high,
        optimal_cut     = cutpt$optimal_cutpoint,
        sensitivity     = cutpt$sensitivity,
        specificity     = cutpt$specificity,
        auc             = cutpt$AUC,
        n_low           = sum(sub$class_bin == 0),
        n_high          = sum(sub$class_bin == 1)
      )
      
      k <- k + 1
    }
    
    # Final table sorted in the correct class order
    results[[reg]] <- do.call(rbind, cut_list)
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


##------------------------------------------------------------------------------
##Run random forest models, looping through each Natural Subregion separately
##------------------------------------------------------------------------------

# Loop through each NSR group
for (nsr_set in nsr_groups) {
  for (nsr_name in unique(nsr_set$NSRNAME)) {
    message("Processing NSR: ", nsr_name)
    
    tryCatch({
      
      ##------------------------------------------------
      ## Stack rasters and extract data to input points
      ##------------------------------------------------
      
      # Subset data for the NSR and create spatial points
      nsr_data <- subset(nsr_set, NSRNAME == nsr_name)
      nsr_shape <- NSR[NSR$NSRNAME == nsr_name]
      
      # Create ID column (stable row index for joining)
      nsr_data$ID <- seq_len(nrow(nsr_data))
      
      nsr_points <- vect(nsr_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")
      
      # Mask and trim rasters to NSR boundary
      masked <- lapply(rasters, function(r) (mask(r, nsr_shape)))
      
      
      # Create raster stack and name layers consistently
      stack_rasters <- rast(masked)
      names(stack_rasters) <- predictors
      
      stack_rasters
      
      # Extract raster values at data points
      raster_values <- terra::extract(stack_rasters, nsr_points)  
      nsr_combined <- left_join(nsr_data, raster_values, by = "ID")
      
      # Assign numeric regime values
      nsr_combined <- assign_regime_values(nsr_combined)
      
      # Drop columns with all NA and convert character columns to factors
      
      nsr_combined <- nsr_combined %>%
        mutate(across(where(is.character), as.factor))
      
      # Skip if essential numeric columns missing
      if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(nsr_combined))) {
        message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
        next
      }
      
      
      # Check missing predictors
      setdiff(predictors, names(nsr_combined))
      
      ##------------------------------------------------
      ## Fit random forest model
      ##------------------------------------------------
      
      model <- rfsrc(rf_formula, data = nsr_combined,
                     importance = "permute", na.action = "na.impute")
      
      point_pred <- predict(model)$regrOutput
      
      nsr_combined$Moist_pred_point    <- point_pred$Moisture_regime$predicted
      nsr_combined$Nutrient_pred_point <- point_pred$Nutrient_regime$predicted
      
      
      # Save point predictions
      write.csv(
        nsr_combined,
        paste0("Ecosite_output/Thresholds/PointPredictions_", gsub(" ", "_", nsr_name), ".csv"),
        row.names = FALSE
      )
      message("Point predictions saved for ", nsr_name)
      
      ##---------------------------------------------
      ## Compute optimal thresholds
      ##---------------------------------------------
      
      cuts <- get_cutpoints(nsr_combined)
      
      moist_cut_table <- fix_cutpoints(cuts$Moisture)
      nut_cut_table   <- fix_cutpoints(cuts$Nutrient)
      
      # Save results
      write.csv(
        moist_cut_table,
        paste0("Ecosite_output/Thresholds/Cutpoints_Moisture_", gsub(" ", "_", nsr_name), ".csv"),
        row.names = FALSE
      )
      
      write.csv(
        nut_cut_table,
        paste0("Ecosite_output/Thresholds/Cutpoints_Nutrient_", gsub(" ", "_", nsr_name), ".csv"),
        row.names = FALSE
      )
      
      message("Cutpoints saved for ", nsr_name)
      
      # Print model summaries
      message("\n--- Moisture_regime Model Summary ---")
      print(model, outcome = "Moisture_regime")
      message("\n--- Nutrient_regime Model Summary ---")
      
      ##------------------------------------------------------------------------------
      ## Create variable importance plots for moisture and nutrient regime separately
      ##------------------------------------------------------------------------------
   
      # Variable importance dataframes
      moist_vi <- data.frame(
        Variable = names(model$regrOutput$Moisture_regime$importance),
        Importance = model$regrOutput$Moisture_regime$importance
      )
      nut_vi <- data.frame(
        Variable = names(model$regrOutput$Nutrient_regime$importance),
        Importance = model$regrOutput$Nutrient_regime$importance
      )
      
      # Reduce text size slightly for side-by-side plotting
      text_size_title <- 12
      text_size_axis_title <- 10
      text_size_axis_text <- 8
      main_title_size <- 16
      
      # Plot variable importance for Moisture regime
      p1 <- ggplot(moist_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        labs(title = paste0("Variable Importance for Moisture regime"),
             x = "Variable", y = "Importance") +
        theme_minimal(base_family = "Arial") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = text_size_title),
              axis.title = element_text(size = text_size_axis_title),
              axis.text = element_text(size = text_size_axis_text))
      
      # Plot variable importance for Nutrient regime
      p2 <- ggplot(nut_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
        geom_col(fill = "darkgreen") +
        coord_flip() +
        labs(title = paste0("Variable Importance for Nutrient regime"),
             x = "Variable", y = "Importance") +
        theme_minimal(base_family = "Arial") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = text_size_title),
              axis.title = element_text(size = text_size_axis_title),
              axis.text = element_text(size = text_size_axis_text))
      
      # Save side-by-side plots 
      png_filename <- paste0("Ecosite_output/Thresholds/VarImp_", gsub(" ", "_", nsr_name), ".png")
      png(png_filename, width = 1600, height = 900, res = 150)  # slightly taller for main title
      grid.arrange(p1, p2, ncol = 2, top = textGrob(
        paste(nsr_name),
        gp = gpar(fontsize = main_title_size, fontface = "bold", fontfamily = "Arial")
      ))
      dev.off()
      
      message("Variable importance plots saved for ", nsr_name)
      
      ##-----------------------------------------------------------------------------------------
      ## Create partial dependence plots of top 3 variables for moisture and nutrient separately
      ##-----------------------------------------------------------------------------------------
      
      # Helper to save one plot.variable output to a PNG file
      save_partial_plot <- function(model, var, outcome, filename) {
        png(filename, width = 700, height = 700)
        
        par(
          mar = c(5, 5, 5, 2) + 0.1,
          cex.lab = 1.6,   # axis labels
          cex.axis = 1.4,  # axis tick labels
          cex.main = 1.8   # title
        )
        
        plot.variable(
          model,
          xvar.names = var,
          outcome = outcome,
          partial = TRUE,
          show.plots = TRUE,
          main = paste0("Partial Dependence of ", outcome, " on ", var),
          xlab = var,
          ylab = paste0("Predicted ", outcome)
        )
        
        dev.off()
      }
      # Get top 3 variables for each outcome
      top3_moist <- moist_vi %>%
        arrange(desc(Importance)) %>%
        slice(1:3) %>%
        pull(Variable)
      
      top3_nut <- nut_vi %>%
        arrange(desc(Importance)) %>%
        slice(1:3) %>%
        pull(Variable)
      
      # Folder for temp partial plots
      tmp_dir <- tempdir()
      
      # Filenames for partial plots
      pd_files <- c()
      
      # Save Moisture regime partial dependence plots
      for (i in seq_along(top3_moist)) {
        fn <- file.path(tmp_dir, paste0("partial_moist_", i, ".png"))
        save_partial_plot(model, top3_moist[i], "Moisture_regime", fn)
        pd_files <- c(pd_files, fn)
      }
      
      # Save Nutrient regime partial dependence plots
      for (i in seq_along(top3_nut)) {
        fn <- file.path(tmp_dir, paste0("partial_nut_", i, ".png"))
        save_partial_plot(model, top3_nut[i], "Nutrient_regime", fn)
        pd_files <- c(pd_files, fn)
      }
      
      # Read images back as grobs
      plot_grobs <- lapply(pd_files, function(f) {
        img <- readPNG(f)
        rasterGrob(img, interpolate = TRUE)
      })
      
      # Output combined 6-panel figure with shared title
      pd_png_filename <- paste0("Ecosite_output/Thresholds/PartialDependence_", gsub(" ", "_", nsr_name), ".png")
      png(pd_png_filename, width = 1800, height = 900)
      
      grid.arrange(
        grobs = plot_grobs,
        nrow = 2, ncol = 3,
        top = grid::textGrob(
          paste0("Partial Dependence Plots of top 3 variables for ", nsr_name),
          gp = grid::gpar(fontsize = 20, fontface = "bold")
        )
      )
      
      dev.off()
      message("Partial dependence plots saved for ", nsr_name)
      
      ##------------------------------------------
      ## Prediction and raster creation 
      ##------------------------------------------
      
      pred_df <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
      pred_result <- predict(model, newdata = pred_df, na.action = "na.omit")
      
      pred_df$Moist_pred <- pred_result$regrOutput$Moisture_regime$predicted
      pred_df$Nut_pred <- pred_result$regrOutput$Nutrient_regime$predicted
      
      moist_rast <- rast(pred_df[, c("x", "y", "Moist_pred")], type = "xyz")
      nut_rast <- rast(pred_df[, c("x", "y", "Nut_pred")], type = "xyz")
      
      writeRaster(moist_rast, paste0("Ecosite_output/Thresholds/Moisture_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      writeRaster(nut_rast, paste0("Ecosite_output/Thresholds/Nutrient_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      
      
      ##------------------------------------------------------------
      ## Build classification breaks from optimal thresholds
      ##------------------------------------------------------------
      
      
      # Nutrient cutpoints → sorted break vector
      nutrient_breaks <- c(-Inf,
                           nut_cut_table$optimal_cut,
                           Inf)
      nutrient_labels <- c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
      
      nut_rast_cat <- classify(
        nut_rast,
        rcl = cbind(nutrient_breaks[-length(nutrient_breaks)],
                    nutrient_breaks[-1],
                    1:length(nutrient_labels))
      )

      
      # Moisture cutpoints → sorted break vector
      moisture_breaks <- c(-Inf,
                           moist_cut_table$optimal_cut,
                           Inf)
      moisture_labels <- c("Xeric", "Subxeric", "Submesic", "Mesic",
                           "Subhygric", "Hygric", "Subhydric", "Hydric")
      
      moist_rast_cat <- classify(
        moist_rast,
        rcl = cbind(moisture_breaks[-length(moisture_breaks)],
                    moisture_breaks[-1],
                    1:length(moisture_labels))
      )
  
      # --- Combine to create Ecosite numeric code (nutrient first) ---
      ecosite_rast <- (nut_rast_cat * 10) + moist_rast_cat
      names(ecosite_rast) <- paste0("Ecosite_", gsub(" ", "_", nsr_name))
      
      writeRaster(ecosite_rast,
                  paste0("Ecosite_output/Thresholds/Ecosite_", nsr_name, ".tif"),
                  overwrite = TRUE)
      
      
      # --- Save lookup table (nutrient first) ---
      ecosite_lookup <- expand.grid(
        Nutrient = nutrient_labels,
        Moisture = moisture_labels
      ) %>%
        mutate(
          Ecosite_Code = (as.numeric(factor(Nutrient, levels = nutrient_labels)) * 10) +
            as.numeric(factor(Moisture, levels = moisture_labels)),
          Ecosite_Class = paste0(Nutrient, "_", Moisture)
        )
      
      write.csv(ecosite_lookup,
                paste0("Ecosite_output/Thresholds/Ecosite_Lookup_", gsub(" ", "_", nsr_name), ".csv"),
                row.names = FALSE)
      
      message("Finished processing ", nsr_name, " — Ecosite rasters and lookup saved.")
      
    }, error = function(e) {
      message("Error processing ", nsr_name, ": ", e$message)
    })
  }
}

##--------------------------------------------------
##Improved memory function for large subregions
##--------------------------------------------------

NSR7 <- subset(data, NSRNAME %in% c("Central Mixedwood"))

NSR8 <- subset(data, NSRNAME %in% c("Dry Mixedwood"))##Do one at a time or memory crashes

for (nsr_name in unique(NSR7$NSRNAME)) {    ## Change to NSR8 for Dry Mixedwood
  message("Processing NSR: ", nsr_name)
  
  tryCatch({
    ##------------------------------------------------
    ## Stack rasters and extract data to input points
    ##------------------------------------------------
    
    # Subset data for the NSR and create spatial points
    nsr_data <- subset(NSR7, NSRNAME == nsr_name)   ## Change to NSR8 for Dry Mixedwood
    nsr_shape <- NSR[NSR$NSRNAME == nsr_name]
    
    # Create ID column (stable row index for joining)
    nsr_data$ID <- seq_len(nrow(nsr_data))
    
    nsr_points <- vect(nsr_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")
    
    # Mask and trim rasters to NSR boundary
    masked <- lapply(rasters, function(r) (mask(r, nsr_shape)))
    
    
    # Create raster stack and name layers consistently
    stack_rasters <- rast(masked)
    names(stack_rasters) <- predictors
    
    stack_rasters
    
    # Extract raster values at data points
    raster_values <- terra::extract(stack_rasters, nsr_points)  
    nsr_combined <- left_join(nsr_data, raster_values, by = "ID")
    
    # Assign numeric regime values
    nsr_combined <- assign_regime_values(nsr_combined)
    
    # Drop columns with all NA and convert character columns to factors
    
    nsr_combined <- nsr_combined %>%
      mutate(across(where(is.character), as.factor))
    
    # Skip if essential numeric columns missing
    if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(nsr_combined))) {
      message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
      next
    }
    
    
    # Check missing predictors
    setdiff(predictors, names(nsr_combined))
    
    ##------------------------------------------------
    ## Fit random forest model
    ##------------------------------------------------
    
    
    model <- rfsrc(rf_formula, data = nsr_combined,
                   importance = "permute", na.action = "na.impute")
    
    point_pred <- predict(model)$regrOutput
    
    nsr_combined$Moist_pred_point    <- point_pred$Moisture_regime$predicted
    nsr_combined$Nutrient_pred_point <- point_pred$Nutrient_regime$predicted
    
    
    # Save point predictions
    write.csv(
      nsr_combined,
      paste0("Ecosite_output/Thresholds/PointPredictions_", gsub(" ", "_", nsr_name), ".csv"),
      row.names = FALSE
    )
    message("Point predictions saved for ", nsr_name)
    
    ##--------------------------------------------- 
    ## Compute optimal thresholds
    ##---------------------------------------------
    
    cuts <- get_cutpoints(nsr_combined)
    
    moist_cut_table <- fix_cutpoints(cuts$Moisture)
    nut_cut_table   <- fix_cutpoints(cuts$Nutrient)
    
    # Save results
    write.csv(
      moist_cut_table,
      paste0("Ecosite_output/Thresholds/Cutpoints_Moisture_", gsub(" ", "_", nsr_name), ".csv"),
      row.names = FALSE
    )
    
    write.csv(
      nut_cut_table,
      paste0("Ecosite_output/Thresholds/Cutpoints_Nutrient_", gsub(" ", "_", nsr_name), ".csv"),
      row.names = FALSE
    )
    
    message("Cutpoints saved for ", nsr_name)
    
    # Print summaries
    message("\n--- Moisture_regime Model Summary ---")
    print(model, outcome = "Moisture_regime")
    message("\n--- Nutrient_regime Model Summary ---")
    print(model, outcome = "Nutrient_regime")
    
    ##------------------------------------------------------------------------------
    ## Create variable importance plots for moisture and nutrient regime separately
    ##------------------------------------------------------------------------------
    
    # Variable importance dataframes
    moist_vi <- data.frame(
      Variable = names(model$regrOutput$Moisture_regime$importance),
      Importance = model$regrOutput$Moisture_regime$importance
    )
    nut_vi <- data.frame(
      Variable = names(model$regrOutput$Nutrient_regime$importance),
      Importance = model$regrOutput$Nutrient_regime$importance
    )
    
    # Reduce text size slightly for side-by-side plotting
    text_size_title <- 12
    text_size_axis_title <- 10
    text_size_axis_text <- 8
    main_title_size <- 16
    
    # Plot variable importance for Moisture regime
    p1 <- ggplot(moist_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(title = paste0("Variable Importance for Moisture regime"),
           x = "Variable", y = "Importance") +
      theme_minimal(base_family = "Arial") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = text_size_title),
            axis.title = element_text(size = text_size_axis_title),
            axis.text = element_text(size = text_size_axis_text))
    
    # Plot variable importance for Nutrient regime
    p2 <- ggplot(nut_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_col(fill = "darkgreen") +
      coord_flip() +
      labs(title = paste0("Variable Importance for Nutrient regime"),
           x = "Variable", y = "Importance") +
      theme_minimal(base_family = "Arial") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = text_size_title),
            axis.title = element_text(size = text_size_axis_title),
            axis.text = element_text(size = text_size_axis_text))
    
    # Save side-by-side plots 
    png_filename <- paste0("Ecosite_output/Thresholds/VarImp_", gsub(" ", "_", nsr_name), ".png")
    png(png_filename, width = 1600, height = 900, res = 150)  # slightly taller for main title
    grid.arrange(p1, p2, ncol = 2, top = textGrob(
      paste(nsr_name),
      gp = gpar(fontsize = main_title_size, fontface = "bold", fontfamily = "Arial")
    ))
    dev.off()
    
    message("Variable importance plots saved for ", nsr_name)
    
    ##-----------------------------------------------------------------------------------------
    ## Create partial dependence plots of top 3 variables for moisture and nutrient separately
    ##-----------------------------------------------------------------------------------------
    
    # Helper to save one plot.variable output to a PNG file
    save_partial_plot <- function(model, var, outcome, filename) {
      png(filename, width = 700, height = 700)
      
      par(
        mar = c(5, 5, 5, 2) + 0.1,
        cex.lab = 1.6,   # axis labels
        cex.axis = 1.4,  # axis tick labels
        cex.main = 1.8   # title
      )
      
      plot.variable(
        model,
        xvar.names = var,
        outcome = outcome,
        partial = TRUE,
        show.plots = TRUE,
        main = paste0("Partial Dependence of ", outcome, " on ", var),
        xlab = var,
        ylab = paste0("Predicted ", outcome)
      )
      
      dev.off()
    }
    # Get top 3 variables for each outcome
    top3_moist <- moist_vi %>%
      arrange(desc(Importance)) %>%
      slice(1:3) %>%
      pull(Variable)
    
    top3_nut <- nut_vi %>%
      arrange(desc(Importance)) %>%
      slice(1:3) %>%
      pull(Variable)
    
    # Folder for temp partial plots
    tmp_dir <- tempdir()
    
    # Filenames for partial plots
    pd_files <- c()
    
    # Save Moist_num partial dependence plots
    for (i in seq_along(top3_moist)) {
      fn <- file.path(tmp_dir, paste0("partial_moist_", i, ".png"))
      save_partial_plot(model, top3_moist[i], "Moisture_regime", fn)
      pd_files <- c(pd_files, fn)
    }
    
    # Save Nut_num partial dependence plots
    for (i in seq_along(top3_nut)) {
      fn <- file.path(tmp_dir, paste0("partial_nut_", i, ".png"))
      save_partial_plot(model, top3_nut[i], "Nutrient_regime", fn)
      pd_files <- c(pd_files, fn)
    }
    
    # Read images back as grobs
    plot_grobs <- lapply(pd_files, function(f) {
      img <- readPNG(f)
      rasterGrob(img, interpolate = TRUE)
    })
    
    # Output combined 6-panel figure with shared title
    pd_png_filename <- paste0("Ecosite_output/Thresholds/PartialDependence_", gsub(" ", "_", nsr_name), ".png")
    png(pd_png_filename, width = 1800, height = 900)
    
    grid.arrange(
      grobs = plot_grobs,
      nrow = 2, ncol = 3,
      top = grid::textGrob(
        paste0("Partial Dependence Plots of top 3 variables for ", nsr_name),
        gp = grid::gpar(fontsize = 20, fontface = "bold")
      )
    )
    
    dev.off()
    message("Partial dependence plots saved for ", nsr_name)
    
    ##------------------------------------------
    ## Prediction and raster creation 
    ##------------------------------------------
    
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
      
      moist_vals[i:i_end] <- pred_chunk$regrOutput$Moisture_regime$predicted
      nut_vals[i:i_end] <- pred_chunk$regrOutput$Nutrient_regime$predicted
    }
    
    # Attach predictions to full df
    pred_df_all$Moist_pred <- moist_vals
    pred_df_all$Nut_pred <- nut_vals
    
    # Convert to rasters
    moist_rast <- rast(pred_df_all[, c("x", "y", "Moist_pred")], type = "xyz")
    nut_rast <- rast(pred_df_all[, c("x", "y", "Nut_pred")], type = "xyz")
    
    writeRaster(moist_rast, paste0("Ecosite_output/Thresholds/Moist_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    writeRaster(nut_rast, paste0("Ecosite_output/Thresholds/Nut_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    
    ##------------------------------------------------------------
    ## Build classification breaks from optimal thresholds
    ##------------------------------------------------------------
  
    # Nutrient cutpoints → sorted break vector
    nutrient_breaks <- c(-Inf,
                         nut_cut_table$optimal_cut,
                         Inf)
    nutrient_labels <- c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
    
    nut_rast_cat <- classify(
      nut_rast,
      rcl = cbind(nutrient_breaks[-length(nutrient_breaks)],
                  nutrient_breaks[-1],
                  1:length(nutrient_labels))
    )
    
    
    # Moisture cutpoints → sorted break vector
    moisture_breaks <- c(-Inf,
                         moist_cut_table$optimal_cut,
                         Inf)
    moisture_labels <- c("Xeric", "Subxeric", "Submesic", "Mesic",
                         "Subhygric", "Hygric", "Subhydric", "Hydric")
    
    moist_rast_cat <- classify(
      moist_rast,
      rcl = cbind(moisture_breaks[-length(moisture_breaks)],
                  moisture_breaks[-1],
                  1:length(moisture_labels))
    )
    
    # --- Combine to create Ecosite numeric code (nutrient first) ---
    ecosite_rast <- (nut_rast_cat * 10) + moist_rast_cat
    names(ecosite_rast) <- paste0("Ecosite_", gsub(" ", "_", nsr_name))
    
    writeRaster(ecosite_rast,
                paste0("Ecosite_output/Thresholds/Ecosite_", nsr_name, ".tif"),
                overwrite = TRUE)
    
    
    # --- Save lookup table (nutrient first) ---
    ecosite_lookup <- expand.grid(
      Nutrient = nutrient_labels,
      Moisture = moisture_labels
    ) %>%
      mutate(
        Ecosite_Code = (as.numeric(factor(Nutrient, levels = nutrient_labels)) * 10) +
          as.numeric(factor(Moisture, levels = moisture_labels)),
        Ecosite_Class = paste0(Nutrient, "_", Moisture)
      )
    
    write.csv(ecosite_lookup,
              paste0("Ecosite_output/Thresholds/Ecosite_Lookup_", gsub(" ", "_", nsr_name), ".csv"),
              row.names = FALSE)
    
    message("Finished processing ", nsr_name, " — Ecosite rasters and lookup saved.")
    
  }, error = function(e) {
    message("Error processing ", nsr_name, ": ", e$message)
  })
}



##Combine Athabasca Plain and Peace-Athabasca Delta (only 9 points in Delta, no model predictions made)

NSR9 <- subset(data, NSRNAME %in% c("Athabasca Plain","Peace-Athabasca Delta"))
View(NSR9)



process_combined_nsr <- function(nsr_name, combined_subregions, NSR9, NSR, rasters,
                                 predictors, rf_formula) {
  tryCatch({
    message("Processing Combined NSR: ", nsr_name)
    
    ##------------------------------------------------
    ## Stack rasters and extract data to input points
    ##------------------------------------------------
    combined_data <- subset(NSR9, NSRNAME %in% combined_subregions)
    if (nrow(combined_data) == 0) {
      message("No data found for ", nsr_name)
      return(NULL)
    }
    
    # Create stable ID
    combined_data$ID <- seq_len(nrow(combined_data))
    
    # Convert to spatial points
    combined_points <- vect(
      combined_data, 
      geom = c("Easting", "Northing"), 
      crs = "EPSG:3400"
    )
    
    #Dissolve multiple polygons into one
    combined_shape <- NSR[NSR$NSRNAME %in% combined_subregions, ]
    combined_shape <- aggregate(combined_shape)   # dissolve
    
    #Mask and stack rasters
    masked <- lapply(rasters, function(r) mask(r, combined_shape))
    
    
    stack_rasters <- rast(masked)
    names(stack_rasters) <- predictors
    
    #Extract raster values
    raster_values <- terra::extract(stack_rasters, combined_points)
    
    # Left join by ID so missing points keep attribute data
    combined_df <- dplyr::left_join(combined_data, raster_values, by = "ID")
    
    #Assign numeric values
    combined_df <- assign_regime_values(combined_df)
    
    # Convert character → factor
    combined_df <- combined_df %>%
      dplyr::mutate(across(where(is.character), as.factor))
    
    #Check required variables
    if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(combined_df))) {
      message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
      return(NULL)
    }
    
    ##------------------------------------------------
    ## Fit random forest model
    ##------------------------------------------------
    model <- rfsrc(
      rf_formula,
      data = combined_df,
      importance = "permute",
      na.action = "na.impute"
    )
    
   
    point_pred <- predict(model)$regrOutput
    
    combined_df$Moist_pred_point    <- point_pred$Moisture_regime$predicted
    combined_df$Nutrient_pred_point <- point_pred$Nutrient_regime$predicted

    out_file <- paste0("Ecosite_output/Thresholds/PointPredictions_",
                       gsub(" ", "_", nsr_name), ".csv")
    write.csv(combined_df, out_file, row.names = FALSE)
    message("Point predictions saved for ", nsr_name)
    
    ##--------------------------------------------- 
    ## Compute optimal thresholds
    ##---------------------------------------------
    
    cuts <- get_cutpoints(nsr_combined)
    
    moist_cut_table <- fix_cutpoints(cuts$Moisture)
    nut_cut_table   <- fix_cutpoints(cuts$Nutrient)
    
    # Save results
    write.csv(
      moist_cut_table,
      paste0("Ecosite_output/Thresholds/Cutpoints_Moisture_", gsub(" ", "_", nsr_name), ".csv"),
      row.names = FALSE
    )
    
    write.csv(
      nut_cut_table,
      paste0("Ecosite_output/Thresholds/Cutpoints_Nutrient_", gsub(" ", "_", nsr_name), ".csv"),
      row.names = FALSE
    )
    
    message("Cutpoints saved for ", nsr_name)
    
    # Print summaries
    message("\n--- Moisture_regime Model Summary ---")
    print(model, outcome = "Moisture_regime")
    message("\n--- Nutrient_regime Model Summary ---")
    print(model, outcome = "Nutrient_regime")
    
    ##------------------------------------------------------------------------------
    ## Create variable importance plots for moisture and nutrient regime separately
    ##------------------------------------------------------------------------------
    
    # Variable importance dataframes
    moist_vi <- data.frame(
      Variable = names(model$regrOutput$Moisture_regime$importance),
      Importance = model$regrOutput$Moisture_regime$importance
    )
    nut_vi <- data.frame(
      Variable = names(model$regrOutput$Nutrient_regime$importance),
      Importance = model$regrOutput$Nutrient_regime$importance
    )
    
    # Reduce text size slightly for side-by-side
    text_size_title <- 12
    text_size_axis_title <- 10
    text_size_axis_text <- 8
    main_title_size <- 16
    
    # Plot variable importance for Moist_num
    p1 <- ggplot(moist_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(title = paste0("Variable Importance for Moisture regime"),
           x = "Variable", y = "Importance") +
      theme_minimal(base_family = "Arial") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = text_size_title),
            axis.title = element_text(size = text_size_axis_title),
            axis.text = element_text(size = text_size_axis_text))
    
    # Plot variable importance for Nut_num
    p2 <- ggplot(nut_vi, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_col(fill = "darkgreen") +
      coord_flip() +
      labs(title = paste0("Variable Importance for Nutrient regime"),
           x = "Variable", y = "Importance") +
      theme_minimal(base_family = "Arial") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold", size = text_size_title),
            axis.title = element_text(size = text_size_axis_title),
            axis.text = element_text(size = text_size_axis_text))
    
    # Save side-by-side plots 
    png_filename <- paste0("Ecosite_output/Thresholds/VarImp_", gsub(" ", "_", nsr_name), ".png")
    png(png_filename, width = 1600, height = 900, res = 150)  # slightly taller for main title
    grid.arrange(p1, p2, ncol = 2, top = textGrob(
      paste(nsr_name),
      gp = gpar(fontsize = main_title_size, fontface = "bold", fontfamily = "Arial")
    ))
    dev.off()
    
    message("Variable importance plots saved for ", nsr_name)
    
    ##-----------------------------------------------------------------------------------------
    ## Create partial dependence plots of top 3 variables for moisture and nutrient separately
    ##-----------------------------------------------------------------------------------------
    
    # Helper to save one plot.variable output to a PNG file
    save_partial_plot <- function(model, var, outcome, filename) {
      png(filename, width = 700, height = 700)
      
      par(
        mar = c(5, 5, 5, 2) + 0.1,
        cex.lab = 1.6,   # axis labels
        cex.axis = 1.4,  # axis tick labels
        cex.main = 1.8   # title
      )
      
      plot.variable(
        model,
        xvar.names = var,
        outcome = outcome,
        partial = TRUE,
        show.plots = TRUE,
        main = paste0("Partial Dependence of ", outcome, " on ", var),
        xlab = var,
        ylab = paste0("Predicted ", outcome)
      )
      
      dev.off()
    }
    # Get top 3 variables for each outcome
    top3_moist <- moist_vi %>%
      arrange(desc(Importance)) %>%
      slice(1:3) %>%
      pull(Variable)
    
    top3_nut <- nut_vi %>%
      arrange(desc(Importance)) %>%
      slice(1:3) %>%
      pull(Variable)
    
    # Folder for temp partial plots
    tmp_dir <- tempdir()
    
    # Filenames for partial plots
    pd_files <- c()
    
    # Save Moist_num partial dependence plots
    for (i in seq_along(top3_moist)) {
      fn <- file.path(tmp_dir, paste0("partial_moist_", i, ".png"))
      save_partial_plot(model, top3_moist[i], "Moisture_regime", fn)
      pd_files <- c(pd_files, fn)
    }
    
    # Save Nut_num partial dependence plots
    for (i in seq_along(top3_nut)) {
      fn <- file.path(tmp_dir, paste0("partial_nut_", i, ".png"))
      save_partial_plot(model, top3_nut[i], "Nutrient_regime", fn)
      pd_files <- c(pd_files, fn)
    }
    
    # Read images back as grobs
    plot_grobs <- lapply(pd_files, function(f) {
      img <- readPNG(f)
      rasterGrob(img, interpolate = TRUE)
    })
    
    # Output combined 6-panel figure with shared title
    pd_png_filename <- paste0("Ecosite_output/Thresholds/PartialDependence_", gsub(" ", "_", nsr_name), ".png")
    png(pd_png_filename, width = 1800, height = 900)
    
    grid.arrange(
      grobs = plot_grobs,
      nrow = 2, ncol = 3,
      top = grid::textGrob(
        paste0("Partial Dependence Plots of top 3 variables for ", nsr_name),
        gp = grid::gpar(fontsize = 20, fontface = "bold")
      )
    )
    
    dev.off()
    message("Partial dependence plots saved for ", nsr_name)
    
    ##------------------------------------------
    ## Prediction and raster creation 
    ##------------------------------------------
    
    pred_df <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
    pred_result <- predict(model, newdata = pred_df, na.action = "na.omit")
    
    pred_df$Moist_pred <- pred_result$regrOutput$Moisture_regime$predicted
    pred_df$Nut_pred <- pred_result$regrOutput$Nutrient_regime$predicted
    
    moist_rast <- rast(pred_df[, c("x", "y", "Moist_pred")], type = "xyz")
    nut_rast <- rast(pred_df[, c("x", "y", "Nut_pred")], type = "xyz")
    
    writeRaster(moist_rast, paste0("Ecosite_output/Thresholds/Moisture_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    writeRaster(nut_rast, paste0("Ecosite_output/Thresholds/Nutrient_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    
    ##------------------------------------------------------------
    ## Build classification breaks from optimal thresholds
    ##------------------------------------------------------------
    
    # Nutrient cutpoints → sorted break vector
    nutrient_breaks <- c(-Inf,
                         nut_cut_table$optimal_cut,
                         Inf)
    nutrient_labels <- c("Very Poor", "Poor", "Medium", "Rich", "Very Rich")
    
    nut_rast_cat <- classify(
      nut_rast,
      rcl = cbind(nutrient_breaks[-length(nutrient_breaks)],
                  nutrient_breaks[-1],
                  1:length(nutrient_labels))
    )
    
    
    # Moisture cutpoints → sorted break vector
    moisture_breaks <- c(-Inf,
                         moist_cut_table$optimal_cut,
                         Inf)
    moisture_labels <- c("Xeric", "Subxeric", "Submesic", "Mesic",
                         "Subhygric", "Hygric", "Subhydric", "Hydric")
    
    moist_rast_cat <- classify(
      moist_rast,
      rcl = cbind(moisture_breaks[-length(moisture_breaks)],
                  moisture_breaks[-1],
                  1:length(moisture_labels))
    )
    
    # --- Combine to create Ecosite numeric code (nutrient first) ---
    ecosite_rast <- (nut_rast_cat * 10) + moist_rast_cat
    names(ecosite_rast) <- paste0("Ecosite_", gsub(" ", "_", nsr_name))
    
    writeRaster(ecosite_rast,
                paste0("Ecosite_output/Thresholds/Ecosite_", nsr_name, ".tif"),
                overwrite = TRUE)
    
    
    # --- Save lookup table (nutrient first) ---
    ecosite_lookup <- expand.grid(
      Nutrient = nutrient_labels,
      Moisture = moisture_labels
    ) %>%
      mutate(
        Ecosite_Code = (as.numeric(factor(Nutrient, levels = nutrient_labels)) * 10) +
          as.numeric(factor(Moisture, levels = moisture_labels)),
        Ecosite_Class = paste0(Nutrient, "_", Moisture)
      )
    
    write.csv(ecosite_lookup,
              paste0("Ecosite_output/Thresholds/Ecosite_Lookup_", gsub(" ", "_", nsr_name), ".csv"),
              row.names = FALSE)
    
    message("Finished processing ", nsr_name, " — Ecosite rasters and lookup saved.")
    
  }, error = function(e) {
    message("Error processing ", nsr_name, ": ", e$message)
  })
}


# Define combined subregion name and members
nsr_name <- "Athabasca_Plain__Peace_Athabasca_Delta"
combined_subregions <- c("Athabasca Plain", "Peace-Athabasca Delta")


##Run function
process_combined_nsr(
  nsr_name = "Athabasca_Plain__Peace_Athabasca_Delta",
  combined_subregions = c("Athabasca Plain", "Peace-Athabasca Delta"),
  NSR9 = NSR9,
  NSR = NSR,
  rasters = rasters,
  predictors = predictors,
  rf_formula = rf_formula
)

