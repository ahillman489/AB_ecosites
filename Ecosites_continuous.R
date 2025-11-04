library(terra)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)


data <- read.csv("C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/All_sites_updated.csv")

data$DataSource <- as.factor(data$DataSource)
data$Nutrient.Regime <- as.factor(data$Nutrient.Regime)
data$Moisture.Regime <- as.factor(data$Moisture.Regime)
data$Ecosite <- as.factor(data$Ecosite)
data$NSRNAME <- as.factor(data$NSRNAME)
data$Site <- as.factor(data$Site)


##Remove Open Water observations, plus swamp, alkali, and not treed for exploratory analyses
data1 <- subset(data, Nutrient.Regime !='Open Water'& Nutrient.Regime!= 'Alkali'& Nutrient.Regime!= 'Swamp'
                & Nutrient.Regime!= 'Not Treed'& Moisture.Regime!="Swamp"& Moisture.Regime!="No Value Provided"&
                  Nutrient.Regime!="No Value Provided" & Nutrient.Regime!=" " & Moisture.Regime !=" ")

# Load spatial and raster data only once
NSR <- vect("W:/EDM/CCVegMod4/Ashley_ecosites/Alberta_Ecosites/GIS/Natural_Sub_Regions_Alberta_2005_shapefile.shp")

# Load raster layers
raster_paths <- list(
  Proportion_sand <- "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_sand_prop.tif",
  Proportion_clay <- "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_clay_prop.tif",
  Proportion_silt <- "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_silt_prop.tif",
  Soil_organic_carbon <- "C:/Users/ahillman/OneDrive - NRCan RNCan/Alberta Ecosites/SG_SOC.tif",
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

# List of NSR names to loop over
nsr_list <- unique(data1$NSRNAME)

# Predictor variable names (must match raster and data)
predictors <- c("Proportion_sand", "Proportion_clay", "Proportion_silt","Soil_organic_carbon", "Surficial_geology",
                "Slope","Elevation","CTI", "TPI", "Profile_curvature", "Tangential_curvature","pH")


# Build formula safely
formula_str <- paste("Multivar(Moisture_regime, Nutrient_regime) ~", paste(predictors, collapse = " + "))
rf_formula <- as.formula(formula_str)
rf_formula

# Create output directory if it doesn't exist
dir.create("Ecosite_output/No_climate", showWarnings = FALSE, recursive = TRUE)

##Split data into chunks for processing memory
#NSR1 <- subset(data1, NSRNAME=="Lower Boreal Highlands")
NSR1 <- subset(data1, NSRNAME %in% c("Upper Boreal Highlands","Boreal Subarctic","Lower Boreal Highlands"))
NSR2 <- subset(data1, NSRNAME %in% c("Northern Mixedwood","Kazan Uplands", "Peace River Parkland"))
NSR3 <- subset(data1, NSRNAME %in% c("Foothills Parkland","Lower Foothills","Upper Foothills"))
NSR4 <- subset(data1, NSRNAME %in% c("Alpine","Subalpine","Montane"))
##NSR5 <- subset(data1, NSRNAME %in% c("Central Parkland","Dry Mixedgrass","Foothills Fescue","Mixedgrass","Northern Fescue"))

# Combine into a list of subsets
nsr_groups <- list(NSR1)
nsr_groups <- list(NSR1,NSR2,NSR3,NSR4)


# Loop through each group
for (nsr_set in nsr_groups) {
  for (nsr_name in unique(nsr_set$NSRNAME)) {
    message("Processing NSR: ", nsr_name)

    tryCatch({
      # Subset data for the NSR and create spatial points
      nsr_data <- subset(nsr_set, NSRNAME == nsr_name)
      nsr_points <- vect(nsr_data, geom = c("Easting", "Northing"), crs = "EPSG:3400")

      # Subset NSR shapefile
      nsr_shape <- NSR[NSR$NSRNAME == nsr_name]

      # Mask and trim rasters to NSR boundary
      masked <- lapply(rasters, function(r) trim(mask(r, nsr_shape)))

      # Align all rasters to the first raster's geometry
      template <- masked[[1]]
      masked_aligned <- lapply(masked, function(r) {
        if (!compareGeom(r, template, stopOnError = FALSE)) {
          resample(r, template, method = "bilinear")
        } else {
          r
        }
      })

      # Create raster stack and name layers consistently
      stack_rasters <- rast(masked_aligned)
      names(stack_rasters) <- predictors

      # Extract raster values at data points
      raster_values <- terra::extract(stack_rasters, nsr_points)[, -1]  # drop ID column
      nsr_combined <- cbind(nsr_data, raster_values)

      # Assign numeric regime values
      nsr_combined <- assign_regime_values(nsr_combined)

      # Drop columns with all NA and convert character columns to factors
      nsr_combined <- nsr_combined[, colSums(is.na(nsr_combined)) < nrow(nsr_combined)]
      nsr_combined <- nsr_combined %>%
        mutate(across(where(is.character), as.factor))

      # Skip if essential numeric columns missing
      if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(nsr_combined))) {
        message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
        next
      }

      # Fit multivariate random forest model
      model <- rfsrc(rf_formula, data = nsr_combined, importance = "permute")

      # Print summaries
      message("\n--- Moisture_regime Model Summary ---")
      print(model, outcome = "Moisture_regime")
      message("\n--- Nutrient_regime Model Summary ---")
      print(model, outcome = "Nutrient_regime")

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

      # Save side-by-side plots with larger width
      # Save side-by-side plots with main title
      png_filename <- paste0("Ecosite_output/No_climate/VarImp_", gsub(" ", "_", nsr_name), ".png")
      png(png_filename, width = 1600, height = 900, res = 150)  # slightly taller for main title
      grid.arrange(p1, p2, ncol = 2, top = textGrob(
        paste(nsr_name),
        gp = gpar(fontsize = main_title_size, fontface = "bold", fontfamily = "Arial")
      ))
      dev.off()

      message("Variable importance plots saved for ", nsr_name)

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
      pd_png_filename <- paste0("Ecosite_output/No_climate/PartialDependence_", gsub(" ", "_", nsr_name), ".png")
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

      # --- Prediction and raster creation ---
      pred_df <- as.data.frame(stack_rasters, xy = TRUE, na.rm = TRUE)
      pred_result <- predict(model, newdata = pred_df, na.action = "na.omit")

      pred_df$Moist_pred <- pred_result$regrOutput$Moisture_regime$predicted
      pred_df$Nut_pred <- pred_result$regrOutput$Nutrient_regime$predicted

      moist_rast <- rast(pred_df[, c("x", "y", "Moist_pred")], type = "xyz")
      nut_rast <- rast(pred_df[, c("x", "y", "Nut_pred")], type = "xyz")

      writeRaster(moist_rast, paste0("Ecosite_output/No_climate/Moisture_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
      writeRaster(nut_rast, paste0("Ecosite_output/No_climate/Nutrient_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)

      # --- Classification ---
      moist_breaks <- quantile(values(moist_rast), probs = seq(0, 1, length.out = 9), na.rm = TRUE)
      moist_cat <- cbind(moist_breaks[-length(moist_breaks)], moist_breaks[-1], 1:8)
      moist_rast_cat <- classify(moist_rast, moist_cat)
      levels(moist_rast_cat) <- data.frame(
        ID = 1:8,
        Category = c("Xeric", "Subxeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
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

      writeRaster(ecosite_rast, paste0("Ecosite_output/No_climate/Ecosite_", nsr_name, ".tif"), overwrite = TRUE)

      message("Finished processing ", nsr_name)

    }, error = function(e) {
      message("Error processing ", nsr_name, ": ", e$message)
    })
  }
}


##Improve memory for large subregions

NSR7 <- subset(data1, NSRNAME %in% c("Central Mixedwood"))
View(NSR7)
NSR8 <- subset(data1, NSRNAME %in% c("Dry Mixedwood"))##Do one at a time or memory crashes

for (nsr_name in unique(NSR8$NSRNAME)) {
  message("Processing NSR: ", nsr_name)

  tryCatch({

    # Subset and convert to spatial
    nsr_data <- subset(NSR8, NSRNAME == nsr_name)
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


    # Skip if essential numeric columns missing
    if (!all(c("Moisture_regime", "Nutrient_regime") %in% names(nsr_combined))) {
      message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
      next
    }

    # Fit multivariate random forest model
    model <- rfsrc(rf_formula, data = nsr_combined, importance = "permute")

    # Print summaries
    message("\n--- Moisture_regime Model Summary ---")
    print(model, outcome = "Moisture_regime")
    message("\n--- Nutrient_regime Model Summary ---")
    print(model, outcome = "Nutrient_regime")

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

    # Save side-by-side plots with larger width
    # Save side-by-side plots with main title
    png_filename <- paste0("Ecosite_output/No_climate/VarImp_", gsub(" ", "_", nsr_name), ".png")
    png(png_filename, width = 1600, height = 900, res = 150)  # slightly taller for main title
    grid.arrange(p1, p2, ncol = 2, top = textGrob(
      paste(nsr_name),
      gp = gpar(fontsize = main_title_size, fontface = "bold", fontfamily = "Arial")
    ))
    dev.off()

    message("Variable importance plots saved for ", nsr_name)


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
    pd_png_filename <- paste0("Ecosite_output/No_climate/PartialDependence_", gsub(" ", "_", nsr_name), ".png")
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

    writeRaster(moist_rast, paste0("Ecosite_output/No_climate/Moist_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)
    writeRaster(nut_rast, paste0("Ecosite_output/No_climate/Nut_", gsub(" ", "_", nsr_name), ".tif"), overwrite = TRUE)

    # Classify
    # --- Classification ---
    moist_breaks <- quantile(values(moist_rast), probs = seq(0, 1, length.out = 9), na.rm = TRUE)
    moist_cat <- cbind(moist_breaks[-length(moist_breaks)], moist_breaks[-1], 1:8)
    moist_rast_cat <- classify(moist_rast, moist_cat)
    levels(moist_rast_cat) <- data.frame(
      ID = 1:8,
      Category = c("Xeric", "Subxeric", "Submesic", "Mesic", "Subhygric", "Hygric", "Subhydric", "Hydric")
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

    writeRaster(ecosite_rast, paste0("Ecosite_output/No_climate/Ecosite_", nsr_name, ".tif"), overwrite = TRUE)

    message("Finished processing ", nsr_name)

  }, error = function(e) {
    message("Error processing ", nsr_name, ": ", e$message)
  })
}



##Combine Athabasca Plain and Peace-Athabasca Delta (only 9 points in Delta, no model predictions made)

NSR9 <- subset(data1, NSRNAME %in% c("Athabasca Plain","Peace-Athabasca Delta"))
View(NSR9)

# Define combined subregion name and members
nsr_name <- "Athabasca_Plain__Peace_Athabasca_Delta"
combined_subregions <- c("Athabasca Plain", "Peace-Athabasca Delta")


# Define combined subregion name and members
nsr_name <- "Athabasca_Plain__Peace_Athabasca_Delta"
combined_subregions <- c("Athabasca Plain", "Peace-Athabasca Delta")

process_combined_nsr(
  nsr_name = "Athabasca_Plain__Peace_Athabasca_Delta",
  combined_subregions = c("Athabasca Plain", "Peace-Athabasca Delta"),
  NSR = NSR9,
  rasters = rasters,
  predictors = predictors,
  rf_formula = rf_formula
)




