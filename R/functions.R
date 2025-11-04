

# Helper function to assign numeric values
assign_regime_values <- function(df) {
  df %>%
    mutate(Moisture_regime = case_when(
      Moisture.Regime == 'Xeric' ~ 0.5,
      Moisture.Regime == 'Subxeric' ~ 1.0,
      Moisture.Regime == 'Submesic' ~ 1.5,
      Moisture.Regime == 'Mesic' ~ 2.0,
      Moisture.Regime == 'Subhygric' ~ 2.5,
      Moisture.Regime == 'Hygric' ~ 3.0,
      Moisture.Regime == 'Subhydric' ~ 3.5,
      Moisture.Regime == 'Hydric' ~ 4.0
    ),
    Nutrient_regime = case_when(
      Nutrient.Regime == 'Very Poor' ~ 0.5,
      Nutrient.Regime == 'Poor' ~ 1.0,
      Nutrient.Regime == 'Medium' ~ 2.0,
      Nutrient.Regime == 'Rich' ~ 3.0,
      Nutrient.Regime == 'Very Rich' ~ 3.5
    ))
}

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


    # Check required vars
    if (!("Moisture_regime" %in% names(combined_df)) || !("Nutrient_regime" %in% names(combined_df))) {
      message("Skipping ", nsr_name, ": Moisture_regime or Nutrient_regime missing.")
      return(NULL)
    }

    # Fit model
    model <- rfsrc(rf_formula, data = combined_df, importance = TRUE)

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
