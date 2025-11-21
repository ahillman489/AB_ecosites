# Build models


# Clear Workspace ---------------------------------------------------------
rm(list=ls())
gc()


# Load Libraries ----------------------------------------------------------

library(terra)
library(dplyr)
library(randomForestSRC)
library(ggplot2)
library(gridExtra)
library(grid)
library(png)

# Load Functions ----------------------------------------------------------

source("R/functions.R", local=TRUE)

# Load Data ---------------------------------------------------------------

nsr_combined <- vect('cache/nsr_combined.gpkg')


predictors <- c("Proportion_sand",
                "Proportion_clay",
                "Soil_organic_carbon",
                "Surficial_geology",
                "Slope",
                "Elevation",
                #"CTI",
                "TPI",
                #"Profile_curvature",
                #"Tangential_curvature",
                "pH",
                "roughness",
                "soil_depth")

if (!all(predictors %in% names(nsr_combined))) {stop("Predictors missing from data.")}

# Build formula safely
formula_str <- paste("Multivar(Moisture_regime, Nutrient_regime) ~", paste(predictors, collapse = " + "))
rf_formula <- as.formula(formula_str)
rf_formula


# Create output directory if it doesn't exist
dir.create("Ecosite_output/all_province", showWarnings = FALSE, recursive = TRUE)

df <- assign_regime_values(as.data.frame(nsr_combined))

# Fit multivariate random forest model
model <- rfsrc(rf_formula, data = df, importance = "permute")

# Print summaries
message("\n--- Moisture_regime Model Summary ---")
print(model, outcome = "Moisture_regime")
message("\n--- Nutrient_regime Model Summary ---")
print(model, outcome = "Nutrient_regime")

# Balance input data
# sample from nsr_combined such that there will be 60 rows from each ecosite. Replacement is okay.
# use sample() to build the new dataframe
set.seed(123)  # For reproducibility
balanced_data <- as.data.frame(nsr_combined) %>%
  group_by(Ecosite_type) %>%
  slice_sample(n = 20, replace = TRUE) %>%
  ungroup()
balanced_data <- assign_regime_values(as.data.frame(balanced_data))

model1 <- rfsrc(rf_formula, data = balanced_data, importance = "permute")



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

# -------------------------------------------------------------------------
# Hannah Meyer's Spatial Variable Selection Method
# -------------------------------------------------------------------------

library(purrr)
library(tidyr)

# Helper function for spatial cross-validation
spatial_cv_performance <- function(vars, data, formula_target, nfolds = 5, seed = 42) {
  set.seed(seed)
  folds <- sample(rep(1:nfolds, length.out = nrow(data)))
  perf <- numeric(nfolds)
  for (i in 1:nfolds) {
    train <- data[folds != i, ]
    test  <- data[folds == i, ]
    formula_str <- paste(formula_target, "~", paste(vars, collapse = " + "))
    rf_formula <- as.formula(formula_str)
    model <- rfsrc(rf_formula, data = train)
    preds <- predict(model, newdata = test)$predicted
    # Use mean squared error for regression, adjust if classification
    perf[i] <- mean((test[[formula_target]] - preds)^2, na.rm = TRUE)
  }
  mean(perf)
}

# 1. Find best 2-variable combination
response_vars <- c("Moisture_regime", "Nutrient_regime") # adjust if needed
predictor_combos <- combn(predictors, 2, simplify = FALSE)
results_2var <- map_dfr(predictor_combos, function(vars) {
  tibble(
    Variables = list(vars),
    Moisture_MSE = spatial_cv_performance(vars, df, "Moisture_regime"),
    Nutrient_MSE = spatial_cv_performance(vars, df, "Nutrient_regime")
  )
})
best_2var_moist <- results_2var %>% arrange(Moisture_MSE) %>% slice(1)
best_2var_nutr  <- results_2var %>% arrange(Nutrient_MSE) %>% slice(1)

# 2. Iteratively add variables
screened_vars_moist <- unlist(best_2var_moist$Variables)
screened_vars_nutr  <- unlist(best_2var_nutr$Variables)
remaining_vars_moist <- setdiff(predictors, screened_vars_moist)
remaining_vars_nutr  <- setdiff(predictors, screened_vars_nutr)

improve <- TRUE
while (improve && length(remaining_vars_moist) > 0) {
  candidate_results <- map_dfr(remaining_vars_moist, function(v) {
    vars <- c(screened_vars_moist, v)
    tibble(
      Variable = v,
      MSE = spatial_cv_performance(vars, df, "Moisture_regime")
    )
  })
  best_candidate <- candidate_results %>% arrange(MSE) %>% slice(1)
  if (best_candidate$MSE < spatial_cv_performance(screened_vars_moist, df, "Moisture_regime")) {
    screened_vars_moist <- c(screened_vars_moist, best_candidate$Variable)
    remaining_vars_moist <- setdiff(remaining_vars_moist, best_candidate$Variable)
  } else {
    improve <- FALSE
  }
}

improve <- TRUE
while (improve && length(remaining_vars_nutr) > 0) {
  candidate_results <- map_dfr(remaining_vars_nutr, function(v) {
    vars <- c(screened_vars_nutr, v)
    tibble(
      Variable = v,
      MSE = spatial_cv_performance(vars, df, "Nutrient_regime")
    )
  })
  best_candidate <- candidate_results %>% arrange(MSE) %>% slice(1)
  if (best_candidate$MSE < spatial_cv_performance(screened_vars_nutr, df, "Nutrient_regime")) {
    screened_vars_nutr <- c(screened_vars_nutr, best_candidate$Variable)
    remaining_vars_nutr <- setdiff(remaining_vars_nutr, best_candidate$Variable)
  } else {
    improve <- FALSE
  }
}

# 3. Fit final screened models
formula_moist <- as.formula(paste("Moisture_regime ~", paste(screened_vars_moist, collapse = " + ")))
formula_nutr  <- as.formula(paste("Nutrient_regime ~", paste(screened_vars_nutr, collapse = " + ")))
model_screened_moist <- rfsrc(formula_moist, data = df)
model_screened_nutr  <- rfsrc(formula_nutr, data = df)

# 4. Compare model performances
message("\n--- Screened Moisture_regime Model Summary ---")
print(model_screened_moist)
message("\n--- Screened Nutrient_regime Model Summary ---")
print(model_screened_nutr)

# Save variable selection results
write.csv(data.frame(Screened_Moisture_Vars = screened_vars_moist), "Ecosite_output/all_province/screened_moisture_vars.csv", row.names = FALSE)
write.csv(data.frame(Screened_Nutrient_Vars = screened_vars_nutr), "Ecosite_output/all_province/screened_nutrient_vars.csv", row.names = FALSE)
