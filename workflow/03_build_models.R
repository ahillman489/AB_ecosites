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
balanced_data <- nsr_combined %>%
  group_by(Ecosite_type) %>%
  slice_sample(n = 20, replace = TRUE) %>%
  ungroup()

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
