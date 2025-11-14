# Build prediction rasters


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


# Load raster layers
raster_paths <- list(
  Proportion_sand <- "data/spatial/mean_Sand_0_60.tif",
  Proportion_clay <- "data/spatial/mean_Clay_0_60.tif",
  Soil_organic_carbon <- "data/spatial/mean_SOC_0_60.tif",
  Surficial_geology = "data/spatial/GreatGroups_10TM.tif",
  Slope = "data/spatial/Slope.tif",
  Elevation = "data/spatial/DEM_100m.tif",
  CTI = "data/spatial/AB_cti.tif", # topographic wetness index
  TPI = "data/spatial/AB_tpi.tif", # topographic position index
  Profile_curvature = "data/spatial/AB_pcurv.tif",
  Tangential_curvature = "data/spatial/AB_tcurv.tif",
  pH = "data/spatial/mean_pH_0_60.tif",
  roughness = "data/spatial/10m/roughness.tif",
  soil_depth = "data/spatial/soil_depth_ab.tif"
)

rasters <- lapply(raster_paths, rast)

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
