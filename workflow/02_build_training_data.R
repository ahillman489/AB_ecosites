# Build training data for models

# This script is where the spatial data is read in and compared with the point
# data to produce dataframes that can be read into the random forest model,
# or any other model we choose to explore.

# Clear Workspace ---------------------------------------------------------

rm(list = ls())
gc()


# Import Libraries --------------------------------------------------------

library(terra)
library(dplyr)
library(qs)


# Import Functions --------------------------------------------------------

source("R/functions.R", local=TRUE)


# Import Data -------------------------------------------------------------

data <- read.csv("data/Updated_Sites_Public.csv")

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
NSR <- vect("data/spatial/Natural_Sub_Regions_Alberta_2005_shapefile.shp")

southern_nsrs <- c("Mixedgrass", "Dry Mixedgrass",
                   "Foothills Fescue", "Foothills Parkland",
                   "Northern Fescue", "Central Parkland")

# Load raster layers
raster_paths <- list(
  Proportion_sand = "data/spatial/mean_Sand_0_60.tif",
  Proportion_clay = "data/spatial/mean_Clay_0_60.tif",
  Soil_organic_carbon = "data/spatial/mean_SOC_0_60.tif",
  Surficial_geology = "data/spatial/GreatGroups_10TM.tif",
  Slope = "data/spatial/10m/slope.tif",
  Elevation = "data/spatial/dem_10.tif",
  CTI = "data/spatial/10m/wetness.tif", # topographic wetness index
  TPI = "data/spatial/10m/tpi.tif", # topographic position index
 # Profile_curvature = "data/spatial/AB_pcurv.tif",
  #Tangential_curvature = "data/spatial/AB_tcurv.tif",
  pH = "data/spatial/mean_pH_0_60.tif",
  roughness = "data/spatial/10m/roughness.tif",
  soil_depth = "data/spatial/soil_depth_ab.tif"
)

rasters <- lapply(raster_paths, rast)

# Predictor variable names (must match raster and data)
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


# Subset data for the NSR and create spatial points

nsr_points <- vect(data1,
                   geom = c("Easting", "Northing"),
                   crs = "EPSG:3400")

# Subset NSR shapefile
nsr_shape <- NSR[!(NSR$NSRNAME %in% (southern_nsrs))]
nsr_points <- nsr_points[!(nsr_points$NSRNAME %in% (southern_nsrs))]

nsr_combined <- terra::deepcopy(nsr_points)

for (p in predictors){
  r <- rasters[[p]]
  raster_values <- terra::extract(r, nsr_points)[, -1]
  nsr_combined[[p]] <- raster_values
}

terra::writeVector(nsr_combined, 'cache/nsr_combined.gpkg', filetype="GPKG")

