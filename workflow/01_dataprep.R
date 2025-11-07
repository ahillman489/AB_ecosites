# Data Preparation for Modeling Ecosites
# C Mallon


# Clear Workspace ---------------------------------------------------------

rm(list = ls())
gc()


# Import Libraries --------------------------------------------------------

library(terra)
library(reproducible)
library(fasterRaster)


# Parameters --------------------------------------------------------------

focal_nsr <- c("Lower Boreal Highlands") # Set NA for no filter
first_production <- TRUE

# Functions

# delete object cleanly
rmObject <- function(x) {
  rm(x)
  gc()
}

# Setup Directories -------------------------------------------------------

if (!dir.exists('cache')) {
  dir.create('cache')
}

if (!dir.exists("data/spatial/10m")) {
  dir.create("data/spatial/10m")
}

# Load Data ---------------------------------------------------------------

NSR <- vect("data/spatial/Natural_Sub_Regions_Alberta_2005_shapefile.shp")
if (!is.na(focal_nsr)){
  NSR <- NSR[NSR$NSRNAME %in% focal_nsr]
}

dem <- rast('data/spatial/dem_10.tif')

# Calculate slope from dem
if (!file.exists("data/spatial/10m/slope.tif")){
  slope <- terrain(x=dem, v="slope", neighbors=8, unit="degrees", filename="data/spatial/10m/slope.tif")
  } else {
    slope <- rast("data/spatial/10m/slope.tif")
  }
if (first_production) {
  rm(slope)
  gc()
}

# Calculate aspect from dem
if (!file.exists("data/spatial/10m/aspect.tif")){
  aspect <- terrain(x=dem, v="aspect", neighbors=8, unit="degrees", filename="data/spatial/10m/aspect.tif")
} else {
  aspect <- rast("data/spatial/10m/aspect.tif")
}
if (first_production) {
  rm(aspect)
  gc()
}

# Calcualte TPI
if (!file.exists("data/spatial/10m/tpi.tif")){
  tpi <- terrain(x=dem, v="TPI", neighbors=8, filename="data/spatial/10m/tpi.tif")
} else {
  tpi <- rast("data/spatial/10m/tpi.tif")
}
if (first_production) {
  rm(tpi)
  gc()
}

# Calculate TRI
if (!file.exists("data/spatial/10m/tri.tif")){
  tri <- terrain(x=dem, v="TRI", neighbors=8, filename="data/spatial/10m/tri.tif")
} else {
  tri <- rast("data/spatial/10m/tri.tif")
}
if (first_production) {
  rm(tri)
  gc()
}

# Calculate roughness
if (!file.exists("data/spatial/10m/roughness.tif")){
  roughness <- terrain(x=dem, v="roughness", neighbors=8, filename="data/spatial/10m/roughness.tif")
} else {
  roughness <- rast("data/spatial/10m/roughness.tif")
}
if (first_production) {
  rm(roughness)
  gc()
}

# Calculate flowdir
if (!file.exists("data/spatial/10m/flowdir.tif")){
  flowdir <- terrain(x=dem, v="flowdir", filename="data/spatial/10m/flowdir.tif")
} else {
  flowdir <- rast("data/spatial/10m/flowdir.tif")
}
if (first_production) {
  rm(flowdir)
  gc()
}

# Calculate topographic wetness index
if (!file.exists('data/spatial/10m/wetness.tif')){
  g_elev <- fast('data/spatial/dem_10.tif')
  twi <- wetness(g_elev)
  fasterRaster::writeRaster(twi, 'data/spatial/10m/wetness.tif')
  twi <- rast('data/spatial/dem_10.tif')
} else {
  twi <- rast('data/spatial/dem_10.tif')
}


# Soil Depth
sd <- rast("data/spatial/depth_to_bedrock_canada_100m.tif")

nsr <- vect("data/spatial/Natural_Sub_Regions_Alberta_2005_shapefile.shp")

# dissolve all polygons in nsr into a single polygon
nsr_dissolved <- aggregate(nsr)
#Project nsr_dissolved into the crs of sd
nsr_dissolved_proj <- project(nsr_dissolved, crs(sd))
#Clip raster to nsr_dissolved_proj
sd_clipped <- mask(crop(sd, nsr_dissolved_proj), nsr_dissolved_proj)

sd_p <- project(sd_clipped, crs(nsr_dissolved))

terra::writeRaster(sd_p, "data/spatial/soil_depth_ab.tif")
