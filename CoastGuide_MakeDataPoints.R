
library(climr)
library(terra)
library(data.table)
library(ccissr)

## PRISM DEM
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
dem.bc <- rast(paste(dir, "PRISM_dem/PRISM_dem.asc", sep=""))

bec <- vect("C:/Users/CMAHONY/OneDrive - Government of BC/Data/BEC13_Draft.gdb", layer="BEC13_v1")
coast <- aggregate(bec)
coast <- simplifyGeom(coast, tol = 100)
coast <- project(coast, dem.bc)
# plot(coast)

dem.noram <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/DEM/northamerica_elevation_cec_2023.tif") # option to use a local copy for faster processing.
# dem.noram <- project(dem.noram, dem, method="near") #project 250m source dem to the study area grid. method="near" to preserve elevation variance
dem.noram <- project(dem.noram, coast) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
dem.coast <- crop(dem.noram, coast) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
dem.coast <- mask(dem.coast, coast) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
# plot(dem.coast)

# Reproject BEC to match DEM's CRS
bec_wgs84 <- project(bec, crs(dem.coast))

# Extract grid points (coordinates and elevation) from DEM
coords <- as.data.table(as.data.frame(dem.coast, xy = TRUE))
setnames(coords, c("x", "y", "northamerica_elevation_cec_2023"), c("lon", "lat", "elevation"))

# Convert grid points to a SpatVector
points <- vect(coords, geom = c("lon", "lat"), crs = crs(dem.coast))

# Create a raster of MAP_LABEL from bec
bec_raster <- rasterize(bec_wgs84, dem.coast, field = "MAP_LABEL", touches = TRUE)

# Extract MAP_LABEL values for each DEM point
coords[, MAP_LABEL := extract(bec_raster, as.matrix(.SD)), .SDcols = c("lon", "lat")]

# Remove points with no label if needed
coords <- coords[!is.na(MAP_LABEL)]

# Add unique ID
coords[, id := .I]

# Rearrange columns
pts <- coords[, .(id, MAP_LABEL, lon, lat, elevation)]
colnames(pts) <- c("id", "BGC", "lon", "lat", "elev") # rename column names to what climr expects

# Remove whitespace from the BGC column
pts[, BGC := factor(gsub("\\s+", "", as.character(BGC)))]

# Rename certain units
pts[BGC %in% c("MHws", "MHwsp"), BGC := factor(ifelse(BGC == "MHws", "MHms", "MHmsp"))]

# subset units to lmh 
bgcs <- as.vector(read.csv("inputs/units_LMH77.csv")$BGC)
pts <- pts[BGC %in% bgcs]
pts$BGC <- droplevels(pts$BGC)
table(pts$BGC)

sampled_pts <- pts[, .SD[sample(.N, min(.N, min(table(pts$BGC))))], by = BGC]
table(sampled_pts$BGC)

# Export
write.csv(sampled_pts, "inputs/pts_lmh77.csv", row.names = F)

