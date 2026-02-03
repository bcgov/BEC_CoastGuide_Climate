
library(climr)
library(terra)
library(data.table)
library(ccissr)

## PRISM DEM
# dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/PRISM_BC/"
dem.bc <- rast(paste(dir, "PRISM_dem/PRISM_dem.asc", sep=""))

units <- fread("C:/Users/CMAHONY/GitHub/misc_temp/BEC_CoastGuide_Climate/inputs/units_SouthArea.csv")

bec <- vect("C:/Users/CMAHONY/OneDrive - Government of BC/Data/BECv14_Draft_30Jan2026/BECv14_Draft_30Jan2026.shp")
bec$BGC_LABEL <- gsub(" ", "", bec$BGC_LABEL)
bec <- subset(bec, bec$BGC_LABEL %in% units$BGC)
as.vector(unique(bec$BGC_LABEL))
# missing <- units$BGC[-which(units$BGC %in% as.vector(unique(bec$BGC_LABEL)))]
# plot(bec)

bec.latlon <- simplifyGeom(bec, tol = 100)
bec.latlon <- project(bec.latlon, dem.bc)
# plot(studyarea)

dem.noram <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/DEM/northamerica_elevation_cec_2023.tif") # option to use a local copy for faster processing.
# dem.noram <- project(dem.noram, dem, method="near") #project 250m source dem to the study area grid. method="near" to preserve elevation variance
dem.noram <- project(dem.noram, dem.bc) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
dem.studyarea <- crop(dem.noram, bec.latlon) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
dem.studyarea <- mask(dem.studyarea, bec.latlon) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
# plot(dem.studyarea)

# Extract grid points (coordinates and elevation) from DEM
coords <- as.data.table(as.data.frame(dem.studyarea, xy = TRUE))
setnames(coords, c("x", "y", "northamerica_elevation_cec_2023"), c("lon", "lat", "elevation"))

# Convert grid points to a SpatVector
points <- vect(coords, geom = c("lon", "lat"), crs = crs(dem.studyarea))

# Create a raster of BGC_LABEL from bec
bec_raster <- rasterize(bec.latlon, dem.studyarea, field = "BGC_LABEL", touches = TRUE)

# Extract BGC_LABEL values for each DEM point
coords[, BGC_LABEL := extract(bec_raster, as.matrix(.SD)), .SDcols = c("lon", "lat")]

# Remove points with no label if needed
coords <- coords[!is.na(BGC_LABEL)]

# Add unique ID
coords[, id := .I]

# Rearrange columns
pts <- coords[, .(id, BGC_LABEL, lon, lat, elevation)]
colnames(pts) <- c("id", "BGC", "lon", "lat", "elev") # rename column names to what climr expects

# Remove whitespace from the BGC column
pts[, BGC := factor(gsub("\\s+", "", as.character(BGC)))]

# Rename certain units
# pts[BGC %in% c("MHws", "MHwsp"), BGC := factor(ifelse(BGC == "MHws", "MHms", "MHmsp"))]

# # subset units to lmh 
# bgcs <- as.vector(read.csv("inputs/units_LMH77.csv")$BGC)
# pts <- pts[BGC %in% bgcs]
# pts$BGC <- droplevels(pts$BGC)
# table(pts$BGC)

sampled_pts <- pts[, .SD[sample(.N, min(.N, min(table(pts$BGC))))], by = BGC]
table(sampled_pts$BGC)

# Export
write.csv(sampled_pts, "inputs/pts_SouthArea.csv", row.names = F)

