
# remotes::install_github("bcgov/climr@devl")

library(terra)
library(data.table)
library(RColorBrewer)
library(scales)
library(rnaturalearth)
library(bcmaps)
library(climr)

monthcodes <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
elements <- c("Tmin", "Tmax", "Pr")
element.names <- c("Tmin (\u00b0C)", "Tmax (\u00b0C)", "precipitation (mm)")

studyarea.bc <- ext(c(-132, -120.3, 48, 55))

###########################
## land mask
###########################

# Define the bounding box for the region of interest
lat_min <- 46
lat_max <- 60
lon_min <- -145
lon_max <- -120

# Load high-resolution land polygons
land <- ne_download(scale = "large", type = "land", category = "physical", returnclass = "sf")

# Convert land polygons to terra vectors (this step is necessary to work with terra)
land_terra <- vect(land)

# Define the bounding box as a terra vector (polygon)
bbox <- ext(lon_min, lon_max, lat_min, lat_max)
bbox_poly <- vect(bbox, crs = "EPSG:4326")  # EPSG:4326 for WGS 84 (latitude/longitude)

# Crop the land polygons to the bounding box area
land_cropped <- crop(land_terra, bbox_poly)

# Check if land_cropped is valid (should not be empty)
print(land_cropped)

# Create the land mask by merging (union) the land polygons in the cropped area
land_mask <- union(land_cropped)

# Check if land_mask is valid (should not be empty)
print(land_mask)

# Ensure both geometries are in the same CRS
if (crs(land_mask) != crs(bbox_poly)) {
  land_mask <- project(land_mask, crs(bbox_poly))
}

# Check the CRS of both geometries after projection
print(crs(land_mask))
print(crs(bbox_poly))

# Create the ocean mask by subtracting the land mask from the bounding box using setdiff
ocean_mask <- erase(bbox_poly, land_mask)

# Check the ocean_mask result
plot(ocean_mask)

###########################
## PRISM data
###########################

## import PRISM stations
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
dem.bc <- rast(paste(dir, "PRISM_dem/PRISM_dem.asc", sep=""))
dem.bc <- crop(dem.bc, studyarea.bc)

# Ocean mask
# dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif") #250m dem downloaded from http://www.cec.org/north-american-environmental-atlas/elevation-2023/
dem.noram <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/DEM/northamerica_elevation_cec_2023.tif") # option to use a local copy for faster processing.
# dem.noram <- project(dem.noram, dem, method="near") #project 250m source dem to the study area grid. method="near" to preserve elevation variance
dem.noram <- project(dem.noram, "EPSG:4326") #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
ocean.bc <- crop(dem.noram, studyarea.bc)
values(ocean.bc)[!is.finite(values(ocean.bc))] <- -99
values(ocean.bc)[values(ocean.bc)> -99] <- NA

## climr query points
grid <- as.data.frame(dem.bc, cells = TRUE, xy = TRUE)
colnames(grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

## climr data
clim <- downscale(grid, which_refmap = "refmap_climr", vars = list_vars())

#######################
## Key Map
#######################

var <- "MAT"

X <-  dem.bc
X[clim[, id]] <- clim[, get(var)] 

lim_upper <- quantile(values(X, na.rm=T), 0.99)
lim_lower <- quantile(values(X, na.rm=T), 0.01)
values(X)[which(values(X)>lim_upper)] <- lim_upper
values(X)[which(values(X)<lim_lower)] <- lim_lower

inc=0.05
breaks=seq(min(values(X), na.rm = T)-inc, max(values(X), na.rm = T)+inc, inc)
ColScheme <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks)-1)
ColScheme <- colorRampPalette(brewer.pal(11, "Greys"))(length(breaks)-1)


png(filename=paste("plots/CoastGuide.map", var, "greyscale.png",sep="."), type="cairo", units="in", width=6.5, height=6.25, pointsize=10, res=600)
# pdf(file=paste("plots/CoastGuide.map", var, "greyscale.pdf",sep="."), width=6.5, height=6.25, pointsize=10)
par(mfrow=c(1,1), mar=c(0,0,0,0))
image(X, col=ColScheme, breaks=breaks, axes=F)
plot(ocean_mask, add=T, col="white", border=F)
rect(-140, 51.5, -130.9, 54.2, col="white", border = NA)
legend_ramp(X, title = "Mean annual Temperature (\u00b0C)", ColScheme = ColScheme, breaks = breaks, pos=c(0.08, 0.1, 0.325, 0.65), log = NULL, horizontal = FALSE, title.height = 6)
box()

regions <- c("NCBC", "SCBC")
lon1 <- -130.7
region1 <- ext(c(lon1, lon1+4.5, 53.2, 53.3))
lon1 <- -128
region2 <- ext(c(lon1, lon1+7, 49.7, 49.8))
for(i in 1:2){
  region <- get(paste("region", i, sep=""))
  plot(region, add=T, lwd=1.5)
  text(ext(region)[c(2,1)][i], mean(ext(region)[3:4]), c("North transect", "South transect")[i], font=2, pos=c(4,2)[i], offset=0.2)
}
box()


# -----------------------
# North Arrow

usr <- ext(X)

x0 <- usr[1] + 0.065 * diff(usr[1:2])
y0 <- usr[3] + 0.8 * diff(usr[3:4])

arrows(
  x0, y0,
  x0, y0 + 0.095 * diff(usr[3:4]),
  lwd = 2,
  length = 0.1
)

text(
  x0, y0 + 0.115 * diff(usr[3:4]),
  "N",
  cex = 1.2,
  font = 2
)

# -----------------------
# Scale Bar

# km per degree latitude
km_per_deg_lat <- 111.32

# scale bar length (km)
km_len <- 100
deg_len <- km_len / km_per_deg_lat

usr <- ext(X)

# top-left placement
x0 <- usr[1] + 0.035 * diff(usr[1:2])
y0 <- usr[3] + 0.8 * diff(usr[3:4])

# main bar
segments(
  x0, y0,
  x0, y0 + deg_len,
  lwd = 2
)

# end ticks
tick_width <- 0.01 * diff(usr[1:2])

segments(x0 - tick_width, y0, x0 + tick_width, y0)
segments(x0 - tick_width, y0 + deg_len, x0 + tick_width, y0 + deg_len)

# label
text(x0 - 1.8 * tick_width, y0 + deg_len/2, paste0(km_len, " km"), cex = 0.9, srt=90)

#======================================
#======================================
# profile transects
#======================================
#======================================

for(i in 1:2){
  transect <- get(paste("region", i, sep=""))  
  caseStudy <- regions[i]
  
  ## DEM
  dem <- crop(dem.bc, transect)

  #################################
  ## Precipitation plot
  #################################
  row <- 6
  rowlat <- yFromRow(dem, row)

  row_height <- res(dem)[2]  # Resolution in the Y direction
  
  # Calculate the latitude for the specific row
  rowlat <- ymax(dem) - (row - 0.5) * row_height
  
  row_extent <- ext(dem)
  row_extent <- ext(row_extent$xmin, row_extent$xmax, rowlat - row_height / 2, rowlat + row_height / 2)
  
  # Crop the raster to the row extent
  row_raster <- crop(dem, row_extent)
  
  # Extract values from the cropped raster
  y.dem <- values(row_raster)
  y.dem <- y.dem+50
  y.dem[is.na(y.dem)] <- 0
  xvals <- 1:length(y.dem)
 
  var1 <- "Tave_01"
  var2 <- "Tave_07"

    x0 <- crop(X, transect)
  X1 <- X2 <- dem.bc
  X1[clim[, id]] <- clim[, get(var1)] 
  X2[clim[, id]] <- clim[, get(var2)] 
  x1 <- crop(X1, transect)
  x2 <- crop(X2, transect)
  
  ## latitudinal cross section
  if(i==1){
    par(fig = c(0.39, 0.99, 0.8, 0.99), new = TRUE) 
  } else {
    par(fig = c(0.01, 0.61, 0.01, 0.2), new = TRUE) 
  }
  par(mar=c(0,0,0,0), mgp=c(2, 0.2, 0))
  ylim.lower=c(0,0)[which(regions==caseStudy)]
  ylim.upper=c(5000, 5000)[which(regions==caseStudy)]
  shiftFactor <- -18
  scaleFactor <- 130
  plot(xvals, y.dem, col="white", xaxs="i", yaxs="i", ylim=c(ylim.lower,ylim.upper), ylab="", xlab="Degrees longitude", yaxt="n", xaxt="n")
  rect(-999, -999, 9999, 9999, col="white", border = NA)
  box()
  # axis(1, at=xvals, labels=round(2^seq(1,16)), tck=0, las=2)
  latseq <- seq(ceiling(ext(dem)[1]), floor(ext(dem)[2]),2)
  # axis(1, at=colFromX(dem, latseq), labels=latseq, tck=0)

  polygon(c(xvals, length(y.dem), 1), c(y.dem, 0, 0), col="gray", border=F)
  yvals1 <- (values(crop(x1, row_extent))-shiftFactor)*scaleFactor
  lines(xvals, yvals1, col="grey30", lwd=1)
  text(min(which(is.finite(yvals1)))+10, yvals1[min(which(is.finite(yvals1)))], "January", pos=2, cex=0.8)
  yvals2 <- (values(crop(x2, row_extent))-shiftFactor)*scaleFactor
  lines(xvals, yvals2, col="grey50", lwd=2)
  text(min(which(is.finite(yvals2)))+10, yvals2[min(which(is.finite(yvals2)))], "July", pos=2, cex=0.8)
  
  axis(2, at=c(1000,2000), labels=rep("",2), tck=0.025, las=2, cex.axis=0.6)
  text(c(-5,-5), c(400, 1000,2000), c("Elevation", paste0(c(1000,2000), "m")), font=c(2,1,1), pos=4, cex=0.8)
  
  axis.pos <- seq(-10, 10, 10)
  axis(4, at=(axis.pos-shiftFactor)*scaleFactor, labels=rep("", length(axis.pos)), tck=0.025, las=2, cex.axis=0.6)
  text(rep(max(xvals)+5, length(axis.pos)), (axis.pos-shiftFactor)*scaleFactor, paste0(axis.pos, "\u00b0C"), font=c(rep(1, length(axis.pos))), pos=2, cex=0.8)
  if(i==1) mtext(side = 3, line = -1.0, "Temperature", cex=0.8, font=2, adj=0.99)
  
  box()
  
}

dev.off()


