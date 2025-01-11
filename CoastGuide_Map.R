
remotes::install_github("bcgov/climr@devl")

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
plot(studyarea.bc, add=T)

###########################
## PRISM data
###########################

dir <- "C:/Users/CMAHONY/OneDrive - Government of BC/Data/climr_composite/"

# Precipitation rasters
var <- "Pr"
e <- which(elements==var)
files <- list.files(dir, pattern=paste(var,"*.tif", sep="."))
prism.monthly <- rast(paste0(dir, files))
prism.monthly <- crop(prism.monthly, studyarea.bc)
prism.annual <- sum(prism.monthly)

## import PRISM stations
dir <- "//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/"
stn.info <- fread(paste(dir, "Stations/",c("Tmin", "Tmax", "Pr")[e],"_uscdn_8110.csv", sep="")) #read in
for (i in which(names(stn.info)%in%c(month.abb, "Annual"))) stn.info[get(names(stn.info)[i])==c(-9999), (i):=NA, ] # replace -9999 with NA
stn.info <- stn.info[-which(El_Flag=="@"),]
stn.ppt.ann <- stn.info$Annual
stn.ppt.ann[stn.ppt.ann<0] <- NA

## DEM
dem.bc <- rast(paste(dir, "PRISM_dem/PRISM_dem.asc", sep=""))
dem.bc <- project(dem.bc, prism.annual)


# Ocean mask
# dem.noram <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif") #250m dem downloaded from http://www.cec.org/north-american-environmental-atlas/elevation-2023/
dem.noram <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/DEM/northamerica_elevation_cec_2023.tif") # option to use a local copy for faster processing.
# dem.noram <- project(dem.noram, dem, method="near") #project 250m source dem to the study area grid. method="near" to preserve elevation variance
dem.noram <- project(dem.noram, crs(dem.bc)) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
ocean.bc <- crop(dem.noram, studyarea.bc)
values(ocean.bc)[!is.finite(values(ocean.bc))] <- -99
values(ocean.bc)[values(ocean.bc)> -99] <- NA

# Calculate hillshade
slope.bc <- terrain(dem.bc)
aspect.bc <- terrain(dem.bc)
hill.bc <- shade(slope.bc, aspect.bc, angle = 45, direction = 315)
hill.bc <- mask(hill.bc, land.bc)

# coastline
bdy.bc <- vect(bc_bound_hres())
bdy.bc <- project(bdy.bc, dem.bc)
bdy.bc <- simplifyGeom(bdy.bc, tol = 0.001)

#######################
## Key Map
#######################

X <- prism.annual
values(X) <- log2(values(X))
values(X)[!is.finite(values(X))] <- NA
X <- crop(X, studyarea.bc)
X.all <- mask(X, land.bc)
lim_upper <- quantile(values(X.all, na.rm=T), 0.99)
lim_lower <- quantile(values(X.all, na.rm=T), 0.05)
values(X)[which(values(X)>lim_upper)] <- lim_upper
values(X)[which(values(X)<lim_lower)] <- lim_lower

inc=0.05
breaks=seq(min(values(X), na.rm = T)-inc, max(values(X), na.rm = T)+inc, inc)
ColScheme <- colorRampPalette(brewer.pal(9,"YlGnBu"))(length(breaks)-1)


png(filename=paste("plots/CoastGuide.map", var, "png",sep="."), type="cairo", units="in", width=6.5, height=6.25, pointsize=10, res=600)
par(mfrow=c(1,1), mar=c(0,0,0,0))
image(X, col=ColScheme, breaks=breaks, axes=F)
plot(ocean.bc, add=T, col="white")
rect(-140, 51.5, -130.9, 54.2, col="white", border = NA)
legend_ramp(X, title = "Annual Precipitation (mm)", ColScheme = ColScheme, breaks = breaks, pos=c(0.08, 0.1, 0.325, 0.65), log = 2, horizontal = FALSE, title.height = 1)
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
# dev.off()

#======================================
#======================================
# Precipitation transects
#======================================
#======================================

for(i in 1:2){
  slopeData <- "slope.Agg"
  studyarea <- get(paste("region", i, sep=""))  
  caseStudy <- regions[i]
  
  ## DEM
  
  land <- crop(land.bc, studyarea)
  dem <- crop(dem.bc, studyarea)
  slope =  crop(slope.bc, studyarea)
  aspect =  crop(aspect.bc, studyarea)
  hill =  crop(hill.bc, studyarea)
  prism <- crop(prism.annual, studyarea)
  bdy <- crop(bdy.bc, studyarea)
  
  dem <- mask(dem, land)
  
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
 
  x0 <- crop(prism, studyarea)

  ## latitudinal cross section
  if(i==1){
    par(fig = c(0.39, 0.99, 0.8, 0.99), new = TRUE) 
  } else {
    par(fig = c(0.01, 0.61, 0.01, 0.2), new = TRUE) 
  }
  par(mar=c(0,0,0,0), mgp=c(2, 0.2, 0))
  ylim.lower=c(0,0,0)[which(regions==caseStudy)]
  ylim.upper=c(5000, 5000, 5000)[which(regions==caseStudy)]
  shiftFactor <- c(9,8,8.8)[which(regions==caseStudy)]
  scaleFactor <- c(1400,1000,1200)[which(regions==caseStudy)]
  shiftFactor <- 8.7
  scaleFactor <- 1200
  plot(xvals, y.dem, col="white", xaxs="i", yaxs="i", ylim=c(ylim.lower,ylim.upper), ylab="", xlab="Degrees longitude", yaxt="n", xaxt="n")
  rect(-999, -999, 9999, 9999, col="white", border = NA)
  box()
  # axis(1, at=xvals, labels=round(2^seq(1,16)), tck=0, las=2)
  latseq <- seq(ceiling(ext(dem)[1]), floor(ext(dem)[2]),2)
  # axis(1, at=colFromX(dem, latseq), labels=latseq, tck=0)

  polygon(c(xvals, length(y.dem), 1), c(y.dem, 0, 0), col="gray", border=F)
  yvals <- (log2(values(crop(x0, row_extent)))-shiftFactor)*scaleFactor
  lines(xvals, yvals, col="grey50", lwd=3)
  text(min(which(is.finite(yvals)))+10, yvals[min(which(is.finite(yvals)))], "Precip.", pos=2, cex=0.8)
  
  axis(2, at=c(1000,2000), labels=rep("",2), tck=0.025, las=2, cex.axis=0.6)
  text(c(-5,-5), c(400, 1000,2000), c("Elevation", paste0(c(1000,2000), "m")), font=c(2,1,1), pos=4, cex=0.8)
  
  axis.pos <- 2^seq(1,16)
  axis(4, at=(log2(axis.pos)-shiftFactor)*scaleFactor, labels=rep("", length(axis.pos)), tck=0.025, las=2, cex.axis=0.6)
  text(rep(max(xvals)+5, length(axis.pos)), (log2(axis.pos)-shiftFactor)*scaleFactor, round(axis.pos), font=c(rep(1, length(axis.pos))), pos=2, cex=0.8)
  mtext(side = 3, line = -1.0, "Precipitation (mm)", cex=0.8, font=2, adj=0.99)
  
  box()
  
}

dev.off()


