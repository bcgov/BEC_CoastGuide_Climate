
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

###########################
## PRISM data
###########################

dir <- paste("//objectstore2.nrs.bcgov/ffec/Climatologies/PRISM_BC/", sep="")

# Precipitation rasters
var <- "pr"
e <- which(c("tmin", "tmax", "pr")==var)
files <- list.files(dir, pattern=paste(var,"*.tif", sep="."))
prism.monthly <- rast(paste0(dir, files))
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
dem.noram <- project(dem.noram, dem.bc) #project 250m source dem to the study area grid. ended up using bilinear interpolation because method="near" produces underestimation of lapse rates later in the workflow.
land.bc <- dem.noram
values(land.bc)[!is.finite(values(land.bc))] <- NA
values(land.bc)[is.finite(values(land.bc))] <- 1

# Calculate hillshade
slope.bc <- terrain(dem.bc)
aspect.bc <- terrain(dem.bc)
hill.bc <- shade(slope.bc, aspect.bc, angle = 45, direction = 315)
hill.bc <- mask(hill.bc, land.bc)

# coastline
bdy.bc <- vect(bc_bound_hres())
bdy.bc <- project(bdy.bc, dem.bc)

#######################
## Key Map
#######################



png(filename=paste("plots\\CoastGuide.studyareas.png",sep="."), type="cairo", units="in", width=6.5, height=5.4, pointsize=10, res=600)
par(mfrow=c(1,1), mar=c(0,0,0,0))
X <- mask(dem.bc, land.bc)
X <- crop(X, ext(c(-134, -119, 48, 55)))
lim <- quantile(values(X, na.rm=T), 0.99)
values(X)[which(values(X)>lim)] <- lim
# plot(hill, col=alpha(grey(0:100/100), 1), maxpixels=ncell(hill), legend=F)
plot(X, col=terrain.colors(99), xaxt="n", yaxt="n", legend=F)

regions <- c("NCBC", "CCBC", "SCBC")
lon1 <- -133.5
region1 <- ext(c(lon1, lon1+7.125, 53, 53.5))
lon1 <- -128.5
region2 <- ext(c(lon1, lon1+6.7, 50, 50.5))
lon1 <- -127
region3 <- ext(c(lon1, lon1+6.55, 49.25, 49.75))
for(i in 1:3){
  region <- get(paste("region", i, sep=""))
  plot(region, add=T)
  text(ext(region)[1], ext(region)[4]-0.25, regions[i], font=2, pos=4, offset=0.1)
}
box()
dev.off()

#======================================
#======================================
# Precipitation transects
#======================================
#======================================

i=1
par(mfrow=c(6,1), xpd=F)

png(filename=paste("plots\\CoastGuide.PrecipProfile.png",sep="."), type="cairo", units="in", width=6.5, height=5.4, pointsize=10, res=600)
mat <- matrix(c(1,2,3,4,5,6), 6); 
layout(mat, heights = rep(c(0.6,1), times=3))

for(i in 1:3){
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
  row <- c(30,30,30)[which(regions==caseStudy)]
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
  
  
  months <- c(1, 7)
  
  x0 <- crop(prism, studyarea)
  x1 <- crop(prism.monthly[[months[1]]], studyarea)
  x2 <- crop(prism.monthly[[months[2]]], studyarea)
  X <- dem
  
  # x1 <- mask(x1, land)
  # x2 <- mask(x2, land)
  # X <- mask(X, land)
  # hill <- mask(hill, land)
  
  values(prism) <- log2(values(prism))
  values(prism)[!is.finite(values(prism))] <- NA
  
  inc=0.05
  breaks=seq(min(values(prism), na.rm = T)-inc, max(values(prism), na.rm = T)+inc, inc)
  ColScheme <- colorRampPalette(brewer.pal(9,"YlGnBu"))(length(breaks)-1)
  
  
  
  # ## elevation plot
  # 
  # legend.args=list(text='Elevation (m)', side=2, font=2, line=0.5, cex=0.8)
  # lim <- quantile(values(X), 0.99)
  # values(X)[which(values(X)>lim)] <- lim
  # # plot(hill, col=alpha(grey(0:100/100), 1), maxpixels=ncell(hill), legend=F)
  # plot(X, col=terrain.colors(99), legend=F, xaxt="n", yaxt="n")
  # plot(hill, add=T, col=alpha(grey(0:100/100), 0.5), maxpixels=ncell(hill), legend=F)
  # lines(c(-180, extent(X)[2]), c(rowlat,rowlat), col="black", lty=2)
  # text(xFromCol(dem, 1), rowlat+0.05, "Cross-section", pos=4, cex=1.)
  # mtext(paste("(",letters[1],")", sep=""), line=-1.5, adj=0.005, side=3, cex=1)
  
  ## precipitation plots
  
  par(mar=c(0,4,0.2,5), mgp=c(2, 0.25, 0))
  image(prism, col=ColScheme, breaks=breaks, axes=F, legend=F)
  legend_ramp(prism, title = "", ColScheme = ColScheme, breaks = breaks, pos=c(1.005, 1.025, 0, 1), log = 2, horizontal = FALSE)
  axis(2, at=rowlat, labels=paste0(round(rowlat, 1), "\u00b0N"), tck=0, las=2, cex=0.8)
  lines(c(-180, 0), c(rowlat,rowlat), col="black", lty=2)
  plot(bdy, add=T)
  mtext(paste("(",letters[i],")", sep=""), line=-1.5, adj=-0.05, side=3, cex=1)
  # mtext("Annual Precipitation", side=1, line=-1.5, adj=0.01, font=2, cex=0.8)
  # mtext(paste("(",letters[2:26][which(datasets==dataset)],")", sep=""), line=-1.5, adj=0.005, side=3, cex=1)
  #   points(stn.info[, 6:7], pch=21, bg="white", cex=1.25, lwd=0.5)
  #   points(stn.info[, 6:7], pch=21, bg=ColScheme[cut(log2(stn.ppt.ann), breaks=breaks)], cex=1.25, lwd=0.5)
  box()
  
  ## latitudinal cross section
  
  par(mar=c(1.8,4,0,5), mgp=c(2, 0.25, 0))
  ylim.lower=c(0,0,0)[which(regions==caseStudy)]
  ylim.upper=c(5000, 5000, 5000)[which(regions==caseStudy)]
  shiftFactor <- c(9,8,8.8)[which(regions==caseStudy)]
  scaleFactor <- c(1400,1000,1200)[which(regions==caseStudy)]
  shiftFactor <- 8.7
  scaleFactor <- 1200
  plot(xvals, y.dem, col="white", xaxs="i", yaxs="i", ylim=c(ylim.lower,ylim.upper), ylab="", xlab="Degrees longitude", yaxt="n", xaxt="n")
  # axis(1, at=xvals, labels=round(2^seq(1,16)), tck=0, las=2)
  axis(2, at=c(0,1000,2000, 2800), labels=c(paste0(c(0,1000,2000), "m"), "Elev."), tck=0, las=2, cex=0.8)
  # mtext(side = 2, line = 2, "Elev.", cex=0.8, font=2, adj=0, las=2)
  latseq <- seq(ceiling(ext(dem)[1]), floor(ext(dem)[2]),2)
  # axis(1, at=colFromX(dem, latseq), labels=latseq, tck=0)
  par(mgp=c(2, 0.1, 0))
  
  axis(4, at=(log2(2^seq(1,16))-shiftFactor)*scaleFactor, labels=round(2^seq(1,16)), tck=0, las=2)
  mtext(side = 4, line = 2.5, "Precip. (mm)", cex=0.8, font=2)
  polygon(c(xvals, length(y.dem), 1), c(y.dem, 0, 0), col="gray", border=F)
  lines(xvals, (log2(values(crop(x0, row_extent)))-shiftFactor)*scaleFactor, col="grey50", lwd=4)
  # lines(xvals, (log2(values(x1, row, nrows=1))-shiftFactor)*scaleFactor, col="grey50", lwd=4)
  # lines(xvals, (log2(values(x2, row, nrows=1))-shiftFactor)*scaleFactor, col="black", lwd=2)
  # mtext("Cross-section", side=1, line=-1.5, adj=0.01, font=2, cex=0.8)
  # for(i in 1:length(col)){
  #   left <- col[i]-(oddRound(ngb/latFactor)-1)/2
  #   right <- col[i]+(oddRound(ngb/latFactor)-1)/2
  #   arrows(left,300, right, 300, length=0.01, angle=90, code=3)
  #   text(right,  300, letters[6:26][i], pos=4, cex=1, font=2)
  # }
  # legend(c(700,250,710,710,710)[which(regions==caseStudy)], ylim.upper+100, legend=datasets.names, col=c("black", "grey50"), lty=c(1,1), lwd=c(2,4) , bty="n", y.intersp = 0.8, cex=1)
  box()
  
  
}
dev.off()


