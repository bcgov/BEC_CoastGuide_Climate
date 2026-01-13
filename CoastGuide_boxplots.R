
library(climr)
library(terra)
library(data.table)
library(ccissr)

#historical climate for training points
pts <- fread("inputs/pts_lmh77.csv")

elements_log = c("AHM", "DD", "Eref", "FFP", "NFFD", "PAS", "PPT", "SHM", "MSP")

clim.pts <- downscale(xyz = pts,
                      vars = list_vars(), 
                      which_refmap = "refmap_prism", 
                      gcms = list_gcms()[c(1,4,5,6,7,10,11,12)],
                      gcm_periods = list_gcm_periods()[2], 
                      ssps = list_ssps()[2],
                      max_run = 0, 
                      obs_periods = list_obs_periods()[1]
                      )
addVars(clim.pts)
clim.pts.log <- logVars(clim.pts, elements = elements_log, base=2)

vars <- c("elev", "PPT_MJ", "PPT_JAS", "DD5", "DDsub0", "SHM") # alternative variable set (for reviewers)
vars <- c("elev", "MAT", "PPT", "PAS", "CMD", "TD")

data("variables")
var.names <- c("Elevation (m)", "Mean annual\ntemperature (\u00b0C)", "Annual\nprecipitation (mm)", "Precipitation\nas snow (mm)", "Climatic moisture\ndeficit (mm)", "Continentality (\u00b0C)")

# extract out the reference period
clim.pts.ref <- clim.pts.log[PERIOD=="1961_1990", ]
clim.pts.ref <- merge(pts, clim.pts.ref, by = "id") # Merge the two datasets by 'id'
clim.pts.ref <- clim.pts.ref[BGC != "CMAun"] # Remove the "CMAun" class
# bgc_order <- clim.pts.ref[, .(mean_MAT = mean(MAT, na.rm = TRUE)), by = BGC][order(-mean_MAT)]$BGC # Calculate the mean MAT for each BGC group and reorder the factor levels
bgc_order <- sort(unique(clim.pts.ref$BGC)) # alphabetical order
clim.pts.ref$BGC <- factor(clim.pts.ref$BGC, levels = bgc_order)

# Calculate the median climate of the 1961_1990 period
clim.pts.ref.median <- clim.pts.ref[, lapply(.SD, median, na.rm = TRUE), by = BGC, .SDcols = vars]
clim.pts.ref.median[, BGC := factor(BGC, levels = bgc_order)] # Sort clim.pts.ref.median by the reordered BGC
setorder(clim.pts.ref.median, BGC) # Sort clim.pts.ref.median by the reordered BGC

# Calculate the median climate of the 2001-2020 period
clim.pts.obs <- clim.pts.log[PERIOD=="2001_2020", ]
clim.pts.obs <- merge(pts, clim.pts.obs, by = "id") # Merge the two datasets by 'id'
clim.pts.obs <- clim.pts.obs[BGC != "CMAun"] # Remove the "CMAun" class
clim.pts.obs.median <- clim.pts.obs[, lapply(.SD, median, na.rm = TRUE), by = BGC, .SDcols = vars]
clim.pts.obs.median[, BGC := factor(BGC, levels = bgc_order)] # Sort clim.pts.obs.median by the reordered BGC
setorder(clim.pts.obs.median, BGC) # Sort clim.pts.obs.median by the reordered BGC

# Calculate the median climate of the 2021-2040 period
clim.pts.2021 <- clim.pts.log[PERIOD=="2021_2040", ]
clim.pts.2021 <- clim.pts.2021[, lapply(.SD, mean, na.rm = TRUE), by = id, .SDcols = -c(1,2,3,4,5)] #ensemble mean for each point
clim.pts.2021 <- merge(pts, clim.pts.2021, by = "id") # Merge the two datasets by 'id'
clim.pts.2021 <- clim.pts.2021[BGC != "CMAun"] # Remove the "CMAun" class
clim.pts.2021.median <- clim.pts.2021[, lapply(.SD, median, na.rm = TRUE), by = BGC, .SDcols = vars]
clim.pts.2021.median[, BGC := factor(BGC, levels = bgc_order)] # Sort clim.pts.2021.median by the reordered BGC
setorder(clim.pts.2021.median, BGC) # Sort clim.pts.2021.median by the reordered BGC

# # Calculate the median climate of the 2041-2060 period
# clim.pts.2041 <- clim.pts.log[PERIOD=="2041_2060", ]
# clim.pts.2041 <- clim.pts.2041[, lapply(.SD, mean, na.rm = TRUE), by = id, .SDcols = -c(1,2,3,4,5)] #ensemble mean for each point
# clim.pts.2041 <- merge(pts, clim.pts.2041, by = "id") # Merge the two datasets by 'id'
# clim.pts.2041 <- clim.pts.2041[BGC != "CMAun"] # Remove the "CMAun" class
# clim.pts.2041.median <- clim.pts.2041[, lapply(.SD, median, na.rm = TRUE), by = BGC, .SDcols = vars]
# clim.pts.2041.median[, BGC := factor(BGC, levels = bgc_order)] # Sort clim.pts.2041.median by the reordered BGC
# setorder(clim.pts.2041.median, BGC) # Sort clim.pts.2041.median by the reordered BGC


# color scheme
subzones_colours_ref <- fread("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv")
color_map <- setNames(subzones_colours_ref$RGB, subzones_colours_ref$BGC )
bgc_colors <- color_map[levels(clim.pts.ref$BGC)]

# Create a boxplot of AHM grouped by BGC with reordered groups
png(filename=paste("plots\\CoastGuide.boxplots.png",sep="."), type="cairo", units="in", width=6.5, height=7, pointsize=12, res=600)
pdf(file=paste("plots/CoastGuide.boxplots.pdf",sep="."), width=6.5, height=7, pointsize=12)
mat <- matrix(c(7,1,2,3,4,5,6,8), 8) 
layout(mat, heights = c(0.5, rep(1, times=6), 0.5))

par(mar=c(0.5, 5, 0.5, 1), mgp=c(2.75, 0.25, 0), tck=-0.01)

var="MAT"
for(var in vars){
  boxplot(get(var) ~ BGC, data = clim.pts.ref, 
          main = "", 
          range = 0,
          yaxt = if(length(grep(var, elements_log))>0) "n" else "s",
          xaxt = if(var==vars[length(vars)]) "s" else "n",
          xlab = "", 
          ylab = var.names[which(vars==var)],
          ylim = range(c(clim.pts.ref[,get(var)], clim.pts.2021.median[,get(var)]), na.rm = T),
          las = 2,           # Rotate x-axis labels for better visibility
          col = bgc_colors # Set box color
  ) 
  
  # Add grid lines for better readability
  grid(nx = NA, ny = NULL, col = "gray", lty = "dotted")
  
  # Redraw boxplots on top of grid lines
  boxplot(get(var) ~ BGC, data = clim.pts.ref, 
          add=T, 
          main = "", 
          range = 0,
          yaxt = if(length(grep(var, elements_log))>0) "n" else "s",
          xaxt = if(var==vars[length(vars)]) "s" else "n",
          xlab = "", 
          ylab = var.names[which(vars==var)],
          ylim = range(c(clim.pts.ref[,get(var)], clim.pts.2021.median[,get(var)]), na.rm = T),
          las = 2,           # Rotate x-axis labels for better visibility
          border = "gray40", # Set border color
          col = "gray" # Set box color
  ) 
  
  if(var==vars[1]) axis(3, at = 1:length(levels(clim.pts.ref$BGC)), labels = levels(clim.pts.ref$BGC), las=2)
  if(length(grep(var, elements_log))>0){ 
    sequence <- if(var=="PAS") seq(2,16, 2) else seq(1,16)
    axis(2, at=log2(2^sequence), labels=2^sequence, las=2)
  }
  
  if(var!="elev"){
    # points(clim.pts.obs.median[,get(var)])
    points(clim.pts.2021.median[,get(var)], pch=16)
  }
  
}

dev.off()
