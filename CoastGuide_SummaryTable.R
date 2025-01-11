
library(climr)
library(ccissr)

#historical climate for training points
pts <- fread("inputs/pts_lmh77.csv")

clim.pts <- downscale(xyz = pts,
                      which_refmap = "refmap_prism", 
                      vars=list_vars(set = "Annual")
                      )

# combine the tables
clim.pts.ref <- merge(pts, clim.pts, by = "id") # Merge the two datasets by 'id'

# Define the climate variables to calculate summary stats for
climate_vars <- list_vars(set = "Annual") # Replace with actual climate variable names

# Calculate the mean, sd, min, and max for each climate variable grouped by BGC
summary_stats <- clim.pts.ref[, lapply(.SD, function(x) {
  list(mean = mean(x, na.rm = TRUE),
       sd = sd(x, na.rm = TRUE),
       min = min(x, na.rm = TRUE),
       max = max(x, na.rm = TRUE))
}), by = BGC, .SDcols = climate_vars]

# Unnest the list columns into separate columns
summary_table <- summary_stats[, unlist(.SD, recursive = FALSE), by = BGC]

# Rename the columns to reflect the variable and the statistic
new_colnames <- c("BGC", 
                  unlist(lapply(climate_vars, function(var) {
                    c(paste0(var, "_mean"), paste0(var, "_sd"),
                      paste0(var, "_min"), paste0(var, "_max"))
                  })))
setnames(summary_table, old = names(summary_table), new = new_colnames)

# Sort the summary table by BGC alphabetically
setorder(summary_table, BGC)

# View the resulting summary table
write.csv(summary_table, "outputs/climate_summary_LMH77.csv", row.names = F)
