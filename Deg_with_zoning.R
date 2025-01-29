library(nlme)
library(fields)
library(data.table)
library(foreach)
library(doParallel)
#install.packages("doParallel")

setwd(this.path::here())
#source("fdr.R")
memory.limit(size = 25)
# Read in the data
rpms_data <- as.data.table(read.table("bpsgt100ku_rpms_stack_cut_filtered_mean_filled.csv", header = TRUE, sep = ","))
zone_data <- as.data.table(read.table("bpsgt100ku_cut.csv", header = TRUE, sep = ","))
# Merge the RPMs data with zone data
rpms_data <- merge( zone_data, rpms_data, by = c("X", "Y"))
print(rpms_data)

drop_zero_rows <- function(data) {
  
  value_band_cols <- grep("^value_band_", names(data), value = TRUE)
  print(value_band_cols)
  # Create a logical vector where TRUE represents rows where all value bands (1 to 39) are zero
  non_zero_rows <- apply(data[, value_band_cols, with = FALSE], 1, function(row) {
    any(row != 0)
  })
  
  # Subset the data to keep only rows where not all value bands are zero
  filtered_data <- data[non_zero_rows, ]
  
  return(filtered_data)
}

# Applying the function to your data
rpms_data <- drop_zero_rows(rpms_data)

# Print the filtered data
#print(rpms_data_filtered)

# Replace zero values with the row mean
rpms_data[, (4:ncol(rpms_data)) := lapply(.SD, function(x) {
  row_mean <- as.integer(rowMeans(rpms_data[, .SD, .SDcols = 4:ncol(rpms_data)], na.rm = TRUE))
  ifelse(x == 0, row_mean, x)
}), .SDcols = 4:ncol(rpms_data)]

nyears <- 39
year <- 1:nyears

# Detect and register parallel backend
n.cores <- detectCores() - 4  # Use available cores minus 4
cl <- makeCluster(n.cores)
print(n.cores)
registerDoParallel(cl)


# Start time for performance monitoring
start.time <- proc.time()

# Perform analysis for each zone separately
results_list <- list()
zones <- unique(rpms_data$Zone) # Replace "Zone" with the actual column name for zones
print(zones)
# Using the cluster for parallel foreach
results_list <- foreach(zone = zones, .combine = "rbind", .packages = c("nlme", "fields", "data.table")) %do% {
  
  # Subset data for current zone
  subset_zone_data <- rpms_data[Zone == zone]
  print(subset_zone_data)
  p1 <- cbind(data.table::melt(subset_zone_data, id.vars = c("X", "Y", "Zone")), as.matrix((rep(1:nyears, each = nrow(subset_zone_data)))))
  colnames(p1) <- c("X", "Y", "Zone", "description", "y", "x")
  #print(p1)
  # Mean and trend models
  p1.mean <- lm(y ~ 1, data = p1)
  p1.trend <- lm(y ~ x, data = p1)
  
  # Sequential inner loop with for loop
  zone_results <- foreach(i = 1:nrow(subset_zone_data), .combine = rbind, .packages = c('nlme', 'data.table')) %dopar% {
    # Extract current row y values
    y <- as.vector(t(subset_zone_data[i, 4:(nyears + 3)])) # Adjust column indices based on your data structure
    
    if (all(y == 0)) {
      return(NULL)  # Skip the current iteration if all values are zero
    }
    
    # Field mean and trend models
    subset_zone_data.mean <- gls(y ~ 1, correlation = corARMA(p = 1), method = "REML")
    subset_zone_data.trend <- gls(y ~ year, correlation = corARMA(p = 1), method = "REML")
    
    degf.trend <- subset_zone_data.trend$dims$N - subset_zone_data.trend$dims$p
    degf.mean <- subset_zone_data.mean$dims$N - subset_zone_data.mean$dims$p
    
    # Create results dataframe for the current row
    z1 <- data.frame(matrix(nrow = 1, ncol = 17)) # Adjust column number to include Zone column
    
    z1[1, 1] <- subset_zone_data[i,1]  # X
    z1[1, 2] <- subset_zone_data[i, 2]  # Y
    z1[1, 3] <- subset_zone_data[i, 3]  # Zone
    z1[1, 4] <- mean(y)  # mean (field)
    z1[1, 5] <- mean(p1$y)  # mean (population) -- make sure p1 is defined outside of foreach
    z1[1, 6] <- z1[1, 4] - z1[1, 5]  # mean difference
    z1[1, 7] <- sqrt(diag(subset_zone_data.mean$varBeta))[1]  # standard error of mean
    z1[1, 8] <- z1[1, 6] / z1[1, 7]  # t-statistic of difference in means
    z1[1, 9] <- 2 * pt(abs(z1[1, 8]), df = degf.mean, lower.tail = FALSE)  # p-value of difference in means
    z1[1, 10] <- NA  # reserved for adjusted mean difference p-value
    
    z1[1, 11] <- subset_zone_data.trend$coefficients[2]  # slope (field)
    z1[1, 12] <- p1.trend$coefficients[2]  # slope (population) -- make sure p1.trend is defined outside of foreach
    z1[1, 13] <- z1[1, 11] - z1[1, 12]  # slope difference
    z1[1, 14] <- sqrt(diag(subset_zone_data.trend$varBeta))[2]  # standard error of slope coefficient
    z1[1, 15] <- z1[1, 13] / z1[1, 14]  # t-statistic of difference in slope coefficients
    z1[1, 16] <- 2 * pt(abs(z1[1, 15]), df = degf.trend, lower.tail = FALSE)  # p-value of difference in slope coefficients
    z1[1, 17] <- NA  # reserved for adjusted slope difference p-value
    
    return(z1)  # Return results for the current row
  }
  rm(p1)
  rm(p1.mean)
  rm(p1.trend)
  return(zone_results) # Return the zone-specific results

}
# Stop the parallel backend
stopCluster(cl)
print(proc.time() - start.time)

# Combine results for all zones
print(results_list)



# Define a function to apply FDR adjustment for each zone
#apply_fdr_by_zone <- function(subset) {
  
 # subset <- subset[, -1, with = FALSE]
  # check to see if subset is null, skip if so
 # if (is.null(subset) || nrow(subset) == 0) {
 #   return(NULL)  # Return NULL if the subset is empty
 # }
  #value = as.double(unlist(subset[, 9]))
 # a1 <- fdr(value, method = "general", adjustment.method = "mean", qlevel = 0.2)
 # b1 <- rep(1, length(subset[, 9]))
 # b1[a1] <- 0
 # new_value = as.double(unlist(subset[, 16]))
 # a2 <- fdr(new_value, method = "general", adjustment.method = "mean", qlevel = 0.2)
  
 # b2 <- rep(1, length(subset[, 16]))
 # b2[a2] <- 0
  #print(b1)
 # print(b2)
  #subset[, 10] <- as.logical(b1)
 # subset[, 17] <- as.logical(b2)
  
 # return(subset)
#}

# Apply the function to each zone
#z1 <- results_list[, apply_fdr_by_zone(.SD), by = X3, .SDcols = c("X3", names(results_list))]

# Drop extra zone column added to front from func apply above
#z1 <- z1[, -1, with = FALSE]

colnames(results_list) <- c("X", "Y", "Zone", "field(mean)", "pop(mean)", "diff(mean)", "SE(mean)", "t(mean)", "p(|t|<T)(mean)", "Adj. Sig. (mean)",
                  "field(slope)", "pop(slope)", "diff(slope)", "SE(slope)", "t(slope)", "p(|t|<T)(slope)", "Adj. Sig. (slope)")

# Save results to CSV
write.table(results_list, file = "comparison_analysis_zones_rpms_full_with_good_zone_res.csv", sep = ",", row.names = FALSE)
