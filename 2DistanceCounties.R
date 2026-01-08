
##################################################
###   Calculating distances between counties   ###
##################################################

# INPUT DATASETS:
#  - data.csv : county-level dataset with GeoFIPS identifiers 
#  - US Census county shapefiles (tigris::counties)

# OUTPUT DATASETS:
#  - county_distance_matrix_haversine.csv : normalized pairwise county distance matrix
#  - county_weights_matrix_haversine.csv  : spatial weight matrix based on exponential decay
#  - missing_counties_dist.csv : mapping of missing counties to nearest available counties

# NOTES:
#  - Distances are computed using county centroids and the Haversine formula
#  - Distance rows are normalized to sum to one
#  - Spatial weights are computed as exp(-b * distance), with b fixed
#  - Non-continental states and territories are excluded

# Clear the workspace
rm(list=ls())

# Load spatial and plotting libraries
suppressWarnings(
    suppressMessages({
      library(ggplot2)
      library(tigris)
      library(sf)
      library(dplyr)
      library(scales)
}))  

# Read the final analysis dataset and format GeoFIPS as GEOID
data <- read.csv("data.csv")
data$GEOID <- sprintf("%05d", data$GeoFIPS) 
data <- data %>% arrange(GEOID)

GEOID <- sort(unique(data$GEOID))

# Load counties shapefile
counties <- counties(cb = TRUE)
# Load US county geometries and exclude non-continental areas
counties <- counties[!counties$STUSPS %in% c('AK', 'AS', 'MP', 'DC', 
                                             'PR', 'RI', 'HI', 'VI', 'GU'),]
# Counties-to-Planning Regions Approximation
counties <- counties[!(counties$STATEFP == '09' & counties$GEOID == '09120'), ]
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09110'] <- '09013'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09170'] <- '09009'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09180'] <- '09011'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09130'] <- '09007'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09190'] <- '09001'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09140'] <- '09003'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09160'] <- '09005'
counties$GEOID[counties$STATEFP=='09' & counties$GEOID=='09150'] <- '09015' 

# Keep only counties present in the analysis dataset
countiesData <- counties[counties$GEOID %in%GEOID, ]
# if countiesData = counties, comment the previous line and : 
# countiesData <- counties

# Convert spatial data to sf format and ensure consistent ordering
countiesData_sf <- st_as_sf(countiesData)
countiesData_sf <- countiesData_sf %>% arrange(GEOID)

# Get the centroids (coordinates) of each county
centroids <- st_centroid(countiesData_sf)
coords <- st_coordinates(centroids)

# Convert degrees to radians
deg2rad <- function(deg) {
  deg * (pi / 180)
}

# Calculate pairwise distances between counties using the Haversine formula
haversine <- function(lat1, lon1, lat2, lon2) {
  to_radians <- function(deg) deg * pi / 180
  lat1 <- to_radians(lat1)
  lon1 <- to_radians(lon1)
  lat2 <- to_radians(lat2)
  lon2 <- to_radians(lon2)
  
  R <- 6371
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R * c
}

# Approximate Euclidean distance in kilometers between two geographic points
euclidean <- function(lat1, lon1, lat2, lon2) {
  
  lat1 <- as.numeric(lat1)
  lon1 <- as.numeric(lon1)
  lat2 <- as.numeric(lat2)
  lon2 <- as.numeric(lon2)
  
  sqrt((lat2 - lat1)^2 + (lon2 - lon1)^2)
}


# Convert coordinates to radians
coords_rad <- deg2rad(coords)

# Vectorized calculation of all pairwise distances using outer() function
distance_matrix <- outer(1:nrow(coords_rad), 1:nrow(coords_rad), 
             Vectorize(function(i, j) haversine(coords_rad[i, 1], 
             coords_rad[i, 2], coords_rad[j, 1], coords_rad[j, 2])))

normalize_rows <- function(mat) {  mat / rowSums(mat) }
distance_matrix <- normalize_rows(distance_matrix)

# Convert the matrix to a data frame and add county names
county_names <- countiesData_sf$GEOID
distance.df <- as.data.frame(distance_matrix)
rownames(distance.df) <- colnames(distance.df) <- county_names

# Save the distance matrix to a CSV file
write.csv(distance.df, "county_distance_matrix_haversine.csv")
# write.csv(distance.df, "county_ALL_distance_matrix_haversine.csv")

# Define the exponential distance-decay weighting function
weights <- function(x, b){ exp( -b*x) }
weights.df <- apply(distance.df, 1, weights,  b=36.7)
rownames(weights.df) <- colnames(weights.df) <- county_names

apply(apply(weights.df, 2, summary), 1, summary)

write.csv(weights.df, "county_weights_matrix_haversine.csv")


#####################################################
###     Variogram and empirical semivariances     ###
##################################################### 

vg <- list()
iter <- 0
for(year.selec in sort(unique(data$Year))){
  iter <- iter +1
  
  # Subset the data to the year of interest
  data.year <- data[data$Year==year.selec, ]
  
  # Identify the GEOIDs common to both the data and the distance matrix
  common_ids <- intersect(rownames(weights.df), data.year$GEOID)
  
  # Restrict the distance matrix to the observed areas only
  distance.sub <- distance.df[common_ids, common_ids]
  
  # Extract the response variable (Median AQI) and align it with the distance matrix
  y.year <- data.year[match(rownames(distance.sub), data.year$GEOID), "Median.AQI"]
  
  # Extract pairwise distances (upper triangular part only)
  w_ij <- as.matrix(distance.sub)[upper.tri(distance.sub)]
  
  # Compute squared differences between observations for all area pairs
  dy2 <- (outer(y.year, y.year, "-")^2)[upper.tri(distance.sub)]
  
  # Standardize the response variable (optional, for scale invariance)
  y.year <- as.numeric(scale(y.year))
  
  # Define distance classes (bins) for the empirical variogram
  nbins  <- 20
  breaks <- seq(0, max(w_ij), length.out = nbins + 1)
  bins   <- cut(w_ij, breaks = breaks, include.lowest = TRUE)
  
  # Compute the empirical semivariance for each distance class
  vg[[iter]] <- aggregate(dy2 / 2, by = list(bins), mean)
  colnames(vg[[iter]]) <- c("distance_class", "semivariance")
  
}

cols <- gray.colors(8, start = 0.8, end = 0.1) 
quartz(width = 6, height = 6) 
ggplot() +
  geom_line(data = bind_rows(lapply(seq_along(vg), function(i) {
    data.frame(
      DistanceClass = seq_along(vg[[i]]$semivariance),
      Semivariance = vg[[i]]$semivariance,
      Year = factor(sort(unique(data$Year))[i])
    ) })), 
    aes(x = DistanceClass, y = Semivariance, group = Year, color = Year), size = 1.2) +
  geom_point(data = bind_rows(lapply(seq_along(vg), function(i) {
    data.frame(DistanceClass = seq_along(vg[[i]]$semivariance),
               Semivariance = vg[[i]]$semivariance,
               Year = factor(sort(unique(data$Year))[i]) )
  })), aes(x = DistanceClass, y = Semivariance, color = Year), size = 3) +
  scale_color_manual(values = cols) +
  labs(title = "", x = "Distance class", y = "Semivariance") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = 'None',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )


####################################
###    Plotting countiesData     ###
####################################

# Load libraries
suppressWarnings(suppressMessages({
  library(ggplot2)
  library(ggspatial)
  library(viridis)
}))  

countiesData <- counties(cb = TRUE)
countiesData_sf <- st_as_sf(countiesData)
countiesData_sf <- countiesData_sf %>% arrange(GEOID)

# county_names[154]: Palm Beach (Florida) (GEOID 12099)
# county names[58] : Los Angeles (California) (GEOID 06037)

index <-  154

countiesData_sf <- merge(countiesData_sf, data.frame('distance' = distance.df[, index], 
                         'weights' = weights.df[,index], 'GEOID' = county_names), 
                         by='GEOID', all.x=T)

# Plot the map with colors based on the GeoFIPS data
quartz(width = 9, height = 6) 
ggplot(data = countiesData_sf) +    
  geom_sf(aes(fill = distance), color = NA, size = 0) +    
  scale_fill_viridis_c(option = "plasma", name = "Scale", na.value = "grey",
                       labels = scales::label_number(accuracy = 0.0001)) +     
  theme_minimal() +    
  labs(title = " ",
       subtitle = " ",
       fill = "Variable") +    
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2) 

# ggsave('PalmBeach.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)
# ggsave('LosAngeles.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)

quartz(width = 9, height = 6) 
ggplot(data = countiesData_sf) +    
  geom_sf(aes(fill = weights), color = NA, size = 0) +    
  scale_fill_viridis_c(option = "plasma", name = "Weight Sc.", na.value = "grey",
        direction = -1, labels = scales::label_number(accuracy = 0.001)) +     
  theme_minimal() +    
  labs(title = " ",
       subtitle = " ",
       fill = "Variable") +    
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2) 


####################################
###   Reading ALL distances     ###
####################################

# Read full county distance matrix including missing counties
county_ALL_distance <- read.csv('county_ALL_distance_matrix_haversine.csv', sep=',', row.names=1)
colnames(county_ALL_distance)  <- gsub("^X", "", colnames(county_ALL_distance ))
rownames(county_ALL_distance) <- sprintf("%05d", as.integer(rownames(county_ALL_distance)))

# Identify counties missing from the analysis dataset
colnamesST <- substr(colnames(county_ALL_distance), 1, 2)

# Match missing counties to the closest available county, prioritizing same-state matches
positions <- match(GEOID, rownames(county_ALL_distance))

missing_positions <- setdiff(1:dim(county_ALL_distance)[1], positions)
missing_positions.name <- colnames(county_ALL_distance)[missing_positions]
missing_positions.nameST <- substr(missing_positions.name, 1, 2)

positions[is.na(positions)] <- round(runif(length(positions[is.na(positions)]), 308, 319))
positions.name <- colnames(county_ALL_distance)[positions]

missing_positions.min <- rep(NA, length(missing_positions))
counter <- 0
for(i in missing_positions){
  counter=counter+1
  aux.mini <- (county_ALL_distance[positions,i])
  aux.mini[i] <- aux.mini[i] + 10
  aux.mini[ colnamesST[positions] !=  missing_positions.nameST[counter] ] <- 
    aux.mini[ colnamesST[positions] !=  missing_positions.nameST[counter] ]  + 10
  missing_positions.min[counter] <- which.min(aux.mini)
} 
positions.name.min <- positions.name[missing_positions.min]

missing.counties.dist <- data.frame(missing_positions.name, positions.name.min )
missing.counties.dist$missing_positions.name.ST <- substr(missing_positions.name, 1, 2)
missing.counties.dist$positions.name.min.ST <- substr(positions.name.min, 1, 2)
missing.counties.dist$id <- 0
missing.counties.dist$id[missing.counties.dist$missing_positions.name.ST!=
                           missing.counties.dist$positions.name.min.ST] <- 1

# Save mapping of missing counties to nearest substitutes
write.csv(missing.counties.dist, "missing_counties_dist.csv")

