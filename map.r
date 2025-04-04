# Load required libraries
library(ggplot2)
library(maps)
library(dplyr)

# Define station IDs and their coordinates
station_ids <- c("OVE", "FAT", "FCH", "BFL", "EDW", "MHV", "HJO", "MAE", "MCE", 
                 "SAC", "SMF", "SCK", "RDD", "BAB", "RBL", "PTV", "VIS", "MYV")

df_coords <- data.frame(
  station = station_ids,
  lat = c(39.49, 36.78, 36.7321, 35.4344, 34.90542, 35.059, 36.31139, 36.98486, 
          37.28603, 38.5069, 38.69542, 37.89417, 40.509, 39.13609, 40.1519, 
          36.02732, 36.31867, 39.10203),
  lon = c(-121.62, -119.7194, -119.8203, -119.0542, -117.8837, -118.152, 
          -119.6232, -120.1107, -120.5179, -121.495, -121.5908, -121.2383, 
          -122.2934, -121.4366, -122.2536, -119.0629, -119.3929, -121.5688)
)

# Get map data for California (state and counties)
california_map <- map_data("state", region = "california")
california_counties <- map_data("county", region = "california")

# Define Central Valley counties
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa", 
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin", 
                             "stanislaus", "merced", "madera", "fresno", "kings", 
                             "tulare", "kern")

# Filter county data for the Central Valley
central_valley_map <- california_counties %>%
  filter(subregion %in% central_valley_counties)

# Plot all stations on one map
ggplot() +
  # Outline of California
    # Outline of the Central Valley counties
  geom_polygon(data = central_valley_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  # Plot all stations
  geom_point(data = df_coords, aes(x = lon, y = lat),
             color = "red", size = 3) +
  coord_fixed(1.3) +
  theme_void() +
  ggtitle("")









library(ggplot2)
library(maps)
library(dplyr)

# Define station IDs and their coordinates
station_ids <- c("OVE", "FAT", "FCH", "BFL", "EDW", "MHV", "HJO", "MAE", "MCE", 
                 "SAC", "SMF", "SCK", "RDD", "BAB", "RBL", "PTV", "VIS", "MYV")

df_coords <- data.frame(
  station = station_ids,
  lat = c(39.49, 36.78, 36.7321, 35.4344, 34.90542, 35.059, 36.31139, 36.98486, 
          37.28603, 38.5069, 38.69542, 37.89417, 40.509, 39.13609, 40.1519, 
          36.02732, 36.31867, 39.10203),
  lon = c(-121.62, -119.7194, -119.8203, -119.0542, -117.8837, -118.152, 
          -119.6232, -120.1107, -120.5179, -121.495, -121.5908, -121.2383, 
          -122.2934, -121.4366, -122.2536, -119.0629, -119.3929, -121.5688)
)

# Get map data for the Central Valley
california_counties <- map_data("county", region = "california")
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa", 
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin", 
                             "stanislaus", "merced", "madera", "fresno", "kings", 
                             "tulare", "kern")

central_valley_map <- california_counties %>%
  filter(subregion %in% central_valley_counties)

# Plot stations on the Central Valley map
ggplot() +
  # Outline of the Central Valley counties
  geom_polygon(data = central_valley_map, 
               aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  # Plot all station points
  geom_point(data = df_coords, aes(x = lon, y = lat),
             color = "red", size = 2) +
  # Add station ID labels near the points
  geom_text(data = df_coords, 
            aes(x = lon, y = lat, label = station),
            vjust = -1,            # moves text slightly above each point
            color = "black", 
            size = 2) +
  coord_fixed(1.3) +
  theme_void() +
  ggtitle("ASOS Stations")



library(ggplot2)
library(maps)
library(dplyr)

# Read the CSV file and extract unique station entries
pm_data <- read.csv("pm.csv", stringsAsFactors = FALSE)

# Select unique station entries using columns: ID, Latitude, Longitude
df_coords <- unique(pm_data[, c("ID", "Latitude", "Longitude")])

# Rename columns to match plotting code
df_coords <- df_coords %>%
  rename(station = ID,
         lat = Latitude,
         lon = Longitude)

# Get map data for the Central Valley counties in California
california_counties <- map_data("county", region = "california")
central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa", 
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin", 
                             "stanislaus", "merced", "madera", "fresno", "kings", 
                             "tulare", "kern")

central_valley_map <- california_counties %>%
  filter(subregion %in% central_valley_counties)

# Plot all stations on the Central Valley map with station IDs as labels
ggplot() +
  # Outline of the Central Valley counties
  geom_polygon(data = central_valley_map, 
               aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  # Plot all station points
  geom_point(data = df_coords, aes(x = lon, y = lat),
             color = "green",size = 2) +
  # Add station ID labels near the points
  geom_text(data = df_coords, 
            aes(x = lon, y = lat, label = station),
            vjust = -1,   # moves text slightly above each point
            color = "black", 
            size = 0) +
  coord_fixed(1.3) +
  theme_void() +
  ggtitle("PM Stations")

