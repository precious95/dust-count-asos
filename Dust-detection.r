library(dplyr)
library(ggplot2)
library(lubridate)
library(riem)
library(gridExtra)
library(cowplot)
library(tidyr)

# Set the working directory and load required libraries
setwd('/home/precious')
# Station IDs for analysis
station_ids <- c("OVE", "FAT", "FCH", "BFL", "EDW", "MHV", "HJO", "MAE", "MCE", "SAC", "SMF", "SCK", "RDD", "BAB", "RBL", "PTV", "VIS", "MYV")

# Station coordinates
df_coords <- data.frame(
  station = station_ids,
  lat = c(39.49, 36.78, 36.7321, 35.4344, 34.90542, 35.059, 36.31139, 36.98486, 37.28603, 38.5069, 38.69542, 37.89417, 40.509, 39.13609, 40.1519, 36.02732, 36.31867, 39.10203),
  lon = c(-121.62, -119.7194, -119.8203, -119.0542, -117.8837, -118.152, -119.6232, -120.1107, -120.5179, -121.495, -121.5908, -121.2383, -122.2934, -121.4366, -122.2536, -119.0629, -119.3929, -121.5688)
)

# Initialize an empty list to store data for each station
stations_data <- list()

# Loop through the station IDs and download the data
for (station_id in station_ids) {
  station_data <- riem_measures(station = station_id, date_start = "2005-01-01", date_end = "2023-12-31")
  # Add the data to the list
  stations_data[[station_id]] <- station_data
}

stations_data = readRDS('stations_data.rds')

# Combine the data for all stations into one dataframe
asos_data <- bind_rows(stations_data)

# Convert wind speed from knots to meters per second (1 knot = 0.514444 m/s)
asos_data <- asos_data %>%
  mutate(wind_speed_mps = sknt * 0.514444)

# Convert visibility from miles to kilometers (1 mile = 1.60934 km)
asos_data <- asos_data %>%
  mutate(visibility_km = vsby * 1.60934)

dust_events <- asos_data %>%
  filter(
    # Standalone dust codes (DU, BLDU, DS)
    (wind_speed_mps >= 6 & visibility_km <= 10 & 
       grepl("\\b(DU|BLDU|DS)\\b", wxcodes, ignore.case = TRUE)) |
      (wind_speed_mps >= 6 & visibility_km <= 10 & 
         grepl("\\b(DU|BLDU|DS)\\b", metar, ignore.case = TRUE)) |
      
      # Standalone HZ without FU
      (wind_speed_mps >= 6 & visibility_km <= 10 & relh < 70 & 
         grepl("\\b(HZ)\\b", wxcodes, ignore.case = TRUE) & 
         !grepl("\\b(FU)\\b", wxcodes, ignore.case = TRUE)) |
      (wind_speed_mps >= 6 & visibility_km <= 10 & relh < 70 & 
         grepl("\\b(HZ)\\b", metar, ignore.case = TRUE) & 
         !grepl("\\b(FU)\\b", metar, ignore.case = TRUE)) |
      
      # HZ with Dust Codes (DU, BLDU, DS), even if FU is present
      (wind_speed_mps >= 6 & visibility_km <= 10 & 
         grepl("\\b(HZ)\\b", wxcodes, ignore.case = TRUE) & 
         grepl("\\b(DU|BLDU|DS)\\b", wxcodes, ignore.case = TRUE)) |
      (wind_speed_mps >= 6 & visibility_km <= 10 & 
         grepl("\\b(HZ)\\b", metar, ignore.case = TRUE) & 
         grepl("\\b(DU|BLDU|DS)\\b", metar, ignore.case = TRUE)) |
      
      # HZ FU with Dust Codes (retain these cases)
      (wind_speed_mps >= 6 & visibility_km <= 10 & 
         grepl("\\b(HZ)\\b", wxcodes, ignore.case = TRUE) & 
         grepl("\\b(FU)\\b", wxcodes, ignore.case = TRUE) & 
         grepl("\\b(DU|BLDU|DS)\\b", wxcodes, ignore.case = TRUE)) |
      (wind_speed_mps >= 6 & visibility_km <= 10 & 
         grepl("\\b(HZ)\\b", metar, ignore.case = TRUE) & 
         grepl("\\b(FU)\\b", metar, ignore.case = TRUE) & 
         grepl("\\b(DU|BLDU|DS)\\b", metar, ignore.case = TRUE))
  )



# Add additional columns for analysis
dust_events <- dust_events %>%
  mutate(valid = as.POSIXct(valid),
         date = as.Date(valid),
         year = year(valid),
         month = month(valid, label = TRUE, abbr = TRUE),
         hour = hour(valid))

# Sort by station and valid time to ensure proper sequencing
dust_events <- dust_events %>% arrange(station, valid)

# Create lag columns for previous event time and date
dust_events <- dust_events %>%
  group_by(station) %>%
  mutate(previous_event_time = lag(valid),
         previous_event_date = lag(date),
         time_diff_hours = as.numeric(difftime(valid, previous_event_time, units = "hours")))

# Identify events near midnight (between 11 PM and 1 AM)
dust_events <- dust_events %>%
  mutate(near_midnight = hour(valid) >= 23 | hour(valid) <= 1)

# Identify if the event is a continuation from the previous day
dust_events <- dust_events %>%
  mutate(is_continuation = near_midnight & (previous_event_date == date - 1))

# Identify events in the morning (before noon) and evening (after noon)
dust_events <- dust_events %>%
  mutate(morning_event = hour(valid) < 12,
         evening_event = hour(valid) >= 12)

# Identify if there are two separate events on the same day (one morning and one evening)
dust_events <- dust_events %>%
  group_by(station, date) %>%
  mutate(is_double_event = sum(morning_event) > 0 & sum(evening_event) > 0) %>%
  ungroup()

# Define new events
dust_events <- dust_events %>%
  arrange(station, valid) %>%
  group_by(station) %>%
  mutate(
    new_event = case_when(
      !is_continuation ~ TRUE,
      is_double_event & evening_event & lag(morning_event, default = FALSE) ~ TRUE,
      TRUE ~ FALSE
    ),
    event_id = cumsum(new_event)
  ) %>%
  ungroup()

# Assign the month to each event based on the event's start time
event_summary <- dust_events %>%
  group_by(station, event_id) %>%
  summarise(event_start = min(valid)) %>%
  ungroup() %>%
  mutate(month = month(event_start, label = TRUE, abbr = TRUE))

# Ensure 'date' and 'month' columns are present
dust_events <- dust_events %>%
  mutate(
    date = as.Date(valid),
    month = month(valid, label = TRUE, abbr = TRUE)
  )

# Count unique dates with dust events per month
monthly_summary_by_station <- dust_events %>%
  distinct(station, date, month) %>%
  group_by(station, month) %>%
  summarise(count = n(), .groups = 'drop')

# Calculate the average monthly dust events across all stations
monthly_average_summary <- monthly_summary_by_station %>%
  group_by(month) %>%
  summarise(avg_count = mean(count, na.rm = TRUE), .groups = 'drop')

# View the monthly average summary
print(monthly_average_summary)


# Loop for California map and station coordinates
california_map <- map_data("state", region = "california")

central_valley_counties <- c("shasta", "tehama", "glenn", "butte", "colusa", 
                             "sutter", "yuba", "yolo", "sacramento", "san joaquin", 
                             "stanislaus", "merced", "madera", "fresno", "kings", 
                             "tulare", "kern")
california_counties <- map_data("county", region = "california")

# Filter for Central Valley counties
central_valley_map <- california_counties %>%
  filter(subregion %in% central_valley_counties)

ca_map_plots <- list()
for (i in 1:nrow(df_coords)) {
  station_coords <- df_coords[i, ]
  ca_map_plot <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "blue") +
    geom_polygon(data = central_valley_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_point(data = station_coords, aes(x = lon, y = lat), color = "red", size = 3) +
    coord_fixed(1.3) +
    theme_void()
  ca_map_plots[[station_coords$station]] <- ca_map_plot
}


# Step 1: Monthly Dust Events
monthly_plots <- list()
for (station_id in station_ids) {
  station_monthly_summary <- monthly_summary_by_station %>%
    filter(station == station_id)
  
  monthly_plot <- ggplot(station_monthly_summary, aes(x = month, y = count)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(title = paste("Monthly Dust Events for Station", station_id), x = "Month", y = "Count")
  
  combined_plot <- ggdraw() +
    draw_plot(monthly_plot)
  
  monthly_plots[[station_id]] <- combined_plot
}
# Step 2: Yearly Dust Events
yearly_avg_summary <- dust_events %>%
  group_by(station, year) %>%
  summarise(count = n_distinct(date), .groups = 'drop')

yearly_plots <- list()
for (station_id in station_ids) {
  station_yearly_summary <- yearly_avg_summary %>%
    filter(station == station_id)
  
  yearly_plot <- ggplot(station_yearly_summary, aes(x = year, y = count)) +
    geom_bar(stat = "identity", fill = "black", alpha = 0.7) +
    geom_line(color = "blue", size = 0.5, linetype = "dotted") +
    geom_point(color = "blue", size = 0.5) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    labs(title = paste("Yearly Dust Events for Station", station_id), x = "Year", y = "Count") +
    theme_minimal() +
    theme(text = element_text(size = 17),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  combined_plot <- ggdraw() +
    draw_plot(yearly_plot)
  
  yearly_plots[[station_id]] <- combined_plot
}

library(tidyr)
# Step 3: Duration of Dust Events
dust_event_duration <- dust_events %>%
  group_by(station, event_date = floor_date(valid, "day")) %>%
  summarise(start_time = min(valid), end_time = max(valid), .groups = 'drop') %>%
  mutate(duration_hours = as.numeric(difftime(end_time, start_time, units = "hours"))) %>%
  mutate(duration_group = case_when(
    duration_hours >= 0 & duration_hours <= 1 ~ "1",
    duration_hours > 1 & duration_hours <= 2 ~ "2",
    duration_hours > 2 & duration_hours <= 3 ~ "3",
    duration_hours > 3 & duration_hours <= 4 ~ "4",
    duration_hours > 4 & duration_hours <= 5 ~ "5",
    duration_hours > 5 & duration_hours <= 6 ~ "6",
    duration_hours > 6 & duration_hours <= 7 ~ "7",
    duration_hours > 7 & duration_hours <= 8 ~ "8",
    duration_hours > 8 & duration_hours <= 9 ~ "9",
    duration_hours > 9 ~ ">10"
  )) %>%
  mutate(duration_group = factor(duration_group, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", ">10")))

duration_plots <- list()
for (station_id in station_ids) {
  station_duration_summary <- dust_event_duration %>%
    filter(station == station_id) %>%
    group_by(duration_group) %>%
    summarise(count = n(), .groups = 'drop')
  
  duration_plot <- ggplot(station_duration_summary, aes(x = duration_group, y = count)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(title = paste("Duration of Dust Events for Station", station_id), x = "Duration (Hours)", y = "Count") +
    theme_minimal() +
    theme(text = element_text(size = 17),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  combined_plot <- ggdraw() +
    draw_plot(duration_plot)+
    draw_plot(ca_map_plots[[station_id]], x = 0.7, y = 0.6, width = 0.3, height = 0.3)
  
  duration_plots[[station_id]] <- combined_plot
}
# Step 4: Diurnal Distribution of Dust Events
# Step 1: Check for duplicate entries based on the 'valid' timestamp
duplicates <- sum(duplicated(dust_events$valid))
print(paste("Number of duplicate entries:", duplicates))

# Step 2: Check for missing or irregular timestamps
missing_times <- sum(is.na(dust_events$valid))
print(paste("Number of missing or irregular timestamps:", missing_times))

# Step 3: Ensure every unique event is counted across stations
dust_events_clean <- dust_events %>%
  mutate(local_time = with_tz(valid, tzone = "America/Los_Angeles")) %>%
  mutate(hour = hour(local_time))  # Extract hour directly

# Step 4: Summarize the number of dust events starting at each hour for each station
diurnal_distribution_by_station <- dust_events_clean %>%
  group_by(station, hour) %>%
  summarise(count = n(), .groups = 'drop')  # Count all unique events per hour per station

# Step 5: Ensure all hours from 0 to 23 are present for each station
diurnal_distribution_by_station <- diurnal_distribution_by_station %>%
  complete(station, hour = 0:23, fill = list(count = 0))

# Step 6: Summarize the average across all stations
diurnal_avg_distribution <- diurnal_distribution_by_station %>%
  group_by(hour) %>%
  summarise(avg_count = mean(count, na.rm = TRUE), .groups = 'drop')

# View the diurnal distribution averaged across stations
print(diurnal_avg_distribution)

# Diurnal plots
diurnal_plots <- list()
for (station_id in station_ids) {
  station_diurnal_summary <- diurnal_distribution_by_station %>%
    filter(station == station_id)
  
  diurnal_plot <- ggplot(station_diurnal_summary, aes(x = factor(hour), y = count)) +
    geom_bar(stat = "identity", fill = "black") +
    labs(title = paste("Diurnal Duration of Dust Events for Station", station_id), x = "Local Time (Hours)", y = "Count") +
    theme_minimal() +
    theme(text = element_text(size = 17),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  combined_plot <- ggdraw() +
    draw_plot(diurnal_plot)
  
  diurnal_plots[[station_id]] <- combined_plot
}

# Plot as grid matrix for each 3 by 6
png("dt3.png", width = 3000, height = 2000)
grid.arrange(grobs = c(monthly_plots), nrow = 3, ncol = 6)
grid.arrange(grobs = c(yearly_plots), nrow = 3, ncol = 6)
grid.arrange(grobs = c(duration_plots), nrow = 3, ncol = 6)
grid.arrange(grobs = c(diurnal_plots), nrow = 3, ncol = 6)
dev.off()
