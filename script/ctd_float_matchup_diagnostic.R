library(tidyverse)
library(arrow)
library(ncdf4)
library(stringi)

dy180_wmo <- c(6990636, 3901581, 4903659)
jc269_wmo <- c(3902681, 3901586, 1902695)

dat <- read_parquet("data/argo_pq/biocarbon_floats_table.parquet") |> select(PLATFORM_NUMBER, JULD, LONGITUDE, LATITUDE) |> 
  mutate(wmo = as.numeric(PLATFORM_NUMBER)) |> 
  filter(wmo %in% c(dy180_wmo, jc269_wmo))

unique(dat$PLATFORM_NUMBER)


# open ctd and make a summary ---------------------------------------------
jc269_files <- list.files("data/CTD/JC269/csv_files", pattern = "*.csv", full.names = TRUE)
jc269_loc <- data.frame()
for(file in jc269_files){
  dat <- readLines(file)
  
  lon <- as.numeric(stringi::stri_extract_first_regex(dat[grep("LONGITUDE = ", dat)], "-?\\d+\\.?\\d*"))
  lat <- as.numeric(stringi::stri_extract_first_regex(dat[grep("LATITUDE = ", dat)], "\\d+\\.?\\d*"))
  
  datetime <- lubridate::ymd_hm(
    paste(
      stringi::stri_extract_first_regex(dat[grep("DATE = ", dat)], "\\d{8}"),
      stringi::stri_extract_first_regex(dat[grep("TIME = ", dat)], "\\d{4}")
    )
  )                  
  profile <- as.numeric(stringi::stri_extract_first_regex(dat[grep("STNNBR = ", dat)], "\\d+"))
  ctd_loc <- data.frame("profile" = profile,
                        "datetime" = datetime,
                        "lon" = lon,
                        "lat" = lat)
  jc269_loc <- bind_rows(jc269_loc, ctd_loc)
}

dy180_files <- list.files("data/CTD/DY180/74EQ20240522_ct1", pattern = "*.csv", full.names = TRUE)
dy180_loc <- data.frame()

for(file in dy180_files){
  dat <- readLines(file)
  
  lon <- as.numeric(stringi::stri_extract_first_regex(dat[grep("LONGITUDE = ", dat)], "-?\\d+\\.?\\d*"))
  lat <- as.numeric(stringi::stri_extract_first_regex(dat[grep("LATITUDE = ", dat)], "\\d+\\.?\\d*"))
  
  datetime <- lubridate::ymd_hm(
    paste(
      stringi::stri_extract_first_regex(dat[grep("DATE = ", dat)], "\\d{8}"),
      stringi::stri_extract_first_regex(dat[grep("TIME = ", dat)], "\\d{4}")
    )
  )                  
  profile <- as.numeric(stringi::stri_extract_first_regex(dat[grep("STNNBR = ", dat)], "\\d+"))
  ctd_loc <- data.frame("profile" = profile,
                        "datetime" = datetime,
                        "lon" = lon,
                        "lat" = lat)
  dy180_loc <- bind_rows(dy180_loc, ctd_loc)
}

jc269_loc$cruise <- "JC269"
dy180_loc$cruise <- "DY180"
ctd_localisation <- bind_rows(jc269_loc, dy180_loc)

ggplot(ctd_localisation)+
  geom_point(aes(x = datetime, y = cruise))

dat <- dat |> distinct()

ggplot(dat)+
  geom_point(aes(x = JULD, y = PLATFORM_NUMBER))+
  geom_vline(aes(xintercept = datetime, color = cruise), data = ctd_localisation)+
  theme_minimal()+
  ylab("Float WMO")+
  xlab("Date")+
  scale_color_brewer(palette = "Set1", name = "CTD cast")

#ggsave("output/plots/ctd_float_time_summary_non_biocarbon.png", dpi = 300, width = 22, height = 15, units = "cm")

float_localisation <- dat |> 
  mutate(type = "Float", profile_id = paste(wmo, substr(JULD, 1, 10), sep = "_")) |> 
  rename("datetime" = JULD, "lon" = LONGITUDE, "lat" = LATITUDE, "mission" = wmo) |> 
  select(-PLATFORM_NUMBER) |> 
  mutate(mission = as.character(mission))

ctd_localisation <- ctd_localisation |> 
  mutate(profile_id = paste(cruise, profile, sep = "_"),
         type = "CTD") |> 
  rename("mission" = cruise) |> 
  select(, - profile)

all_position <- bind_rows(ctd_localisation, float_localisation)

library(geosphere)  # for distance calculation

# Example: all_position dataframe already loaded

# Separate CTD and Float profiles, rename to avoid clashes, and drop 'type'
ctd <- all_position %>%
  filter(type == "CTD") %>%
  select(-type) %>%
  rename_with(~ paste0(., "_ctd"))

floats <- all_position %>%
  filter(type == "Float") %>%
  select(-type) %>%
  rename_with(~ paste0(., "_float"))

# Cross CTD with floats
matches <- expand_grid(ctd, floats) %>%
  mutate(time_diff = abs(difftime(datetime_ctd, datetime_float, units = "hours"))) %>%
  filter(time_diff <= 72) %>%
  group_by(profile_id_ctd, mission_float) %>%
  slice_min(order_by = time_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    spatial_dist_km = distHaversine(
      cbind(lon_ctd, lat_ctd),
      cbind(lon_float, lat_float)
    ) / 1000
  )

ggplot(matches, aes(x = as.numeric(time_diff), y = spatial_dist_km)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Time vs Spatial Distance (per float mission)",
    x = "Time Difference (hours)",
    y = "Spatial Distance (km)",
    color = "Float Mission (WMO)"
  ) +
  theme_minimal()

matches_filtered <- matches |> filter(spatial_dist_km < 100)

ggplot(matches_filtered, aes(x = as.numeric(time_diff), y = spatial_dist_km)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Time vs Spatial Distance (per float mission)",
    x = "Time Difference (hours)",
    y = "Spatial Distance (km)",
    color = "Float Mission (WMO)"
  ) +
  theme_minimal()


ggplot(matches_filtered, aes(x = as.numeric(time_diff))) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Time Differences (CTD vs Float)",
    x = "Time Difference (hours)",
    y = "Count"
  ) +
  theme_minimal()

position_to_share <- float_localisation |> filter(mission %in% c(4903659, 6990636))

#plot position to share with coastline
ggplot()+
  borders("world", colour = "gray85", fill = "gray80")+
  geom_path(data = position_to_share, aes(x = lon, y = lat, group = mission))+
  geom_point(data = position_to_share, aes(x = lon, y = lat, color = mission), size = 2)+
  coord_quickmap(xlim = c(-33, -10), ylim = c(58, 64))+
  theme_minimal()+
  labs(title = "Float Positions", x = "Longitude", y = "Latitude", color = "Float WMO")

ggsave("output/plots/float_positions_to_share.png", dpi = 300, width = 20, height = 15, units = "cm")

write_csv(position_to_share, "output/float_positions_4903659_and_6990636.csv")

# Prepare matches with time in hours# Prepafilter_()re matches with time in hours
matches_binned <- matches_filtered %>%
  mutate(
    time_hours = as.numeric(time_diff, units = "hours"),
    time_bin = cut(
      time_hours,
      breaks = c(0, 3, 6, 12, 24, 48, 72),
      labels = c("3", "6", "12", "24", "48", "<72"),
      right = TRUE
    ),
    dist_bin = cut(
      spatial_dist_km,
      breaks = c(0, 1, 5, 10, 20, 50, 100),
      labels = c("1", "5", "10", "20", "50", "<100"),
      right = TRUE
    )
  )

heatmap_data <- matches_binned %>%
  count(dist_bin, time_bin)

# Plot heatmap with counts in cells
ggplot(heatmap_data, aes(x = time_bin, y = dist_bin, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(
    title = "CTD–Float Matchup Density",
    x = expression(Delta~time~(hours)),
    y = "distance (km)",
    fill = "Count"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

summary_table <- matches %>%
  transmute(
    profile_id_ctd,
    cruise = mission_ctd,
    datetime_ctd,
    wmo = mission_float,
    profile_id_float,
    datetime_float,
    time_diff_hours = round(as.numeric(time_diff, units = "hours"), 2),
    dist_km = round(spatial_dist_km, 1)
  ) |> 
  select(-profile_id_float)

write_excel_csv(summary_table, "output/CTD_float_matchup_summary.csv")



# time and space interpolation --------------------------------------------
library(zoo)

floats <- all_position %>%
  filter(type == "Float") %>%
  mutate(datetime = floor_date(datetime, "hour")) %>%  # round down to hour
  arrange(mission, datetime)

# Build hourly grid per mission
float_grid <- floats %>%
  group_by(mission) %>%
  summarise(start = min(datetime), end = max(datetime), .groups = "drop") %>%
  rowwise() %>%
  mutate(datetime = list(seq(from = start, to = end, by = "1 hour"))) %>%
  unnest(datetime)

# Interpolate lon/lat onto hourly grid
floats_interp <- floats %>%
  select(mission, datetime, lon, lat) %>%
  right_join(float_grid, by = c("mission", "datetime")) %>%
  arrange(mission, datetime) %>%
  group_by(mission) %>%
  mutate(
    lon = na.approx(lon, x = as.numeric(datetime), na.rm = FALSE),
    lat = na.approx(lat, x = as.numeric(datetime), na.rm = FALSE)
  ) %>%
  ungroup()

ggplot(arrange(floats_interp, datetime))+
  geom_point(aes(x = lon, y = lat, color = mission))+
  coord_quickmap()

floats_interp_to_match <- floats_interp |> 
  rename_with(~ paste0(., "_float"))
matches <- expand_grid(ctd, floats_interp_to_match) %>%
  mutate(time_diff = abs(difftime(datetime_ctd, datetime_float, units = "hours"))) %>%
  filter(time_diff <= 12) %>%
  group_by(profile_id_ctd, mission_float) %>%
  slice_min(order_by = time_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    spatial_dist_km = distHaversine(
      cbind(lon_ctd, lat_ctd),
      cbind(lon_float, lat_float)
    ) / 1000
  )


matches_filter <- filter(matches, spatial_dist_km < 1000)

ggplot(matches)+
  geom_line(aes(x = datetime_ctd, y = spatial_dist_km, colour = mission_float))+
  geom_point(aes(x = datetime_ctd, y = spatial_dist_km))+
  facet_wrap(.~ mission_ctd, scales = "free", ncol = 1)+
  theme_bw()+
  ggtitle("Distance of the CTD from the theoritical float position")+
  ylab("Distance from the float (km)")

ggsave("output/plots/ctd_float_distance_non_biocarbon.png", dpi = 300, width = 20, height = 25, units = "cm")

stn_id <- read_csv("data/station_ctd_id.csv")

stn_id <- stn_id |> mutate(profile_id_ctd = paste(Cruise, ctd, sep = "_"))

matches <- left_join(matches, stn_id)

matches_boxplot <- matches |> 
  filter(!is.na(Location) & mission_float != 3901586)

ggplot(matches_boxplot)+
  geom_boxplot(aes(x = stn, y = spatial_dist_km, fill = Location))+
  theme_bw()+
  xlab("Station name")
ggsave("output/plots/station_to_floats_distance.png", dpi = 300, width = 30, height = 25, units = "cm")


ggplot(matches_boxplot)+
  geom_boxplot(aes(x = stn, y = spatial_dist_km, fill = mission_float), position = "dodge")+
  theme_bw()+
  xlab("Station name")

ggplot(matches_boxplot)+
  geom_boxplot(aes(x = Location, y = spatial_dist_km, fill = Location))+
  theme_bw()+
  xlab("Station name")+
  facet_wrap(.~mission_float)


ggplot(matches_boxplot)+
  geom_point(aes(x = stn, y = spatial_dist_km, color = Location), size = 2)+
  theme_bw()+
  xlab("Station name")+
  facet_wrap(.~ mission_float, ncol = 1)
ggsave("output/plots/floats_distance_to_stations.png", dpi = 300, width = 30, height = 25, units = "cm")

