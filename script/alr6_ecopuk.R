library(tidyverse)
library(fuzzyjoin)

# Load libraries
library(data.table)

# 1️⃣ Read and combine navigation files
alr4_files <- list.files("data/ALR/nav_alr6/", pattern = "science.csv", full.names = TRUE)

# fread is faster than read_csv
dat_list <- lapply(alr4_files, function(file) {
  temp <- fread(file)
  temp[, .(timestamp, 
           lon = position_lon, 
           lat = position_lat, 
           depth)]
})

dat <- rbindlist(dat_list, use.names = TRUE, fill = TRUE)

# 2️⃣ Summarize data by rounded lon/lat
dat_to_plot <- dat[, .(
  timestamp = mean(timestamp, na.rm = TRUE),
  depth = mean(depth, na.rm = TRUE)
), by = .(lon = round(lon, 2), lat = round(lat, 2))]

# 3️⃣ Plot
ggplot(dat_to_plot) +
  geom_point(aes(x = lon, y = lat, colour = -depth)) +
  scale_color_viridis_c(name = "Depth") +
  theme_minimal()+
  coord_quickmap()

# 4️⃣ Read fluorometer log
fluo_files <- list.files(
  "data/ALR/eco_puk_alr6/",
  pattern = "fluorometer.*",
  full.names = TRUE
)

# 2️⃣ Read each file with fread()
fluo_list <- lapply(fluo_files, function(f) {
  fread(f, sep = "\t", header = FALSE, skip = 1, fill = TRUE)
})

# 3️⃣ Combine into one big data.table
fluo <- rbindlist(fluo_list, use.names = TRUE, fill = TRUE)

# 5️⃣ Separate timestamp vs science rows
fluo[, flag := ifelse(!is.na(V3), "science", "timestamp")]

fluo_science <- fluo[flag == "science"][c(1:779061)]
fluo_timestamp <- fluo[flag == "timestamp"]

fluo_timestamp[, V1_clean := as.numeric(gsub("[^0-9\\.]", "", V1))]
fluo_science[, timestamp := fluo_timestamp$V1_clean]
fluo_science[, flag := NULL]

# Optional: Rename for clarity
setnames(fluo_science, old = c("V2","V3","V4"), new = c("param1","param2","param3"), skip_absent = TRUE)

setnames(dat, "timestamp", "nav_timestamp")
setkey(dat, nav_timestamp)
setkey(fluo_science, timestamp)

fluo_joined <- dat[fluo_science, roll = "nearest"]

fluo_to_plot <- fluo_joined |> 
  mutate(lon = round(lon, 2), lat = round(lat, 2)) |> 
  group_by(lon, lat) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

ggplot(dat_to_plot) +
  geom_point(aes(x = lon, y = lat, colour = "Track")) +
  geom_point(aes(x = lon, y = lat, colour = "EcoPuk"), data = fluo_to_plot) +
  theme_minimal()+
  coord_quickmap()
ggsave("ALR6_traj.jpg", width = 30, height = 20, units = "cm", dpi = 300)

fluo_joined <- fluo_joined |> mutate(datetime = as.POSIXct(nav_timestamp, origin = "1970-01-01", tz = "UTC"))

ggplot(fluo_joined)+
  geom_line(aes(x = datetime, y = -depth))

fluo_joined <- fluo_joined |> 
  mutate("auv" = "ALR6") |> 
  mutate(fluo = (param3 - 49) * 0.0073,
         beta = (V6 - 49) * 4.513e-7)|> 
  select(nav_timestamp, datetime, lon, lat, depth, fluo, beta)

ggplot(filter(fluo_joined, depth < 100))+
  geom_path(aes(x = datetime, y = -depth, color = fluo))+
  scale_color_viridis_c(name = "Fluorescence (mg m-3)")

write_csv(fluo_joined, "output/ALR6_fluorometer_data.csv")

ggsave("ALR6_transect.jpg", width = 30, height = 20, units = "cm", dpi = 300)