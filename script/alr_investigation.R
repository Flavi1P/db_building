library(tidyverse)

alr4 <- read_csv("output/ALR4_fluorometer_data.csv")
alr6 <- read_csv("output/ALR6_fluorometer_data.csv")

dat <- bind_rows(alr4, alr6)

dat_to_plot <- dat |> 
  mutate(lon = round(lon, 2), lat = round(lat, 2)) |> 
  group_by(lon, lat, auv) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

ggplot(dat_to_plot)+
  geom_point(aes(x = lon, y = lat, color = auv))+
  coord_quickmap()

ggplot(alr4)+
  geom_line(aes(x = datetime, y = -depth))+
  ggtitle("ALR4")

ggplot(alr6)+
  geom_line(aes(x = datetime, y = -depth))+
  ggtitle("ALR6")


library(rnaturalearth)
library(rnaturalearthdata)

iceland <- ne_countries(scale = "medium", country = "iceland", returnclass = "sf")

shelf_boundary <- data.frame(
  lon = c(-25, -13, -13, -25, -25),  # Example coordinates
  lat = c(63, 63, 64.5, 64.5, 63)
)

ggplot() +
  geom_sf(data = iceland, fill = "lightgray") +
  geom_polygon(data = shelf_boundary, aes(x = lon, y = lat), 
               fill = "lightblue", alpha = 0.3, color = "blue") +
  geom_point(data = dat_to_plot, aes(x = lon, y = lat, color = auv)) +
  coord_sf(xlim = c(-26, -12), ylim = c(60, 65))

ggsave("ALR_traj_shelf.jpg", width = 30, height = 20, units = "cm", dpi = 300)
datetime_shelf <- dat |>
  filter(lat >= 63) |> 
  group_by(auv) |> 
  summarise(datetime = max(datetime)) |> 
  ungroup()

# Then, create the plot separately
ggplot(dat) +
  geom_line(aes(x = datetime, y = -depth)) +
  geom_vline(data = datetime_shelf, aes(xintercept = datetime), 
             color = "red", linetype = "dashed") +
  facet_wrap(~auv)

ggsave("ALR_transect_shelf.jpg", width = 30, height = 20, units = "cm", dpi = 300)
