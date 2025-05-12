library(tidyverse)
library(janitor)
library(castr)
library(oce)
library(sf)
library(gganimate)
bathymetry <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_J_1000.shp") |> 
  st_cast("MULTILINESTRING")

bathymetry_2000 <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_I_2000.shp") |> 
  st_cast("MULTILINESTRING")

db <- read_csv("data/argo_pp_estimations.csv") |> 
  clean_names()

db <- db |> mutate(prof_index = paste0(platform_number, juld))

loc <- db |> select(platform_number, juld, longitude, latitude) |>
  mutate(platform_number = as.character(platform_number),
         date = as.Date(juld)) |>
  unique()

ggplot(loc)+
  geom_point(aes(x = longitude, y = latitude, colour = platform_number))+
  coord_quickmap()+
  scale_color_brewer(palette = 'Set1')
# 
stations <- data.frame('stations' = c(60.00016 ,-24.0004, 60.01929, -23.7195, 59.43984,-22.6396, 59.24316,-22.0786,60.23638,-21.2035, 56.21601, -26.8302, 57.12447,-26.1376, 57.2301,-28.3399,58.99212,-24.6729,60.12108, -18.9494, 61.17, -21.89, 62.15, -20.07),
                       'coord' = rep(c('lat', 'lon'), 12),
                       'names' = c(paste0('SS_', sort(rep(c(1:5), 2))), paste0('R_', sort(rep(c(1:5), 2))), 'F1', 'F1', 'F2', 'F2'))
stations <- pivot_wider(stations, names_from = 'coord', values_from = 'stations')

stations <- stations |>
  mutate(names = case_when(names == 'SS_1' ~ 'CIB',
                           names != 'SS_1' ~ names))

last_pos <- loc |>
  group_by(platform_number) |>
  filter(juld == max(juld)) |>
  mutate(date = as.Date(juld)) |> 
  ungroup()

ggplot(loc) +
  geom_point(data = last_pos, aes(x = longitude, y = latitude, color = platform_number), size = 5, alpha = 0.8) +
  geom_path(aes(x = longitude, y = latitude, group = platform_number, color = platform_number), linewidth = 2) +
  geom_point(data = stations, aes(x = lon, y = lat), shape = 4, size = 5)+
  geom_text(data = stations, aes(x = lon, y = lat + 0.3, label = names), size = 6)+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +
  scale_color_brewer(palette = 'Set1', name = 'WMO')+
  theme_bw(base_size = 16) +
  xlab('Longitude')+
  ylab('Latitude')+
  coord_sf(xlim = c(-27, -14), ylim = c(57, 64))

ggsave("output/plots/pp_floats/floats_map.png", dpi = 300, width = 20, height = 15, units = "cm")


library(ggplot2)
library(scales)

ggplot(loc) +
  geom_point(data = last_pos, aes(x = longitude, y = latitude), size = 5, alpha = 0.8) +
  geom_path(aes(x = longitude, y = latitude, group = platform_number, color = date), linewidth = 2) +
  geom_point(data = stations, aes(x = lon, y = lat), shape = 4, size = 5) +
  geom_text(data = stations, aes(x = lon, y = lat + 0.3, label = names), size = 6) +
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed") +
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed") +
  borders("world", colour = "black", fill = "gray80") +
  scale_color_viridis_c(name = 'Date' ) +
  theme_bw(base_size = 16) +
  xlab('Longitude') +
  ylab('Latitude') +
  coord_sf(xlim = c(-27, -14), ylim = c(57, 64))
ggsave("output/plots/pp_floats/floats_map_date.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(loc) +
  geom_point(data = last_pos, aes(x = longitude, y = latitude), size = 5, alpha = 0.8) +
  geom_path(aes(x = longitude, y = latitude, group = platform_number, color = date), linewidth = 2) +
  geom_point(data = stations, aes(x = lon, y = lat), shape = 4, size = 5) +
  geom_text(data = stations, aes(x = lon, y = lat + 0.3, label = names), size = 6) +
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed") +
  geom_rect(aes(xmin = -26, xmax = -15.5, ymin = 59, ymax = 62.5), alpha = 0, color = "Red")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed") +
  borders("world", colour = "black", fill = "gray80") +
  scale_color_viridis_c(name = 'Date' ) +
  theme_bw(base_size = 16) +
  xlab('Longitude') +
  ylab('Latitude') +
  coord_sf(xlim = c(-27, -14), ylim = c(57, 64))

ggsave("output/plots/pp_floats/floats_map_date_rectangle.png", dpi = 300, width = 20, height = 15, units = "cm")

# gif <- p + transition_reveal(juld)+
#   labs(subtitle = "Date: {frame_along}")
# 
# animate(gif, height = 1200, width = 1600, fps = 30, duration = 15,
#         end_pause = 60, res = 100)
# anim_save("floats_icb_2024__60s.gif")

#p

db <- db |> 
  group_by(prof_index) |> 
  mutate(ct = gsw_CT_from_t(psal, temp, pres_rounded),
         sigma0 = gsw_sigma0(psal, ct),
         chla_smoothed = smooth(chla_adjusted, k = 7, n = 1),
         bbp470_smoothed = smooth(bbp470, k = 7, n = 1),
         carbon_smoothed = smooth(carbon, k = 7),
         pp_smoothed = smooth(pp, k = 7))

stats <- db %>% group_by(platform_number, juld) %>%
  summarise(
    thermocline = clined(ct, pres_rounded, n.smooth=2, k=2),
    pycnocline = clined(psal, pres_rounded),
    DCM = maxd(chla_adjusted, pres_rounded, n.smooth=2, k=3),
    MLD = mld(sigma0, pres_rounded, ref.depths=0:5, criteria = 0.03, default.depth=200),
    chla_mld_stock = integrate(chla_adjusted, pres_rounded, from=0, to=200),
    pp_stock = integrate(pp_smoothed, pres_rounded, from=0, to=200),
    carbon_stock = integrate(carbon_smoothed, pres_rounded, from = 0, to = 200)
  )

db <- left_join(db, stats)

db_icb <- db |> filter(longitude < -15.5 & longitude > -26) |> 
  filter(latitude < 62.5 & latitude > 59) |>
  filter(pres_rounded <= 200) |> 
  ungroup() |> 
  mutate(date = as.Date(juld))

ggplot(db_icb)+
  geom_path(aes(x = chla_smoothed, y = -pres_rounded, group = prof_index))+
  facet_wrap(.~ platform_number)+
  xlab("Chla Adjusted")+
  ylab("Pres")+
  theme_bw()

#ggsave("output/plots/pp_floats/chla_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(filter(db_icb, carbon_smoothed < 150))+
  geom_path(aes(x = carbon_smoothed, y = -pres_rounded, group = prof_index))+
  facet_wrap(.~ platform_number)+
  xlab("Carbon phyto")+
  ylab("Pres")+
  theme_bw()

#ggsave("output/plots/pp_floats/carbon_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(filter(db_icb, carbon_smoothed < 150))+
  geom_path(aes(x = carbon_smoothed, y = -pres_rounded, group = prof_index))+
  xlab("Carbon phyto")+
  ylab("Pres")+
  theme_bw()

ggplot(db)+
  geom_point(aes(x = juld, y = -zeu, color = 'zeu'))+
  geom_point(aes(x = juld, y = -MLD, color = 'mld'))+
  scale_color_brewer(palette = 'Set1')

ggplot(filter(db_icb, pres_rounded < 50))+
  geom_path(aes(x = juld, y = pp, color = pres_rounded, group = pres_rounded))+
  facet_wrap(.~ platform_number)+
  ylab("NPP (mg C m-3 d-1)")+
  xlab("Date")+
  scale_color_viridis_c(nam = "Depth")+
  theme_bw()
#ggsave("output/plots/pp_floats/NPP_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(filter(db_icb, pres_rounded < 50))+
  geom_path(aes(x = juld, y = mu, color = pres_rounded, group = pres_rounded))+
  facet_wrap(.~ platform_number)+
  ylab("mu (d-1)")+
  xlab("Date")+
  scale_color_viridis_c(nam = "Depth")+
  theme_bw()
#ggsave("output/plots/pp_floats/mu_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")


db_icb = db_icb |> mutate(unique_ind = paste(platform_number, pres_rounded))
ggplot(filter(db_icb, pres_rounded < 50))+
  geom_path(aes(x = juld, y = pp, color = pres_rounded, group = unique_ind))+
  ylab("NPP (mg C m-3 d-1)")+
  xlab("Date")+
  scale_color_viridis_c(nam = "Depth")+
  theme_bw()
#ggsave("output/plots/pp_floats/NPP_from_all_floats.png", dpi = 300, width = 20, height = 15, units = "cm")



ggplot(db_icb)+
  geom_path(aes(x = pp, y = -pres_rounded, group = prof_index))+
  facet_wrap(.~ platform_number)+
  xlab("NPP (mg C m-3 d-1)")+
  ylab("Date")+
  theme_bw()
#ggsave("output/plots/pp_floats/NPP_profiles_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(db_icb)+
  geom_path(aes(x = mu, y = -pres_rounded, group = prof_index))+
  facet_wrap(.~ platform_number)+
  xlab("mu (d-1)")+
  ylab("Date")+
  theme_bw()
#ggsave("output/plots/pp_floats/mu_profiles_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")


tet = db |> filter(platform_number == 4903532)

ggplot(tet)+
  geom_path(aes(x = chla_adjusted, y = -pres_rounded, group = juld))
yr <- filter(db, platform_number == 1902637) |> 
  mutate(yield = pp/mu)

ggplot(yr)+
  geom_tile(aes(x = juld, y = - pres_rounded, fill = yield))+
  scale_color_viridis_c()

ggplot(db)+
  geom_point(aes(x = mu, y = pp, color = pres_rounded))+
  scale_color_viridis_c()

yr2 <- filter(db, platform_number == 7902223) |> 
  mutate(yield = pp/mu)

ggplot(filter(yr2, juld < as.Date('2024-09-15')))+
  geom_tile(aes(x = juld, y = - pres_rounded, fill = pp))+
  scale_fill_viridis_c()

db_map <- db |> group_by(platform_number, juld, longitude, latitude) |> 
  select(pp, mu) |> 
  summarise_all(sum, na.rm = TRUE) |> 
  ungroup()

ggplot(filter(db_map, latitude > 58))+
  geom_point(aes(x = longitude, y = latitude, color = pp))+
  scale_color_viridis_c()+
  coord_quickmap()

seasons = function(x){
  if(x %in% 3:5) return("Spring")
  if(x %in% 6:8) return("Summer")
  if(x %in% 9:11) return("Fall")
  if(x %in% c(12,1,2)) return("Winter")
  
}

db$season = sapply(month(db$juld), seasons)

season_avg <- db |> select(season, pres_rounded, pp, mu, carbon, chla_adjusted) |> 
  group_by(season, pres_rounded) |> 
  summarise_all(list(mean, sd), na.rm = TRUE) |> 
  ungroup()

ggplot(season_avg)+
  geom_path(aes(x = pp_fn1, y = - pres_rounded))+
  geom_line(aes(x = pp_fn1 - pp_fn2, y = - pres_rounded), linetype = "dotted")+
  geom_line(aes(x = pp_fn1 + pp_fn2, y = - pres_rounded), linetype = "dotted")+
  facet_wrap(.~season)+
  ylim(-100, 0)+
  theme_bw()+
  xlab("NPP (mg C m-3 d-1)")+
  ylab("Depth (m)")
#ggsave("output/plots/pp_floats/NPP_ season_profiles_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(season_avg)+
  geom_path(aes(x = mu_fn1, y = - pres_rounded))+
  geom_line(aes(x = mu_fn1 - mu_fn2, y = - pres_rounded), linetype = "dotted")+
  geom_line(aes(x = mu_fn1 + mu_fn2, y = - pres_rounded), linetype = "dotted")+
  facet_wrap(.~season)+
  ylim(-100, 0)+
  theme_bw()+
  xlab("mu (d-1)")+
  ylab("Depth (m)")

#ggsave("output/plots/pp_floats/mu_season_profiles_from_floats.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(db)+
  geom_point(aes(x = juld, y = satellite_par))

winter <- filter(db, season == "Winter")

winter <- left_join(winter, stats)
ggplot(filter(winter, MLD != 300))+
  geom_point(aes(x = temp, y = - pres_rounded))+
  geom_hline(aes(yintercept = -MLD))


# net growth or loss ------------------------------------------------------

ggplot(db)+
  geom_point(aes(x = juld, y = -pres_rounded, color = pp))+
  facet_wrap(.~ platform_number)+
  scale_color_viridis_c()

ggplot(filter(db, pp_stock < 2500))+
  geom_line(aes(x = juld, y = pp_stock, color = as.factor(platform_number)))

#Let's start first with 1902637

flt <- filter(db, platform_number == 4903659)

ggplot(flt)+
  geom_point(aes(x = juld, y = chla_mld_stock))
ggplot(flt)+
  geom_point(aes(x = juld, y = -MLD, color = 'mld'))+
  geom_point(aes(x = juld, y = -zeu, color = 'zeu'))+
  scale_color_brewer(palette = 'Set1')+
  ylab('Depth (m)')

ggplot(flt)+
  geom_point(aes(x = juld, y = chla_mld_stock, color = satellite_par))+
  scale_color_viridis_c()

ggplot(flt)+
  geom_path(aes(x = chla_adjusted, y = -pres_rounded, group = juld))

surf <- filter(db, pres_rounded < 1 & downwelling_par > 1) |>
  group_by(juld) |> 
  select(downwelling_par, satellite_par) |> 
  summarise_all(mean) |> 
  na.omit()

ggplot(surf)+
  geom_point(aes(y = satellite_par, x = downwelling_par))+
  xlab("Argo PAR")+
  ylab("MODIS PAR")



# transect ----------------------------------------------------------------


db_icb <- db |> filter(longitude < -17 & longitude > -26) |> 
  filter(latitude < 62.5 & latitude > 59) |>
  filter(pres_rounded <= 200) |> 
  ungroup() |> 
  mutate(date = as.Date(juld))

depth_grid <- seq(0, 200, by = 1)  # Regular grid for depth (every 1 meter)
time_grid <- seq(min(db_icb$date), max(db_icb$date), by = "1 week")  # Regular daily time grid

# Create a data frame for the regular grid
grid <- expand.grid(
  juld = time_grid,
  pres_rounded = depth_grid
)

db_icb <- db_icb %>%
  mutate(date_numeric = as.numeric(date))

grid <- grid %>%
  mutate(date_numeric = as.numeric(juld))

# # Use the `interp` function from the akima package for 2D interpolation
# interp_result <- with(db_icb, 
#                       akima::interp(
#                         x = date_numeric,
#                         y = pres_rounded,
#                         z = pp,
#                         xo = grid$date_numeric,
#                         yo = grid$pres_rounded,
#                         duplicate = "mean"
#                       )
# )
# 
# # Convert the result back into a tidy data frame
# interpolated_data <- as_tibble(expand.grid(
#   date = interp_result$x,
#   pres_rounded = interp_result$y
# )) %>%
#   mutate(pp = c(interp_result$z))  # Flatten the interpolated matrix

db_icb <- db_icb |> mutate(zeu_mld = MLD-zeu)

ggplot(db_icb)+
  geom_point(aes(x = date, y = chla_mld_stock, color = as.factor(platform_number)))+
  scale_color_brewer(palette = "Set1") +
  scale_x_date(
    date_breaks = "1 month",    # Set breaks at the start of every month
    date_labels = "%b"      # Format labels as "Jan 2024", "Feb 2024", etc.
  ) +
  theme_minimal() +             # Use a minimal theme (optional)
  labs(
    x = "Date",
    y = "Chlorophyll-a Stock (200m)",
    color = "Float WMO",
    title = "Chlorophyll-a Stock vs Date"
  )

#ggsave("output/plots/pp_floats/chla_stock_ts.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(filter(db_icb, pp_stock < 2500))+
  geom_point(aes(x = date, y = pp_stock, color = as.factor(platform_number)))+
  scale_color_brewer(palette = "Set1") +
  scale_x_date(
    date_breaks = "1 month",    # Set breaks at the start of every month
    date_labels = "%b"      # Format labels as "Jan 2024", "Feb 2024", etc.
  ) +
  theme_minimal() +             # Use a minimal theme (optional)
  labs(
    x = "Date",
    y = "Integrated NPP (200m) (mg C m-2 d-1)",
    color = "Float WMO"
  )

#ggsave("output/plots/pp_floats/npp_integrated_ts.png", dpi = 300, width = 20, height = 15, units = "cm")

t <- filter(db_icb, pp_stock < 2500)
ggplot(db_icb)+
  geom_point(aes(x = date, y = chla_mld_stock/carbon_stock, color = as.factor(platform_number)))+
  scale_color_brewer(palette = "Set1") +
  scale_x_date(
    date_breaks = "1 month",    # Set breaks at the start of every month
    date_labels = "%b"      # Format labels as "Jan 2024", "Feb 2024", etc.
  ) +
  theme_minimal() +             # Use a minimal theme (optional)
  labs(
    x = "Date",
    y = "Integrated NPP (200m) (mg C m-2 d-1)",
    color = "Float WMO"
  )
