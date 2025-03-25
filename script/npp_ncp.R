library(tidyverse)
library(arrow)
library(zoo)
library(patchwork)
library(sf)
library(castr)
library(gsw)

bathymetry <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_J_1000.shp") |> 
  st_cast("MULTILINESTRING")

bathymetry_2000 <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_I_2000.shp") |> 
  st_cast("MULTILINESTRING")

argo <- read_parquet("data/argo_pq/biocarbon_floats_table.parquet") |> mutate(depth = round(PRES))



new_ref <- select(argo, PLATFORM_NUMBER, JULD) |> unique()|> dplyr::filter(PLATFORM_NUMBER != "5904183 ")
depth <- tibble("depth" = c(0:200))

new_ref <- new_ref |> crossing(depth) |> left_join(argo)

new_ref <- new_ref |> mutate(JULD_date = lubridate::date(JULD),
                             prof_id = paste0(PLATFORM_NUMBER, JULD_date)) |>
  filter(JULD_date > lubridate::date("2023-12-31")) |> 
  filter(!prof_id %in% c("4903659 2024-06-29", "4903659 2024-06-23")) |> group_by(prof_id) %>%  # Interpolate within each profile
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.approx(.x, depth, na.rm = FALSE))) %>%
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  mutate(across(c(NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.locf(.x, na.rm = FALSE))) %>%
  ungroup()


write_parquet(new_ref, "data/argo_pq/float_cleaned.parquet")

npp_df <- read_csv("data/argo_pp_estimations_float4903532.csv") |> 
  mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth),
         sigma0 = gsw_sigma0(PSAL, TEMP))


prof_dat <- npp_df |> 
  group_by(PLATFORM_NUMBER, JULD) |> 
  mutate(nitrate_smoothed = smooth(NITRATE_ADJUSTED, k = 5, n = 2)) |> 
  summarise(MLD = mld(sigma0, depth, ref.depths=0:5, criteria = 0.03, default.depth=200),
            nitracline = clined(nitrate_smoothed, depth))

npp_df <- left_join(npp_df, prof_dat)

ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = CHLA_ADJUSTED))+
  geom_line(aes(x = JULD, y = - MLD), colour = "white")+
  scale_fill_viridis_c()

ggsave("output/plots/pp_floats/chla_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = pp))+
  geom_line(aes(x = JULD, y = - MLD), colour = "white")+
  scale_fill_viridis_c(name = "NPP (mg C m-3)")

ggsave("output/plots/pp_floats/npp_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = NITRATE_ADJUSTED))+
  geom_line(aes(x = JULD, y = - MLD), colour = "white")+
  scale_fill_viridis_c(name = "Nitrate (micromol kg-1)")

ggsave("output/plots/pp_floats/nitrate_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")

location <- select(argo2, TIME, LONGITUDE, LATITUDE) |> unique()

integrated_df <- npp_df |> select(PLATFORM_NUMBER, JULD, pp, NITRATE_ADJUSTED) |>
  group_by(PLATFORM_NUMBER, JULD) |> 
  summarise_all(sum) |> 
  ungroup() |> 
  left_join(prof_dat) |> 
  left_join(location) |> 
  mutate(nitrate_smoothed = smooth.spline(JULD, NITRATE_ADJUSTED, spar = 0.6)$y)
  

integrated_df$date <- lubridate::date(integrated_df$JULD)

ggplot(integrated_df)+
  geom_point(aes(x = date, y = NITRATE_ADJUSTED))+
  geom_path(aes(x = date, y = nitrate_smoothed))

ggsave("output/plots/pp_floats/nitrate_ts.png", dpi = 300, width = 20, height = 15, units = "cm")



ggplot(integrated_df)+
  geom_point(aes(x = date, y = -nitracline))


ggplot(integrated_df)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE, colour = pp/12))+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +  
  xlab('Longitude')+
  ylab('Latitude')+
  scale_color_viridis_c(name = "Intregrated NPP mmol C .m-2 .d-1")+
  coord_sf(xlim = c(-32, -12), ylim = c(57, 64))+
  theme_minimal()

ggsave("output/plots/pp_floats/npp4903532_map.png", dpi = 300, width = 20, height = 15, units = "cm")

integrated_df$nitrate_dd <- c(NA, diff(integrated_df$nitrate_smoothed))
integrated_df$diff_date <- c(NA, diff(integrated_df$date))

integrated_df <- integrated_df |> mutate(dd_smoothed = rollmean(nitrate_dd, k = 3, fill = "extend"))

ggplot(integrated_df)+
  geom_point(aes(x = date, y = dd_smoothed))

integrated_df <- integrated_df |> mutate(ncp = -(dd_smoothed*6.6),
                                         pp = pp/10,
                                         budget = pp - ncp)


ggplot(integrated_df)+
  geom_path(aes(x = date, y = pp))+
  ylab("Integrated NPP")

ggsave("output/plots/pp_floats/npp4903532_line.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(integrated_df)+
  geom_path(aes(x = date, y = ncp))+
  ylab("Integrated NCP")

ggsave("output/plots/pp_floats/ncp4903532_line.png", dpi = 300, width = 20, height = 15, units = "cm")

integrated_df <- integrated_df |> mutate(budget = pp - ncp,
                                         period = case_when(date > lubridate::date("2024-01-01") & date < lubridate::date("2024-04-20") ~ "Pre bloom",
                                                            date >= lubridate::date("2024-05-20") & date <= lubridate::date("2024-07-15") ~ "Main bloom",
                                                            date >= lubridate::date("2024-07-15") & date < lubridate::date("2024-10-01") ~ "Declining productivity",
                                                            date < lubridate::date("2025-03-01") ~ "Post bloom"))

stats_table <- integrated_df |> select(pp, ncp, period) |> group_by(period) |> summarise_all(mean, na.rm = TRUE) |> 
  pivot_longer(cols = c(pp, ncp), names_to = "productivity", values_to = "val")

uncertainty <- integrated_df |> select(pp, ncp, period) |> group_by(period) |> summarise_all(sd, na.rm = TRUE) |> 
  pivot_longer(cols = c(pp, ncp), names_to = "productivity", values_to = "sd")

stats_table <- stats_table |> left_join(uncertainty)

stats_table$period <- factor(stats_table$period, levels = c("Pre bloom", "Main bloom", "Declining productivity", "Post bloom"))

p1 <- ggplot(stats_table)+
  geom_col(aes(x = period, y = val, fill = productivity), position="dodge")+
  geom_errorbar(aes(x = period, ymin = val - sd, ymax = val+sd, group = productivity), position = "dodge")+
  theme_bw()+
  scale_fill_brewer(palette = "Set1", name = "")+
  ylab("Productivity (mmol C m-2 d-1)")+
  xlab("")

p2 <- ggplot(integrated_df)+
  geom_path(aes(x = date, y = pp, colour = "NPP"))+
  geom_path(aes(x = date, y = ncp, colour = "NCP"))+
  geom_vline(aes(xintercept = lubridate::date("2024-04-20")), linetype = "dashed")+
  geom_vline(aes(xintercept = lubridate::date("2024-07-15")), linetype = "dashed")+
  geom_vline(aes(xintercept = lubridate::date("2024-11-01")), linetype = "dashed")+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  ylab("Integrated productivity (mmol C m-2 d-1)")

  
p1/p2  +  plot_layout(height = c(2, 1))

ggsave("output/plots/pp_floats/npp_vs_ncp.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(integrated_df)+
  geom_path(aes(x = pp, y = ncp, color = date), linewidth = 2)+
  scale_color_viridis_c()

ggsave("output/plots/pp_floats/npp_vs_ncp_ts.png", dpi = 300, width = 20, height = 15, units = "cm")

