library(tidyverse)
library(arrow)
library(zoo)
library(patchwork)
library(sf)
library(castr)
library(gsw)
library(broom)
bathymetry <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_J_1000.shp") |> 
  st_cast("MULTILINESTRING")

bathymetry_2000 <- st_read("C:/Users/flapet/OneDrive - NOC/Documents/NRT_viz/biocarbon_nrt_data_viz/Data/ne_10m_bathymetry_all/ne_10m_bathymetry_I_2000.shp") |> 
  st_cast("MULTILINESTRING")


# formatting dataset (to run only once) -----------------------------------

argo <- read_parquet("data/argo_pq/biocarbon_floats_table.parquet") |> mutate(depth = round(PRES))



new_ref <- select(argo, PLATFORM_NUMBER, JULD) |> unique()|> dplyr::filter(PLATFORM_NUMBER != "5904183 ")
depth <- tibble("depth" = c(0:200))

new_ref <- new_ref |> crossing(depth) |> left_join(argo)

new_ref <- new_ref |> mutate(JULD_date = lubridate::date(JULD),
                             prof_id = paste0(PLATFORM_NUMBER, JULD_date)) |>
  filter(JULD_date > lubridate::date("2023-12-31")) |>
  filter(!prof_id %in% c("4903659 2024-06-29", "4903659 2024-06-23", "6904185 2024-01-12", "1902304 2024-01-19")) |> group_by(prof_id) %>%  # Interpolate within each profile
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.approx(.x, depth, na.rm = FALSE))) %>%
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  mutate(across(c(NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.locf(.x, na.rm = FALSE))) %>%
  ungroup()


#write_parquet(new_ref, "data/argo_pq/float_cleaned.parquet")


# KD490 part --------------------------------------------------------------

test <- read_parquet("data/argo_pq/float_cleaned.parquet")

kd490 <- filter(test, !is.na(DOWN_IRRADIANCE490))

chl_kd <- kd490 %>% group_by(prof_id) %>%
  filter(PRES < 40 & PRES > 10) |>
  mutate(DOWN_IRRADIANCE490 = case_when(DOWN_IRRADIANCE490 < 0 ~ 4e-09,
                                        TRUE ~ DOWN_IRRADIANCE490),
         lm_dwn490 = log(DOWN_IRRADIANCE490)) |> 
  do(model = lm(lm_dwn490 ~ PRES, data = .)) %>% ungroup() %>% 
  transmute(prof_id, coef = map(model, tidy)) %>% 
  unnest(coef) %>% 
  filter(term == "PRES") 

chl_kd <- chl_kd |> 
  mutate(chl_estimate = ((-estimate - 0.0166) / (0.0773))^(1/0.6715))

chl_mean <- kd490 %>% group_by(prof_id) %>%
  filter(PRES < 40 & PRES > 10) |> select(prof_id, CHLA_ADJUSTED) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  rename("mean_chla" = CHLA_ADJUSTED)

kd490_df <- kd490 |> left_join(select(chl_kd, prof_id, chl_estimate)) |> left_join(chl_mean)

plot_kd <- kd490_df |> 
  select(JULD, chl_estimate, mean_chla, LONGITUDE, LATITUDE, PLATFORM_NUMBER) |> 
  unique() |> 
  mutate(f_chl = (mean_chla *2)/ chl_estimate)

ggplot(plot_kd)+
  geom_point(aes(x = JULD, y = chl_estimate, color = "Kd490"))+
  geom_point(aes(x = JULD, y = mean_chla *2, color = "Fluorescence"))+
  ylim(0, 2.5)+
  scale_color_brewer(palette = "Set1", name = "Chl estimate")+
  ylab("Chl")+
  xlab("Date")+
  theme_bw()

ggsave("output/plots/pp_floats/chla_estimates_comparison.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(plot_kd)+
  geom_point(aes(x = JULD, y = f_chl))+
  geom_smooth(aes(x = JULD, y = f_chl), method = "lm", se = FALSE)+
  ylim(0, 3)+
  ylab("Fluo/Chl")+
  xlab("Date")+
  theme_bw()
ggsave("output/plots/pp_floats/f_chl_timeseries.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(plot_kd)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE, color = PLATFORM_NUMBER))+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +  
  xlab('Longitude')+
  ylab('Latitude')+
  coord_sf(xlim = c(-32, -12), ylim = c(57, 64))+
  theme_minimal()+
  ggtitle("KD490 dataset")

ggsave("output/plots/pp_floats/kd490_dataset.png", dpi = 300, width = 20, height = 15, units = "cm")

kd_profile <- filter(kd490, prof_id == "4903659 2024-07-24")

ggplot(kd_profile)+
  geom_point(aes(x = DOWN_IRRADIANCE490, y = -PRES))+
  geom_hline((aes(yintercept = -40)))+
  geom_hline((aes(yintercept = -10)))+
  ylim(-100, 0)
ggsave("output/plots/pp_floats/kd490_profile.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(kd_profile)+
  geom_point(aes(x = CHLA_ADJUSTED, y = -PRES))+
  geom_hline((aes(yintercept = -40)))+
  geom_hline((aes(yintercept = -10)))+
  ylim(-100, 0)

ggplot(filter(kd_profile, PRES < 40 & PRES > 10))+
  geom_point(aes(x = PRES, y = log(DOWN_IRRADIANCE490)))
ggsave("output/plots/pp_floats/kd490_regression.png", dpi = 300, width = 20, height = 15, units = "cm")


# NCP computation ---------------------------------------------------------

npp_df <- read_csv("data/argo_pp_estimations_floats.csv") |> 
  mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth),
         sigma0 = gsw_sigma0(PSAL, TEMP)) |> 
  filter(PLATFORM_NUMBER == "1902304" | PLATFORM_NUMBER == "4903532")


prof_dat <- npp_df |> 
  mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth),
         sigma0 = gsw_sigma0(PSAL, TEMP)) |> 
  filter(PLATFORM_NUMBER == "1902304" | PLATFORM_NUMBER == "4903532") |> 
  group_by(PLATFORM_NUMBER, JULD) |> 
  mutate(nitrate_smoothed = smooth(NITRATE_ADJUSTED, k = 5, n = 2)) |> 
  summarise(MLD = mld(sigma0, depth, ref.depths=0:5, criteria = 0.03, default.depth=200),
            nitracline = clined(nitrate_smoothed, depth),
            chl_cline = clined(CHLA_ADJUSTED, depth)) |>
  ungroup() |> 
  group_by(PLATFORM_NUMBER) |> 
  mutate(delta_MLD = MLD - lag(MLD)) |> 
  replace_na(list(delta_MLD = 0))

npp_df <- left_join(npp_df, prof_dat)

ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = NITRATE_ADJUSTED))+
  geom_line(aes(x = JULD, y = - MLD))+
  scale_fill_viridis_c()+
  ylim(-300, 0)+
  facet_wrap(.~PLATFORM_NUMBER)

#ggsave("output/plots/pp_floats/chla_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = pp))+
  geom_line(aes(x = JULD, y = - MLD), colour = "white")+
  scale_fill_viridis_c(name = "NPP (mg C m-3)")

#ggsave("output/plots/pp_floats/npp_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")


ggplot(npp_df)+
  geom_tile(aes(y = -depth, x = JULD, fill = NITRATE_ADJUSTED))+
  scale_fill_viridis_c(name = "Nitrate (micromol kg-1)")

#ggsave("output/plots/pp_floats/nitrate_transect_4903532.png", dpi = 300, width = 20, height = 15, units = "cm")

location <- select(npp_df, JULD, LONGITUDE, LATITUDE) |> unique()

#Making a 200:1000 depth argo table

argo <- read_parquet("data/argo_pq/biocarbon_floats_table.parquet") |> mutate(depth = round(PRES))

new_ref <- select(argo, PLATFORM_NUMBER, JULD) |> unique()|> dplyr::filter(PLATFORM_NUMBER %in% c("1902304 ", "4903532 "))
depth <- tibble("depth" = c(0:1000))

new_ref <- new_ref |> crossing(depth) |> left_join(argo)

new_ref <- new_ref |> mutate(JULD_date = lubridate::date(JULD),
                             prof_id = paste0(PLATFORM_NUMBER, JULD_date)) |>
  filter(JULD_date > lubridate::date("2023-12-31")) |>
  filter(!prof_id %in% c("1902304 2024-01-19")) |> group_by(prof_id) %>%  # Interpolate within each profile
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.approx(.x, depth, na.rm = FALSE))) %>%
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  mutate(across(c(NITRATE_ADJUSTED, DOXY_ADJUSTED),
                ~ na.locf(.x, na.rm = FALSE))) %>%
  ungroup()

new_ref$PLATFORM_NUMBER <- as.numeric(new_ref$PLATFORM_NUMBER)


integrated_df <- npp_df |>
  mutate(layer = case_when(depth <= MLD ~ "MLD",
                           depth > MLD & depth <= max(prof_dat$MLD) ~ "SUB MLD",
                           depth > max(prof_dat$MLD)~"DEEP")) |> 
  select(PLATFORM_NUMBER, JULD, NITRATE_ADJUSTED, layer, MLD, depth) |>
  group_by(PLATFORM_NUMBER, JULD, layer) |> 
  summarise("integrated_nitrate" = sum(NITRATE_ADJUSTED), "mean_nitrate" = mean(NITRATE_ADJUSTED), layer_top = min(depth), layer_bottom = max(depth)) |> 
  ungroup() |> 
  left_join(prof_dat) |> 
  group_by(PLATFORM_NUMBER, layer) |> 
  mutate(mean_nitrate_smoothed = smooth.spline(JULD, mean_nitrate, spar = 0.6)$y,
         integrated_nitrate_smoothed = smooth.spline(JULD, integrated_nitrate, spar = 0.6)$y)
  

integrated_df$date <- lubridate::date(integrated_df$JULD)

ggplot(integrated_df)+
  geom_point(aes(x = date, y = mean_nitrate, color = as.factor(PLATFORM_NUMBER), shape = layer))+
  geom_path(aes(x = date, y = mean_nitrate_smoothed, color = as.factor(PLATFORM_NUMBER), linetype = layer))+
  scale_color_brewer(palette = "Set1", name = "wmo")


ggplot(integrated_df)+
  geom_point(aes(x = date, y = integrated_nitrate, color = as.factor(PLATFORM_NUMBER), shape = layer))+
  geom_path(aes(x = date, y = integrated_nitrate_smoothed, color = as.factor(PLATFORM_NUMBER), linetype = layer))+
  scale_color_brewer(palette = "Set1", name = "wmo")

#ggsave("output/plots/pp_floats/nitrate_ts.png", dpi = 300, width = 20, height = 15, units = "cm")

ncp_ready <- left_join(npp_df, integrated_df)

ggplot(integrated_df)+
  geom_point(aes(x = date, y = -nitracline))

ggplot(integrated_df)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE, colour = as.factor(PLATFORM_NUMBER)))+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +  
  xlab('Longitude')+
  ylab('Latitude')+
  scale_color_brewer(name = "WMO", palette = "Set1")+
  coord_sf(xlim = c(-32, -12), ylim = c(57, 64))+
  theme_minimal()

ggsave("output/plots/pp_floats/comparison_map.png", dpi = 300, width = 20, height = 15, units = "cm")

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

#ggsave("output/plots/pp_floats/npp4903532_map.png", dpi = 300, width = 20, height = 15, units = "cm")

integrated_df$nitrate_dd <- c(NA, diff(integrated_df$nitrate_smoothed))
integrated_df$diff_date <- c(NA, diff(integrated_df$date))

integrated_df <- integrated_df |> mutate(dd_smoothed = rollmean(nitrate_dd, k = 3, fill = "extend"))

ggplot(integrated_df)+
  geom_point(aes(x = date, y = dd_smoothed, color = as.factor(PLATFORM_NUMBER)))

integrated_df <- integrated_df |>
  group_by(PLATFORM_NUMBER) |> 
  mutate(ncp = -(dd_smoothed*6.6),
                                         pp = pp/diff_date,
                                         budget = pp - ncp) |> 
  mutate(pp = na.approx(pp, na.rm = FALSE))


ggplot(integrated_df)+
  geom_line(aes(x = date, y = pp, color = as.factor(PLATFORM_NUMBER)))+
  ylab("Integrated NPP")+
  scale_color_brewer(palette = "Set1", name = "WMO")

#ggsave("output/plots/pp_floats/nppcomparison_line.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(integrated_df)+
  geom_line(aes(x = date, y = ncp, color = as.factor(PLATFORM_NUMBER)))+
  ylab("Integrated NPP")+
  scale_color_brewer(palette = "Set1", name = "WMO")

#ggsave("output/plots/pp_floats/ncp_only_comparison_line.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(integrated_df)+
  geom_path(aes(x = date, y = ncp, color = as.factor(PLATFORM_NUMBER)))+
  geom_path(aes(x = date, y = pp*12, color = as.factor(PLATFORM_NUMBER)), linetype = "dashed")+
  ylab("Integrated NCP")+
  scale_color_brewer(palette = "Set1", name = "WMO")
ggsave("output/plots/pp_floats/ncpcomparison_line.png", dpi = 300, width = 20, height = 15, units = "cm")

#ggsave("output/plots/pp_floats/ncp4903532_line.png", dpi = 300, width = 20, height = 15, units = "cm")

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

plot_int <- integrated_df |> na.omit()
ggplot(integrated_df)+
  geom_line(aes(x = date, y = pp * 12, colour = "NPP"))+
  geom_path(aes(x = date, y = ncp, colour = "NCP"))+
  theme_bw()+
  scale_color_brewer(palette = "Set1")+
  ylab("Integrated productivity (mmol C m-2 d-1)")+
  ylim(0, 350)
#ggsave("output/plots/pp_floats/npp_vs_ncp.png", dpi = 300, width = 20, height = 15, units = "cm")

#ggsave("output/plots/pp_floats/positive_npp_vs_ncp.png", dpi = 300, width = 20, height = 15, units = "cm")

ggplot(integrated_df)+
  geom_path(aes(x = pp, y = ncp, color = date), linewidth = 2)+
  scale_color_viridis_c()

#ggsave("output/plots/pp_floats/npp_vs_ncp_ts.png", dpi = 300, width = 20, height = 15, units = "cm")

