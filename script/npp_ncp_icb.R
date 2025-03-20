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

argo <- read_parquet("data/argo_pq/biocarbon_nitrate_floats_table.parquet") |> mutate(depth = round(PRES))
argo2 <- read_parquet("data/argo_pq/icb_floats_table.parquet") |> mutate(depth = round(PRES))

argo2 <- argo2 |> filter(PRES < 202 & LATITUDE > 55)

new_ref <- select(argo, PLATFORM_NUMBER, JULD) |> unique()
depth <- tibble("depth" = c(0:200))

new_ref <- new_ref |> crossing(depth) |>
  left_join(argo) |>
  mutate(JULD_date = lubridate::date(JULD),
         prof_id = paste0(PLATFORM_NUMBER, JULD_date))

#Check for profiles with too many NaNs

variables_to_check <- c("TEMP", "PSAL", "CHLA_ADJUSTED", "BBP700_ADJUSTED", "NITRATE_ADJUSTED")

missing_profiles <- new_ref |> 
  pivot_longer(cols = all_of(variables_to_check), names_to = "variable", values_to = "value") |> 
  group_by(prof_id, variable) |> 
  summarise(nona_count = sum(!is.na(value))) |> 
  filter(nona_count <= 4) |> 
  pull(prof_id) |> 
  unique()

single_date <- new_ref |> select(prof_id, JULD) |> unique()
no_data <- new_ref |> select(prof_id, JULD) |> unique() |> filter(prof_id %in% missing_profiles)

interp_data <- new_ref |>
  filter(!prof_id %in% missing_profiles) |> group_by(prof_id) %>%  # Interpolate within each profile
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.approx(.x, depth, na.rm = FALSE))) %>%
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  mutate(across(c(NITRATE_ADJUSTED, DOXY_ADJUSTED), 
                ~ na.locf(.x, na.rm = FALSE))) %>%
  ungroup()



ggplot()+
  geom_point(data = single_date, aes(y = 1, x = JULD, colour = "full data"))+
  geom_point(data = no_data, aes(y = 2, x = JULD, colour = "missing data"))+
  ylim(0,3)

write_parquet(interp_data, "data/argo_pq/float_icb_nitrate_cleaned.parquet")

location_fluo <- argo2 |>
  filter(CHLA_ADJUSTED_QC <= 3) |> select(TIME, LONGITUDE, LATITUDE) |> unique() |> 
  mutate("var" = "Fluo",
         prof_id = paste0(TIME, LONGITUDE, sep = "_"))

location_nitrate <- argo2 |>
  filter(NITRATE_ADJUSTED_QC <= 2) |> select(TIME, LONGITUDE, LATITUDE) |> 
  unique() |> 
  mutate("var" = "Nitrate",
         prof_id = paste0(TIME, LONGITUDE, sep = "_"))

location_sum <- bind_rows(location_fluo, location_nitrate)


ggplot(location_fluo)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE, colour = TIME))+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +  
  xlab('Longitude')+
  ylab('Latitude')+
  scale_color_viridis_c(name = "Date")+
  coord_sf(xlim = c(-45, -15), ylim = c(55, 66))+
  theme_minimal()+
  ggtitle("Fluorescence")+


ggplot(location_nitrate)+
  geom_point(aes(x = LONGITUDE, y = LATITUDE, colour = TIME))+
  geom_sf(data = bathymetry, color = "Black", linetype = "dashed")+
  geom_sf(data = bathymetry_2000, color = "Grey", linetype = "dashed")+
  borders("world", colour = "black", fill = "gray80") +  
  xlab('Longitude')+
  ylab('Latitude')+
  scale_color_viridis_c(name = "Date")+
  coord_sf(xlim = c(-45, -15), ylim = c(55, 66))+
  theme_minimal()+
  ggtitle("Nitrate")

ggplot(location_sum)+
  geom_point(aes(x = TIME, y = var))

matched_profiles <- location_fluo$prof_id[location_fluo$prof_id %in% location_nitrate$prof_id]

fluo_nitrate_df <- argo2 |> 
  mutate(prof_id = paste0(TIME, LONGITUDE, sep = "_")) |> 
  filter(prof_id %in% matched_profiles) |> 
  select(prof_id, TIME, LONGITUDE, LATITUDE, depth, PSAL, PSAL_QC, TEMP, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED) |> 
  filter(PSAL_QC <= 3) |> 
  group_by(prof_id, TIME, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()
  
test <- filter(fluo_nitrate_df, !is.na(CHLA_ADJUSTED))

nona <- argo2 |> filter(!is.na(CHLA_ADJUSTED) | !is.na(NITRATE_ADJUSTED)) |> 
  mutate(prof_id = paste0(TIME, LONGITUDE, sep = "_")) 

new_ref <- select(nona, prof_id) |> unique()
depth <- tibble("depth" = c(0:200))

new_ref <- new_ref |> crossing(depth) |>
  left_join(nona) |>
  mutate(JULD_date = lubridate::date(TIME))

interp_data <- new_ref |> group_by(prof_id) %>%
  filter(!prof_id %in% c("2024-01-12 12:20:00-42.5394833333333_", "2024-06-29 14:33:01-17.9631536666667_")) |> 
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED), 
                ~ na.approx(.x, depth, na.rm = FALSE))) %>%
  mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED), 
                ~ na.locf(.x, na.rm = FALSE, fromLast = TRUE))) %>%
  ungroup() |> 
  group_by(prof_id, depth) |> 
  summarise_all(mean) |> 
  ungroup()

table(is.na(interp_data$NITRATE_ADJUSTED))

test <- filter(interp_data, prof_id == unique(nona$prof_id)[120])
ggplot(test)+
  geom_point(aes(x = CHLA_ADJUSTED, y = -depth))

fluorescence_icb <- filter(interp_data, !is.na(CHLA_ADJUSTED) & !is.na(BBP700_ADJUSTED))

fluorescence_icb <- fluorescence_icb |> mutate(mth = lubridate::month(TIME))
test <- filter(fluorescence_icb, prof_id == unique(fluorescence_icb$prof_id)[12])
ggplot(fluorescence_icb)+
  geom_point(aes(x = CHLA_ADJUSTED, y = BBP700_ADJUSTED, colour = mth))

plt_df <- filter(fluorescence_icb, depth < 40) |> 
  mutate(ratio = CHLA_ADJUSTED/BBP700_ADJUSTED) |> 
  filter(ratio < 1000 & ratio > 0) |> 
  group_by(prof_id, mth, TIME, LONGITUDE, LATITUDE) |> 
  summarise(ratio = mean(ratio)) |> 
  ungroup()

ggplot(na.omit(plt_df))+
  geom_boxplot(aes(y = ratio, x = as.factor(mth)))+
  xlab("Month")+
  ylab("[Chla]/bbp")

dates_df <- interp_data |> select(prof_id, TIME) |> unique() |> na.omit() |> rename("JULD" = "TIME")

interp_data <- interp_data |> left_join(dates_df) |> select(-TIME)

#write_parquet(interp_data, "data/argo_pq/float_icb_fluorescence_cleaned.parquet")

npp_est <- vroom::vroom('data/argo_icb_pp_estimations_floats.csv')

table(is.na(npp_est$pp))

library(zoo)

plt_df <- npp_est |>
  filter(pp > 0) |>
  mutate(mth = lubridate::month(JULD),
         year = lubridate::year(JULD),
         doy = lubridate::yday(JULD)) |> 
  group_by(prof_id, mth, JULD, year, doy) |> 
  summarise(integrated_pp = sum(pp),
            integrated_nitrate = sum(NITRATE_ADJUSTED)) |> 
  ungroup() |> 
  select(-prof_id) |> 
  group_by(doy, year) |> 
  summarise_all(mean) |> 
  ungroup() |> 
  arrange(JULD) |> 
  filter(!is.na(integrated_pp))

ggplot(plt_df)+
  geom_path(aes(x = doy, y = integrated_pp, color = as.factor(year), group = as.factor(year)))+
  scale_color_brewer(palette = 'Paired', name = "Year")

ggplot(plt_df)+
  geom_path(aes(x = doy, y = integrated_pp, color = as.factor(year), group = as.factor(year)))+
  scale_color_brewer(palette = 'Paired', name = "Year")+
  facet_wrap(.~year, ncol = 3)

ggplot(na.omit(plt_df))+
  geom_boxplot(aes(y = integrated_pp, x = as.factor(mth), fill = as.factor(year)))

ts_npp <- npp_est |> 
  filter(pp>0) |> 
  group_by(prof_id, JULD) |> 
  summarise(integrated_pp = sum(pp), integrated_chl = sum(CHLA_ADJUSTED), mu = mean(mu)) |> 
  ungroup() |> 
  select(integrated_pp, integrated_chl, mu, JULD) |>
  na.omit() |> 
  mutate(date = lubridate::date(JULD)) |>
  group_by(date) |> 
  summarise_all(mean) |> 
  ungroup() |> 
  select(-JULD) |> 
  filter(integrated_chl <= 100)

ggplot(ts_npp)+
  geom_point(aes(x = integrated_chl, y = integrated_pp, colour = date))+
  scale_color_viridis_c()

ggplot(ts_npp)+
  geom_path(aes(x = date, y = integrated_pp))

ggplot(filter(ts_npp, mu >= 0))+
  geom_path(aes(x = date, y = mu))
