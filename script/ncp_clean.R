library(tidyverse)
library(lubridate)
library(castr)
library(zoo)

dat <- read_csv("output/float_nitrate_data_corrected.csv")


# MLD and Zeu time series -------------------------------------------------

#Create a dataset with MLD and Zeu for each profile 
dat_prof <- dat |> select(float_wmo, prof_number, date, MLD, zeu) |> 
  distinct() |> 
  arrange(date)

#Smooth the MLD and Zeu time series using a spline
# Remove NAs before smoothing
dat_smooth <- dat_prof %>%
  filter(!is.na(date), !is.na(MLD))

# Fit spline
fit_mld <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$MLD, spar = 0.6)
fit_zeu <- smooth.spline(as.numeric(dat_smooth$date), dat_smooth$zeu, spar = 0.6)
# Predict smoothed values for all dates
pred_mld <- predict(fit_mld, as.numeric(dat_prof$date))$y
pred_zeu <- predict(fit_zeu, as.numeric(dat_prof$date))$y


# Add back to main data frame
dat_prof <- dat_prof %>%
  mutate(mld_smooth = pred_mld,
         zeu_smooth = pred_zeu)


ggplot(dat_prof)+
  geom_line(aes(x = date, y = -mld_smooth, color = "MLD"))+
  geom_point(aes(x = date, y = -mld_smooth), shape = 1)+
  geom_line(aes(x = date, y = -zeu_smooth, color = "Zeu"))+
  labs(y = "Depth (m)", color = "Parameter")+
  theme_minimal()
  
dat <- left_join(dat, select(dat_prof, prof_number, float_wmo, mld_smooth, zeu_smooth), by = c("float_wmo", "prof_number"))


# Integration depth computation -------------------------------------------


#Define Integration depth (i.e. deepest value of MLD or Zeu between t and t-1)
integration_depth_data <- dat |>
  group_by(float_wmo, prof_number, date) |>
  summarise(
    mld_smooth = unique(mld_smooth, na.rm = TRUE),
    zeu_smooth = unique(zeu_smooth, na.rm = TRUE)) |> 
  ungroup()

#Regrid the data on a 10day interval timeseries
time_grid <- tibble(date_10day = seq(min(integration_depth_data$date),
                                     max(integration_depth_data$date),
                                     by = "12 days"))

prof_dat_smoothed <- integration_depth_data %>%
  mutate(date_10day = as.Date(cut(date, breaks = "12 days"))) %>%
  group_by(date_10day) %>%
  summarise(mld = mean(mld_smooth, na.rm = TRUE),
            zeu = mean(zeu_smooth, na.rm = TRUE)) |> 
  mutate(
    prev_MLD = dplyr::lag(mld),
    prev_zeu = dplyr::lag(zeu),
    NCP_integration_depth = pmax(mld, prev_MLD, zeu, prev_zeu, na.rm = TRUE),
    next_integration_depth = lead(NCP_integration_depth)
  ) |> 
  full_join(time_grid, by = "date_10day") |>
  arrange(date_10day) %>%
  mutate(
    NCP_integration_depth = zoo::na.approx(NCP_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    next_integration_depth = zoo::na.approx(next_integration_depth, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)
  ) %>%
  filter(date_10day %in% time_grid$date_10day) |> 
  mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
         NCP_integration_depth_smooth = rollapply(NCP_integration_depth, width = 3, FUN = mean, align = "center", fill = NA),
         next_integration_depth_smooth = rollapply(next_integration_depth, width = 3, FUN = mean, align = "center", fill = NA)
  )

ggplot(prof_dat_smoothed)+
  geom_line(aes(x = date_10day, y = -NCP_integration_depth_smooth))+
  geom_line(aes(x = date_10day, y = -next_integration_depth_smooth), linetype = "dashed")

prof_dat_smoothed <- na.omit(prof_dat_smoothed)


# Interpolation of data  --------------------------------------------------

#Create a dataset with the median value of nitrate integrals every 10 days for smoothing
dat_smoothed <- dat %>%
  mutate(date_10day = as.Date(cut(date, breaks = "12 days"))) %>%
  group_by(date_10day, depth) %>%
  summarise(
    nitrate = mean(nitrate_corrected, na.rm = TRUE),
    npp = mean(cbpm_npp, na.rm = TRUE)) |> 
  ungroup() |> 
  full_join(time_grid, by = "date_10day") %>%     # ensure regular 10-day spacing
  arrange(date_10day) %>%
  group_by(depth) |> 
  mutate(
    nitrate = na.approx(nitrate, x = as.numeric(date_10day), na.rm = FALSE, rule = 2),
    npp = na.approx(npp, x = as.numeric(date_10day), na.rm = FALSE, rule = 2)) %>%
  ungroup() |> 
  filter(date_10day %in% time_grid$date_10day)

ggplot(dat_smoothed)+
  geom_point(aes(x = nitrate, y = - depth, color = date_10day))

dat_smoothed <- dat_smoothed |> 
  left_join(prof_dat_smoothed, by = "date_10day") |> 
  na.omit()


# Integration of nitrate over integration depth ---------------------------

dat_final <- dat_smoothed |>
  group_by(date_10day, mld, zeu) |>
  filter(!is.na(next_integration_depth)) |> 
  summarise(
    int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(NCP_integration_depth)),
    next_int_N_mmol_m2 = integrate(nitrate, depth, from = 0, to = unique(next_integration_depth)),
    npp = integrate(npp, depth, from = 0, to = 200),
    mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) - 20, to = unique(NCP_integration_depth) - 10) / 10,
    sub_mld_concentration = integrate(nitrate, depth, from = unique(NCP_integration_depth) + 20, to = unique(NCP_integration_depth) + 30) / 10
  ) |> 
  ungroup()

ggplot(dat_final)+
  geom_line(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_point(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_line(aes(x = date_10day, y = next_int_N_mmol_m2), linetype = "dashed")

#Smoothing the integration
dat_final_smoothed <- dat_final %>%
  mutate(
    int_N_smooth = rollapply(int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA),
    next_int_N_smooth = rollapply(next_int_N_mmol_m2, width = 3, FUN = mean, align = "center", fill = NA)
  )

ggplot(dat_final_smoothed)+
  geom_line(aes(x = date_10day, y = int_N_mmol_m2))+
  geom_line(aes(x = date_10day, y = int_N_smooth), linetype = "dashed")


# Computing NCP -----------------------------------------------------------


ncp_results <- dat_final_smoothed |> 
  arrange(date_10day) |> 
  mutate(dt = as.numeric(date_10day - dplyr::lag(date_10day)),
         diff_mld = mld - lag(mld),
         we = pmax(0, diff_mld),
         delta = sub_mld_concentration - mld_concentration,
         diff = delta * we,
         nitrate_consumption =  dplyr::lag(next_int_N_smooth) - int_N_smooth,
         c_consumption = nitrate_consumption * 6.625 / dt,
         NCP = c_consumption + diff) |> 
  ungroup()

ncp_results <- ncp_results |> 
  mutate(ncp_smooth = rollapply(NCP, width = 6, FUN = mean, align = "center", fill = NA),
         ncp_sd = rollapply(NCP, width = 6, FUN = sd, align = "center", fill = NA),
         npp_smooth = rollapply(npp, width = 6, FUN = median, align = "center", fill = NA))

ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = NCP))+
  geom_line(aes(x = date_10day, y = npp_smooth/12))+
  geom_point(aes(x = date_10day, y = NCP, color = we))+
  scale_color_viridis_c()


ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = npp_smooth, color = "NPP"))+
  geom_ribbon(aes(x = date_10day, ymin = ncp_smooth*12 - ncp_sd, ymax = ncp_smooth*12 + ncp_sd), alpha = 0.3)+
  geom_line(aes(x = date_10day, y = ncp_smooth*12, color = "NCP"))+
  labs(y = "mmol C m-2 d-1")+
  theme_bw()


ggplot(ncp_results)+
  geom_line(aes(x = date_10day, y = -mld))+
  geom_point(aes(x = date_10day, y = -mld, color = diff))+
  scale_color_viridis_c()

ggplot(ncp_results)+
  geom_point(aes(x = date_10day, y = - mld, color = NCP))+
  scale_color_distiller(palette = "RdBu")+
  theme_dark()
# draft -------------------------------------------------------------------

# # Ensure data are ordered by profile and date
# dat <- dat |> arrange(prof_number, date, depth)
# 
# # Step 1: Compute integrated nitrate down to each possible depth limit
# integrate_to_depth <- function(df, depth_limit) {
#   sub <- df |> filter(depth <= depth_limit)
#   if (nrow(sub) < 2) return(NA)  # Need at least 2 points for integration
#   trapz(sub$depth, sub$nitrate)
# }
# 
# # Step 2: Compute integrated nitrate to 200 m for sanity check
# dat_integrated_200 <- dat |>
#   group_by(prof_number, date, MLD) |>
#   summarise(int_N_200_mmol_m2 = integrate_to_depth(cur_data(), 200), .groups = "drop")
# 
# # Step 3: Compute NCP integration depth dynamically
# dat_integrated <- dat_integrated_200 |>
#   arrange(date) |>
#   mutate(
#     prev_MLD = lag(MLD),
#     NCP_integration_depth = pmax(MLD, prev_MLD, 40, na.rm = TRUE),
#     next_integration_depth = lead(NCP_integration_depth)
#   )
# 
# dat_joined <- left_join(dat, dat_integrated)
# 
# dat_final <- dat_joined |>
#   group_by(prof_number, date, MLD, NCP_integration_depth, next_integration_depth) |>
#   summarise(
#     int_N_mmol_m2 = integrate_to_depth(pick(nitrate, depth), unique(NCP_integration_depth)),
#     next_int_N_mmol_m2 = integrate_to_depth(pick(nitrate, depth), unique(next_integration_depth)),
#     .groups = "drop"
#   )
# 
# dat_final <- dat_final |> 
#   mutate(nitrate_consumption = int_N_mmol_m2 - lag(next_int_N_mmol_m2),
#          c_consumption = nitrate_consumption * 6.625 / 10,
#          NCP = - c_consumption * 12)  # Convert to mmol C m-2
# 
# ggplot(filter(dat_final, NCP > -4000))+
#   geom_line(aes(x = date, y = NCP))
