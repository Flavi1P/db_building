library(arrow)
library(tidyverse)
library(readxl)
library(zoo)
library(imputeTS)

dat <- read_parquet("data/glider/all_gliders_npp.parquet")

npp_dat <- dat |>
  group_by(glider, date, hour, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

daily_zeu <- dat |> select(datetime, zeu, depth) |> 
  na.omit() |> 
  mutate(date = date(datetime)) |> 
  group_by(date) |> 
  summarise_all(mean) |> 
  ungroup()

ggplot(daily_zeu)+
  geom_path(aes(x = date, y = zeu))

# Looking at PAR comparison glider vs modis -------------------------------


datetime_par <- dat |> select(date, hour, par_insitu, par, depth) |>
  filter(depth == 0 & par_insitu > 0) |> 
  group_by(date, hour) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

all_datetime_par <- datetime_par |>
  expand(date, hour) |>
  left_join(datetime_par) |> 
  mutate(datetime = lubridate::as_datetime(paste(date, hour, sep = " "), format = "%Y-%m-%d %H"))

day1 <- lubridate::as_datetime("2024-07-01 01", format = "%Y-%m-%d %H")
day2 <- lubridate::as_datetime("2024-10-02 00", format = "%Y-%m-%d %H")
one_week <- all_datetime_par |> filter(datetime > day1 & datetime < day2)

ggplot_na_distribution(all_datetime_par$par_insitu)

imp <- na_interpolation(all_datetime_par$par_insitu, option = "spline")
ggplot_na_imputations(all_datetime_par$par_insitu, imp)

all_datetime_par <- all_datetime_par |> 
  mutate(par_insitu = imp,
         par_insitu = case_when(par_insitu < 0 ~ 0,
                                TRUE ~ par_insitu))

daily_par <- all_datetime_par |> 
  select(datetime, par_insitu) |>
  mutate(date = date(datetime)) |> 
  select(-datetime) |> 
  group_by(date) |> 
  summarise_all(mean) |> 
  ungroup()

daily_modis <- dat |> select(datetime, par, depth) |> 
  filter(depth < 3) |> 
  na.omit() |> 
  mutate(date = date(datetime)) |>
  select(-datetime, -depth) |> 
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

daily_par <- left_join(daily_par, daily_modis)

ggplot(daily_par)+
  geom_point(aes(x = date, y = par_insitu), color = "steelblue")+
  geom_path(aes(x = date, y = par_insitu), color = "steelblue1")+
  geom_point(aes(x = date, y = par), color = "darksalmon")+
  geom_path(aes(x = date, y = par), color = "darksalmon")+
  theme_bw()

ggplot(daily_par)+
  geom_point(aes(x = par, y = par_insitu), color = "steelblue")+
  geom_smooth(aes(x = par, y = par_insitu), method = "lm")+
  theme_bw()


all_datetime <- seq(from = min(datetime_par$datetime), to = max(datetime_par$datetime), by = "hour")

npp_dat <- npp_dat |> 
  mutate(datetime = lubridate::as_datetime(paste(date, hour, sep = " "), format = "%Y-%m-%d %H"))

plot_dat <- npp_dat |> filter(npp < 50 & depth < 50)


int_npp <- npp_dat |> 
  group_by(datetime, glider) |> 
  summarise(npp_int = sum(npp),
            par = mean(par)) |>
  ungroup()

full_time <- tibble ("datetime" = all_datetime)

ggplot(int_npp)+
  geom_path(aes(x = datetime, y = npp_int, colour = glider))

# all_datetime <- seq(from = min(int_npp$datetime), to = max(int_npp$datetime), by = "hour")
# 
# full_time <- tibble(datetime = rep(all_datetime, 4),
#                     glider = c(rep("cabot", length(all_datetime)),
#                                rep("churchill", length(all_datetime)),
#                                rep("doombar", length(all_datetime)),
#                                rep("nelson", length(all_datetime))))
# 
# int_npp_filled <- full_time |> 
#   left_join(int_npp, by = c("glider", "datetime")) |> 
#   mutate(datatype = case_when(
#     is.na(npp_int) ~ "interpolated",
#     TRUE ~ "measured" 
#   ))
# 
# 
# test <- filter(int_npp_filled, glider == "doombar") |> select(-glider, -datatype)
# test$datetime <- as.POSIXct(test$datetime)
# ts_data <- zoo(test[, -1], order.by = test$datetime)
# 
# imputed_data <- na_kalman(ts_data)
# df_imputed <- data.frame(
#   datetime = index(imputed_data),
#   coredata(imputed_data)) |>
#   mutate(date = date(datetime)) |> 
#   group_by(date) |>
#   select(-datetime) |> 
#   summarise_all(list(mean, sd))
# 
# ggplot(df_imputed)+
#   geom_path(aes(x = date, y = npp_int_fn1))+
#   geom_path(aes(x = date, y = npp_int_fn1+npp_int_fn2), linetype = "dashed")+
#   geom_path(aes(x = date, y = npp_int_fn1-npp_int_fn2), linetype = "dashed")
# 
# 
# int_npp_filled <- int_npp_filled |> 
#   group_by(glider) |> 
#   mutate(npp_int = zoo::na.approx(npp_int, na.rm = FALSE)) |> 
#   ungroup() |> 
# 
# 
# 
# ggplot(int_npp_filled) +
#   geom_path(aes(x = datetime, y = npp_int, color = glider)) +
#   labs(title = "Interpolated NPP over Time", y = "NPP Int", x = "Datetime") +
#   theme_minimal()

daily_npp <- int_npp |> 
  mutate(date = date(datetime)) |> 
  group_by(date, glider) |> 
  summarise(daily_npp = mean(npp_int, na.rm = TRUE),
            var_npp = sd(npp_int, na.rm = TRUE)) |> 
  ungroup()

ggplot(daily_npp)+
  geom_path(aes(x = date, y = daily_npp, colour = glider))+
  geom_path(aes(x = date, y = daily_npp - var_npp, colour = glider), linetype = "dashed")+
  geom_path(aes(x = date, y = daily_npp + var_npp, colour = glider), linetype = "dashed")+
  theme_bw()
  

gpp <- read_excel("~/Downloads/Median_GPP_data.xlsx", 
                  col_types = c("date", "numeric"))

names(gpp) <- c("date", "gpp")

gpp <- gpp |> mutate(date = date(date))

synth_pp <- daily_npp |> left_join(gpp)

daily_par <- daily_par |> 
  mutate(new_par = (par_insitu/max(par_insitu)) * 1000)

synth_pp <- left_join(synth_pp, select(daily_par, date, new_par))

ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp, colour = "NPP from CbPM", group = glider))+
  geom_path(aes(x = date, y = gpp, colour = "GPP"))+
  geom_path(aes(x = date, y = new_par), linetype = "dashed", color = 'Grey')+
  theme_bw()+
  ggtitle("NPP estimation from CbPM compared to GPP estimation for Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")




ggplot(tot_dat)+
  geom_path(aes(x = date, y = npp))+
  geom_path(aes(x = date, y = new_par), linetype = "dashed")

ggplot(filter(tot_dat, par_insitu > 0))+
  geom_point(aes(x = par_insitu, y = par))+
  geom_smooth(aes(x = par_insitu, y = par), method = "lm")


# VgPM --------------------------------------------------------------------


calculate_Pbmax <- function(temp) {
  if (temp <= 20) {
    return(1.3 + (6.6 - 1.3) * (temp / 20))
  } else {
    return(6.6 * exp(-0.05 * (temp - 20)))
  }
}

alpha_B <- 3.1e-2

data <- data %>%
  mutate(
    temperature = 18
    Pbmax = map_dbl(temperature, calculate_Pbmax),  # Calculate Pbmax per row
    PP_inst = Pbmax * chl * (1 - exp(-PAR * (alpha_B / Pbmax)))
  )