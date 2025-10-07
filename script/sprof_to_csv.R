library(ncdf4)
library(dplyr)
library(tidyr)
library(lubridate)

# open the Sprof NetCDF
nc <- nc_open("data/argo_nc/7902223/7902223_Sprof.nc")

# extract metadata
lon   <- ncvar_get(nc, "LONGITUDE")
lat   <- ncvar_get(nc, "LATITUDE")
juld  <- ncvar_get(nc, "JULD")  # days since 1950-01-01 00:00:00
date  <- as.Date(juld, origin = "1950-01-01")

# extract depth and variables (dimensions usually [levels x profiles])
depth <- ncvar_get(nc, "PRES")  
chla  <- ncvar_get(nc, "CHLA_ADJUSTED")
bbp   <- ncvar_get(nc, "BBP700_ADJUSTED")   # sometimes "BBP700_ADJUSTED"
ipar  <- ncvar_get(nc, "DOWNWELLING_PAR")
temp  <- ncvar_get(nc, "TEMP")
sal   <- ncvar_get(nc, "PSAL")
nitrate <- ncvar_get(nc, "NITRATE")

# close connection
nc_close(nc)

# get dimensions
n_levels   <- dim(depth)[1]
n_profiles <- dim(depth)[2]

# repeat metadata across depth levels
lon_rep  <- rep(lon, each = n_levels)
lat_rep  <- rep(lat, each = n_levels)
date_rep <- rep(date, each = n_levels)

# flatten arrays into vectors
df <- tibble(
  lon   = lon_rep,
  lat   = lat_rep,
  date  = date_rep,
  depth = as.vector(depth),
  chla  = as.vector(chla),
  bbp   = as.vector(bbp),
  iPAR  = as.vector(ipar),
  temp  = as.vector(temp),
  sal   = as.vector(sal),
  nitrate = as.vector(nitrate)
)

df_nona <- filter(df, !is.na(temp))
df_nona <- df_nona |> mutate(depth = round(depth))

df_interp <- df_nona %>%
  group_by(lon, lat, date) %>%
  arrange(depth, .by_group = TRUE) %>%
  reframe(
    depth = seq(0,
                1000, 
                by = 1)
  ) %>%
  left_join(df_nona, by = c("lon", "lat", "date", "depth")) %>%
  arrange(lon, lat, date, depth) %>%
  group_by(lon, lat, date) %>%
  mutate(across(
    c(chla, bbp, iPAR, temp, sal, nitrate),
    ~ if (sum(!is.na(.x)) > 1) {
      na.approx(.x, x = depth, na.rm = FALSE, rule = 2)
    } else {
      rep(NA_real_, length(.x))  # if no interpolation possible
    }
  )) %>%
  ungroup() |> 
  group_by(lon, lat, date, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

df_interp <- df_interp %>%
  arrange(date) %>% 
  mutate(prof_number = dense_rank(date))


# Plotting raw data -------------------------------------------------------

# create folders
dir.create("output/plots/temp_sal", recursive = TRUE, showWarnings = FALSE)
dir.create("output/plots/chla_bbp", recursive = TRUE, showWarnings = FALSE)

for (p in unique(df_interp$prof_number)) {
  
  df_prof <- df_interp %>% filter(prof_number == p)
  profile_date <- unique(df_prof$date)
  
  ## ---- Plot 1: Temperature + Salinity ----
  scale_factor_ts <- max(df_prof$temp, na.rm = TRUE) / 
    max(df_prof$sal,  na.rm = TRUE)
  
  p1 <- ggplot(df_prof, aes(y = depth)) +
    geom_path(aes(x = temp, colour = "Temperature")) +
    geom_path(aes(x = sal * scale_factor_ts, colour = "Salinity")) +
    scale_y_reverse(expand = c(0,0)) +
    scale_x_continuous(
      name = "Temperature (°C)",
      sec.axis = sec_axis(~ . / scale_factor_ts, name = "Salinity (PSU)")
    ) +
    labs(
      title = paste("Profile", p, "-", profile_date),
      y = "Depth (m)",
      colour = ""
    ) +
    theme_minimal()
  
  ## ---- Plot 2: Chla + BBP (dual x axis, top 100 m) ----
  df_prof100 <- df_prof %>% filter(depth <= 300)
  
  # scaling factor between Chla and BBP ranges
  scale_factor <- max(df_prof100$chla, na.rm = TRUE) / 
    max(df_prof100$bbp,  na.rm = TRUE)
  
  p2 <- ggplot(df_prof100, aes(y = depth)) +
    geom_path(aes(x = chla, colour = "Chla")) +
    geom_path(aes(x = bbp * scale_factor, colour = "BBP")) +
    scale_y_reverse(expand = c(0,0)) +
    scale_x_continuous(
      name = "Chla (mg/m³)",
      sec.axis = sec_axis(~ . / scale_factor, name = "BBP (1/m)")
    ) +
    labs(
      title = paste("Profile", p, "-", profile_date, "(0–300 m)"),
      y = "Depth (m)",
      colour = ""
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path("output/plots/chla_bbp", paste0("profile_", p, "_chla_bbp.png")),
    plot = p2,
    width = 6, height = 8
  )
}


# save the csv ------------------------------------------------------------


write_csv(df_interp, "argo_7902223_interp.csv")
