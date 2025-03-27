library(arrow)
library(tidyverse)
library(readxl)
library(zoo)
library(imputeTS)
library(ggrepel)
library(castr)

# color palette creation --------------------------------------------------
gbh <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

# dataloading -------------------------------------------------------------


dat <- read_parquet("data/glider/all_gliders_npp.parquet")

#data averaged per glider and date-hour
npp_dat <- dat |>
  group_by(glider, date, hour, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

#Zeu value from cbpm algorithm
daily_zeu <- dat |> select(datetime, zeu, depth) |> 
  na.omit() |> 
  mutate(date = date(datetime)) |> 
  group_by(date) |> 
  summarise_all(mean) |> 
  ungroup()

# daily_parze <- npp_dat |>
#   filter(glider %in% c("churchill", "nelson")) |> 
#   select(hour, date, par_insitu, depth) |> 
#   filter(hour > 10 & hour < 16) |> 
#   group_by(date, hour, depth) |>
#   summarise_all(mean, na.rm = TRUE) |> 
#   ungroup() |> 
#   mutate(prof_id = paste(date, hour, sep = "_")) |> 
#   filter(!prof_id %in% c("2024-06-04_13", "2024-06-14_16", "2024-07-03_13")) |> 
#   group_by(date, hour) |> 
#   mutate(zeu_par = ze_from_par(par_insitu, depth)) |>
#   ungroup() |> 
#   select(date, zeu_par) |> 
#   group_by(date) |> 
#   summarise_all(mean, na.rm = TRUE) |> 
#   ungroup()
#   
# daily_zeu <- daily_zeu |> left_join(daily_parze)

ggplot(daily_zeu)+
  geom_path(aes(x = date, y = zeu, color = "Ze from Chla"))+
  theme_bw()

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

day1 <- lubridate::as_datetime("2024-07-15 01", format = "%Y-%m-%d %H")
day2 <- lubridate::as_datetime("2024-07-19 00", format = "%Y-%m-%d %H")
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
  mutate(par_insitu = par_insitu *3600) |> 
  summarise_all(sum) |>
  ungroup() |> 
  mutate(par_insitu = par_insitu * 10e-6)

daily_modis <- dat |> select(datetime, par, depth) |> 
  filter(depth < 3) |> 
  na.omit() |> 
  mutate(date = date(datetime)) |>
  select(-datetime, -depth) |> 
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

daily_par <- left_join(daily_par, daily_modis) |> 
  arrange(date) |> 
  mutate(new_par = (par_insitu/max(par_insitu)) * 1000,
         weekly_par = rollapply(par_insitu, width =4, FUN = sum, align = "right", fill = NA)
  )


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

npp_dat <- npp_dat |> 
  mutate(datetime = lubridate::as_datetime(paste(date, hour, sep = " "), format = "%Y-%m-%d %H"))

plot_dat <- npp_dat |> filter(npp < 50 & depth < 50)

# NPP models comparison ---------------------------------------------------



int_npp <- npp_dat |> 
  group_by(datetime, glider) |> 
  summarise(npp_int = sum(npp, na.rm = TRUE),
            vgpm = mean(npp_vgpm),
            par = mean(par),
            zeu = mean(zeu, na.rm = TRUE)) |>
  ungroup() 

library(forecast)
int_npp <- int_npp |> 
  mutate(outlier = row_number() %in% tsoutliers(vgpm)$index,
         vgpm = tsclean(vgpm))

ggplot(int_npp)+
  geom_path(aes(x = datetime, y = npp_int, colour = glider))

ggplot(int_npp)+
  geom_path(aes(x = datetime, y = vgpm, colour = glider))


daily_npp <- int_npp |> 
  mutate(date = date(datetime)) |> 
  group_by(date) |> 
  summarise(daily_npp = mean(npp_int, na.rm = TRUE),
            daily_vgpm = mean(vgpm, na.rm = TRUE),
            var_npp = sd(npp_int, na.rm = TRUE),
            var_vgpm = sd(vgpm, na.rm = TRUE)) |> 
  ungroup() |> 
  left_join(daily_par)
  

gpp <- read_excel("~/Downloads/Median_GPP_data.xlsx", 
                  col_types = c("date", "numeric"))

names(gpp) <- c("date", "gpp")

gpp <- gpp |> mutate(date = date(date))

synth_pp <- daily_npp |> left_join(gpp)


ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp))+
  geom_path(aes(x = date, y = daily_npp - var_npp), linetype = "dashed")+
  geom_path(aes(x = date, y = daily_npp + var_npp), linetype = "dashed")+
  geom_path(aes(x = date, y = weekly_par/2), color = "grey")+
  theme_bw()+
  ylab("NPP from CbPM (mg C m-2 d-1)")

ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_vgpm/2))+
  geom_path(aes(x = date, y = daily_vgpm/2 - var_vgpm/2), linetype = "dashed")+
  geom_path(aes(x = date, y = daily_vgpm/2 + var_vgpm/2), linetype = "dashed")+
  geom_path(aes(x = date, y = weekly_par*2), color = "grey")+
  theme_bw()+
  ylab("NPP from VGPM (mg C m-2 d-1)")


ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp, colour = "NPP from CbPM"))+
  geom_path(aes(x = date, y = daily_vgpm/2, colour = "NPP from VGPM"))+
  geom_path(aes(x = date, y = gpp, colour = "GPP"))+
  geom_path(aes(x = date, y = weekly_par * 2), color = "lightgrey")+
  theme_bw()+
  ggtitle("NPP estimation from CbPM compared to GPP estimation for Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")

npp_df <- read_csv("data/argo_pp_estimations_floats.csv") |> 
  mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth),
         sigma0 = gsw_sigma0(PSAL, TEMP)) |> 
  mutate(date = date(JULD))
  

vgpm_floats <- read_csv("data/argo_pp_estimations_floats_2.csv") |> 
  mutate(date = date(JULD)) |> 
  select(date, npp_vgpm) |> 
  na.omit() |> 
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

npp_df <- left_join(npp_df, vgpm_floats)

synth_float <- npp_df |>
  group_by(JULD) |> 
  summarise(npp = sum(pp, na.rm = TRUE)) |> 
  ungroup() |> 
  mutate(date = date(JULD)) |>
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  left_join(vgpm_floats) |> 
  arrange(date) |> 
  na.omit() |> 
  mutate(cbpm_smoothed = rollapply(npp, width =5, FUN = mean, align = "right", fill = NA),
         vgpm_smoothed = rollapply(npp_vgpm, width =5, FUN = mean, align = "right", fill = NA))

ggplot(synth_pp)+
  geom_path(aes(x = date, y = cbpm_smoothed, color = "CbPM applied to floats"), data = synth_float, linewidth = 2)+
  geom_path(aes(x = date, y = vgpm_smoothed, color = "VGPM applied to floats"), data = synth_float, linewidth = 2)+
  theme_bw(base_size = 14)+
  ggtitle("NPP estimation from CbPM and VGPM compared to GPP estimation from Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")+
  xlim(c(date("2024-04-01"), date("2024-12-01")))+
  scale_color_manual(values = c("CbPM applied to gliders" = "steelblue1",
                                "CbPM applied to floats" = "steelblue4",
                                "VGPM applied to gliders" = "tomato",
                                "VGPM applied to floats" = "tomato4",
                                "GPP estimation" = "olivedrab4",
                                "Light history" = "grey"),
                     name = "")

ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp, colour = "CbPM applied to gliders"), linewidth = 2)+
  geom_path(aes(x = date, y = daily_vgpm/2, colour = "VGPM applied to gliders"),  linewidth = 2)+
  geom_path(aes(x = date, y = gpp, colour = "GPP estimation"),  linewidth = 2)+
  geom_path(aes(x = date, y = cbpm_smoothed, color = "CbPM applied to floats"), data = synth_float, linewidth = 2)+
  geom_path(aes(x = date, y = vgpm_smoothed, color = "VGPM applied to floats"), data = synth_float, linewidth = 2)+
  theme_bw(base_size = 14)+
  ggtitle("NPP estimation from CbPM and VGPM compared to GPP estimation from Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")+
  xlim(c(date("2024-04-01"), date("2024-12-01")))+
  scale_color_manual(values = c("CbPM applied to gliders" = "steelblue1",
                                "CbPM applied to floats" = "steelblue4",
                                "VGPM applied to gliders" = "tomato",
                                "VGPM applied to floats" = "tomato4",
                                "GPP estimation" = "olivedrab4",
                                "Light history" = "grey"),
                     name = "")

ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp, colour = "CbPM applied to gliders"), linewidth = 2)+
  geom_path(aes(x = date, y = daily_vgpm/2, colour = "VGPM applied to gliders"),  linewidth = 2)+
  geom_path(aes(x = date, y = gpp, colour = "GPP estimation"),  linewidth = 2)+
  geom_path(aes(x = date, y = cbpm_smoothed, color = "CbPM applied to floats"), data = synth_float, linewidth = 2)+
  geom_path(aes(x = date, y = vgpm_smoothed, color = "VGPM applied to floats"), data = synth_float, linewidth = 2)+
  theme_bw(base_size = 14)+
  ggtitle("NPP estimation from CbPM and VGPM compared to GPP estimation from Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")+
  xlim(c(date("2024-06-01"), date("2024-10-01")))+
  scale_color_manual(values = c("CbPM applied to gliders" = "steelblue1",
                                "CbPM applied to floats" = "steelblue4",
                                "VGPM applied to gliders" = "tomato",
                                "VGPM applied to floats" = "tomato4",
                                "GPP estimation" = "olivedrab4",
                                "Light history" = "grey"),
                     name = "")

ggsave("output/plots/pp_floats/pp_algorithms_comparison.png", dpi = 300, width = 30, height = 18, units = "cm")



tot_dat <- left_join(select(synth_pp, -par), select(npp_dat, -par_insitu))


# Perform PCA analysis of the whole datatset, avergaed per day and --------


daily_surf <- tot_dat |> 
  filter(depth < 5) |> 
  select(-glider) |>
  mutate(par_insitu_normalized = par_insitu/max(par_insitu, na.rm = TRUE),
         par_weekly_normalized = weekly_par/max(weekly_par, na.rm = TRUE),
         par_normalized = par/max(par)) |> 
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup()

ggplot(daily_surf)+
  geom_smooth(aes(x = par_insitu_normalized, y = npp, color = "daily"), method = "lm", se = FALSE)+
  geom_point(aes(x = par_insitu_normalized, y = npp, color = "daily"), alpha = 0.5)+
  geom_smooth(aes(x = par_weekly_normalized, y = npp, color = "4 days avg"), method = "lm", se = FALSE)+
  geom_point(aes(x = par_weekly_normalized, y = npp, color = "4 days avg"), alpha = 0.5)+
  geom_smooth(aes(x = par_normalized, y = npp, color = "modis"), method = "lm", se = FALSE)+
  geom_point(aes(x = par_normalized, y = npp, color = "modis"), alpha = 0.5)+
  theme_bw()+
  xlab("Integrated PAR (normalized)")+
  ylab("NPP from CbPM (mg C m-2 d-1)")+
  scale_color_discrete(name = "")
  
pca_df <- select(daily_surf, -datetime, new_par, -weekly_par, -hour, -depth, -profile_index, -par_insitu_normalized, -par_normalized, -par_weekly_normalized, -new_par, -var_npp, -daily_npp, -var_vgpm, -daily_vgpm) |> 
  na.omit() |> 
  mutate(fluo = fluo/2)

float_pca_df <- npp_df |> 
  mutate(date = date(JULD),
         platform = "float") |> 
  filter(depth < 50) |> 
  select(date,
         longitude = LONGITUDE,
         latitude = LATITUDE,
         temperature = TEMP,
         salinity = PSAL,
         fluo = CHLA_ADJUSTED,
         bbp700 = bbp700_cleaned,
         par = satellite_par,
         mu,
         npp = pp,
         zeu,
         npp_vgpm,
         platform) |>
  group_by(date, platform) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  na.omit()

pca_df <- pca_df |> mutate(platform = "glider") |> 
  select(-par_insitu,
         -gpp) |> 
  bind_rows(float_pca_df) |> 
  select(-longitude, - latitude) |> na.omit()

ggplot(pca_df)+
  geom_point(aes(x = date, y = fluo, colour = platform))

ggplot(pca_df)+
  geom_point(aes(x = date, y = bbp700, colour = platform))


# Select only numeric columns for PCA
numeric_cols <- pca_df %>% select(where(is.numeric))

# Perform PCA
pca_result <- prcomp(numeric_cols, center = TRUE, scale. = TRUE)


pca_loadings <- as.data.frame(pca_result$rotation)
loadings_scale <- 6  # Adjust this if needed to control arrow length
pca_loadings <- pca_loadings %>%
  rownames_to_column("variable") %>%
  mutate(PC1 = PC1 * loadings_scale,
         PC2 = PC2 * loadings_scale)



# Create a dataframe for visualization
pca_scores <- as.data.frame(pca_result$x)
pca_scores$date <- pca_df$date  # Keep the dates for labeling if needed
pca_scores$platform <- pca_df$platform  # Keep the dates for labeling if needed

# Scree plot to visualize variance explained by each PC
ggplot(data = data.frame(PC = 1:length(pca_result$sdev),
                                       Variance = (pca_result$sdev^2) / sum(pca_result$sdev^2))) +
  geom_bar(aes(x = PC, y = Variance), stat = "identity", fill = "#69b3a2") +
  geom_line(aes(x = PC, y = Variance), color = "blue") +
  labs(title = "Scree Plot", x = "Principal Components", y = "Variance Explained") +
  theme_minimal()

# PCA biplot visualization
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = date), size = 3) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "steelblue3", linewidth = 1) +
  geom_text_repel(data = pca_loadings, aes(x = PC1, y = PC2, label = variable), size = 6)+
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  scale_color_gradientn(
    colors = gbh,
    name = "Date",
    breaks = pretty(pca_scores$date, n = 4),
    labels = scales::date_format("%b")  
  )+
  theme_minimal()

ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = platform), size = 3) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "steelblue3", linewidth = 1) +
  geom_text_repel(data = pca_loadings, aes(x = PC1, y = PC2, label = variable), size = 6)+
  labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
  theme_minimal()

ggplot(pca_scores, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = date), size = 3) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "steelblue3", linewidth = 1) +
  geom_text_repel(data = pca_loadings, aes(x = PC2, y = PC3, label = variable), size = 6)+
  labs(title = "PCA Biplot", x = "PC2", y = "PC3") +
  scale_color_gradientn(
    colors = gbh,
    name = "Date",
    breaks = pretty(pca_scores$date, n = 4),
    labels = scales::date_format("%b")  
  )+
  theme_minimal()

location <- npp_dat |> 
  select(date, latitude, longitude) |>
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  left_join(pca_scores)

ggplot(location)+
  geom_point(aes(x = longitude, y = latitude, color = PC2))



# some plot of gliders transect -------------------------------------------


doombar <- filter(dat, glider == "doombar")

ggplot(doombar)+
  geom_tile(aes(y = -depth, x = datetime, color = npp))
