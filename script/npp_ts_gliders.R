library(arrow)
library(tidyverse)
library(readxl)
library(zoo)
library(imputeTS)
library(ggrepel)
library(castr)
library(gsw)

# color palette creation --------------------------------------------------
gbh <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

# dataloading -------------------------------------------------------------


dat <- read_parquet("data/glider/all_gliders_npp.parquet") |> 
  filter(npp < 150)

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

day1 <- lubridate::as_datetime("2024-06-15 01", format = "%Y-%m-%d %H")
day2 <- lubridate::as_datetime("2024-06-19 00", format = "%Y-%m-%d %H")
one_week <- all_datetime_par |> filter(datetime > day1 & datetime < day2)

ggplot_na_distribution(one_week$par_insitu)
ggsave("output/plots/pp_floats/par_ts_filled.png", dpi = 300, width = 16, height = 9, units = "cm")

imp <- na_interpolation(all_datetime_par$par_insitu, option = "spline")
ggplot_na_imputations(all_datetime_par$par_insitu, imp)

all_datetime_par <- all_datetime_par |> 
  mutate(par_insitu = imp,
         par_insitu = case_when(par_insitu < 0 ~ 0,
                                TRUE ~ par_insitu),
         par_insitu = case_when(hour < 6 | hour > 20 ~ 0,
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
  ) |> 
  mutate(par_normalized = par_insitu/max(par_insitu),
         modis_normalized = par/max(par))


ggplot(daily_par)+
  geom_point(aes(x = date, y = par_normalized), color = "steelblue")+
  geom_path(aes(x = date, y = par_normalized), color = "steelblue1")+
  geom_point(aes(x = date, y = modis_normalized), color = "darksalmon")+
  geom_path(aes(x = date, y = modis_normalized), color = "darksalmon")+
  theme_bw()


ggsave("output/plots/pp_floats/par_ts_comparison.png", dpi = 300, width = 16, height = 9, units = "cm")
ggplot(daily_par)+
  geom_point(aes(x = par_normalized, y = modis_normalized), color = "steelblue")+
  geom_smooth(aes(x = par_normalized, y = modis_normalized), method = "lm")+
  theme_bw()
ggsave("output/plots/pp_floats/par_correlation_comparison.png", dpi = 300, width = 16, height = 9, units = "cm")

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
  geom_path(aes(x = date, y = daily_npp - var_npp), alpha = 0.2, linetype = "dashed")+
  geom_path(aes(x = date, y = daily_npp + var_npp), alpha = 0.2, linetype = "dashed")+
  geom_path(aes(x = date, y = weekly_par), color = "grey")+
  theme_bw()+
  ylab("NPP from CbPM (mg C m-2 d-1)")

ggsave("output/plots/pp_floats/pglider_cbpm.png", dpi = 300, width = 16, height = 9, units = "cm")

ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_vgpm))+
  geom_path(aes(x = date, y = daily_vgpm - var_vgpm), alpha = 0.2, linetype = "dashed")+
  geom_path(aes(x = date, y = daily_vgpm + var_vgpm), alpha = 0.2, linetype = "dashed")+
  theme_bw(base_size = 14)+
  ylab(expression("NPP estimated from VGPM" ~ (mg ~ C ~ m^{-2} ~ d^{-1})))+
  xlab("Date")
ggsave("output/plots/pp_floats/glider_vgpm.png", dpi = 300, width = 22, height = 13, units = "cm")


ggplot(synth_pp)+
  geom_path(aes(x = date, y = daily_npp, colour = "NPP from CbPM"))+
  geom_path(aes(x = date, y = daily_vgpm, colour = "NPP from VGPM"))+
  geom_path(aes(x = date, y = gpp, colour = "GPP"))+
  geom_path(aes(x = date, y = weekly_par * 2), color = "lightgrey")+
  theme_bw()+
  ggtitle("NPP estimation from CbPM compared to GPP estimation for Oxygen dial variability")+
  ylab("Producitivity (mg C m-2 d-1)")
ggsave("output/plots/pp_floats/glider_vgpm_cbpm.png", dpi = 300, width = 16, height = 9, units = "cm")

npp_df <- read_csv("data/argo_pp_estimations_floats_2.csv") |> 
  mutate(ct = gsw_CT_from_t(PSAL, TEMP, PRES_rounded),
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

no_pp <- npp_df |> select(JULD, pp) |> group_by(JULD) |> summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  filter(is.na(pp)) |> 
  pull(JULD)


synth_float <- npp_df |>
  filter(!JULD %in% no_pp) |> 
  group_by(JULD) |> 
  mutate(pp = na_interpolation(pp)) |> 
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
  geom_path(aes(x = date, y = daily_vgpm, colour = "VGPM applied to gliders"),  linewidth = 2)+
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
  geom_path(aes(x = date, y = daily_vgpm, colour = "VGPM applied to gliders"),  linewidth = 2)+
  geom_path(aes(x = date, y = gpp, colour = "GPP estimation"),  linewidth = 2)+
  geom_path(aes(x = date, y = cbpm_smoothed, color = "CbPM applied to floats"), data = synth_float, linewidth = 2)+
  geom_path(aes(x = date, y = vgpm_smoothed, color = "VGPM applied to floats"), data = synth_float, linewidth = 2)+
  theme_bw(base_size = 16)+
  ylab("Depth integrated NPP (mg C m-2 d-1)")+
  xlab("Date")+
  xlim(c(date("2024-06-01"), date("2024-10-01")))+
  scale_color_manual(values = c("CbPM applied to gliders" = "steelblue1",
                                "CbPM applied to floats" = "steelblue4",
                                "VGPM applied to gliders" = "tomato",
                                "VGPM applied to floats" = "tomato4",
                                "GPP estimation" = "olivedrab4",
                                "Light history" = "grey"),
                     name = "")

ggsave("output/plots/pp_floats/pp_algorithms_comparison2.png", dpi = 300, width = 28, height = 12, units = "cm")


corr_df <- synth_pp |> 
  left_join(select(synth_float, date, vgpm_smoothed, cbpm_smoothed)) |> 
  select(daily_npp, daily_vgpm, cbpm_smoothed, vgpm_smoothed, gpp) |> 
  pivot_longer(daily_npp:vgpm_smoothed, names_to = "method", values_to = "npp") |> 
  mutate(method = case_when(method == "cbpm_smoothed" ~ "CbPM Argo",
                            method == "daily_npp" ~ "CbPM Glider",
                            method == "vgpm_smoothed" ~ "VGPM Argo",
                            method == "daily_vgpm" ~ "VGPM Glider"))
  
plot2 <- ggplot(corr_df)+
  geom_smooth(aes(x = gpp, y = npp, color = method), method = "lm", se = FALSE, linewidth = 2)+
  geom_point(aes(x = gpp, y = npp, color = method), alpha = 0.8)+
  geom_line(aes(x = gpp, y = gpp, color = "GPP"), linewidth = 2)+
  theme_bw(base_size = 16)+
  scale_color_manual(values = c("CbPM Glider" = "steelblue1",
                                "CbPM Argo" = "steelblue4",
                                "VGPM Glider" = "tomato",
                                "VGPM Argo" = "tomato4",
                                "GPP" = "olivedrab4"),
                     name = "Method")+
  xlab("GPP from dial oxygen variability (mg C m-2 d-1)")+
  ylab("Depth integrated NPP (mg C m-2 d-1)")

ggsave("output/plots/pp_floats/pcorrelation_between_algorithms.png", dpi = 300, width = 40, height = 15, units = "cm")

library(patchwork)


plot1+plot2+plot_layout(guides = 'collect')

ggsave("output/plots/pp_floats/comparison_plot.png", dpi = 300, width = 35, height = 20, units = "cm")

floats_daily_pp <- synth_float |> select(date, cbpm_smoothed, vgpm_smoothed)

monthly_pp <- synth_pp |>
  left_join(floats_daily_pp) |> 
  select(date, daily_npp, daily_vgpm, gpp, cbpm_smoothed, vgpm_smoothed) |> 
  filter(date >= date("2024-06-01") & date < date("2024-10-01")) |> 
  mutate(month = month(date)) |> 
  pivot_longer(daily_npp:vgpm_smoothed, names_to = "method", values_to = "pp")

ggplot(monthly_pp)+
  geom_boxplot(aes(x = as.factor(month), y = pp, color = method))+
  xlab("Month")+
  ylab("Primary productivity (mg C m-2 d-1)")+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c("daily_npp" = "steelblue1",
                                "cbpm_smoothed" = "steelblue4",
                                "daily_vgpm" = "tomato",
                                "vgpm_smoothed" = "tomato4",
                                "gpp" = "olivedrab4"),
                     labels = c("Float CbPM", "Glider CbPM", "Glider VGPM", "GPP", "Float VGPM"),
                     name = "")

ggsave("output/plots/pp_floats/pp_algorithms_comparison_boxplot.png", dpi = 300, width = 30, height = 18, units = "cm")

tot_dat <- left_join(select(synth_pp, -par), select(npp_dat, -par_insitu))

floats_august <- npp_df |> 
  filter(date >= date("2024-08-01") & date < date("2024-08-10"))|> 
  select(date,
         longitude = LONGITUDE,
         latitude = LATITUDE,
         depth = PRES_rounded,
         bbp700 = bbp700_cleaned,
         chla = CHLA_ADJUSTED,
         pp,
         npp_vgpm) |> 
  mutate(platform = "float")

gliders_august <-  tot_dat |> 
  filter(date >= date("2024-08-01") & date < date("2024-08-10"))|>
  select(date,
         longitude,
         latitude,
         depth,
         bbp700,
         chla = fluo,
         npp_cbpm = npp,
         npp_vgpm) |> 
  mutate(platform = "glider") |> 
  bind_rows(floats_august)

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

# ggplot(daily_surf)+
#   geom_smooth(aes(x = par_insitu_normalized, y = npp, color = "daily"), method = "lm", se = FALSE)+
#   geom_point(aes(x = par_insitu_normalized, y = npp, color = "daily"), alpha = 0.5)+
#   geom_smooth(aes(x = par_weekly_normalized, y = npp, color = "4 days avg"), method = "lm", se = FALSE)+
#   geom_point(aes(x = par_weekly_normalized, y = npp, color = "4 days avg"), alpha = 0.5)+
#   geom_smooth(aes(x = par_normalized, y = npp, color = "modis"), method = "lm", se = FALSE)+
#   geom_point(aes(x = par_normalized, y = npp, color = "modis"), alpha = 0.5)+
#   theme_bw()+
#   xlab("Integrated PAR (normalized)")+
#   ylab("NPP from CbPM (mg C m-2 d-1)")+
#   scale_color_discrete(name = "")
  
pca_df <- select(daily_surf, -datetime, -new_par, -weekly_par, -hour, -depth, -profile_index, -par_insitu_normalized, -par_normalized, -par_weekly_normalized, -new_par, -var_npp, -daily_npp, -var_vgpm, -daily_vgpm) |> 
  na.omit() |> 
  mutate(chla_bbp = fluo/bbp700)

float_pca_df <- npp_df |> 
  mutate(date = date(JULD),
         platform = "float") |> 
  filter(PRES_rounded < 5) |> 
  left_join(gpp) |> 
  filter(date >= date("2024-06-01") & date < date("2024-10-01")) |> 
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
         gpp,
         platform) |>
  group_by(date, platform) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() 
  na.omit() |> 
  mutate(chla_bbp = fluo/bbp700)

pca_df <- pca_df |> mutate(platform = "glider") |>
  select(-par_insitu,
         -modis_normalized) |> 
  bind_rows(float_pca_df) |> 
  select(-longitude, - latitude) |> na.omit()

ggplot(pca_df)+
  geom_point(aes(x = date, y = fluo, colour = platform))+
  theme_bw()
ggsave("output/plots/pp_floats/fluo_comparison.png", dpi = 300, width = 16, height = 9, units = "cm")

ggplot(pca_df)+
  geom_point(aes(x = date, y = bbp700, colour = platform))+
  theme_bw()
ggsave("output/plots/pp_floats/bbp_comparison.png", dpi = 300, width = 16, height = 9, units = "cm")

ggplot(pca_df)+
  geom_point(aes(x = date, y = chla_bbp, colour = platform))+
  theme_bw()
ggsave("output/plots/pp_floats/ratio_comparison.png", dpi = 300, width = 16, height = 9, units = "cm")

ggplot(pca_df)+
  geom_point(aes(x = date, y = par, colour = platform))+
  theme_bw()

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

pca_loadings <- pca_loadings |> 
  mutate(variable = case_when (variable == "npp" ~ "CbPM",
                               variable == "npp_vgpm" ~ "VGPM",
                               TRUE ~ variable))

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
    breaks = pretty(pca_scores$date),
    labels = scales::date_format("%b")  
  )+
  theme_minimal(base_size = 14)+
  xlab("PC1 (40%)")+
  ylab("PC2 (24%)")

ggsave("output/plots/pp_floats/pca_biplot.png", dpi = 300, width = 30, height = 18, units = "cm")

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

ggplot(arrange(pca_scores, date))+
  geom_path(aes(x = date, y = PC2))

location <- npp_dat |> 
  select(date, latitude, longitude) |>
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  left_join(pca_scores)

ggplot(location)+
  geom_point(aes(x = longitude, y = latitude, color = PC2))+
  scale_color_distiller(palette = "RdBu")



# some plot of gliders transect -------------------------------------------


plot_df <- dat |> 
  filter(npp < 150 & fluo < 4) |> 
  mutate(date = date(datetime)) |> 
  select(glider, date, depth, latitude, longitude, temperature, salinity, fluo, bbp700, npp, npp_vgpm) |> 
  group_by(glider, date, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  left_join(select(daily_zeu, date, zeu))

ggplot(plot_df)+
  geom_tile(aes(y = -depth, x = date, fill = fluo))+
  geom_path(aes(x = date, y = -zeu))+
  scale_fill_distiller(palette = "Greens", direction = 1)+
  ylim(-100,0)+
  theme_minimal()+
  facet_wrap(.~glider, ncol = 2)

ggplot(plot_df)+
  geom_tile(aes(y = -depth, x = date, fill = npp))+
  geom_path(aes(x = date, y = -zeu))+
  scale_fill_distiller(palette = "YlGnBu", direction = -1)+
  ylim(-100,0)+
  theme_minimal()+
  facet_wrap(.~glider, ncol = 2)

ggplot(plot_df)+
  geom_tile(aes(y = -depth, x = date, fill = npp))+
  geom_path(aes(x = date, y = -zeu))+
  scale_fill_distiller(palette = "Spectral", direction = -1)+
  ylim(-100,0)+
  theme_minimal()+
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(.~glider, ncol = 2)

plot_df_error <- plot_df |> 
  select(-glider) |> 
  group_by(date, depth) |> 
  summarise_all(sd, na.rm = TRUE) |> 
  ungroup()

ggplot(plot_df_error)+
  geom_tile(aes(y = -depth, x = date, fill = npp))+
  scale_fill_distiller(palette = "RdBu")+
  ylim(-100,0)+
  theme_minimal()

mean_npp <- plot_df |> 
  select(date, depth, npp) |> 
  group_by(date, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  select(date, depth, mean_npp = npp)

plot_df_diff <- plot_df |> left_join(mean_npp) |> 
  mutate(npp_diff = npp - mean_npp)

ggplot(plot_df_diff)+
  geom_tile(aes(y = -depth, x = date, fill = npp_diff))+
  scale_fill_distiller(palette = "RdBu")+
  ylim(-50,0)+
  theme_minimal()+
  facet_wrap(.~glider, ncol = 2)

plot_df_loc <- plot_df |> select(glider, date, npp, latitude, longitude) |> 
  group_by(date, glider) |> 
  summarise(npp_int = sum(npp, na.rm = TRUE),
                longitude = mean(longitude),
                latitude = mean(latitude)) |> 
  ungroup()

ggplot(plot_df_loc)+
  geom_point(aes(x = longitude, y = latitude, color = npp_int, shape = glider), size = 2)+
  scale_color_viridis_c()+
  coord_quickmap()+
  theme_minimal()


mean_npp <- plot_df_loc |> 
  select(date, npp_int) |> 
  group_by(date) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  select(date, mean_npp_int = npp_int)

mean_npp_transect <- plot_df |> 
  select(date, depth, npp, zeu, temperature, salinity, longitude, latitude) |>
  group_by(date, depth) |> 
  summarise_all(mean, na.rm = TRUE) |> 
  ungroup() |> 
  mutate(ct = gsw_CT_from_t(salinity, temperature, depth),
         rho = gsw_sigma0(salinity, ct)) |> 
  select(date, depth, npp, zeu, rho)

mld_date <- mean_npp_transect |> 
  group_by(date) |> 
  summarise(mld_depth = mld(rho, depth)) |> 
  ungroup()

mean_npp_transect <- left_join(mean_npp_transect, mld_date)

ggplot(mean_npp_transect) +
  geom_tile(aes(y = -depth, x = date, fill = npp, color = npp)) +
  geom_path(aes(x = date, y = -zeu), linetype = "dashed") +
  geom_text(aes(x = date, y = -zeu + 3, label = "Zeu"), data = filter(mean_npp_transect, date == max(date) - 5)) +
  geom_path(aes(x = date, y = -mld_depth), color = "white") +
  geom_text(aes(x = date, y = -(mld_depth + 15), label = "MLD"), color = "white", 
            data = filter(mean_npp_transect, date == max(date) - 5)) +
  scale_color_distiller(palette = "Spectral", direction = -1, guide = "none") +
  scale_fill_distiller(
    palette = "Spectral", direction = -1, 
    name = expression(atop("NPP from CbPM", "(mg C m"^{-3} ~ "d"^{-1} * ")"))
  ) +
  ylim(-100, 0) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing = unit(0.1, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(0.4, "cm"),
    
    # Reduce space between ticks and the plot
    axis.ticks.length = unit(-0.5, "cm"),  # Moves ticks closer to the plot
    axis.text.x = element_text(margin = margin(t = 2)),  # Adjust space for x-axis labels
    axis.text.y = element_text(margin = margin(r = 2))   # Adjust space for y-axis labels
  ) +
  ylab("Depth (m)") +
  xlab("Date")


ggsave("output/plots/pp_floats/npp_transect.png", dpi = 300, width = 22, height = 13, units = "cm")

loc_df_diff <- plot_df_loc |> left_join(mean_npp) |> 
  mutate(npp_diff_int = npp_int - mean_npp_int)

ggplot(loc_df_diff)+
  geom_point(aes(x = longitude, y = latitude, color = npp_diff_int, shape = glider), size = 2)+
  scale_color_distiller(palette = "RdBu")+
  coord_quickmap()+
  theme_minimal()


