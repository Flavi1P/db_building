library(tidyverse)
library(vroom)
# Load dataset
dat <- vroom("data/ALR/alr6_ramses_timed.csv")

smaller_trios <- filter(dat, depth > 0 & depth < 60)


samples <- smaller_trios |> select(nearest_trios_time, depth) |> unique()

test <- smaller_trios |> filter(nearest_trios_time %in% samples$nearest_trios_time[c(2307:2339)])

#test <- smaller_trios |> filter(nearest_trios_time %in% samples$nearest_trios_time[c(1052:1200)])

ggplot(test)+
  geom_path(aes(x = Wavelength, y = intesity, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()

ggsave("ALR6_trios_spectra_sample_before_correction.jpg", width = 30, height = 20, units = "cm", dpi = 300)

test_image <- test |> mutate(round_depth = round(depth)) |> 
  group_by(round_depth, Wavelength) |> 
  summarize(intesity = mean(intesity, na.rm = TRUE)) |> 
  mutate(depth_normalised_intensity = intesity / max(intesity)) |>
  ungroup()

ggplot(test_image)+
  geom_tile(aes(x = Wavelength, y = -round_depth, fill = intesity))+
  scale_fill_viridis_c()

ggplot(test_image)+
  geom_path(aes(x = Wavelength, y = intesity, color = round_depth, group = round_depth))+
  scale_color_viridis_c()

#Let's correct the intensity on the first observation considering an integration time of 4096ms

#Dark computation

test_sample <- test |> filter(depth == min(depth)) 

black_dat <- read_table("data/ALR/Radiometer_TRIOS_RAMSES/Calibration/Pre-cruise/Back_SAM_8820.dat", 
           col_names = FALSE, skip = 39)
colnames(black_dat) <- c("pixel", "b0", "b1", "status")

wavelnegth_lookup <- unique(test |> select(Wavelength)) |> mutate(n_wl = row_number())

black_dat <- black_dat |> 
  mutate(pixel = as.integer(pixel)) |> 
  inner_join(wavelnegth_lookup, by = c("pixel" = "n_wl")) |> 
  select(Wavelength, b0, b1) |> 
  na.omit() |> 
  mutate(b0 = as.double(b0),
         b1 = as.double(b1))



test_corrected <- test_sample |> 
  left_join(black_dat, by = "Wavelength") |> 
  mutate(Bn = b0 - 2048/8192 * b1,
         Cn = intesity - Bn)

ggplot(test_corrected)+
  geom_path(aes(x = Wavelength, y = Cn, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()

#Reminder of dark

Dark_pixel_start = 237
Dark_pixel_end = 254

dark_wl = wavelnegth_lookup %>%
  filter(n_wl >= Dark_pixel_start & n_wl <= Dark_pixel_end) |> 
  pull(Wavelength)

A = mean(test_corrected |> 
           filter(Wavelength %in% dark_wl) |> 
           pull(Cn), na.rm = TRUE)
test_corrected <- test_corrected |>
  mutate(Dn = Cn - A)

ggplot(test_corrected)+
  geom_path(aes(x = Wavelength, y = Dn, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()

#Normalised per max integration depth

test_corrected <- test_corrected |> 
  mutate(En = Dn * 8192/4096)

ggplot(test_corrected)+
  geom_path(aes(x = Wavelength, y = En, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()

#Normalised to physical unit

cal_aq <- read_table("data/ALR/Radiometer_TRIOS_RAMSES/Calibration/Pre-cruise/CalAQ_SAM_8820.dat", 
                     col_names = FALSE, col_types = cols(X2 = col_double()), 
                     skip = 39)
colnames(cal_aq) <- c("pixel", "Sn", "XX", "status")

cal_aq <- cal_aq |> 
  mutate(pixel = as.integer(pixel)) |> 
  inner_join(wavelnegth_lookup, by = c("pixel" = "n_wl")) |> 
  select(Wavelength, Sn) |> 
  na.omit()

test_corrected <- test_corrected |> 
  left_join(cal_aq, by = "Wavelength") |> 
  mutate(Fn = En / Sn)

ggplot(test_corrected)+
  geom_path(aes(x = Wavelength, y = Fn, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()


# Vhat GPT doing its thing ------------------------------------------------


calibrate_trios <- function(test, int_time_ms = 256,
                            black_file = "data/ALR/Radiometer_TRIOS_RAMSES/Calibration/Pre-cruise/Back_SAM_8820.dat",
                            calaq_file  = "data/ALR/Radiometer_TRIOS_RAMSES/Calibration/Pre-cruise/CalAQ_SAM_8820.dat",
                            plot = TRUE) {
  
  # --- 2. Read black (dark) calibration file ---
  black_dat <- read_table(black_file, col_names = FALSE, skip = 39)
  colnames(black_dat) <- c("pixel", "b0", "b1", "status")
  
  wavelength_lookup <- unique(test |> select(Wavelength)) |> 
    mutate(n_wl = row_number())
  
  black_dat <- black_dat |> 
    mutate(pixel = as.integer(pixel)) |> 
    inner_join(wavelength_lookup, by = c("pixel" = "n_wl")) |> 
    select(Wavelength, b0, b1) |> 
    mutate(across(c(b0, b1), as.double)) |> 
    drop_na()
  
  # --- 3. Correct for dark current ---
  test_corrected <- test |> 
    left_join(black_dat, by = "Wavelength") |> 
    mutate(
      Bn = b0 - int_time_ms/8192 * b1,
      Cn = intesity - Bn
    )
  
  # --- 4. Compute dark reference pixels ---
  dark_pixel_start <- 237
  dark_pixel_end   <- 254
  
  dark_wl <- wavelength_lookup %>%
    filter(n_wl >= dark_pixel_start & n_wl <= dark_pixel_end) |> 
    pull(Wavelength)
  
  A <- mean(
    test_corrected |> 
      filter(Wavelength %in% dark_wl) |> 
      pull(Cn),
    na.rm = TRUE
  )
  
  test_corrected <- test_corrected |> 
    mutate(Dn = Cn - A)
  
  # --- 5. Normalize by integration time (simulate response scaling) ---
  # The instrument has binary levels: 4, 8, 16, …, 4096 ms
  # We normalize relative to maximum integration (4096 ms)
  test_corrected <- test_corrected |> 
    mutate(En = Dn * 8192 / int_time_ms)
  
  # --- 6. Convert to physical units using CalAQ file ---
  cal_aq <- read_table(calaq_file, col_names = FALSE, 
                       col_types = cols(X2 = col_double()), skip = 39)
  colnames(cal_aq) <- c("pixel", "Sn", "XX", "status")
  
  cal_aq <- cal_aq |> 
    mutate(pixel = as.integer(pixel)) |> 
    inner_join(wavelength_lookup, by = c("pixel" = "n_wl")) |> 
    select(Wavelength, Sn) |> 
    drop_na()
  
  test_corrected <- test_corrected |> 
    left_join(cal_aq, by = "Wavelength") |> 
    mutate(Fn = En / Sn)
  
  # --- 7. Optional visualization ---
  if (plot) {
    p <- ggplot(test_corrected) +
      geom_path(aes(x = Wavelength, y = Fn, color = depth,
                    group = nearest_trios_time)) +
      scale_color_viridis_c() +
      labs(
        title = paste("TRIOS Spectra (Integration =", int_time_ms, "ms)"),
        y = "Corrected Radiance (Fn)",
        x = "Wavelength (nm)"
      ) +
      theme_minimal()
    print(p)
  }
  
  return(test_corrected)
}

# Example usage of the function

test_sample <- filter(test, depth == min(depth))

integration_levels <- c(4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)

# Apply the function for multiple integration times and visualize
results <- map_dfr(integration_levels, ~{
  tmp <- calibrate_trios(test, int_time_ms = .x, plot = FALSE)
  tmp |> mutate(int_time_ms = .x)
})

# Plot how spectral shape changes with integration time
ggplot(filter(results, int_time_ms == 512)) +
  geom_path(aes(x = Wavelength, y = Fn, color = factor(int_time_ms),
                group = int_time_ms)) +
  scale_color_viridis_d(name = "Integration (ms)") +
  labs(title = "Effect of Integration Time on TRIOS Spectra",
       x = "Wavelength (nm)",
       y = "Corrected Radiance (Fn)") +
  theme_minimal()

ggsave("ALR6_trios_spectra_integration_effect_fixed.jpg", width = 30, height = 20, units = "cm", dpi = 300)

#Apply 512 as default integration depth
int512 <- calibrate_trios(test, int_time_ms = 512, plot = TRUE)

dat489.757 <- results |> filter(Wavelength == 489.757)

ggplot(filter(dat489.757, int_time_ms > 128))+
  geom_point(aes(x = Fn, y = -depth, color = as.factor(int_time_ms)))+
  scale_color_viridis_d()+
  theme_minimal()+
  ggtitle("Ed490 from surf to 30m at 2pm UTC Integration time as color")+
  xlab("Ed 490")

ggsave("ALR6_trios_Ed490_integration_effect.jpg", width = 30, height = 20, units = "cm", dpi = 300)

surf_value = dat489.757 |> 
  filter(depth < 10 & int_time_ms == 256) |> 
  pull(Fn) |> 
  mean(na.rm = TRUE)

surf_value = 13

kd_fit <- dat489.757 |> 
  select(depth) |> 
  unique() |> 
  mutate(Ed_synth = surf_value * exp(-0.08 * depth))

full_dat <- dat489.757 |> 
  left_join(kd_fit, by = "depth")

ggplot(filter(full_dat, int_time_ms > 128))+
  geom_point(aes(x = Fn, y = -depth, color = as.factor(int_time_ms)))+
  geom_point(aes(y = -depth, x = Ed_synth), color = "red", data = kd_fit)+
  scale_color_viridis_d()+
  theme_minimal()+
  ggtitle("Ed490 from surf to 30m at 2pm UTC Integration time as color")+
  xlab("Ed 490")

ggsave("ALR6_trios_Ed490_integration_effect_with_guessed_kd.jpg", width = 30, height = 20, units = "cm", dpi = 300)

full_dat <- full_dat |> 
  mutate(diff_synth = Fn - Ed_synth) |> 
  group_by(depth) |> 
  filter(diff_synth == min(abs(diff_synth), na.rm = TRUE)) |> 
  ungroup()

ggplot(filter(full_dat, log(Fn) ))+
  geom_point(aes(x = Fn, y = -depth, color = as.factor(int_time_ms)))+
  geom_path(aes(y = -depth, x = Ed_synth), color = "red", data = kd_fit)+
  scale_color_viridis_d()+
  theme_minimal()+
  ggtitle("Ed490 from surf to 30m at 2pm UTC Integration time as color")+
  xlab("Ed 490")

ggsave("ALR6_trios_Ed490_integration_effect_with_guessed_kd_best_fit.jpg", width = 30, height = 20, units = "cm", dpi = 300)

int_time_results <- full_dat |> 
  select(depth, int_time_ms) |> 
  rename("true_integration" = int_time_ms)

results <- results|> 
  left_join(int_time_results) |> 
  filter(int_time_ms == true_integration)



ggplot(filter(results, Wavelength < 800 & Wavelength > 340))+
  geom_path(aes(x = Wavelength, y = Fn, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()+
  theme_minimal()+
  ylab("Ed")

ggsave("ALR6_trios_spectra_after_correction.jpg", width = 30, height = 20, units = "cm", dpi = 300)
#compute KD lambda on the results dataframe
ggplot(filter(results, Wavelength < 700))+
  geom_path(aes(x = log(Fn), y = - depth, color = Wavelength, group = Wavelength))+
  scale_color_viridis_c()+
  theme_minimal()+
  ylim(-15, 0)


ggplot(test)+
  geom_point(aes(x = nearest_trios_time, y = -depth))
