# R script to compute phytoplankton productivity, net population growth rate, and division rate
# Based on methods from https://www.nature.com/articles/s41467-017-02143-6

# Required libraries
library(tidyverse)
library(lubridate)
library(zoo)
library(arrow)
library(oce)  # For seawater density calculations

argo2 <- read_parquet("data/argo_pq/icb_floats_table.parquet") |> mutate(depth = round(PRES))

# Function to compute mixed layer depth (MLD) based on potential density threshold
compute_mld <- function(depth, temp, sal, threshold = 0.03) {
  # Calculate potential density using seawater package
  density <- swSigmaTheta(sal, temp, depth, referencePressure = 10)
  # Find depth where density difference from 10m exceeds threshold
  ref_density <- density[which.min(abs(depth - 10))]
  delta_sigma <- abs(density - ref_density)
  mld_idx <- which(delta_sigma > threshold)[1]
  if (is.na(mld_idx)) return(max(depth))  # Default to max depth if no threshold met
  return(depth[mld_idx])
}

# Function to compute euphotic layer depth (He) based on PAR threshold
compute_euphotic_depth <- function(depth, ipar, par_threshold = 0.1) {
  # Fit exponential decay to iPAR to get attenuation coefficient (K)
  noon_ipar <- ipar  # Assume iPAR is at local noon
  surface_ipar <- noon_ipar[which.min(depth)]
  if (surface_ipar == 0) return(max(depth))  # Avoid division by zero
  ipar_ratio <- noon_ipar / surface_ipar
  fit <- glm(ipar ~ depth, subset = depth <= 50, family = poisson)
  K <- -coef(fit)[2]  # Diffuse attenuation coefficient
  # Daily-averaged PAR at surface (mol photons m^-2 d^-1)
  par_surface <- mean(noon_ipar) * 86400 / 1e6  # Convert to daily mol photons
  if (par_surface == 0) return(max(depth))
  # Euphotic depth where PAR drops to threshold
  He <- -1 / K * log(par_threshold / par_surface)
  return(min(He, max(depth)))  # Cap at max depth
}

# Function to convert b_bp to phytoplankton carbon biomass (P)
compute_phytoplankton_carbon <- function(b_bp, depth, conversion_factor = 11260, background_depth = 900) {
  # Estimate background b_bp below 900m
  background_bbp <- min(b_bp[depth > background_depth], na.rm = TRUE)
  # Convert b_bp to phytoplankton carbon (mgC m^-3)
  P <- (b_bp - background_bbp) * conversion_factor
  P[P < 0] <- 0  # Ensure non-negative biomass
  return(P)
}

# Function to compute phytoplankton division rate (mu) at depth z and time t
compute_mu <- function(chl, P, ipar, temp, depth, alpha_B = 6.4e-6, P_max_B_20 = 1.3e-3) {
  # Calculate P_max_B as a function of temperature
  ml_temp <- mean(temp[depth <= max(depth)/4], na.rm = TRUE)  # Approx mixed layer temp
  P_max_B <- P_max_B_20 * 1.065^(ml_temp - 20)
  # Compute mu(z,t) using the photosynthesis model
  mu <- (chl / P) * P_max_B * (1 - exp(-alpha_B * ipar / P_max_B))
  mu[is.na(mu) | P == 0] <- 0  # Handle division by zero
  return(mu)
}

# Function to compute depth-integrated primary production and division rate
compute_production <- function(P, mu, depth, mld, He) {
  # Integrate over the productive layer (max of MLD or euphotic depth)
  L <- max(mld, He, na.rm = TRUE)
  depth_subset <- depth <= 0 & depth >= -L
  if (sum(depth_subset) < 2) return(list(PP = 0, mu_p = 0))
  # Depth-integrated primary production (mgC m^-2 d^-1)
  PP <- trapz(depth[depth_subset], mu[depth_subset] * P[depth_subset]) * 86400
  # Depth-integrated phytoplankton biomass
  P_int <- trapz(depth[depth_subset], P[depth_subset])
  # Population division rate
  mu_p <- PP / P_int
  return(list(PP = PP, mu_p = mu_p))
}

# Function to compute net population growth rate (r_p)
compute_rp <- function(P, depth, mld, time, He) {
  # Interpolate to 3-day intervals
  time_seq <- seq(min(time, na.rm = TRUE), max(time, na.rm = TRUE), by = 3*86400)
  P_ml <- numeric(length(time_seq))
  P_int <- numeric(length(time_seq))
  H <- numeric(length(time_seq))
  
  for (i in 1:length(time_seq)) {
    idx <- which.min(abs(time - time_seq[i]))
    profile_data <- data[data$time == time[idx], ]
    depth_profile <- profile_data$depth
    P_profile <- P[data$time == time[idx]]
    depth_subset <- depth_profile <= 0 & depth_profile >= -mld[idx] & !is.na(P_profile)
    
    if (sum(depth_subset) < 2) {
      P_ml[i] <- NA
      P_int[i] <- NA
    } else {
      P_ml[i] <- mean(P_profile[depth_subset], na.rm = TRUE)
      P_int[i] <- trapz(depth_profile[depth_subset], P_profile[depth_subset])
    }
    H[i] <- mld[idx]
  }
  # Apply 24-day running average
  P_ml_smooth <- zoo::rollmean(P_ml, k = 8, fill = NA, align = "center")
  P_int_smooth <- zoo::rollmean(P_int, k = 8, fill = NA, align = "center")
  H_smooth <- zoo::rollmean(H, k = 8, fill = NA, align = "center")
  # Compute derivatives
  dP_ml_dt <- diff(P_ml_smooth) / (3 * 86400)
  dP_int_dt <- diff(P_int_smooth) / (3 * 86400)
  dH_dt <- diff(H_smooth) / (3 * 86400)
  # Compute r_p based on mixed layer dynamics
  rp <- numeric(length(dP_ml_dt))
  for (i in 1:length(rp)) {
    if (is.na(H_smooth[i]) || is.na(P_ml_smooth[i]) || is.na(P_int_smooth[i])) {
      rp[i] <- NA
    } else if (dH_dt[i] > 0) {  # Mixed layer deepening
      rp[i] <- dP_int_dt[i] / P_int_smooth[i]
    } else {  # Mixed layer shoaling
      rp[i] <- dP_ml_dt[i] / P_ml_smooth[i]
    }
  }
  return(list(rp = rp, time_seq = time_seq[-1]))
}

# Main function to process BGC-Argo float data
process_bgc_argo <- function(data) {
  # Assume data is a data frame with columns: time, depth, temp, sal, chl, b_bp, ipar
  # Time in POSIXct, depth in meters (negative downwards), others in appropriate units
  
  # Initialize output
  results <- list()
  
  # Compute mixed layer depth and euphotic depth for each profile
  profiles <- unique(data$time)
  mld <- numeric(length(profiles))
  He <- numeric(length(profiles))
  for (i in 1:length(profiles)) {
    idx <- data$time == profiles[i]
    mld[i] <- compute_mld(data$depth[idx], data$temp[idx], data$sal[idx])
    He[i] <- compute_euphotic_depth(data$depth[idx], data$ipar[idx])
  }
  
  # Convert b_bp to phytoplankton carbon
  P <- compute_phytoplankton_carbon(data$b_bp, data$depth)
  
  # Compute division rate mu(z,t)
  mu <- compute_mu(data$chl, P, data$ipar, data$temp, data$depth)
  
  # Compute primary production and population division rate
  PP <- numeric(length(profiles))
  mu_p <- numeric(length(profiles))
  for (i in 1:length(profiles)) {
    idx <- data$time == profiles[i]
    prod <- compute_production(P[idx], mu[idx], data$depth[idx], mld[i], He[i])
    PP[i] <- prod$PP
    mu_p[i] <- prod$mu_p
  }
  
  # Compute net population growth rate
  rp_results <- compute_rp(P, data$depth, mld, profiles, He)
  
  # Compile results
  results <- data.frame(
    time = profiles,
    mld = mld,
    He = He,
    PP = PP,
    mu_p = mu_p
  ) %>%
    left_join(data.frame(time = rp_results$time_seq, rp = rp_results$rp), by = "time")
  
  return(results)
}

df <-  vroom::vroom('data/argo_pp_north_atlantic_estimations_floats.csv') |> 
  select(time = JULD, depth = PRES_rounded, temp = TEMP, sal = PSAL, chl = CHLA_ADJUSTED, b_bp = BBP700_ADJUSTED, ipar = DOWNWELLING_PAR) |> 
  na.omit() |> 
  mutate(ipar = case_when(ipar <= 0 ~ 0,
                          TRUE ~ ipar))

results <- process_bgc_argo(df)
