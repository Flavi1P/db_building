library(tidyverse)
library(zoo)
library(gsw)
library(castr)
library(slider)
library(vroom)
library(lubridate)

source("script/cbpm_r/cbpm_argo.r")
source("script/cbpm_r/despike.r")

# --- Load and preprocess float data -----------------------------------------
dat <- read_csv("argo_3901586_interp.csv")

dat <- dat |> 
  mutate(sa = gsw_SA_from_SP(sal, depth, lon ,lat),
         ct = gsw_CT_from_t(sa, temp, depth),
         sigma0 = gsw_sigma0(sa, ct))

# Calculate Zeu
zeu_calc <- function(par, depth){
  if(all(is.na(par)) | all(par < 1)){
    return(NA)
  } else {
    par_surf = par[which(depth == 1)]
    threshold = par_surf * 0.01
    return(depth[which(par <= threshold)[1]])
  }
}

prof_dat <- dat |>
  group_by(prof_number, date) |>
  arrange(depth) |> 
  summarise(
    MLD = mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03, default.depth = 300),
    chl_cline = clined(chla, depth),
    zeu = zeu_calc(iPAR, depth),
    .groups = "drop"
  ) |>
  filter(!is.na(zeu)) |> 
  mutate(
    mld_s = pmax(smooth.spline(date, MLD, spar = 0.7)$y, 2),
    zeu_s = pmax(smooth.spline(date, zeu, spar = 0.6)$y, 2)
  )

dat_combined <- left_join(dat, prof_dat, by = c("prof_number","date"))

# --- NPP (CbPM) -------------------------------------------------------------
npp <- vroom("data/argo_pp_estimations_floats_09_2025.csv") |> 
  filter(PLATFORM_NUMBER == 3901586) |> 
  mutate(date = date(JULD)) |> 
  janitor::clean_names() |> 
  select(date, day, month, year, satellite_par) |> 
  unique()

dat_all_pp <- left_join(dat_combined, npp)

dat_all_pp <- dat_all_pp |>
  group_by(prof_number) |> 
  mutate(
    bbp_baseline = despike(bbp, 10)$baseline,
    bbp_baseline = ifelse(is.na(bbp_baseline), bbp, bbp_baseline),
    carbon = bbp_baseline / (470/440) * 12128 + 0.59,
    cbpm_npp = cbpm_argo(chla, carbon, unique(satellite_par),
                         unique(na.omit(year)), unique(na.omit(month)), 
                         unique(na.omit(day)), unique(na.omit(lat)))$pp,   # mg C m-3 d-1
    cbpm_mu  = cbpm_argo(chla, carbon, unique(satellite_par),
                         unique(na.omit(year)), unique(na.omit(month)), 
                         unique(na.omit(day)), unique(na.omit(lat)))$mu_z
  ) |> ungroup()

# --- Integrate profiles -----------------------------------------------------
integrated_val <- dat_all_pp |> 
  filter(!is.na(zeu_s) & !is.na(mld_s)) |> 
  group_by(prof_number, date) |> 
  summarise(
    chl_mld     = integrate(chla, depth, from = 0, to = abs(unique(MLD))),
    chl_zeu     = integrate(chla, depth, from = 0, to = abs(unique(zeu))),
    npp_mld     = integrate(cbpm_npp, depth, from = 0, to = abs(unique(MLD))),
    npp_zeu     = integrate(cbpm_npp, depth, from = 0, to = abs(unique(zeu))),
    carbon_mld  = integrate(carbon, depth, from = 0, to = abs(unique(MLD))),
    carbon_zeu  = integrate(carbon, depth, from = 0, to = abs(unique(zeu))),
    nitrate_mld = integrate(nitrate, depth, from = 0, to = abs(unique(MLD))),
    nitrate_zeu = integrate(nitrate, depth, from = 0, to = abs(unique(zeu))),
    .groups = "drop"
  ) |> 
  arrange(date) |>
  mutate(
    across(c(chl_mld, chl_zeu, npp_mld, npp_zeu, carbon_mld, carbon_zeu,
             nitrate_mld, nitrate_zeu),
           list(s = ~slide_index_dbl(.x, date, mean, .before=15, .after=15, na.rm=TRUE),
                std = ~slide_index_dbl(.x, date, sd,   .before=15, .after=15, na.rm=TRUE)),
           .names = "{.col}_{.fn}")
  )

dat_combined <- left_join(dat_combined, integrated_val, by = c("prof_number","date"))

# --- Net growth rates and NCP -----------------------------------------------
prof_dat1 <- prof_dat |>
  left_join(integrated_val, by = c("prof_number","date")) |> 
  mutate(yday = yday(date)) |> 
  arrange(date) |> 
  mutate(
    p         = if_else(mld_s > zeu_s, chl_mld_s, chl_zeu_s),
    nitrate_p = if_else(mld_s > zeu_s, nitrate_mld_s, nitrate_zeu_s),
    p_diff        = p - lag(p),
    nitrate_p_diff = nitrate_p - lag(nitrate_p),
    date_diff     = yday - lag(yday),
    rp  = (1/p) * (p_diff/date_diff),  # relative phytoplankton growth
    ncp = 0.0796 * (nitrate_p_diff/date_diff)  # mg C m-2 d-1
  ) |> 
  mutate(
    rp_s  = slide_index_dbl(rp, date, mean, .before=10, .after=10, na.rm=TRUE),
    rp_std= slide_index_dbl(rp, date, sd,   .before=10, .after=10, na.rm=TRUE),
    ncp_s = slide_index_dbl(ncp, date, mean, .before=10, .after=10, na.rm=TRUE),
    ncp_std=slide_index_dbl(ncp, date, sd,   .before=10, .after=10, na.rm=TRUE)
  )

# --- Comparison plot NPP vs NCP --------------------------------------------
ggplot(prof_dat1) +
  geom_line(aes(x=date, y=npp_mld_s, color="NPP MLD")) +
  geom_line(aes(x=date, y=npp_zeu_s, color="NPP Zeu")) +
  geom_line(aes(x=date, y=-ncp_s*50, color="NCP (Redfield)")) +
  ylab("mg C m-2 d-1") +
  theme_minimal()

ggplot(prof_dat1)+
  geom_line(aes(x = date, y = -ncp_s*12))
