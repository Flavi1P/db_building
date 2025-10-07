library(tidyverse)
library(zoo)
library(gsw)
library(castr)
library(slider)
library(patchwork)
source("script/cbpm_r/cbpm_argo.r")
source("script/cbpm_r/despike.r")

dat <- read_csv("argo_7902223_interp.csv")

dat <- dat |> mutate(sa = gsw_SA_from_SP(sal, depth, lon ,lat), ct = gsw_CT_from_t(sa, temp, depth), sigma0 = gsw_sigma0(sa, ct))

# NPP data ---------------------------------------------------------------

library(vroom)

npp <- vroom("data/argo_pp_estimations_floats_09_2025.csv") |> filter(PLATFORM_NUMBER == 7902223)

npp <- npp |> mutate(date = lubridate::date(JULD)) |> 
  janitor::clean_names() |> 
  select(date, day, month, year, satellite_par) |> 
  unique()


dat_all_pp <- left_join(dat, npp)


dat_all_pp <- dat_all_pp |>
  filter(!is.na(satellite_par)) |> 
  group_by(prof_number) |> 
  mutate(
    bbp_baseline = despike(bbp, 10)$baseline,
    bbp_baseline = case_when(is.na(bbp_baseline) ~ bbp,
                             TRUE ~ bbp_baseline),
    bbp_dark = mean(bbp_baseline[which(depth >= 900)], na.rm = TRUE),
    bb_baseline = bbp_baseline - bbp_dark,
    carbon = bbp_baseline / (470/440) * 12128 + 0.59,
    izeu_cbpm = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$mzeu,
    cbpm_npp = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$pp,
    cbpm_mu = cbpm_argo(chla, carbon, unique(satellite_par), unique(na.omit(year)), unique(na.omit(month)), unique(na.omit(day)), unique(na.omit(lat)))$mu_z)


# Prof dat ---------------------------------------------------------------


zeu_calc <- function(par, depth){
  # find the depth where par first exceeds the threshold
  if(all(is.na(par)) | all(par < 1)){
    return(NA)
  } else {
    par_surf = par[which(depth == 1)]
    threshold = par_surf * 0.01
    return(depth[which(par <= threshold)[1]])
  }
}


prod_depth <- function(npp, depth){
  if(all(is.na(npp))){
    return(NA)
  } else {
    return(depth[which(npp < 0.01)[1]])
  }
}
prof_dat <- dat_all_pp |>
  group_by(prof_number, date) |>
  arrange(depth) |> 
  summarise(
    MLD = mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03, default.depth = 300),
    chl_cline = clined(chla, depth),
    zeu_cbpm = mean(izeu_cbpm, na.rm = TRUE),
    zeu = zeu_calc(iPAR, depth),
    z_prod = prod_depth(cbpm_npp, depth)
  ) |>
  ungroup() |>
  filter(!is.na(zeu) & !is.na(MLD)) |> 
  mutate(
    mld_s = pmax(smooth.spline(date, MLD, spar = 0.7)$y, 2),
    zeu_s = pmax(smooth.spline(date, zeu, spar = 0.6)$y, 2),
    z_prod_s = pmax(smooth.spline(date, z_prod, spar = 0.6)$y, 2)
  ) |>
  mutate(
    delta_MLD = mld_s - lag(mld_s),
    delta_zeu = zeu_s - lag(zeu_s)
  ) |>
  replace_na(list(delta_MLD = 0)) |> 
  mutate(prev_mld = lag(mld_s),
         next_mld = lead(mld_s),
         prev_zeu = lag(zeu),
         next_zeu = lead(zeu))


p0 <- ggplot(prof_dat)+
  geom_point(aes(x = date, y = -MLD, colour = "MLD"))+
  geom_line(aes(x = date, y = -mld_s, colour = "MLD"))+
  geom_point(aes(x = date, y = -zeu, colour = "Zeu"))+
  geom_line(aes(x = date, y = -zeu_s, colour = "Zeu"))+
  geom_line(aes(x = date, y = -z_prod_s, color = "Prod layer (> 0.001 mg.m2.d-1)"))+
  ylim(-800,0)

p0
dat_all_pp <- left_join(dat_all_pp, prof_dat)




# integrating ------------------------------------------------------------

integrated_val <- dat_all_pp |> 
  filter(!is.na(zeu_s) & !is.na(mld_s)) |> 
  group_by(prof_number, date) |> 
  summarise(chl_mld = integrate(chla, depth, from =0, to = abs(unique(MLD))),
            chl_mld_deepest = integrate(chla, depth, from =0, to = abs(unique(MLD))),
            chl_zeu = integrate(chla, depth, from =0, to = abs(unique(zeu))),
            npp_mld = integrate(cbpm_npp, depth, from = 0, to = abs(unique(MLD))),
            npp_zeu = integrate(cbpm_npp, depth, from = 0, to = abs(unique(zeu))),
            carbon_mld = integrate(carbon, depth, from = 0, to = abs(unique(MLD))),
            carbon_zeu = integrate(carbon, depth, from = 0, to = abs(unique(zeu))),
            nitrate_mld = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(MLD))),
            nitrate_prev = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(prev_mld))),
            nitrate_next = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(next_mld))),
            nitrate_zeu = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(zeu))),
            nitrate_prev_zeu = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(prev_zeu))),
            nitrate_next_zeu = integrate(nitrate * 1000, depth, from = 0, to = abs(unique(next_zeu))),
            nitrate_prod = integrate(nitrate * 1000, depth, from = 0, to = 200),
            carbon_prev = integrate(carbon, depth, from = 0, to = abs(unique(prev_mld))),
            carbon_next = integrate(carbon, depth, from = 0, to = abs(unique(next_mld))),
            carbon_prev_zeu = integrate(carbon, depth, from = 0, to = abs(unique(prev_zeu))),
            carbon_next_zeu = integrate(carbon, depth, from = 0, to = abs(unique(next_zeu))))

integrated_val <- integrated_val |>
  ungroup() |> 
  arrange(date) |>
  mutate(
    chl_mld_s = slide_index_dbl(chl_mld, date, mean, 
                                .before = 15, .after = 15, na.rm = TRUE),
    chl_zeu_s = slide_index_dbl(chl_zeu, date, mean, 
                                .before = 15, .after = 15, na.rm = TRUE),
    chl_mld_std = slide_index_dbl(chl_mld, date, sd, 
                                  .before = 15, .after = 15, na.rm = TRUE),
    chl_zeu_std = slide_index_dbl(chl_zeu, date, sd, 
                                  .before = 15, .after = 15, na.rm = TRUE),
    npp_zeu_s = slide_index_dbl(npp_zeu, date, mean, 
                                .before = 15, .after = 15, na.rm = TRUE),
    npp_mld_s = slide_index_dbl(npp_mld, date, mean, 
                                .before = 15, .after = 15, na.rm = TRUE),
    npp_mld_std = slide_index_dbl(npp_mld, date, sd, 
                                  .before = 15, .after = 15, na.rm = TRUE),
    npp_zeu_std = slide_index_dbl(npp_zeu, date, sd, 
                                  .before = 15, .after = 15, na.rm = TRUE),
    carbon_zeu_s = slide_index_dbl(carbon_zeu, date, mean, 
                                   .before = 15, .after = 15, na.rm = TRUE),
    carbon_mld_s = slide_index_dbl(carbon_mld, date, mean, 
                                   .before = 15, .after = 15, na.rm = TRUE),
    carbon_zeu_std = slide_index_dbl(carbon_zeu, date, sd, 
                                     .before = 15, .after = 15, na.rm = TRUE),
    carbon_mld_std = slide_index_dbl(carbon_mld, date, sd, 
                                     .before = 15, .after = 15, na.rm = TRUE),
    nitrate_zeu_s = slide_index_dbl(nitrate_zeu, date, mean, 
                                    .before = 15, .after = 15, na.rm = TRUE),
    nitrate_mld_s = slide_index_dbl(nitrate_mld, date, mean, 
                                    .before = 15, .after = 15, na.rm = TRUE),
    nitrate_zeu_std = slide_index_dbl(nitrate_zeu, date, sd, 
                                      .before = 15, .after = 15, na.rm = TRUE),
    nitrate_mld_std = slide_index_dbl(nitrate_mld, date, sd, 
                                      .before = 15, .after = 15, na.rm = TRUE),
    nitrate_prod_s = slide_index_dbl(nitrate_prod, date, mean, 
                                    .before = 15, .after = 15, na.rm = TRUE),
    nitrate_prod_std = slide_index_dbl(nitrate_prod, date, sd, 
                                       .before = 15, .after = 15, na.rm = TRUE),
    nitrate_prev_s = slide_index_dbl(nitrate_prev, date, mean, 
                                    .before = 15, .after = 15, na.rm = TRUE),
    nitrate_next_s = slide_index_dbl(nitrate_next, date, mean,
                                     .before = 15, .after = 15, na.rm = TRUE),
    nitrate_prev_zeu_s = slide_index_dbl(nitrate_prev_zeu, date, mean, 
                                        .before = 15, .after = 15, na.rm = TRUE),
    nitrate_next_zeu_s = slide_index_dbl(nitrate_next_zeu, date, mean,
                                         .before = 15, .after = 15, na.rm = TRUE),
    carbon_prev_s = slide_index_dbl(carbon_prev, date, mean, 
                                    .before = 15, .after = 15, na.rm = TRUE),
    carbon_next_s = slide_index_dbl(carbon_next, date, mean,
                                     .before = 15, .after = 15, na.rm = TRUE),
    carbon_prev_zeu_s = slide_index_dbl(carbon_prev_zeu, date, mean, 
                                        .before = 15, .after = 15, na.rm = TRUE),
    carbon_next_zeu_s = slide_index_dbl(carbon_next_zeu, date, mean,
                                         .before = 15, .after = 15, na.rm = TRUE)
  )



p1 <- ggplot(na.omit(integrated_val))+
  geom_line(aes(x = date, y = chl_mld_s,  color = "MLD"))+
  geom_errorbar(aes(x = date, ymin = chl_mld_s - chl_mld_std, ymax = chl_mld_s + chl_mld_std, color = "MLD"), width = 0.5)+
  geom_line(aes(x = date, y = chl_zeu_s, color = "Zeu"))+
  geom_errorbar(aes(x = date, ymin = chl_zeu_s - chl_zeu_std, ymax = chl_zeu_s + chl_zeu_std, color = "Zeu"), width = 0.5)+
  ylab("Integrated Chl (mg.m-2)")+
  scale_color_brewer(palette = "Set2")

p2 <- ggplot(na.omit(integrated_val))+
  geom_line(aes(x = date, y = nitrate_mld_s, color = "MLD"))+
  geom_errorbar(aes(x = date, ymin = nitrate_mld_s - nitrate_mld_std, ymax = nitrate_mld_s + nitrate_mld_std, color = "MLD"), width = 0.5)+
  geom_line(aes(x = date, y = nitrate_zeu_s, color = "Zeu"))+
  geom_errorbar(aes(x = date, ymin = nitrate_zeu_s - nitrate_zeu_std, ymax = nitrate_zeu_s + nitrate_zeu_std, color = "Zeu"))+
  ylab("Integrated nitrates (umol.m-2)")+
  scale_color_brewer(palette = "Set2")

p3 <- ggplot(na.omit(integrated_val))+
  geom_line(aes(x = date, y = carbon_mld_s, color = "MLD"))+
  geom_errorbar(aes(x = date, ymin = carbon_mld_s - npp_mld_std, ymax = carbon_mld_s + carbon_mld_std, color = "MLD"), width = 0.5)+
  geom_line(aes(x = date, y = npp_zeu_s, color = "Zeu"))+
  geom_errorbar(aes(x = date, ymin = npp_zeu_s - carbon_zeu_std, ymax = npp_zeu_s + carbon_zeu_std, color = "Zeu"), width = 0.5)+
  ylab("Phyto Carbon stock (mg.m-2)")+
  scale_color_brewer(palette = "Set2")

p4 <- ggplot(na.omit(integrated_val))+
  geom_line(aes(x = date, y = npp_mld_s, color = "MLD"))+
  geom_errorbar(aes(x = date, ymin = npp_mld_s - npp_mld_std, ymax = npp_mld_s + npp_mld_std, , color = "MLD"), width = 0.5)+
  geom_line(aes(x = date, y = npp_zeu_s, , color = "Zeu"))+
  geom_errorbar(aes(x = date, ymin = npp_zeu_s - npp_zeu_std, ymax = npp_zeu_s + npp_zeu_std, color = "Zeu"), width = 0.5)+
  ylab("Integrated NPP mg C m-2 d-1 (CbPM)")+
  scale_color_brewer(palette = "Set2")


p1 + p2 + p3 + p4
#ggsave("output/plots/pp_floats/float_3901586_integrated_vals.png", dpi = 300, width = 30, height = 20, units = 'cm')

dat_combined <- left_join(dat_all_pp, integrated_val)


# net growth rate calculation ---------------------------------------------

prof_dat1 <- prof_dat |>
  left_join(integrated_val) |> 
  mutate(yday = lubridate::yday(date)) |> 
  arrange(date)

prof_dat1 <- prof_dat1 |> 
  filter(!is.na(zeu)) |> 
  mutate(shallow = case_when(mld_s < zeu ~  "MLD",
                             mld_s > zeu ~ "zeu"),
         MLD_diff = mld_s - lag(mld_s),
         zeu_diff = zeu - lag(zeu),
         chl_mld_diff = chl_mld_s - lag(chl_mld_s),
         chl_zeu_diff = chl_zeu_s - lag(chl_zeu_s),
         #carbon_mld_diff = carbon_mld_s - lag(carbon_mld_s),
         carbon_zeu_diff = carbon_zeu_s - lag(carbon_zeu_s),
         date_diff = yday - lag(yday),
         p = case_when(MLD > zeu ~ carbon_mld_s,
                       MLD < zeu ~ carbon_zeu_s),
         p_diff = p - lag(p),
         nitrate_mld_diff = case_when(MLD > zeu & prev_mld > mld_s ~ nitrate_prev_s - lag(nitrate_mld_s),
                                      MLD > zeu & prev_mld < mld_s ~ nitrate_mld_s - lag(nitrate_next_s),
                                      MLD < zeu & prev_zeu > zeu ~ nitrate_prev_zeu_s - lag(nitrate_zeu_s),
                                      MLD < zeu & prev_zeu < zeu ~ nitrate_zeu_s - lag(nitrate_next_zeu_s),
                                      is.na(prev_mld)~NA),
         carbon_mld_diff = case_when(MLD > zeu & prev_mld > mld_s ~ carbon_prev_s - lag(carbon_mld_s),
                                      MLD > zeu & prev_mld < mld_s ~ carbon_mld_s - lag(carbon_next_s),
                                      MLD < zeu & prev_zeu > zeu ~ carbon_zeu_s - lag(carbon_zeu_s),
                                      MLD < zeu & prev_zeu < zeu ~ carbon_zeu_s - lag(carbon_next_zeu_s),
                                      is.na(prev_mld)~NA),
         deepest_limit = case_when(MLD > zeu ~ "MLD",
                                       MLD < zeu ~ "Zeu"),
         nitrate_zeu_diff = nitrate_zeu_s - lag(nitrate_zeu_s),
         nitrate_prod_diff = nitrate_prod_s - lag(nitrate_prod_s),
         nitrate_mld_diff_s = slide_index_dbl(nitrate_mld_diff, date, mean, .before = 15, .after = 15, na.rm = TRUE),
         nitrate_prod_diff_s = slide_index_dbl(nitrate_prod_diff, date, mean, .before = 15, .after = 15, na.rm = TRUE),
         nitrate_p = case_when(mld_s > zeu ~ nitrate_mld_s,
                               mld_s < zeu ~ nitrate_zeu_s),
         nitrate_p_diff = nitrate_p - lag(nitrate_p),
         mol_c_diff = nitrate_mld_diff_s * 10e-6 * 6.6,
         mg_c_diff = mol_c_diff * 12 * 1000)

prof_dat1 <- prof_dat1 |> mutate(rate_mld = (1/mld_s) * (MLD_diff/date_diff))


ggplot(prof_dat1)+
  geom_point(aes(x= date, y = rate_mld, color = MLD_diff))+
  scale_color_distiller(palette = "RdBu")+
  theme_dark()

prof_dat1 <- prof_dat1 |> 
  mutate(rp = case_when(MLD > zeu & rate_mld > 0 ~ (1/p) * (p_diff/date_diff),
                        MLD > zeu & rate_mld < 0 ~ (1/carbon_mld_s) * (chl_mld_diff/date_diff),
                        MLD < zeu ~ (1/p) * (p_diff / date_diff)),
         rp_simple = (1/carbon_mld_diff) * (carbon_mld_diff/date_diff), 
         rp_simple_s = slide_index_dbl(rp_simple, date, mean, .before = 10, .after = 10, na.rm = TRUE),
         rp_simple_std = slide_index_dbl(rp_simple, date, sd, .before = 10, .after = 10, na.rm = TRUE),
         ncp = case_when(mld_s > zeu & rate_mld > 0 ~ 0.0796 * (nitrate_p_diff/date_diff),
                         mld_s > zeu & rate_mld < 0 ~ 0.0796 * (nitrate_mld_diff/date_diff),
                         mld_s < zeu ~ 0.0796 * (nitrate_p_diff / date_diff)),
         ncp_mld = mol_c_diff/date_diff,
         ncp_deepest = case_when(mld_s < zeu ~ 0.0796 * (nitrate_zeu_diff/date_diff),
                                 mld_s > zeu ~ 0.0796 * (nitrate_mld_diff/date_diff)),
         rp_s = slide_index_dbl(rp, date, mean, 
                                .before = 10, .after = 10, na.rm = TRUE),
         rp_std = slide_index_dbl(rp, date, sd, 
                                  .before = 10, .after = 10, na.rm = TRUE),
         ncp_s = slide_index_dbl(ncp, date, mean, .before = 10, .after = 10, na.rm = TRUE),
         ncp_std = slide_index_dbl(ncp, date, sd, .before = 10, .after = 10, na.rm = TRUE),
         ncp_mld_s = slide_index_dbl(ncp_mld, date, mean, .before = 10, .after = 10, na.rm = TRUE))


ggplot(prof_dat1)+
  geom_point(aes(x = date, y = nitrate_mld_diff, color = deepest_limit))+
  geom_line(aes(x = date, y = nitrate_mld_diff_s))


ggplot(filter(prof_dat1))+
  geom_line(aes(x = date, y = rp_s))
  geom_errorbar(aes(x = date, ymin = rp - rp_std, ymax = rp + rp_std, color = deepest_limit), width = 0.5)

#ggsave("output/plots/pp_floats/net_growth_rate.png", dpi = 300, width = 30, height = 20, units = 'cm')

ggplot(prof_dat1)+
  geom_line(aes(x = date, y = -ncp_mld_s, color = "mld"))

ggplot(prof_dat1)+
  geom_line(aes(x = date, y = npp_mld_s, color = "NPP MLD"))+
  geom_line(aes(x = date, y = npp_zeu_s, color = "NPP Zeu"))+
  geom_line(aes(x = date, y = -ncp_mld_s*100, color = "NCP"))+
  theme_bw(base_size = 14)+
  ylab("productivity (mg C m-2 d-1)")

#ggsave("output/plots/pp_floats/productivity_tot.png", dpi = 300, width = 30, height = 20, units = 'cm')

# regularisation ----------------------------------------------------------

# sequence of every 3rd day of the year
yday_list <- seq(3, 363, 3)
date_list <- lubridate::parse_date_time(yday_list, order = "j") - (3600 * 24 * 366)


MLD_interp = approx(prof_dat1$yday, prof_dat1$MLD, xout = yday_list, rule = 2, na.rm = FALSE)$y
zeu_interp = approx(prof_dat1$yday, prof_dat1$zeu, xout = yday_list, rule = 2, na.rm = FALSE)$y
chl_mld_interp = approx(prof_dat1$yday, prof_dat1$chl_mld, xout = yday_list, rule = 2, na.rm = FALSE)$y
chl_zeu_interp = approx(prof_dat1$yday, prof_dat1$chl_zeu, xout = yday_list, rule = 2, na.rm = FALSE)$y

prof_dat_extend <- tibble(
  yday = yday_list,
  date = date_list,
  type = "final",
  MLD = MLD_interp,
  zeu = zeu_interp,
  chl_mld = chl_mld_interp,
  chl_zeu = chl_zeu_interp
)

prof_dat_extend <- prof_dat_extend |> 
  mutate(mld_s = rollmean(MLD, k = 24, fill = NA, align = "center"),
         zeu_s = rollmean(zeu, k = 24, fill = NA, align = "center"))

ggplot(prof_dat_extend)+
  geom_line(aes(x = date, y = - mld_s))+
  geom_point(aes(x = date, y = - MLD), size = 2)+
  geom_point(aes(x = date, y = - zeu, color = "Zeu"))




prof_dat_extend <- prof_dat_extend |> 
  filter(!is.na(mld_s)) |> 
  mutate(MLD_diff = mld_s - lag(mld_s),
         zeu_diff = zeu_s - lag(zeu_s),
         chl_mld_diff = chl_mld - lag(chl_mld),
         chl_zeu_diff = chl_zeu - lag(chl_zeu),
         date_diff = yday - lag(yday),
         rate_mld = (1/mld_s) * (MLD_diff/date_diff),
         p = case_when(mld_s > zeu ~ chl_mld,
                       mld_s < zeu ~ chl_zeu),
         p_diff = p - lag(p))

ggplot(prof_dat_extend)+
  geom_point(aes(x= date, y = rate_mld, color = MLD_diff))+
  scale_color_distiller(palette = "RdBu")+
  theme_dark()

prof_dat1 <- prof_dat_extend |> 
  mutate(rp = case_when(mld_s > zeu & rate_mld > 0 ~ (1/p) * (p_diff/date_diff),
                        mld_s > zeu & rate_mld < 0 ~ (1/chl_mld) * (chl_mld_diff/date_diff),
                        mld_s < zeu ~ (1/p) * (p_diff / date_diff)),
         rp_s = slide_index_dbl(rp, date, mean, 
                                .before = 10, .after = 10, na.rm = TRUE))

ggplot(prof_dat1)+
  geom_line(aes(x = date, y = rp_s))




# profile diagnostic ------------------------------------------------------


prof_number_to_select <- 66

data_to_plot <- filter(dat_all_pp, prof_number == prof_number_to_select) |> filter(depth < 200)

ggplot(data_to_plot)+
  geom_path(aes(x = cbpm_npp, y = - depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)+
  xlim(0,0.1)

ggplot(data_to_plot)+
  geom_point(aes(x = sigma0, y = -depth))+
  geom_hline(aes(yintercept = -mld_s), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)

ggplot(data_to_plot)+
  geom_point(aes(x = chla, y = -depth))
geom_hline(aes(yintercept = -mld_s), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)

ggplot(data_to_plot)+
  geom_point(aes(x = iPAR, y = -depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)


ggplot(data_to_plot)+
  geom_point(aes(x = nitrate, y = -depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)

