library(tidyverse)
library(zoo)
library(gsw)
library(castr)
library(slider)

dat <- read_csv("argo_3901586_interp.csv")

dat <- dat |> 
  mutate(sa = gsw_SA_from_SP(sal, depth, lon ,lat),
         ct = gsw_CT_from_t(sa, temp, depth),
         sigma0 = gsw_sigma0(sa, ct))

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
prof_dat <- dat |>
  group_by(prof_number, date) |>
  arrange(depth) |> 
  summarise(
    MLD = mld(sigma0, depth, ref.depths = 0:2, criteria = 0.03, default.depth = 300),
    chl_cline = clined(chla, depth),
    zeu = zeu_calc(iPAR, depth)
  ) |>
  ungroup() |>
  filter(!is.na(zeu)) |> 
  mutate(
    mld_s = pmax(smooth.spline(date, MLD, spar = 0.7)$y, 2),
    zeu_s = pmax(smooth.spline(date, zeu, spar = 0.6)$y, 2)
  ) |>
  mutate(
    delta_MLD = mld_s - lag(mld_s),
    delta_zeu = zeu_s - lag(zeu_s)
  ) |>
  replace_na(list(delta_MLD = 0))


ggplot(na.omit(prof_dat))+
  geom_point(aes(x = date, y = -MLD, colour = "MLD"))+
  geom_line(aes(x = date, y = -mld_s, colour = "MLD"))+
  geom_point(aes(x = date, y = -zeu, colour = "Zeu"))+
  geom_line(aes(x = date, y = -zeu_s, colour = "Zeu"))


dat_combined <- left_join(dat, prof_dat)

integrated_val <- dat_combined |> 
  filter(!is.na(zeu_s) & !is.na(mld_s)) |> 
  group_by(prof_number, date) |> 
  summarise(chl_mld = integrate(chla, depth, from =0, to = abs(unique(MLD))),
            chl_zeu = integrate(chla, depth, from =0, to = abs(unique(zeu))))

integrated_val <- integrated_val |>
  ungroup() |> 
  arrange(date) |>
  mutate(
    chl_mld_s = slide_index_dbl(chl_mld, date, mean, 
                                .before = 10, .after = 10, na.rm = TRUE),
    chl_zeu_s = slide_index_dbl(chl_zeu, date, mean, 
                                .before = 10, .after = 10, na.rm = TRUE)
  )
ggplot(na.omit(integrated_val))+
  geom_point(aes(x = date, y = chl_mld), color = "red")+
  geom_line(aes(x = date, y = chl_mld_s), color = "red")+
  geom_point(aes(x = date, y = chl_zeu), color = "green")+
  geom_line(aes(x = date, y = chl_zeu_s), color = "green")+
  ylab("Chl stock")

dat_combined <- left_join(dat_combined, integrated_val)

prof_number_to_select <- 70

data_to_plot <- filter(dat_combined, prof_number == prof_number_to_select)

ggplot(data_to_plot)+
  geom_point(aes(x = sigma0, y = -depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)

ggplot(data_to_plot)+
  geom_point(aes(x = chla, y = -depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)

ggplot(data_to_plot)+
  geom_point(aes(x = iPAR, y = -depth))+
  geom_hline(aes(yintercept = -MLD), color = "red")+
  geom_hline(aes(yintercept = -zeu), color = "blue")+
  ylim(-200, 0)


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
         date_diff = yday - lag(yday),
         p = case_when(mld_s > zeu ~ chl_mld_s,
                       mld_s < zeu ~ chl_zeu_s),
         p_diff = p - lag(p))

prof_dat1 <- prof_dat1 |> mutate(rate_mld = (1/mld_s) * (MLD_diff/date_diff))

ggplot(prof_dat1)+
  geom_point(aes(x= date, y = rate_mld, color = MLD_diff))+
  geom_text(aes(x = date, y = rate_mld, label = prof_number))+
  scale_color_distiller(palette = "RdBu")+
  theme_dark()

prof_dat1 <- prof_dat1 |> 
  mutate(rp = case_when(mld_s > zeu & rate_mld > 0 ~ (1/p) * (p_diff/date_diff),
                        mld_s > zeu & rate_mld < 0 ~ (1/chl_mld_s) * (chl_mld_diff/date_diff),
                        mld_s < zeu ~ (1/p) * (p_diff / date_diff)),
         rp_s = slide_index_dbl(rp, date, mean, 
                                .before = 10, .after = 10, na.rm = TRUE))


ggplot(prof_dat1)+
  geom_point(aes(x = date, y = rp_s))



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
         p = case_when(mld_s > zeu ~ chl_mld_s,
                       mld_s < zeu ~ chl_zeu_s),
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

