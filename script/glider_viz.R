library(tidyverse)
library(vroom)
library(gsw)
library(castr)

cabot <- read_csv("data/glider/csv/Cabot_raw_profile.csv")
churchill <- vroom("data/glider/csv/Churchill_raw_profile.csv")
nelson <- vroom("data/glider/csv/Nelson_raw_profile.csv")
doombar <- vroom("data/glider/csv/Doombar_raw_profile.csv")


gliders <- bind_rows(cabot, churchill, nelson, doombar)

gliders_asc <- gliders |> filter(direction == 'asc')

write_csv(gliders_asc, 'data/glider/csv/compiled_asc.csv')

check <- data.frame(table(cabot$profile))

dat_plot <- filter(cabot, profile %in% c(1030:1064))
ggplot(dat_plot)+
  geom_point(aes(x = profile, y = -pres, color = direction))

test <- filter(cabot, profile == 1005)

ggplot(test)+
  geom_point(aes(x = TEMP, y = - pres))

test <- test |>
  arrange(pres) |> 
  mutate(sp = gsw_SP_from_C(CNDC * 10, TEMP, pres),
                       sa = gsw_SA_from_SP(sp, pres, -24, 60),
                       ct = gsw_CT_from_pt(sa, TEMP),
                       sigma = gsw_sigma0(sa, ct),
                       mld = mld(sigma, pres))

ggplot(test)+
  geom_point(aes(x = time, y = -pres))+
  geom_hline(aes(yintercept = - mld))

ggplot(gliders_asc)+
  geom_point(aes(x = CNDC, y = TEMP, color = pres))

miss <- filter(gliders, direction == 'NULL')
