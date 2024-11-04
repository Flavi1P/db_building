library(tidyverse)
library(patchwork)
df <- read_csv("data/cidre/data_sensor_1_20240528T050612.csv")

colnames(df) <- paste0('x', seq(1,12))

df_low <- filter(df, x7 < 200)

ggplot(df)+
  geom_point(aes(x = x7, y = -x3))+
  ylim(-120, -80)+
  xlim(45, 70)



ctd4 <- read_csv("data/cidre/LOV_ECO_ctd004s/data_sensor_0_20240527T071543.csv", 
                 col_names = FALSE)

colnames(ctd4) <- paste0('x', seq(1,12))

ggplot(ctd4)+
  geom_point(aes(x = x7, y = -x3))+
  ylim(-450, -400)+
  xlim(45,60)


ctd1 <- read_csv("data/cidre/LOV_ECO_ctd001s/data_sensor_0_20240526T080613.csv", 
                 col_names = FALSE)
ctd3 <- read_csv("data/cidre/LOV_ECO_ctd003s/data_sensor_0_20240527T040338.csv", 
                 col_names = FALSE)

colnames(ctd1) <- paste0('x', seq(1,12))
colnames(ctd3) <- paste0('x', seq(1,12))

ggplot(ctd4)+
  geom_point(aes(x = x7, y = -x3))+
  ylim(-450, -400)+
  xlim(45,60)
ggplot(ctd1)+
  geom_point(aes(x = x7, y = -x3))+
  ylim(-1000, -500)+
  xlim(45,75)+
  ggtitle('Only one Echo')+

ggplot(ctd3)+
  geom_point(aes(x = x7, y = -x3))+
  ylim(-1000, -500)+
  xlim(45,75)+
  ggtitle('Two Echo')

ctd3 <- ctd3 |> mutate(depth = round(x3)) |> 
  group_by(depth) |> 
  summarise(noise = sd(x7)) |> 
  ungroup()

ctd1 <- ctd1 |> mutate(depth = round(x3)) |> 
  group_by(depth) |> 
  summarise(noise = sd(x7)) |> 
  ungroup()

ggplot(ctd3)+
  geom_point(aes(x = noise, y = -depth))+
  ylim(-1000, -500)+
  xlim(0,15)+
  ggtitle('Only one echo')+
ggplot(ctd1)+
  geom_point(aes(x = noise, y = -depth))+
  ylim(-1000, -500)+
  xlim(0,15)+
  ggtitle('two echo')
