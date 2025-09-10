library(tidyverse)
library(zoo)
library(arrow)
library(gsw)
library(castr)
# Making a 200:1000 depth argo table

argo <- read_parquet("data/argo_pq/biocarbon_floats_table.parquet") |>
    mutate(depth = round(PRES))

new_ref <- select(argo, PLATFORM_NUMBER, JULD) |>
    unique() |>
    dplyr::filter(PLATFORM_NUMBER %in% c("1902304 ", "4903532 "))
depth <- tibble(depth = c(0:1000))

new_ref <- new_ref |>
    crossing(depth) |>
    left_join(argo)

new_ref <- new_ref |>
    mutate(JULD_date = lubridate::date(JULD), prof_id = paste0(PLATFORM_NUMBER, JULD_date)) |>
    filter(JULD_date > lubridate::date("2023-12-31")) |>
    filter(!prof_id %in% c("1902304 2024-01-19")) |>
    group_by(prof_id) %>%
    # Interpolate within each profile
mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), ~na.approx(.x, depth, na.rm = FALSE))) %>%
    mutate(across(c(LONGITUDE, LATITUDE, TEMP, PSAL, CHLA_ADJUSTED, BBP700_ADJUSTED, NITRATE_ADJUSTED, DOXY_ADJUSTED), ~na.locf(.x, na.rm = FALSE,
        fromLast = TRUE))) %>%
    mutate(across(c(NITRATE_ADJUSTED, DOXY_ADJUSTED), ~na.locf(.x, na.rm = FALSE))) %>%
    ungroup()

new_ref$PLATFORM_NUMBER <- as.numeric(new_ref$PLATFORM_NUMBER)

ncp_df <- new_ref |>
    filter(PLATFORM_NUMBER == "1902304" | PLATFORM_NUMBER == "4903532") |>
    mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth), sigma0 = gsw_sigma0(PSAL, TEMP))

prof_dat <- ncp_df |>
    group_by(PLATFORM_NUMBER, JULD) |>
    mutate(nitrate_smoothed = smooth(NITRATE_ADJUSTED, k = 5, n = 2)) |>
    summarise(MLD = mld(sigma0, depth, ref.depths = 0:5, criteria = 0.03, default.depth = 200), nitracline = clined(nitrate_smoothed, depth),
        chl_cline = clined(CHLA_ADJUSTED, depth)) |>
    ungroup() |>
    group_by(PLATFORM_NUMBER) |>
    mutate(mld_s = smooth.spline(JULD, MLD, spar = 0.6)$y) |>
    mutate(delta_MLD = mld_s - lag(mld_s)) |>
    replace_na(list(delta_MLD = 0))

ncp_df <- left_join(ncp_df, prof_dat) |>
    mutate(layer = case_when(depth <= mld_s ~ "MLD", depth > mld_s & depth <= max(prof_dat$mld_s) + 20 ~ "SUB MLD", depth > max(prof_dat$mld_s) + 20 ~
        "DEEP"))

ggplot(ncp_df) + geom_tile(aes(y = -depth, x = JULD, fill = NITRATE_ADJUSTED)) + geom_line(aes(x = JULD, y = -MLD)) + geom_line(aes(x = JULD,
    y = -mld_s), linetype = "dashed") + scale_fill_viridis_c() + ylim(-1000, 0) + facet_wrap(. ~ PLATFORM_NUMBER)

location <- select(ncp_df, JULD, LONGITUDE, LATITUDE) |>
    unique()


integrated_df <- ncp_df |>
    select(PLATFORM_NUMBER, JULD, NITRATE_ADJUSTED, layer, mld_s, depth) |>
    group_by(PLATFORM_NUMBER, JULD, layer) |>
    summarise(integrated_nitrate = sum(NITRATE_ADJUSTED), mean_nitrate = mean(NITRATE_ADJUSTED), layer_top = min(depth), layer_bottom = max(depth)) |>
    ungroup() |>
    left_join(prof_dat) |>
    group_by(PLATFORM_NUMBER, layer) |>
    mutate(mean_nitrate_smoothed = smooth.spline(JULD, mean_nitrate, spar = 0.6)$y, integrated_nitrate_smoothed = smooth.spline(JULD, integrated_nitrate,
        spar = 0.6)$y)


integrated_df$date <- lubridate::date(integrated_df$JULD)

ggplot(integrated_df) + geom_point(aes(x = date, y = mean_nitrate, color = layer, shape = as.factor(PLATFORM_NUMBER))) + geom_path(aes(x = date,
    y = mean_nitrate_smoothed, color = layer, linetype = as.factor(PLATFORM_NUMBER))) + scale_color_brewer(palette = "Set1", name = "wmo")

# ggsave('output/plots/pp_floats/multi_layer_mean_nitrate.png', dpi = 300, width = 20, height = 15, units = 'cm')

ggplot(integrated_df) + geom_point(aes(x = date, y = integrated_nitrate, color = layer, shape = as.factor(PLATFORM_NUMBER))) + geom_path(aes(x = date,
    y = integrated_nitrate_smoothed, color = layer, linetype = as.factor(PLATFORM_NUMBER))) + scale_color_brewer(palette = "Set1", name = "wmo")
# ggsave('output/plots/pp_floats/multi_layer_integrated_nitrate.png', dpi = 300, width = 20, height = 15, units = 'cm')

integrated_wide <- integrated_df |>
    select(-integrated_nitrate, -mean_nitrate) |>
    pivot_wider(names_from = layer, values_from = c(integrated_nitrate_smoothed, mean_nitrate_smoothed, layer_bottom, layer_top)) |>
    janitor::clean_names() |>
    mutate(ncp = -((lag(integrated_nitrate_smoothed_mld) + ((delta_mld/(mld + abs(delta_mld))) * lag(integrated_nitrate_smoothed_sub_mld))) -
        integrated_nitrate_smoothed_mld) * 6.6/as.numeric(juld - lag(juld)))

ggplot(integrated_wide) + geom_line(aes(x = date, y = ncp/12, color = as.factor(platform_number)))

# ggsave('output/plots/pp_floats/ncp_multi_layer.png', dpi = 300, width = 20, height = 15, units = 'cm')

npp_df <- read_csv("data/argo_pp_estimations_floats.csv") |>
    mutate(ct = gsw_CT_from_t(PSAL, TEMP, depth), sigma0 = gsw_sigma0(PSAL, TEMP)) |>
    filter(PLATFORM_NUMBER == "1902304" | PLATFORM_NUMBER == "4903532") |>
    janitor::clean_names() |>
    select(platform_number, juld, pp) |>
    group_by(platform_number, juld) |>
    summarise(integrated_npp = sum(pp, na.rm = TRUE)) |>
    ungroup() |>
    mutate(date = lubridate::date(juld)) |>
    select(-juld)

comparison <- left_join(integrated_wide, npp_df) |>
    pivot_longer(c(ncp, integrated_npp), values_to = "pp", names_to = "type") |>
    mutate(type = case_when(type == "integrated_npp" ~ "integrated NPP", type == "ncp" ~ "integrated NCP"))

ggplot(comparison) + geom_line(aes(x = date, y = pp/12, colour = type)) + geom_vline(aes(xintercept = lubridate::as_date("2024-06-15")), linetype = "dashed") +
    facet_wrap(. ~ platform_number) + theme_bw(base_size = 18) + ylab("mmol C m-2 d-1") + scale_color_brewer(palette = "Set1", name = "")
# ggsave('output/plots/pp_floats/ncp_multi_layer_float_comp_v2.png', dpi = 300, width = 20, height = 15, units = 'cm')

ggplot(filter(comparison, platform_number == 1902304)) + geom_line(aes(x = date, y = pp/12, colour = type)) + geom_vline(aes(xintercept = lubridate::as_date("2024-06-15")),
    linetype = "dashed") + theme_bw(base_size = 18) + ylab("mmol C m-2 d-1") + scale_color_brewer(palette = "Set1", name = "")
# ggsave('output/plots/pp_floats/npp_ncp_1902304.png', dpi = 300, width = 20, height = 15, units = 'cm')

plot <- filter(ncp_df, JULD == unique(ncp_df$JULD)[15])

ggplot(plot) + geom_path(aes(x = NITRATE_ADJUSTED, y = -depth)) + geom_hline(aes(yintercept = -mld_s)) + ggtitle(unique(plot$JULD_date))
