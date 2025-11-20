library(tidyverse)

eco_data <- read_csv("output/ALR6_fluorometer_data.csv")

file <- "data/ALR/Radiometer_TRIOS_RAMSES/RamsesG2_2024-06-11_23-52-38_file_1.log"     # <- change to your filename
raw <- read_lines(file)

# ---- 2. Extract wavelengths (first line only) ----
wavelength_line <- raw[grepl("^OCSTime:.*Wavelengths:", raw)]
wavelengths <- wavelength_line %>%
  str_extract("Wavelengths:.*") %>%
  str_remove("Wavelengths:\\s*") %>%
  str_split(",") %>%
  unlist() %>%
  as.numeric() 

wavelengths <- unique(wavelengths)
# ---- 3. Extract all measurement lines (with counts) ----
measure_lines <- raw[grepl("WaveLengthCounts:", raw)]

# Helper function to parse a single line
parse_record <- function(line) {
  list(
    OCSTime  = str_extract(line, "OCSTime:\\s*[^,]+") %>% str_remove("OCSTime:\\s*"),
    DeviceTime = str_extract(line, "DeviceTime:\\s*[^,]+") %>% str_remove("DeviceTime:\\s*"),
    Pressure = str_extract(line, "Pressure:\\s*[^,]+") %>% str_remove("Pressure:\\s*") %>% as.numeric(),
    PreInclinationAngle  = str_extract(line, "PreInclinationAngle:\\s*[^,]+") %>% str_remove("PreInclinationAngle:\\s*") %>% as.numeric(),
    PostInclinationAngle = str_extract(line, "PostInclinationAngle:\\s*[^,]+") %>% str_remove("PostInclinationAngle:\\s*") %>% as.numeric(),
    Counts = line %>%
      str_extract("WaveLengthCounts:.*") %>%
      str_remove("WaveLengthCounts:\\s*") %>%
      str_split(",") %>%
      unlist() %>%
      as.numeric()
  )
}

# ---- 4. Parse all records ----
parsed <- map(measure_lines, parse_record)

# ---- 5. Convert to dataframe ----
df <- map_dfr(parsed, function(rec) {
  tibble(
    OCSTime = rec$OCSTime,
    DeviceTime = rec$DeviceTime,
    Pressure = rec$Pressure,
    PreInclination = rec$PreInclinationAngle,
    PostInclination = rec$PostInclinationAngle,
    Wavelength = wavelengths,
    Counts = rec$Counts
  )
})

df <- df %>%
  mutate(datetime = lubridate::as_datetime(OCSTime, tz="UTC"))
# ---- 6. Example: show first rows ----

#ggplot(df)+
#  geom_point(aes(x = Wavelength, y = Counts, color = datetime))

trios_timestamps <- df |> select(datetime) |> unique()

library(data.table)

# convert to data.table
eco_dt   <- as.data.table(eco_data)[, .(depth, datetime, lon, lat)]
trios_dt <- as.data.table(trios_timestamps)[, .(datetime)]

# make sure both datetime columns are POSIXct
eco_dt[,   datetime := as.POSIXct(datetime, tz = "UTC")]
trios_dt[, datetime := as.POSIXct(datetime, tz = "UTC")]

# sort trios timestamps (required for rolling join)
setkey(trios_dt, datetime)

# find nearest timestamp by rolling join (bidirectional)
eco_dt[
  trios_dt,
  nearest_trios_time := i.datetime,
  on = .(datetime),
  roll = "nearest"
]

eco_data_matched <- eco_dt |> as_tibble() |> na.omit()

eco_data_to_match <- eco_data_matched |> 
  mutate(diff_time = nearest_trios_time - datetime) |> 
  group_by(nearest_trios_time) |> 
  filter(diff_time == min(diff_time)) |> 
  select(nearest_trios_time, lon, lat, depth) |> 
  distinct()

trios_matched <- eco_data_to_match |>
  left_join(
    df |> 
      select(-OCSTime, -DeviceTime, -PreInclination, -PostInclination),
    by = c("nearest_trios_time" = "datetime")
  )

trios_matched <- trios_matched |> 
  mutate(intesity = Counts / 653535)

smaller_trios <- filter(trios_matched, depth > 1 & depth < 30 & Wavelength > 400 & Wavelength < 800)

ggplot(smaller_trios)+
  geom_path(aes(x = Wavelength, y = intesity, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()
ggsave("ALR6_trios_spectra.jpg", width = 30, height = 20, units = "cm", dpi = 300)

ggplot(eco_data_to_match)+
  geom_point(aes(x = nearest_trios_time, y = -depth))

samples <- smaller_trios |> select(nearest_trios_time) |> unique()

test <- smaller_trios |> filter(nearest_trios_time %in% samples$nearest_trios_time[c(1052:1200)])

ggplot(test)+
  geom_path(aes(x = Wavelength, y = intesity, color = depth, group = nearest_trios_time))+
  scale_color_viridis_c()

write_csv(test, "output/ALR6_sample_data.csv")

write_csv(trios_matched, "data/ALR/alr6_ramses_timed.csv")
