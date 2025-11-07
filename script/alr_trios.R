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

ggplot(df)+
  geom_point(aes(x = Wavelength, y = Counts, color = datetime))

trios_timestamps <- df |> select(datetime) |> unique()

eco_data_matched <- eco_data |> 
  select(depth, datetime, lon, lat) |> 
  distinct() |> 
  rowwise() |> 
  mutate(
    datetime = lubridate::as_datetime(datetime, tz="UTC"),
    nearest_trios_time = trios_timestamps$datetime[which.min(abs(trios_timestamps$datetime - datetime))]
  ) |> 
  ungroup()

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

ggplot(trios_matched)+
  geom_point(aes(x = Wavelength, y = Counts, color = depth))+
  scale_color_viridis_c()
