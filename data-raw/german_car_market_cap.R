
library(usethis)
library(tidyverse)
library(readxl)

german_car_market_cap <- read_excel("data-raw/german_car_market_cap.xlsx") %>%
  pivot_longer(cols = 2:4, names_to = "car_manufacturer", values_to = "market_cap") %>%
  mutate(year = as.integer(year)) %>%
  mutate(car_manufacturer = factor(car_manufacturer, levels = unique(car_manufacturer))) %>%
  as.data.frame()

german_car_market_cap

use_data(german_car_market_cap, internal = FALSE, overwrite = TRUE)


