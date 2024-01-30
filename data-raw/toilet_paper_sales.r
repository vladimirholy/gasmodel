
library(tidyverse)
library(readxl)

toilet_paper_sales <- read_xlsx("data-raw/toilet_paper_sales.xlsx") %>%
  mutate(date = as.Date(date)) %>%
  mutate(promo = as.integer(promo)) %>%
  mutate(quantity = as.integer(quantity))  %>%
  mutate(wday = wday(date, week_start = 1)) %>%
  mutate(val = 1L) %>%
  pivot_wider(names_from = "wday", values_from = "val", names_prefix = "day_", values_fill = 0L) %>%
  arrange(date) %>%
  select(date, monday = day_1, tuesday = day_2, wednesday = day_3, thursday = day_4, friday = day_5, saturday = day_6, sunday = day_7, promo, quantity) %>%
  as.data.frame()

use_data(toilet_paper_sales, internal = FALSE, overwrite = TRUE)
