
library(usethis)
library(tidyverse)
library(lubridate)
library(hms)
library(readxl)

bookshop_orders <- read_excel("data-raw/bookshop_orders.xlsx") %>%
  select(-id) %>%
  filter(quantity > 0) %>%
  mutate(time = dmy_hms(paste(date, time))) %>%
  rowid_to_column("id") %>%
  rename(datetime = time) %>%
  mutate(date = as.Date(datetime)) %>%
  mutate(time = as_hms(datetime)) %>%
  mutate(duration = as.integer(as.numeric(datetime - lag(datetime)) / 60)) %>%
  mutate(duration_adj = recode(as.numeric(duration), "0" = 0.5))

model_spl <- smooth.spline(as.vector(bookshop_orders$time[-1]), bookshop_orders$duration_adj[-1], df = 10)

bookshop_orders <- bookshop_orders %>%
  mutate(duration_spl = predict(model_spl, x = as.vector(time))$y) %>%
  mutate(duration_adj = duration_adj / duration_spl) %>%
  select(id, datetime, quantity, duration, duration_adj) %>%
  as.data.frame()

use_data(bookshop_orders, internal = FALSE, overwrite = TRUE)
