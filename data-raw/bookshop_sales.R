
library(usethis)
library(tidyverse)
library(lubridate)
library(readxl)

bookshop_sales <- read_excel("data-raw/bookshop_sales.xlsx") %>%
  mutate(order = as.integer(id - min(id) + 1)) %>%
  mutate(time = dmy_hms(paste(date, time))) %>%
  mutate(quantity = as.integer(volume)) %>%
  filter(time < "2018-12-21") %>%
  select(order, time, quantity) %>%
  as.data.frame()

bookshop_sales

use_data(bookshop_sales, internal = FALSE, overwrite = TRUE)


