
library(usethis)
library(tidyverse)
library(readr)

file_data <- "data-raw/sp500_daily.csv"
file_data_new <- "data-raw/sp500_daily_new.csv"

if (file.exists(file_data_new)) {
  df_old <- read_csv(file_data)
  last_day <- max(as.Date(df_old$Date, format = "%m/%d/%Y"))
  df_new <- read_csv(file_data_new) %>%
    filter(as.Date(Date, format = "%m/%d/%Y") > last_day)
  df_write <- rbind(df_new, df_old)
  write_csv(df_write, file = file_data)
  unlink(file_data_new)
}

sp500_daily <- read_csv(file_data) %>%
  select(date = Date, open = Open, high = High, low = Low, close = "Close/Last") %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y")) %>%
  filter(open > 0)

use_data(sp500_daily, internal = FALSE, overwrite = TRUE)


