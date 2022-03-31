
library(usethis)
library(tidyverse)
library(readxl)

distr_table <- read_excel("data-raw/distr_table.xlsx") %>%
  mutate(type = factor(type, levels = c("binary", "categorical", "ranking", "count", "integer", "directional", "interval", "compositional", "duration", "real"))) %>%
  mutate(dim = factor(dim, levels = c("uni", "multi"))) %>%
  mutate(orthog = as.logical(orthog)) %>%
  mutate(default = as.logical(default)) %>%
  arrange(distr, param) %>%
  as.data.frame()

distr_table

use_data(distr_table, internal = TRUE, overwrite = TRUE)


