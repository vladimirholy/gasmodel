
library(usethis)
library(tidyverse)
library(readxl)

rankings <- read_excel("data-raw/ice_hockey_championships.xlsx", sheet = "rankings") %>%
  column_to_rownames("year") %>%
  as.matrix() %>%
  apply(1, function(x) { if (any(!is.na(x))) { replace(x, is.na(x), Inf) } else { x } }) %>%
  t()

hosts <- read_excel("data-raw/ice_hockey_championships.xlsx", sheet = "hosts") %>%
  column_to_rownames("year") %>%
  as.matrix()

ice_hockey_championships <- list(rankings = rankings, hosts = hosts)

use_data(ice_hockey_championships, internal = FALSE, overwrite = TRUE)


