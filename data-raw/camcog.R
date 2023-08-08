# code to prepare camcog dataset
library(dplyr)
library(forcats)
camcog <- read.csv("data-raw/CAMCOG.csv") %>%
  rename(diag_group = Diagnostico) %>%
  mutate(
    diag_group = as.factor(diag_group),
    diag_group = fct_recode(diag_group, CG = "GC",
                            AD = "DA",
                            MCI = "CCL"))

usethis::use_data(camcog, overwrite = TRUE)
