## code to prepare `meta` dataset goes here
library(dplyr)
meta <- read.csv("data-raw/metanalise.csv", header = T) |>
  select(studlab, te, lowerte, upperte) |>
  filter(lowerte!=0)

usethis::use_data(meta, overwrite = TRUE)
