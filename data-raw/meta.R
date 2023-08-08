## code to prepare `meta` dataset goes here
library(dplyr)
meta <- read.csv("data-raw/metanalise.csv", header = T) %>%
  filter(ne != 0 & nc != 0) %>%
  filter(ncpooled != 614) %>% #removing error
  distinct() %>%
  mutate(dist = ifelse(is.na(evente), "Student-T", "Binomial"),
         study_class = ifelse(dist == "Student-T", 3,
                              ifelse(nepooled == 18224, 1, 2))) %>%
  # removing t-student studies
  filter(study_class != 3) %>%
  mutate(study_class = factor(study_class,
                              levels = c(1,2),
                              labels = c("Follow-up",
                                         "Pharmacotherapy"))) %>%
  select(studlab, study_class, evente, ne, eventc, nc, te, lowerte, upperte)


usethis::use_data(meta, overwrite = TRUE)
