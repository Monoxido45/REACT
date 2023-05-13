library(dplyr)
library(ggplot2)

meta <- read.csv("data-raw/metanalise.csv", header = T)
View(meta)

meta_fil <- meta %>%
  select(studlab, te, lowerte, upperte) %>%
  filter(lowerte!=0)

meta_fil %>%
  mutate(color = ifelse(apply(cbind(meta_fil$lowerte, meta_fil$upperte),
                              1, function(x) findInterval(1, x)) == 1,
         "contains", "does not contain")) %>%
  ggplot(aes(y = studlab, x = te, xmin = lowerte, xmax = upperte, col = color)) +
  geom_point() + geom_errorbarh(height=.1) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_color_manual(values=c("green", "red"))
  scale_y_continuous(name = "", breaks=1:nrow(meta), labels=meta$studlab)

