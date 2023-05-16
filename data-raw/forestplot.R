library(dplyr)
library(ggplot2)
# importing data
meta <- read.csv("data-raw/metanalise.csv", header = T)

meta_fil <- meta %>%
  filter(ne != 0 & nc != 0) %>%
  distinct() %>%
  select(studlab, te, lowerte, upperte)

# replicating forest plot from paper
meta_fil %>%
  mutate(color = ifelse(apply(cbind(meta_fil$lowerte, meta_fil$upperte),
                              1, function(x) findInterval(1, x)) == 1,
         "contains", "does not contain")) %>%
  ggplot(aes(y = studlab, x = te, xmin = lowerte, xmax = upperte, col = color)) +
  geom_point() + geom_errorbarh(height=.1) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_color_manual(values=c("green", "red"))

# testing NNT based test for difference of proportions
obj <- REACT::NNT_indep_test(alpha = 0.05, NNT = 3,
         event_e = 77, n_e = 257, event_c = 8, n_c = 47)

REACT:::plot.simple_REACT(obj)


# graphical analysis
alpha <- .05
# changing meta analysis object
meta_fil <- meta %>%
  filter(ne != 0 & nc != 0) %>%
  filter(ncpooled != 614) %>% #removing error
  distinct() %>%
  mutate(dist = ifelse(is.na(evente), "Student-T", "Binomial"),
         study_class = ifelse(dist == "Student-T", 3,
                              ifelse(nepooled == 18224, 1, 2)))
# getting pooled estimates
pooled <- meta_fil %>%
  mutate(sume = meane*ne, sumc = meanc*nc,
         SQe = sde^2*(ne-1), SQc = sdc^2*(nc-1)) %>%
  group_by(dist, study_class) %>%
  summarise(ne = sum(ne), evente = sum(evente),
            nc = sum(nc), eventc = sum(eventc),
            meane = sum(sume)/sum(ne), meanc = sum(sumc)/sum(nc),
            sde = sqrt(sum(SQe)/(sum(ne))),
            sdc = sqrt(sum(SQc)/(sum(nc)))) %>%
  ungroup() %>%
  mutate(studlab = "Pooled")
# combining with original dataset and obtaining CI
meta_fil <- bind_rows(meta_fil, pooled) %>%
  mutate(meane = ifelse(dist == "Binomial", evente/ne, meane),
         meanc = ifelse(dist == "Binomial", eventc/nc, meanc),
         sde = ifelse(dist == "Binomial", sqrt(abs(meane*(1-meane)/ne)), sde),
         sdc = ifelse(dist == "Binomial", sqrt(abs(meanc*(1-meanc)/nc)), sdc),
         mean = meane - meanc,
         Se = sde^2/ne,
         Sc = sdc^2/nc,
         sd = ifelse(dist == "Binomial",
                     sqrt(sde^2 + sdc^2),
                     sqrt(Se + Sc)),
         df = ifelse(dist == "Binomial", NA,
                     round((Se + Sc)^2/(Se^2/(ne-1) + Sc^2/(nc-1)))),
         CIlower = mean + ifelse(dist=="Binomial",
                                 qnorm(alpha/2)*sd,
                                 qt(alpha/2, df)*sd),
         CIupper = mean + ifelse(dist=="Binomial",
                                 qnorm(1 - alpha/2)*sd,
                                 qt(1 - alpha/2, df)*sd)
  ) %>%
  select(studlab, dist, study_class, ne, meane, sde, nc, meanc, sdc,
         mean, sd, df, CIlower, CIupper)

#Forestplot
NNT <- 10
CI_mat <- as.matrix(meta_fil[c("CIlower", "CIupper")])
g <- REACT::REACT_forestplot(CI_mat, NNT = NNT, study_names = meta_fil$studlab,
                             point_estim = meta_fil$mean)
g$data$studlab <- relevel(g$data$studlab, "Pooled")
g$data <- bind_cols(g$data,
                    study_class = factor(meta_fil$study_class,
                                         levels = c(1,2,3),
                                         labels = c("Follow-up (Binomial)",
                                                    "Pharmacotherapy (Binomial)",
                                                    "Follow-up (T)")))
g + facet_grid(. ~ study_class, scales="free")
