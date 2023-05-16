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
  distinct() %>%
  mutate(dist = ifelse(is.na(evente), "Student-T", "Binomial"),
         meane = ifelse(dist == "Binomial", evente/ne, meane),
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
  select(studlab, dist, ne, meane, sde, nc, meanc, sdc,
         mean, sd, df, CIlower, CIupper,
         lowerte, upperte, te) %>%
  distinct()


NNT <- 3
epsilon <- 1/NNT
idxs <- matrix(findInterval(as.matrix(meta_fil[c("CIlower", "CIupper")]),
               c(-epsilon,epsilon)), ncol = 2)

CI_mat <- as.matrix(meta_fil[c("CIlower", "CIupper")])
REACT::REACT_forestplot(CI_mat, NNT = NNT, study_names = meta_fil$studlab,
                        point_estim = meta_fil$mean)



meta_fil %>%
  mutate(color = factor(ifelse(idxs[,1] == 1 & idxs[,2] == 1, 0,
                               ifelse((idxs[,1] == 0 & idxs[,2] == 0) |
                                        (idxs[,1] == 2 & idxs[,2] == 2), 1,
                                      1/2)),
                        levels = c(0, 1/2, 1),
                        labels = c("accept", "agnostic", "reject"))) %>%
  ggplot(aes(y = studlab, x = mean, xmin = CIlower, xmax = CIupper, col = color)) +
  geom_point(fill = NA) +
  annotate('rect', xmin = -epsilon, xmax = epsilon,
           ymin = 0, ymax = 13, alpha=.5, fill='lightblue') +
  geom_point() +
  geom_errorbarh(height=.1) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(labels = c("accept", "agnostic", "reject"),
                     values=c("green", "yellow", "red"))

