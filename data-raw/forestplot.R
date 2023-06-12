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



# Graphical analysis ------------------------------------------------------
alpha <- .05
# changing meta analysis object
meta_fil <- meta %>%
  filter(ne != 0 & nc != 0) %>%
  filter(ncpooled != 614) %>% #removing error
  distinct() %>%
  mutate(dist = ifelse(is.na(evente), "Student-T", "Binomial"),
         study_class = ifelse(dist == "Student-T", 3,
                              ifelse(nepooled == 18224, 1, 2))) |>
  # removing t-student studies
  filter(study_class != 3)

# getting pooled estimates with heteroscedasticity included
# using {meta} mixed model random effects CI
CI_s <- meta_fil$study_class %>% unique() %>%
  purrr::map_dfr(function(x){
    used_data <- meta_fil |> filter(study_class == x)

    meta_model <- meta::metabin(evente,
            ne,
            eventc,
            nc,
            studlab,
            data= used_data,
            sm= "RD",
            model.glmm = ("CM.AL"),
            method = "Inverse")

    lower <- meta_model %>% purrr::pluck("lower.random")
    upper <- meta_model %>% purrr::pluck("upper.random")
    point_estim <- meta_model %>% purrr::pluck("TE.random")

    list(study_class = x, "CIlower" = lower, "CIupper" = upper, "mean" = point_estim)
  })

# old procedure
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
  filter(dist != "Student-T")
  # changing CI for pooled
meta_fil <- meta_fil %>%
  mutate(
    CIlower = ifelse(studlab == "Pooled" & study_class == 1,
                     CI_s$CIlower[1],
                     ifelse(studlab == "Pooled" & study_class == 2,
                            CI_s$CIlower[2], CIlower)),
    CIupper = ifelse(studlab == "Pooled" & study_class == 1,
                     CI_s$CIupper[1],
                     ifelse(studlab == "Pooled" & study_class == 2,
                            CI_s$CIupper[2], CIupper)),
    mean = ifelse(studlab == "Pooled" & study_class == 1,
                  CI_s$mean[1],
                  ifelse(studlab == "Pooled" & study_class == 2,
                         CI_s$mean[2], mean))
  ) %>%
  # changing point estimation also
  select(studlab, study_class, ne, evente, meane, sde, nc, eventc, meanc, sdc,
         mean, sd, CIlower, CIupper)


#Forestplot
NNT <- c(-0, 6)
CI_mat <- as.matrix(meta_fil[c("CIlower", "CIupper")])
g <- REACT::REACT_forestplot(CI_mat, NNT = NNT, study_names = meta_fil$studlab,
                             point_estim = meta_fil$mean)
g$data$studlab <- relevel(g$data$studlab, "Pooled")
g$data <- bind_cols(g$data,
                    study_class = factor(meta_fil$study_class,
                                         levels = c(1,2),
                                         labels = c("Follow-up",
                                                    "Pharmacotherapy")))

# highlighting pooled
bold.labels <- ifelse(levels(g$data$studlab) == "Pooled", "bold", "plain")

highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue::glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

g + facet_grid(. ~ study_class, scales="free") +
  theme(axis.text.y = element_text(face = bold.labels))+
  scale_y_discrete(labels= function(x) highlight(x, "Pooled", "black"))+
  theme(axis.text.y= ggtext::element_markdown())




# Bayesian tests ----------------------------------------------------------

# NNT forest plot (Bayesian, prior: Beta(0, 0))

library(coda)

my_meta <- meta_fil

alphae <- my_meta$evente
betae <- my_meta$ne - my_meta$evente
alphac <- my_meta$eventc
betac <- my_meta$nc - my_meta$eventc

N <- 1000000
NNT_sim <- matrix(NA, ncol = dim(my_meta)[1], nrow = N)
colnames(NNT_sim) <- my_meta$studlab

set.seed(42)
for(i in 1:dim(my_meta)[1]){
  NNT_sim[,i] <- 1/(rbeta(N, alphac[i], betac[i]) -
                      rbeta(N, alphae[i], betae[i]))
}

NNT_MC <- mcmc(NNT_sim)

my_int <- HPDinterval(NNT_MC)
my_int
rm(NNT_MC)

my_ints <- matrix(NA, ncol = 4, nrow = dim(my_int)[1])
for(i in 1:dim(my_int)[1]){
  where <- findInterval(0, my_int[i,])
  if(where == 2){ #Both negative
    my_ints[i,1] <- my_int[i,1]
    my_ints[i,2] <- my_int[i,2]
  } else if(where == 1){ #Negative and positive
    my_ints[i,1] <- my_int[i,1]
    my_ints[i,2] <- -1
    my_ints[i,3] <- 1
    my_ints[i,4] <- my_int[i,2]
  } else{
    my_ints[i,3] <- my_int[i,1]
    my_ints[i,4] <- my_int[i,2]
  }
}

colnames(my_ints) <- paste("int", 1:4, sep = "")

my_meta <- bind_cols(my_meta, my_ints)

NNT <- c(-10, 6)

my_fun <- function(x){
  if(all(is.na(x))) return(c(NA, NA))
  return(findInterval(NNT, x))
}

c1 <- rowSums(t(apply(my_meta[c("int1", "int2")], 1, my_fun)))/2
c2 <- rowSums(t(apply(my_meta[c("int3", "int4")], 1, my_fun)))/2

cols <- ifelse(!is.na(c1) & !is.na(c2),
               ifelse(c1%%1 != 0 | c2%%1 != 0 |
                        (c1 == 1 & c2 %in% c(0,2)) |
                        (c2 == 1 & c1 %in% c(0,2)), 1/2,
                      ifelse(c1 == 1 & c2 == 1, 1, 0)),
               ifelse(!is.na(c1),
                      ifelse(c1%%1 == 0, ifelse(c1 == 1, 1, 0), 1/2),
                      ifelse(c2%%1 == 0, ifelse(c2 == 1, 1, 0), 1/2)
               )
)
cols <- factor(cols, levels = c(0, 1/2, 1),
               labels = c("Accept", "Agnostic", "Reject"))

meta_data <- data.frame(CIlower1 = my_meta$int1,
                        CIupper1 = my_meta$int2,
                        CIlower2 = my_meta$int3,
                        CIupper2 = my_meta$int4,
                        points = 1/-my_meta$mean,
                        color = cols,
                        study_class = factor(my_meta$study_class,
                                             levels = c(1,2,3),
                                             labels = c("Follow-up",
                                                        "Pharmacotherapy",
                                                        "Follow-up (T)")),
                        studlab = relevel(
                          factor(my_meta$studlab,
                                 levels = rev(sort(unique(my_meta$studlab)))),
                          "Pooled"
                        ))
meta_data %>%
  ggplot2::ggplot(ggplot2::aes(y = studlab, col = color)) +
  ggplot2::geom_errorbarh(aes(xmin = CIlower1, xmax = CIupper1), height=.2, linewdith = 1) +
  ggplot2::geom_errorbarh(aes(xmin = CIlower2, xmax = CIupper2), height=.2, linewdith = 1) +
  ggplot2::geom_point(ggplot2::aes(x = points)) +
  ggplot2::annotate('rect', xmin = -Inf, xmax = NNT[1],
                    ymin = 0, ymax = nrow(meta_data), alpha=.2, fill='dodgerblue3') +
  ggplot2::annotate('rect', xmin = NNT[2], xmax = Inf,
                    ymin = 0, ymax = nrow(meta_data), alpha=.2, fill='dodgerblue3') +
  ggplot2::annotate('rect', xmin = -1, xmax = 1,
                    ymin = 0, ymax = nrow(meta_data), alpha=.2, fill='red') +
  ggplot2::geom_vline(xintercept = NNT[1], linetype = 'dashed') +
  ggplot2::geom_vline(xintercept = NNT[2], linetype = 'dashed') +
  ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                       "Agnostic" = "goldenrod",
                                       "Reject" = "darkred")) +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "NNT-based REACT Forestplot",
                x = "NNT Estimate", y = "Study", color = "Decision") +
  coord_cartesian(xlim=c(-30, 30)) + facet_grid(~study_class, scales="free")
