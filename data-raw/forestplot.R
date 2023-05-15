library(dplyr)
library(ggplot2)

meta <- read.csv("data-raw/metanalise.csv", header = T)
View(meta)

meta_fil <- meta %>%
  filter(ne != 0 & nc != 0) %>%
  distinct() %>%
  select(studlab, te, lowerte, upperte)

View(meta_fil)

meta_fil %>%
  mutate(color = ifelse(apply(cbind(meta_fil$lowerte, meta_fil$upperte),
                              1, function(x) findInterval(1, x)) == 1,
         "contains", "does not contain")) %>%
  ggplot(aes(y = studlab, x = te, xmin = lowerte, xmax = upperte, col = color)) +
  geom_point() + geom_errorbarh(height=.1) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_color_manual(values=c("green", "red"))


NNT_test <- function(confidence, NNT, event_e, n_e, event_c, n_c,
                     plot = F, verbose = T){
  epsilon <- c(-1/NNT, 1/NNT)
  prob_e <- event_e/n_e
  prob_c <- event_c/n_c
  var_e <- prob_e*(1-prob_e)/n_e
  var_c <- prob_c*(1-prob_c)/n_c
  prob_diff <- event_e/n_e - event_c/n_c
  v_diff <- var_e + var_c
  CI <- prob_diff +
    qnorm(c((1-confidence)/2, 1-(1-confidence)/2))*sqrt(v_diff)
  idxs <- findInterval(CI, epsilon)
  test_outcome <- ifelse(all(idxs == 1), "accept",
                         ifelse(all(idxs == 0) | all(idxs == 2),
                                "reject", "remain agnostic"))
  if(verbose){
    cat("REACT results:\n")
    cat("Pragmatic lower bound: ", format(epsilon[1], digits = 3, nsmall = 2, scientific = FALSE))
    cat("\n")
    cat("Pragmatic upper bound: ", format(epsilon[2], digits = 3, nsmall = 2, scientific = FALSE))
    cat("\n")
    cat("Confidence interval:")
    cat("\n")
    cat("lower bound: ", paste0(round(CI[1], digits = 3)),
        "\nupper bound: ",paste0(round(CI[2],digits = 3)), sep = "")
    cat("\n")
    cat("REACT conclusion:\n")
    message("Based on the provided confidence interval we ", test_outcome,
            ifelse(test_outcome == "remain agnostic", ".", " the null hypothesis."))
  }
  if(plot){
    p <- ggplot2::ggplot()+
      ggplot2::coord_cartesian(ylim=c(-0.1, 0.1)) +
      ggplot2::geom_segment(ggplot2::aes(x = CI[1], y = 0,
                                         xend = CI[2], yend = 0), linewidth = 0.75)+
      ggplot2::annotate('rect', xmin = epsilon[1], xmax = epsilon[2],
                        ymin = -0.05, ymax = 0.05, alpha=.3, fill='blue')+
      ggplot2::theme_bw()+
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank())+
      ggplot2::labs(y = "",
                    x = "Parameter values",
                    title = "Hypothesis testing acceptance region and confidence interval")

    show(p)

  }

  return(list(CI = CI,
              pragmatic = epsilon,
              test_outcome = test_outcome))
}

NNT_test(confidence = 0.95, NNT = 3,
         event_e = 77, n_e = 257, event_c = 8, n_c = 47,
         plot = T)

alpha <- .05
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

