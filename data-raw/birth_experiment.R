#Context: Baby born underweight given potentially risky behavior from the mother
#Objective: Test the relative risk of each factor to identify if any of them
#           may not actually be harmful to baby's weight
#Caveats: Extra care required to accept H0: RR = 1, region of equivalence has to
#         be more restrictive (NNT = 100), CER unknown (setting CER = 1 to get a
#         more conservative region).


library(dplyr)
library(ggplot2)
library(lbreg)

# importing birth data
data(Birth)
# transforming all into factor
Birth <- Birth |>
  mutate(alc = as.factor(alc),
         smo = as.factor(smo),
         soc = as.factor(soc)) |>
  mutate(alc = as.factor(case_when(
    alc == 1 | alc == 2 ~ 0,
    alc == 3 ~ 1
  )))

# fixing CER as proposed in literature
CER <- 0.1

# fixing NNT as also proposed
NNT <- 100


# fitting model to BIrth data
mod <- lbreg(lowbw ~ ., data = Birth)

# obtaing relative risk CI's
rr_obj <- relrisk(mod, alpha = 0.05)[-1, ] #Intercept removal

# delimiting tolerance by RRR
rrr <- 1/(NNT*CER)
hyp <- 1
p_int <- c(hyp - rrr, hyp + rrr)

# plotting each covariate results
idxs <- matrix(findInterval(rr_obj[,2:3],
                            p_int), ncol = 2)

cols <- factor(ifelse(idxs[,1] == 1 & idxs[,2] == 1, 0,
                      ifelse((idxs[,1] == 0 & idxs[,2] == 0) |
                               (idxs[,1] == 2 & idxs[,2] == 2), 1,
                             1/2)),
               levels = c(0, 1/2, 1),
               labels = c("Accept", "Agnostic", "Reject"))

mod_data <- data.frame( CIlower = rr_obj[,2],
                        CIupper = rr_obj[,3],
                        points = rr_obj[,1],
                        color = cols,
                        cov = factor(rownames(rr_obj),
                                         levels = rev(sort(unique(rownames(rr_obj))))))

# arrows limits
arr_inf <- min(c(rr_obj[, 2], p_int[1]))
arr_sup <- max(c(rr_obj[, 3], p_int[2]))

# adding annotations
p <- mod_data %>%
  ggplot2::ggplot(ggplot2::aes(y = cov, xmin = CIlower, xmax = CIupper, col = color)) +
  ggplot2::geom_errorbarh(ggplot2::aes(height=.2)) +
  {if(!all(is.na(mod_data$points))) ggplot2::geom_point(ggplot2::aes(x = points))} +
  ggplot2::annotate('rect', xmin = p_int[1], xmax = p_int[2],
                    ymin = 0, ymax = nrow(mod_data), alpha=.2, fill='dodgerblue3') +
  ggplot2::geom_vline(xintercept = hyp, linetype = 'dashed') +
  ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                       "Agnostic" = "goldenrod",
                                       "Reject" = "darkred")) +
  ggplot2::annotate("text", x = p_int[1] - 0.15,
                    y = -0.75, label = "Disfavors Treatment",
                    colour = "orangered")+
  ggplot2::annotate("segment", x = hyp - 0.05, xend = arr_inf,
                    y = -0.5, yend = -0.5, colour = "orangered", size = 1, alpha=0.6,
                    arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
  ggplot2::annotate("text", x = p_int[2] + 0.15,
                    y = -0.75, label = "Favors Treatment",
                    colour = "slateblue4")+
  ggplot2::annotate("segment", x = hyp + 0.05, xend = arr_sup,
                    y = -0.5, yend = -0.5, colour = "slateblue4", size = 1, alpha=0.6,
                    arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
  ggplot2::theme_bw() +
  ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 35)))+
  ggplot2::coord_cartesian(ylim=c(0 + 0.5, nrow(mod_data) - 0.5),
                           clip = "off") +
  ggplot2::labs(title = "RRR-based REACT Procedure",
                x = "Relative Risk", y = "Covariate", color = "Decision")

methods::show(p)

################################################################################

library(lbreg)
data(Birth)

mod <- lbreg(lowbw ~ ., data = Birth)
RR_CI <- relrisk(mod, alpha = 0.05)[-1, ]

cov_names <- c("Alcohol Consumption", "Smoking", "Social Class")

REACT_RRplot(CI_matrix = RR_CI[,2:3], NNT = 100,
             point_estim = RR_CI[,1], covar_names = cov_names)
