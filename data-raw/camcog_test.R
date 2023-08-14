library(dplyr)
library(ggplot2)
library(ggforce)
library(latex2exp)

# importing data (CG - Control, DA - Alzheimer, CCL - MCI?)
camcog <- read.csv("data-raw/CAMCOG.csv")
camcog |> count(Diagnostico)
my_labs <- (camcog |> count(Diagnostico))$Diagnostico

# different counts for each group of diagnostic
# computing means for each group
stats <- camcog |>
  group_by(Diagnostico) |>
  summarise(mu = mean(CAMCOG, na.rm = TRUE),
            sigma = var(CAMCOG, na.rm = TRUE)/n(),
            n = n())

# computing pooled variance
pooled_var <- stats |>
  mutate(new_sigma = sigma*(n - 1)/(sum(n) - 3)) |>
  pull(new_sigma) |> sum()

var_cov <- diag(stats$sigma, nrow = 3, ncol = 3) #1 - CCL, 2 - DA, 3 - GC
pooled_var_cov <- diag(pooled_var, nrow = 3, ncol = 3)
mu <- stats$mu
c <- 1
alpha <- 0.05
delta <- 15


# Testing one ellipsis ----------------------------------------------------
# testing ggplot2
new_cov <- var_cov[c(1, 2), c(1, 2)]
new_mu <- mu[c(1, 2)]
q <- qchisq(1 - alpha, nrow(var_cov))

eig <- Matrix::Schur(new_cov)
eigs = eig$EValues
# maxvalue lambda
idx <- which.max(eig$EValues)
vecs <- eig$Q
main_eigenvec <- eig$Q[,1]

# minor and major axis
a <- sqrt(eigs[1]*q)
b <- sqrt(eigs[2]*q)

angle = atan2(main_eigenvec[2], main_eigenvec[1])

# diferrences between two means to increase xlab and ylab
mu_dif <- abs(new_mu[1] - new_mu[2])

points <- DescTools::DrawEllipse(x = new_mu[1], y = new_mu[2],
                                 radius.x = a, radius.y = b,
                       rot = angle , plot = FALSE, nv = 4000)


uniroot(check_ellipse,
        c(min(points$x) - 10, max(points$x) + 10),
        extendInt = "yes",
        delta = delta, c1 = new_mu[1], c2 = new_mu[2], a = a, b = b,
        tol = 0.0001)



print(check_ellipse(x_cand, x_cand - delta,
                    c1 = new_mu[1], c2 = new_mu[2],
                    a = a, b = b))

print(check_ellipse(x_cand, x_cand + delta,
                    c1 = new_mu[1], c2 = new_mu[2],
                    a = a, b = b))

# points data
points_data <- data.frame(x = points$x,
                          y = points$y)

ggplot(aes(x = x, y = y), data = points_data) +
  geom_polygon(colour="black", fill=NA)+
  geom_abline(intercept = delta, linetype = "dashed")+
  geom_abline(intercept = -delta, linetype = "dashed")+
  coord_cartesian(xlim = c(new_mu[1] - mu_dif - delta, new_mu[1] + mu_dif + delta),
                  ylim = c(new_mu[2] - mu_dif - delta, new_mu[2] + mu_dif + delta))



# Multiple comparisons with ellipse  ----------------------------------------------------
names_list <- c("MCI", "AD", "CG")
# supposing heterocedastic behavior:
p <- REACT::m_comparisons(alpha = 0.05, nrow = 1, ncol = 3,
                     delta = delta, par = mu, f_matrix = var_cov)

# modifying labs
c <- 1
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    tex_xlab = TeX(paste0("$\\mu_{", names_list[i], "}$"))
    tex_ylab = TeX(paste0("$\\mu_{", names_list[j], "}$"))

    p[[c]] <- p[[c]] + labs(x = tex_xlab,
                  y = tex_ylab,
                  title = paste0(names_list[i], " vs ", names_list[j]),
                  colour = "Decision",
                  fill = "Decision")
    c <- c + 1
  }}

ggpubr::ggarrange(plotlist = p, nrow = 1, ncol = 3, common.legend = TRUE)

# supposing homocedastic behavior:
p <- REACT::m_comparisons(alpha = 0.05, nrow = 1, ncol = 3,
                          delta = delta, par = mu, f_matrix = pooled_var_cov)

# modifying labs
c <- 1
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    tex_xlab = TeX(paste0("$\\mu_{", names_list[i], "}$"))
    tex_ylab = TeX(paste0("$\\mu_{", names_list[j], "}$"))

    p[[c]] <- p[[c]] + labs(x = tex_xlab,
                            y = tex_ylab,
                            title = paste0(names_list[i], " vs ", names_list[j]),
                            colour = "Decision",
                            fill = "Decision")
    c <- c + 1
  }}

ggpubr::ggarrange(plotlist = p, nrow = 1, ncol = 3, common.legend = TRUE)






# Multiple comparisons with CI's ------------------------------------------
# function to create the CI graph in article
# auxliary functions
compute_CI <- function(obs1, obs2, alpha){
  sd_tot <- sqrt((var(obs1)/length(obs1)) + (var(obs2)/length(obs2)))
  delta_mu <- mean(obs1) - mean(obs2)
  ci <- c(delta_mu - (qnorm(1 - (alpha/2)) * sd_tot),
          delta_mu + (qnorm(1 - (alpha/2))*sd_tot))
  return(ci)
}

dec_func <- function(CI_fix, delta_used){
  res <- REACT::base_test(CI = CI_fix, tol = delta_used, hyp = 0,verbose = FALSE)
  final <- ifelse(res$test_outcome == "remain agnostic", 1/2,
                  ifelse(res$test_outcome == "accept", 0, 1))
  return(final)
}

plot_CI <- function(camcog_data, names_list = c("MCI", "AD", "CG"),
                    seed = 125, alpha = 0.05, delta = 15){
  # reordering data according to names list
  set.seed(seed)
  camcog_data <- camcog_data |>
    na.omit() |>
    mutate(Diagnostico = factor(Diagnostico)) |>
    mutate(Diagnostico = recode_factor(Diagnostico,
                                       CCL = "MCI", DA = "AD", GC = "CG")) |>
    arrange(match(Diagnostico, names_list))

  # function to compute CI in the camcog context

  plot_list <- list()
  c <- 1
  for(i in 1:(length(names_list) - 1)){
    for(j in (i + 1):length(names_list)){
      name_1 <- names_list[i]
      name_2 <- names_list[j]
      n <- camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        nrow()

      # removing the two first entries from first group
      group_1_obs <- camcog_data |> filter(Diagnostico  == name_1) |>
        slice(1:2) |> pull(CAMCOG)

      # removing the two last entries from second group
      group_2_obs <- camcog_data |> filter(Diagnostico  == name_2) |>
        slice_tail(n = 2) |> pull(CAMCOG)

      # making new data
      new_data <- camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        slice(-c(1,2, n - 1, n))

      # shuffling rows
      new_data <- new_data |> sample_frac()

      # total number of sampling to do
      tot_n <- n - 4

      # vector of CI_upper and lower
      CI_mat <- matrix(nrow = tot_n + 1, ncol = 2)

      # computing variance CI lower and upper by data
      ci <- compute_CI(group_1_obs, group_2_obs, alpha = alpha)
      CI_mat[1,] = ci

      # looping through all data and computing each CI
      decisions <- numeric(nrow(CI_mat))
      decisions[1] <- dec_func(CI_mat[1,], delta)

      for(k in 1:tot_n){
        current_data <- new_data |> slice(k)
        if(current_data$Diagnostico == name_1){
          group_1_obs <- c(group_1_obs, current_data$CAMCOG)
        }else{
          group_2_obs <- c(group_2_obs, current_data$CAMCOG)
        }
        CI_mat[k + 1, ] <-  compute_CI(group_1_obs, group_2_obs, alpha)
        decisions[k + 1] <- dec_func(CI_fix = CI_mat[k + 1, ], delta = delta)
      }
      colnames(CI_mat) <- c("lower", "upper")

      cols <- factor(decisions,
                     levels = c(0, 1/2, 1),
                     labels = c("Accept", "Agnostic", "Reject"))

      data_final <- data.frame(
        x = 4:(nrow(CI_mat) + 3),
        y_min = CI_mat[,1],
        y_max = CI_mat[, 2],
        dec = cols)
      # plotting
      plot_list[[c]] <-  ggplot(
        aes(x = x, ymin = y_min, ymax = y_max, colour = dec, fill = dec),
        data = data_final) +
        geom_errorbar(width = 1,
                      position = position_dodge(0.05))+
        annotate('rect', xmin = -0.5, xmax = 100.5,
                          ymin = -delta, ymax = delta,
                 alpha=.2, fill='dodgerblue3')+
        theme_bw()+
        geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                             "Agnostic" = "goldenrod",
                                             "Reject" = "darkred"),
                                    drop = FALSE)+
      ggplot2::scale_fill_manual(values=c("Accept" = "darkgreen",
                                           "Agnostic" = "goldenrod",
                                           "Reject" = "darkred"),
                                  drop = FALSE)+
      labs(x = "Sample size",
           y = TeX(paste0("$\\mu_{", name_1, "} - \\mu_{", name_2, "}$")),
           title = paste0(name_1, " vs ", name_2),
           fill = "Decision",
           colour = "Decision") +
      coord_cartesian(xlim = c(0, 100),
                      ylim = c(-35, 35))

      c <- c + 1
    }
  }
  ggpubr::ggarrange(plotlist = plot_list, nrow = 1, ncol = 3, common.legend = TRUE)
}

p <- plot_CI(camcog)
