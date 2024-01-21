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
  na.omit() |>
  group_by(Diagnostico) |>
  summarise(mu = mean(CAMCOG, na.rm = TRUE),
            sigma = var(CAMCOG, na.rm = TRUE)/n(),
            sum_squared = sum(CAMCOG^2, na.rm = TRUE),
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
                     tol = delta, par = mu, f_matrix = var_cov)

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

      new_camcog_data <- camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        sample_frac() |>
        mutate(idx = 1:n)

      # removing the two first entries from first group
      group_1_obs_info <- new_camcog_data |> filter(Diagnostico == name_1) |>
        slice(1:2)

      # removing the two last entries from second group
      group_2_obs_info <- new_camcog_data |> filter(Diagnostico  == name_2) |>
        slice(1:2)

      # making new data
      new_data <- new_camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        slice(-c(group_1_obs_info$idx, group_2_obs_info$idx))

      group_1_obs <- group_1_obs_info$CAMCOG
      group_2_obs <- group_2_obs_info$CAMCOG

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

plot_CI(camcog, seed = 785)


plot_multiple_CI <- function(camcog_data, nsim = 5, names_list = c("MCI", "AD", "CG"),
                             seed = 125, alpha = 0.05, delta = 15){
  # reordering data according to names list
  set.seed(seed)
  seeds_list <- sample(1:10^6, nsim, replace = FALSE)
  plot_list <- list()

  # function to compute CI in the camcog context
  c <- 1
  for(i in 1:(length(names_list) - 1)){
    for(j in (i + 1):length(names_list)){
      l_pos <- 1
      data_list <- list()
      # changing shuffling

      for(seed in seeds_list){
        set.seed(seed)
        camcog_data <- camcog_data |>
          na.omit() |>
          mutate(Diagnostico = factor(Diagnostico)) |>
          mutate(Diagnostico = recode_factor(Diagnostico,
                                             CCL = "MCI", DA = "AD", GC = "CG")) |>
          arrange(match(Diagnostico, names_list)) |>
          sample_frac()

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

      # shuffling rows again
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

      data_list[[l_pos]] <- data.frame(
        x = 4:(nrow(CI_mat) + 3),
        y_min = CI_mat[,1],
        y_max = CI_mat[, 2],
        dec = cols,
        sim_id = as.character(l_pos))
      l_pos <- l_pos + 1
      }
      data_final <- do.call(rbind, data_list)
      # plotting
      plot_list[[c]] <-  ggplot(
        aes(x = x, ymin = y_min, ymax = y_max,
            colour = dec, fill = dec, linetype = sim_id),
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
             colour = "Decision",
             linetype = "Simulation") +
        coord_cartesian(xlim = c(0, 100),
                        ylim = c(-35, 35))

      c <- c + 1
    }
  }
  ggpubr::ggarrange(plotlist = plot_list, nrow = 1, ncol = 3, common.legend = TRUE)
}

plot_multiple_CI_separated <- function(camcog_data, nsim = 5, names_list = c("MCI", "AD", "CG"),
                             seed = 125, alpha = 0.05, delta = 15){
  # reordering data according to names list
  set.seed(seed)
  seeds_list <- sample(1:10^6, nsim, replace = FALSE)
  plot_list <- list()

  # function to compute CI in the camcog context
  c <- 1
  l_pos <- 1
  for(seed in seeds_list){
    set.seed(seed)
  for(i in 1:(length(names_list) - 1)){
    for(j in (i + 1):length(names_list)){
      # changing shuffling
      name_1 <- names_list[i]
      name_2 <- names_list[j]
      camcog_data <- camcog_data |>
        na.omit() |>
        mutate(Diagnostico = factor(Diagnostico)) |>
        mutate(Diagnostico = recode_factor(Diagnostico,
                                           CCL = "MCI", DA = "AD", GC = "CG")) |>
        arrange(match(Diagnostico, names_list))

      n <- camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        nrow()

      new_camcog_data <- camcog_data |>
        filter(Diagnostico %in% c(name_1, name_2)) |>
        sample_frac() |>
        mutate(idx = 1:n)

      # removing the two first entries from first group
      group_1_obs_info <- new_camcog_data |> filter(Diagnostico == name_1) |>
        slice(1:2)

      # removing the two last entries from second group
      group_2_obs_info <- new_camcog_data |> filter(Diagnostico  == name_2) |>
        slice(1:2)

        # making new data
        new_data <- new_camcog_data |>
          filter(Diagnostico %in% c(name_1, name_2)) |>
          slice(-c(group_1_obs_info$idx, group_2_obs_info$idx))

        group_1_obs <- group_1_obs_info$CAMCOG
        group_2_obs <- group_2_obs_info$CAMCOG

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

        values <- decisions

        data_final <- data.frame(
          x = 4:(nrow(CI_mat) + 3),
          y_min = CI_mat[,1],
          y_max = CI_mat[, 2],
          dec = cols,
          dec_values = values)

      # plotting
      plot_list[[c]] <-  ggplot(
        aes(x = x, ymin = y_min, ymax = y_max,
            colour = dec, fill = dec),
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
             title = paste0(name_1, " vs ", name_2, " (resampling ", l_pos, ")"),
             fill = "Decision",
             colour = "Decision",
             linetype = "Simulation") +
        coord_cartesian(xlim = c(0, 100),
                        ylim = c(-35, 35))

      c <- c + 1
    }
  }
  l_pos <- l_pos + 1
  }
  ggpubr::ggarrange(plotlist = plot_list, nrow = nsim, ncol = 3, common.legend = TRUE)
}

plot_multiple_CI_separated(camcog, nsim = 3, seed = 750)


# Bayesian CAMCOG ---------------------------------------------------------
# using jeffreys prior, we obtain that for each mu_k:
# mu_k|\bar{X}_k ~ t_{n_k - 1}(\bar{X}_k, S^2/n_k), that is mu_K - \bar{X}/(S/sqrt{n_k}) ~ t_{n_k - 1}
# using stats values
n <- stats$n - 1
sd <- stats$sigma |> sqrt()
xbars <- stats$mu
names_list <- c("MCI", "AD", "CG")
B <- 10^5

# generating several t's
generate_HPD_data <- function(df_s, sd, xbars, B,
                              names_list = c("MCI", "AD", "CG"), alpha = 0.05, seed = 250){
  # generating t samples and computing densities
  set.seed(seed)
  sim_data <- names_list |>
    purrr::map_dfc(B = B, function(x, B){
      idx <- which(names_list == x)
      n_used <- n[idx]
      # generating t's
      sample <- rt(B, df = n_used)
      # computing dt's
      densities <- dt(sample, n_used)
      # joining in dataframe
      colname_1 <- x
      colname_2 <- paste0(x, "_dens")
      data <- data.frame(var1 = sample,
                 var2 = densities)
      colnames(data) <- c(colname_1, colname_2)
      return(data)
    }) |>
    rowwise() |>
    mutate(dens_total = prod(c_across(ends_with("dens")))) |> # generating the total densities for each sample
    ungroup() |>
    filter(dens_total >= quantile(dens_total, alpha)) |> # selecting the 1-alpha highest densities
    select(names_list)

  data_final <- names_list |>
    purrr::map_dfc(function(x){
      idx <- which(names_list == x)
      temp_data <- sim_data |> select(idx)
      (temp_data*sd[idx]) + xbars[idx]
    })
  return(data_final)
}

sim_data <- generate_HPD_data(n, sd, xbars, B)

# multiple comparisons using jeffreys prior
m_comparisons_bayes <- function(df_s, sd, xbars, B, tol,
                                names_list = c("MCI", "AD", "CG"),
                                alpha = 0.05, seed = 250, nrow = 1, ncol = 3){
  # generating HPD data for posterior
  hpd_data <- generate_HPD_data(df_s, sd, xbars, B, names_list, alpha, seed)

  # generating HPD data for prior
  c <- 1
  plot_list <- list()

  for(i in 1:(length(names_list) - 1)){
    for(j in (i + 1):length(names_list)){
    current_data <- hpd_data |> select(c(i, j))

    # points in convex hull
    ch_points <- current_data |>
      slice(chull(current_data |> pull(1),
          current_data |> pull(2)))

    # check convex hull
    check_convex_hull <- function(point_x, point_y, tol){
      abs(point_x - point_y) <= tol
    }
    ch_checking <- check_convex_hull(ch_points[, 1], ch_points[, 2], tol = tol)

    # obtaining convex hull point indexes
    res <- ifelse(all(ch_checking == TRUE), 0,
                  ifelse(all(ch_checking == FALSE), 1, 1/2))

    # transforming res into factor
    cols <- factor(res,
                   levels = c(0, 1/2, 1),
                   labels = c("Accept", "Agnostic", "Reject"))


    # plotting based on res
    # points data
    tex_title = latex2exp::TeX(paste0("$\\mu_{", names_list[i],  "}$ and $\\mu_{", names_list[j], "}$"))
    tex_xlab = latex2exp::TeX(paste0("$\\mu_{", names_list[i], "}$"))
    tex_ylab = latex2exp::TeX(paste0("$\\mu_{", names_list[j], "}$"))

    data_used <- data.frame(x = ch_points[, 1],
                            y = ch_points[, 2],
                            col = cols)

    # building data to construct pragmatic region
    x_min <- min(data_used$x) - tol
    x_max <- max(data_used$x) + tol

    # adding y to each x between x_min and x_max
    y_1 <- seq(x_min - tol, x_max - tol, length.out = 300)
    y_2 <- seq(x_min + tol, x_max + tol, length.out = 300)
    x <- seq(x_min, x_max, length.out = 300)
    prag_data <- data.frame(x = x,
                            y_1 = y_1,
                            y_2 = y_2)

    plot_list[[c]] <- ggplot2::ggplot(
      ggplot2::aes(x = x, y = y, fill = col, colour = col), data = data_used) +
      ggplot2::geom_polygon(alpha = 0.4)+
      ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                           "Agnostic" = "goldenrod",
                                           "Reject" = "darkred"),
                                  drop = FALSE) +
      ggplot2::scale_fill_manual(values=c("Accept" = "darkgreen",
                                          "Agnostic" = "goldenrod",
                                          "Reject" = "darkred"),
                                 drop = FALSE)+
      ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_1, ymax = y_2),
                           data = prag_data, inherit.aes = FALSE,
                           fill = "dodgerblue3", alpha = 0.2)+
      ggplot2::theme_bw()+
      ggplot2::labs(x = tex_xlab,
                    y = tex_ylab,
                    title = tex_title,
                    colour = "Decision",
                    fill = "Decision")
    c <- c + 1
    }
  }
  p <- ggpubr::ggarrange(plotlist = plot_list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  methods::show(p)

  invisible(plot_list)
}

# generating objective bayes camcog image
delta <- 15
m_comparisons_bayes(n, sd, xbars, B = 500, tol = delta)

# using NGI prior
# obtaining posterior parameters
# generating new df_s, sd and xbars based on this prior
generate_post_par <- function(theta, lambda, alpha, beta, stats){
  # generating new df
  new_dfs <- (2*alpha + stats$n)

  # new mean
  new_xbars <- ((stats$n*stats$mu) + (lambda*theta))/(lambda + stats$n)
  delta <- (beta + (stats$sum_squared/2) + (lambda*(theta^2))/2 - (
    (stats$n*stats$mu + (lambda*theta))^2/(2*(lambda + stats$n))))
  print(delta)
  new_sds <- sqrt(2*delta/(
    (lambda + stats$n)*(2*alpha + (stats$n + 1))
  ))

  return(list(dfs = new_dfs, sd = new_sds, xbars = new_xbars))
}

# HPD for prior
generate_prior_HPD_data <- function(theta, lambda, alpha_par, beta,
                                    names_list = c("MCI", "AD", "CG"), alpha = 0.05, seed = 250){
  # generating t samples and computing densities
  set.seed(seed)
  sim_data <- names_list |>
    purrr::map_dfc(B = B, function(x, B){
      # generating t's
      sample <- rt(B, df = 2*alpha_par)
      # computing dt's
      densities <- dt(sample, 2*alpha_par)
      # joining in dataframe
      colname_1 <- x
      colname_2 <- paste0(x, "_dens")
      data <- data.frame(var1 = sample,
                         var2 = densities)
      colnames(data) <- c(colname_1, colname_2)
      return(data)
    }) |>
    rowwise() |>
    mutate(dens_total = prod(c_across(ends_with("dens")))) |> # generating the total densities for each sample
    ungroup() |>
    filter(dens_total >= quantile(dens_total, alpha)) |> # selecting the 1-alpha highest densities
    select(names_list)

  sd <- sqrt((2*beta)/lambda)
  xbar <- theta
  data_final <- names_list |>
    purrr::map_dfc(function(x){
      idx <- which(names_list == x)
      temp_data <- sim_data |> select(idx)
      (temp_data*sd) + xbar
    })

  return(data_final)
}

# multiple comparisons using NGI prior
m_comparisons_bayes_NGI <- function(theta, lambda, alpha_par, beta, stats, B, tol,
                                names_list = c("MCI", "AD", "CG"),
                                alpha = 0.05, seed = 250, nrow = 1, ncol = 3){
  # posterior parameters
  post_par <- generate_post_par(theta, lambda, alpha_par, beta, stats)
  df_s <- post_par$dfs
  sd <- post_par$sd
  xbars <- post_par$xbars

  # generating HPD data for posterior
  hpd_data <- generate_HPD_data(df_s, sd, xbars, B, names_list, alpha, seed)

  # generating HPD data for prior
  df_p <- 2*alpha_par
  hpd_prior_data <- generate_prior_HPD_data(theta, lambda, alpha_par, beta,
                                            names_list, alpha, seed)
  c <- 1
  plot_list <- list()

  for(i in 1:(length(names_list) - 1)){
    for(j in (i + 1):length(names_list)){
      current_data <- hpd_data |> select(c(i, j))
      current_prior_data <- hpd_prior_data |> select(c(i, j))
      plot(current_prior_data |> pull(1), current_prior_data |> pull(2))

      # points in convex hull
      ch_points <- current_data |>
        slice(chull(current_data |> pull(1),
                    current_data |> pull(2)))

      ch_points_prior <- current_prior_data |>
        slice(chull(current_prior_data |> pull(1),
                    current_prior_data |> pull(2)))


      # check convex hull
      check_convex_hull <- function(point_x, point_y, tol){
        abs(point_x - point_y) <= tol
      }
      ch_checking <- check_convex_hull(current_data[, 1], current_data[, 2], tol = tol)
      ch_checking_prior <- check_convex_hull(current_prior_data[, 1], current_prior_data[, 2], tol = tol)

      # obtaining convex hull point indexes
      res <- ifelse(all(ch_checking == TRUE), 0,
                    ifelse(all(ch_checking == FALSE), 1, 1/2))

      res_prior <- ifelse(all(ch_checking_prior == TRUE), 0,
                          ifelse(all(ch_checking_prior == FALSE), 1, 1/2))

      # transforming res into factor
      cols <- factor(res,
                     levels = c(0, 1/2, 1),
                     labels = c("Accept", "Agnostic", "Reject"))

      cols_prior <- factor(res_prior,
                           levels = c(0, 1/2, 1),
                           labels = c("Accept", "Agnostic", "Reject"))


      # plotting based on res
      # points data
      tex_title = latex2exp::TeX(paste0("$\\mu_{", names_list[i],  "}$ and $\\mu_{", names_list[j], "}$"))
      tex_xlab = latex2exp::TeX(paste0("$\\mu_{", names_list[i], "}$"))
      tex_ylab = latex2exp::TeX(paste0("$\\mu_{", names_list[j], "}$"))

      data_used <- data.frame(x = ch_points[, 1],
                              y = ch_points[, 2],
                              col = cols,
                              type = "posterior")

      data_used_prior <- data.frame(x = ch_points_prior[, 1],
                                    y = ch_points_prior[, 2],
                                    col = cols_prior,
                                    type = "prior")
      # joining data to plot
      all_data_used <- bind_rows(data_used, data_used_prior)

      # building data to construct pragmatic region
      x_min <- min(all_data_used$x) - tol
      x_max <- max(all_data_used$x) + tol

      # adding y to each x between x_min and x_max
      y_1 <- seq(x_min - tol, x_max - tol, length.out = 300)
      y_2 <- seq(x_min + tol, x_max + tol, length.out = 300)
      x <- seq(x_min, x_max, length.out = 300)
      prag_data <- data.frame(x = x,
                              y_1 = y_1,
                              y_2 = y_2)

      plot_list[[c]] <- ggplot2::ggplot(
        ggplot2::aes(x = x, y = y, fill = col, colour = col, linetype = type), data = all_data_used) +
        ggplot2::geom_polygon(alpha = 0.4)+
        ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                             "Agnostic" = "goldenrod",
                                             "Reject" = "darkred"),
                                    drop = FALSE) +
        ggplot2::scale_fill_manual(values=c("Accept" = "darkgreen",
                                            "Agnostic" = "goldenrod",
                                            "Reject" = "darkred"),
                                   drop = FALSE)+
        ggplot2::scale_linetype_manual(values=c("prior" = 2,
                                                "posterior" = 1),
                                       guide = guide_legend(override.aes = list(colour = "black",
                                                                                fill = NA)))+
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_1, ymax = y_2),
                             data = prag_data, inherit.aes = FALSE,
                             fill = "dodgerblue3", alpha = 0.2)+
        ggplot2::theme_bw()+
        ggplot2::labs(x = tex_xlab,
                      y = tex_ylab,
                      title = tex_title,
                      colour = "Decision",
                      linetype = "Distribution",
                      fill = "Decision")
      c <- c + 1
    }
  }
  p <- ggpubr::ggarrange(plotlist = plot_list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  methods::show(p)

  invisible(plot_list)
}

# posterior parameters for t distribution
prior_par <- c(80, 1, 2, 2)
delta <- 15
B <- 10^5
m_comparisons_bayes_NGI(prior_par[1], prior_par[2], prior_par[3], prior_par[4],
                    stats, B = B, tol = delta)
# result similar to frequentist and prior is more concentrated

# increasing lambda to obtain an even more concentrated prior
prior_par <- c(80, 2.5, 2, 2)
m_comparisons_bayes_NGI(prior_par[1], prior_par[2], prior_par[3], prior_par[4],
                        stats, B = B, tol = delta)

# increasing alpha and beta and decreasing lambda to obtain more assertive posteriors but wide priors
prior_par <- c(0, 0.1, 50, 50)
m_comparisons_bayes_NGI(prior_par[1], prior_par[2], prior_par[3], prior_par[4],
                        stats, B = B, tol = delta)

# using flat prior but with balance in prior concentration, increasing variance with lambda less than zero
prior_par <- c(80, 0.5, 1, 1)
B <- 10^5
m_comparisons_bayes_NGI(prior_par[1], prior_par[2], prior_par[3], prior_par[4],
                        stats, B = B, tol = delta)

# trying to obtain a smoother region
prior_par <- c(80, 0.25, 3, 3)
B <- 10^5
m_comparisons_bayes_NGI(prior_par[1], prior_par[2], prior_par[3], prior_par[4],
                        stats, B = B, tol = delta)
