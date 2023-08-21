# REACT glm tests and simulations
library(REACT)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

# testing REACT glm with DL95 example -------------------------------------
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

# comparing basic glm with summary and react glm
base_glm_D93 <- glm(counts ~ outcome + treatment, family = poisson())
summary(base_glm_D93)

react_glm_D93 <- REACT.glm(counts ~ outcome + treatment, family = poisson(),
                           tol = 0.5, alpha = 0.05, verbose = TRUE)

# stepAIC
step_mod <- step(base_glm_D93)
summary(step_mod)


# Poisson simulation experiment with pre-estabilished significant covariates  --------
# generating independent covariates from uniform(-1,1)
glm_sim <- function(n, p, k = 5, type = "logistic", seed = 750, eps_beta = 0.005){
  set.seed(seed)
  # covariate matrix
  X <- setNames(1:p, paste0("V",1:p)) |> map_dfc(function(x){runif(n, -1, 1)})

  # generating beta from significant variable from uniform 1 to 2
  beta_sig <- runif(k, min = 2, max = 3)

  # generating beta from non-significant variable from uniform -0.05 to 0.05
  beta_noise <- runif(p - k, min = -eps_beta, max = eps_beta)

  all_beta <- c(beta_sig, beta_noise)

  if(type == "logistic"){
    linear_comp <- as.numeric(as.matrix(X) %*% all_beta)
    p_x <- 1/(1 + exp(- linear_comp))
    return(list(data = cbind(X, y = rbinom(n, 1, p_x)), true_betas = all_beta))
    }
}

n <- 1000
p <- 10
# less noisy scenario
logis_sim <- glm_sim(n, p, k = 5)

# defatul glm
glm_obj <- glm(y ~ ., family = binomial, data = logis_sim$data)

# summary and step AIC
summary(glm_obj)
# step AIC
step(glm_obj)

# REACT
tol <- 0.5
alpha <- 0.05
react_obj <- REACT.glm(y ~ ., family = binomial, data = logis_sim$data,
          tol = tol, alpha = alpha, verbose = TRUE)

# repeating analysis but more noisy
n <- 1000
p <- 10
# less noisy scenario
logis_sim <- glm_sim(n, p, k = 1)

# defatul glm
glm_obj <- glm(y ~ ., family = binomial, data = logis_sim$data)

# summary and step AIC
summary(glm_obj)
# step AIC
step(glm_obj)

# REACT
tol <- 0.5
alpha <- 0.05
react_obj <- REACT.glm(y ~ ., family = binomial, data = logis_sim$data,
                       tol = tol, alpha = alpha, verbose = TRUE)

# assimilating beta = 0 to noise covariates
n <- 750
p <- 10
# less noisy scenario
logis_sim <- glm_sim(n, p, k = 5, eps_beta = 0, seed = 250)

# defatul glm
glm_obj <- glm(y ~ ., family = binomial, data = logis_sim$data)

# summary and step AIC
summary(glm_obj)
# step AIC
step_obj <- step(glm_obj, trace = FALSE)

# REACT
tol <- 0.5
alpha <- 0.05
react_obj <- REACT.glm(y ~ ., family = binomial, data = logis_sim$data,
                       tol = tol, alpha = alpha, verbose = TRUE)

# More intensive monte carlo experiments ----------------------------------
mc_simulation <- function(B, n, p, k = 5, type = "logistic", master_seed = 750, eps_beta = 0.005,
                          alpha = 0.05, tol = 0.5){
  set.seed(master_seed)
  seeds <- sample(1:10^5, size = B, replace = FALSE)
  all_res <- list(anova = as.data.frame(matrix("", nrow = B, ncol = p+1)),
                  stepAIC = as.data.frame(matrix("", nrow = B, ncol = p+1)),
                  react = as.data.frame(matrix("", nrow = B, ncol = p+1)))

  for(i in 1:B){
    sim_obj <- glm_sim(n = n, p = p, k = k, type = type, eps_beta = eps_beta, seed = seeds[i])

    # defatul glm
    glm_obj <- glm(y ~ ., family = binomial, data = sim_obj$data)

    names_list <- coef(glm_obj) |> names()

    # extracting anova decisions by summary
    p_values <- summary(glm_obj)[["coefficients"]][,4]
    all_res$anova[i, which(p_values <= alpha)] <- "significant"
    all_res$anova[i, -which(p_values <= alpha)] <- "non significant"

    # extracting stepwise decisions
    names_step <- step(glm_obj, trace = FALSE) |>
      coef() |>
      names()
    all_res$stepAIC[i, which(names_list %in% names_step)] <- "significant"
    all_res$stepAIC[i, -which(names_list %in% names_step)] <- "non significant"

    # extracting react decisions
    react_obj <- REACT.glm(y ~ ., family = binomial, data = sim_obj$data,
                           tol = tol, alpha = alpha, verbose = FALSE)
    all_res$react[i, which(names_list %in% react_obj$signif_vars)] <- "significant"
    all_res$react[i, which(names_list %in% react_obj$non_signif_vars)] <- "non significant"
    all_res$react[i, which(names_list %in% react_obj$agnostic_vars)] <- "agnostic"
  }
  # obtaining final data
  bind_rows(all_res, .id = "method")
}

# two simulation results
sim_res_tol_05 <- mc_simulation(1000, 1000, 10, eps_beta = 0, tol = 0.5)
sim_res_tol_065 <- mc_simulation(1000, 1000, 10, eps_beta = 0, tol = 0.65)

names(sim_res_tol_05)[-1] <- c("(Intercept)", paste0("V", 1:p))
names(sim_res_tol_065)[-1] <- c("(Intercept)", paste0("V", 1:p))

sim_res_tol_05 <- sim_res_tol_05 |> mutate(tol = "Tolerance = 0.5")
sim_res_tol_065 <- sim_res_tol_065 |> mutate(tol = "Tolerance = 0.65")

sim_res <- rbind(sim_res_tol_05, sim_res_tol_065)

sim_res |>
  pivot_longer(2:(p+2), names_to = "variable") |>
  group_by(method, variable, tol) |>
  count(value) |>
  mutate(variable = factor(variable, levels = c("(Intercept)", paste0("V", 1:p)))) |>
  mutate(n = (n/1000)*100,
         percentage = paste0(round(n, 2), '%')) |>
  ggplot(aes(x = variable, fill = value, y = n)) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.75) +
  geom_text(aes(label=percentage), position=position_stack(vjust = 0.5),
            hjust = 0.5, size = 2.5, col = "black")+
  facet_grid(rows = vars(method), cols = vars(tol))+
  labs(x = "Covariates",
       y = "Percentage",
       fill = "Decision") +
  scale_x_discrete(labels = wrap_format(50)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c("significant" = "darkgreen",
                                      "agnostic" = "goldenrod",
                                      "non significant" = "darkred"),
                             drop = FALSE)






