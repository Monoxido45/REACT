# REACT glm tests and simulations
library(REACT)
library(dplyr)
library(purrr)

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
glm_sim <- function(n, p, k = 5, type = "logistic", seed = 750){
  set.seed(seed)
  # covariate matrix
  X <- setNames(1:p, paste0("V",1:p)) |> map_dfc(function(x){runif(n, -1, 1)})

  # generating beta from significant variable from uniform 1 to 2
  beta_sig <- runif(k, min = 1, max = 2)

  # generating beta from non-significant variable from uniform -0.05 to 0.05
  beta_noise <- runif(p - k, min = -0.05, max = 0.05)

  all_beta <- c(beta_sig, beta_noise)

  if(type == "logistic"){
    linear_comp <- as.numeric(as.matrix(X) %*% all_beta)
    p_x <- 1/(1 + exp(- linear_comp))
    return(list(data = cbind(X, y = rbinom(n, 1, p_x)), true_betas = all_beta))
    }
}

n <- 2500
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
tol <- 0.3
alpha <- 0.05
react_obj <- REACT.glm(y ~ ., family = binomial, data = logis_sim$data,
          tol = tol, alpha = alpha, verbose = TRUE)

# repeating analysis but with less samples
n <- 500
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
tol <- 0.2
alpha <- 0.05
react_obj <- REACT.glm(y ~ ., family = binomial, data = logis_sim$data,
                       tol = tol, alpha = alpha, verbose = TRUE)



