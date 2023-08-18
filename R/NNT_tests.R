#' @title REACT for difference between success proportions of two independent populations
#' @description NNT based REACT for difference between success proportions of two independent populations.
#' @param alpha alpha level. Default is 0.05.
#' @param NNT Number Necessary to Treat used to build pragmatic hypothesis.
#' @param event_e Number of successes in treatment group.
#' @param n_e Sample size of treatment group.
#' @param event_c Number of successes in control (or other) group.
#' @param n_c Sample size of control (or other) group.
#' @param verbose set whether text output should be generated (verbose = TRUE) or not (verbose = FALSE).
#' Default is TRUE.
#' @return "simple REACT" object with the original confidence interval, pragmatic region derived by the chosen NNT and the test outcome
#' @export

NNT_indep_test <- function(NNT, event_e, n_e, event_c, n_c,
                            alpha = 0.05, verbose = TRUE){
  epsilon <- 1/NNT
  prob_e <- event_e/n_e
  prob_c <- event_c/n_c

  var_e <- prob_e*(1-prob_e)/n_e
  var_c <- prob_c*(1-prob_c)/n_c

  prob_diff <- event_e/n_e - event_c/n_c
  v_diff <- var_e + var_c

  CI <- prob_diff +
    stats::qnorm(c(alpha/2, 1 - (alpha)/2))*sqrt(v_diff)

  obj <- REACT::base_test(CI = CI, tol = epsilon, verbose = verbose)

  return(obj)
}
