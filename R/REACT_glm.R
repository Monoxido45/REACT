#' @title REACT GLM covariate selection
#' @description REACT-based feature selection of GLM covariates.
#' @param tol Tolerance to build pragmatic hypothesis. It can be also an asymmetric tolerance (vector of length 2)
#' @param alpha Chosen significance level. Default is 0.05
#' @param verbose Set whether text output should be generated (verbose = TRUE) or not (verbose = FALSE)
#'
#' @return An object of class "simple_REACT" with the original confidence interval, pragmatic region and the test outcome. There is also a plot method available for this kind of object.
#' @export
#'
# TODO: change description and documentation
REACT.glm <- function(formula, family = gaussian, data, weights, subset, tol,
                      na.action, alpha = 0.05, verbose = FALSE, start = NULL, etastart,
                      mustart, offset, control = list(...),model = TRUE, method = "glm.fit",
                      x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...)
{
  # original GLM code with some adaptations
  cal <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  if (identical(method, "glm.fit"))
    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)

  # scaling model.matrix except for the intercept
  X[,-1] <- scale(X[, -1])

  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
                   x = X, y = Y, weights = weights, start = start, etastart = etastart,
                   mustart = mustart, offset = offset, family = family,
                   control = control, intercept = attr(mt, "intercept") >
                     0L, singular.ok = singular.ok))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
                      x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values,
                      weights = weights, offset = offset, family = family,
                      control = control, intercept = TRUE))
    if (!fit2$converged)
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    fit$null.deviance <- fit2$deviance
  }
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  final_obj <- structure(c(fit, list(call = cal, formula = formula, terms = mt,
                        data = data, offset = offset, control = control, method = method,
                        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
                                                                                mf))),
            class = c(fit$class, c("glm", "lm")))

  # obtaining variance covariance matrix of coefficients
  beta_vcov <- stats::vcov(final_obj)

  # changing final obj class
  class(final_obj) <- c(fit$class, c("REACT.glm", "REACT.lm"))

  # adding to final obj the variable selection analysis
  p <- ncol(beta_vcov)
  if(!(length(tol) == 1 || length(tol) == p))
    stop("Tolerance must be a vector of same size than the number of parameters or 1")

  # obtaining Confidence Intervals from ellipsis projection
  CI_s <- cbind(final_obj$coefficients - sqrt(qchisq(1 - alpha, p)*diag(beta_vcov)),
                final_obj$coefficients + sqrt(qchisq(1 - alpha, p)*diag(beta_vcov)))

  # analysing if each CI is inside each region of equivalence or not
  conduct_base_test <- function(CI, tol){
    base_test(CI, tol, verbose = FALSE)$test_outcome
  }
  if(length(tol) == 1)
    tol <- rep(tol, p)

  # TODO: use another loop function instead of for
  # obtaining all results
  all_res <- character(p)
  for(i in 1:p){
    all_res[i] <- base_test(CI_s[i, ], tol[i], verbose = FALSE)$test_outcome
  }
  var_names <- rownames(CI_s)
  signif_variables <- var_names[which(all_res == "reject")]
  n_signif_variables <- var_names[which(all_res == "accept")]
  agnostic_variables <- var_names[which(all_res == "remain agnostic")]

  if(verbose){
    cat("REACT results:\n")
    cat("Fixed tolerances:", tol)
    cat("\n")
    cat("Significant variables: ", signif_variables)
    cat("\n")
    cat("Non significant variables: ", n_signif_variables)
    cat("\n")
    cat("Agnostic variables: ", agnostic_variables)
    cat("\n")
    cat("Some Confidence intervals (matrix form):")
    cat("\n")
    print(if(nrow(CI_s) <= 10){
      CI_s}else{
        CI_s[1:10, ]})
  }

  # adding results to glm object
  final_obj$react_CIs <- CI_s
  final_obj$signif_vars <- signif_variables
  final_obj$non_signif_vars <- n_signif_variables
  final_obj$agnostic_vars <- agnostic_variables
  final_obj$tol <- tol
  final_obj
}

# print function
print.REACT.glm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

# summary function


# TODO: plot function

