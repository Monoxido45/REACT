#' @title REACT GLM covariate selection
#' @description REACT-based general linear model with agnostic and pragmatic feature selection of covariates.
#' It reproduces the same steps as a classical glm with normalized covariates and conducts significance testing
#' of covariates using the REACT procedure, providing three conclusions about each variable: significant,
#' non-significant and agnostic.
#' @param formula An object of class "formula" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param family a description of the error distribution and link function to be used in the model.
#' For glm this can be a character string naming a family function, a family function or the result of a
#' call to a family function.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model.If not found in data, the variables are taken from environment(formula),
#' typically the environment from which glm is called.
#' @param tol Tolerance to build pragmatic hypothesis around each covariate parameter.
#' It can be a vector of different tolerances (the same size as the number of covariates) or a global tolerance
#' @param alpha Chosen significance level to build confidence intervals for testing at level alpha. Default is 0.05.
#' @param verbose Set whether REACT conclusions output should be printed (verbose = TRUE) or not (verbose = FALSE)
#' @param weights An optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param na.action A function which indicates what should happen when the data contain NAs.
#' The default is set by the \link[stats]{na.fail}.
#' @param start Starting values for the parameters in the linear predictor.
#' @param etastart Starting values for the linear predictor.
#' @param mustart Starting values for the vector of means.
#' @param offset This can be used to specify an a priori known component to be included in the linear predictor during fitting.
#' This should be NULL or a numeric vector of length equal to the number of cases.
#' @param control A list of parameters for controlling the fitting process. For glm.fit this is passed to \link[stats]{glm.control}.
#' @param model A logical value indicating whether model frame should be included as a component of the returned value.
#' @param method The method to be used in fitting the model. The default method "glm.fit" uses iteratively reweighted least squares
#' (IWLS): the alternative "model.frame" returns the model frame and does no fitting.
#'
#' User-supplied fitting functions can be supplied either as a function or a character string naming a function, with a function
#' which takes the same arguments as glm.fit. If specified as a character string it is looked up from within the stats namespace.
#' @param x,y For glm: logical values indicating whether the response vector and model matrix used in the fitting process should be returned as components of the returned value.
#'
#' For glm.fit: x is a design matrix of dimension n * p, and y is a vector of observations of length n.
#' @param singular.ok Logical; if FALSE a singular fit is an error.
#' @param constrats An optional list. See the contrasts.arg of model.matrix.default.
#' @param intercept Logical. Should an intercept be included in the null model?
#' @param object An object inheriting from class "glm".
#' @param type Character, partial matching allowed. Type of weights to extract from the fitted model object. Can be abbreviated.
#' @param ... For glm: arguments to be used to form the default control argument if it is not supplied directly.
#'
#' For weights: further arguments passed to or from other methods.
#'
#' @return A personalized object of class "REACT.glm" and "REACT.lm" with the original glm output
#' plus a list of significant, non-significant and agnostic covariates, chosen tolerances and significance levels
#' and parameter confidence intervals.
#' @export
#'
#' @examples
#' # DL95 example
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' react_glm_D93 <- REACT.glm(counts ~ outcome + treatment, family = poisson(),
#' tol = 0.5, alpha = 0.05, verbose = TRUE)
#' react_glm_D93
#'
# TODO: change description and documentation
REACT.glm <- function(formula, family = stats::gaussian, data, weights, subset, tol,
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
  Y <- stats::model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!stats::is.empty.model(mt))
    stats::model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)

  # scaling model.matrix except for the intercept
  X[,-1] <- scale(X[, -1])

  weights <- as.vector(stats::model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(stats::model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  mustart <- stats::model.extract(mf, "mustart")
  etastart <- stats::model.extract(mf, "etastart")
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
                        contrasts = attr(X, "contrasts"), xlevels = stats::.getXlevels(mt,
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
  CI_s <- cbind(final_obj$coefficients - sqrt(stats::qchisq(1 - alpha, p)*diag(beta_vcov)),
                final_obj$coefficients + sqrt(stats::qchisq(1 - alpha, p)*diag(beta_vcov)))

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
  final_obj$alpha <- alpha
  final_obj
}

# print function
print.REACT.glm <- function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(stats::coef(x))) {
    cat("Coefficients:\n")
    print.default(format(stats::coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

# summary function


# TODO: plot function

