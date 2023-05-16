#' @title Uni parametric Agnostic Confidence-based precise test
#' @description REACT function for uni parametric precise hypothesis testing
#' @param CI parameter confidence interval for testing
#' @param tol tolerance to build pragmatic hypothesis. It can be also an assymetric tolerance (vector of length 2)
#' @param hip precise hypothesis of interest - default is 0
#' @param verbose set whether text output should be generated (verbose = TRUE) or not (verbose = FALSE)
#' @export
base_test <- function(CI, tol, hip = 0, verbose = TRUE){
  # designing pragmatic interval
  if(length(tol) == 1){
    p_int <- c(hip - tol, hip + tol)
  }else{
    p_int <- c(hip - tol[1], hip + tol[2])
  }

  # analysing whether CI falls into pragmatic region, or not or not entirely
  idxs <- findInterval(CI, p_int)

  # if both values are equal it means that the CI is entirely inside or outside the pragmatic region
  test_outcome <- ifelse(all(idxs == 1), "accept",
         ifelse(all(idxs == 0) | all(idxs == 2),
         "reject", "remain agnostic"))

  if(verbose){
    cat("REACT results:\n")
    cat("Pragmatic lower bound: ", format(p_int[1], digits = 3, nsmall = 2, scientific = FALSE))
    cat("\n")
    cat("Pragmatic upper bound: ", format(p_int[2], digits = 3, nsmall = 2, scientific = FALSE))
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


  vals <- list(CI = CI,
       pragmatic = p_int,
       test_outcome = test_outcome)

  class(vals) = c("simple_REACT")

  return(vals)
}

# TODO: multiple testing interface (if needed)



#' @title Plot Region of equivalence and confidence interval of precise hypothesis test
#' @param obj simple_REACT object to plot
#' @export
plot.simple_REACT <- function(obj){

  p <- ggplot2::ggplot()+
    ggplot2::coord_cartesian(ylim=c(-0.1, 0.1)) +
    ggplot2::geom_segment(ggplot2::aes(x = obj$CI[1], y = 0,
                                       xend = obj$CI[2], yend = 0), linewidth = 0.75)+
    ggplot2::annotate('rect', xmin = obj$pragmatic[1], xmax = obj$pragmatic[2],
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
  invisible(p)
}


