#' @title REACT forestplot.
#' @description REACT adapted forestplot for agnostic and pragmatic meta analysis data.
#' @param CI_matrix Matrix of lower and upper bounds of confidence intervals related to each study
#' with columns named CIlower and CIupper.
#' @param point_estim Logical or vector of point estimations related to each study parameter. If TRUE,
#' it uses the mean of the confidence intervals as the point estimate. Default is FALSE.
#' @param study_names Vector of study names. Default is NULL.
#' @param hyp precise hypothesis of interest. Default is 0.
#' @param NNT Number Necessary to Treat used to build pragmatic hypothesis.
#' @param fail Boolean specifying whether the probabilities are for failure (TRUE) or success (FALSE).
#' Default is FALSE.
#' @return React forestplot ggplot object
#'
#' @export
REACT_forestplot <- function(CI_matrix,
                             NNT,
                             hyp = 0,
                             point_estim = FALSE,
                             study_names = NULL,
                             fail = FALSE){
  if(fail){
    message("Probability of failure choosed. We will consider P(Treatment) - P(Control) as input. \n")
    CI_matrix <- -rev(CI_matrix)
    meassage("Transforming P(Treatment) - P(Control) into P(Control) - P(Treatment)")
  }else{
    message("Probability of success choosed. We will consider P(Treatment) - P(Control) as input.")
  }
  tol <- 1/NNT
  if(length(tol) == 1){
    p_int <- c(hyp - tol, hyp + tol)
  }else{
    p_int <- c(hyp + tol[1], hyp + tol[2])
  }

  idxs <- matrix(findInterval(CI_matrix,
                              p_int), ncol = 2)

  cols <- factor(ifelse(idxs[,1] == 1 & idxs[,2] == 1, 0,
                        ifelse((idxs[,1] == 0 & idxs[,2] == 0) |
                                 (idxs[,1] == 2 & idxs[,2] == 2), 1,
                               1/2)),
                 levels = c(0, 1/2, 1),
                 labels = c("Accept", "Agnostic", "Reject"))

  if(is.null(study_names)){
    study_names <- as.factor(1:nrow(CI_matrix))
  }
  if(is.logical(point_estim)){
    if(!point_estim){
      point_estim <- NA
    }else{
      point_estim <- rowSums(CI_matrix)/2
    }
  }
  meta_data <- data.frame(CIlower = CI_matrix[,1],
                          CIupper = CI_matrix[,2],
                          points = point_estim,
                          color = cols,
                          studlab = factor(study_names, levels = rev(sort(unique(study_names)))))

  # arrows limits
  #arr_inf <- min(c(CI_matrix[,1], p_int[1])) - 0.05
  #arr_sup <- max(c(CI_matrix[, 2], p_int[2])) + 0.05
  CI_sd <- sqrt(max(var(CI_mat)))
  arr_inf <- min(c(CI_matrix[,1])) - CI_sd
  arr_sup <- max(c(CI_matrix[, 2])) + CI_sd

  # adding annotations
  p <- meta_data %>%
    ggplot2::ggplot(ggplot2::aes(y = studlab, xmin = CIlower, xmax = CIupper, col = color)) +
    ggplot2::geom_errorbarh(ggplot2::aes(height=.2)) +
    {if(!all(is.na(meta_data$points))) ggplot2::geom_point(ggplot2::aes(x = points))} +
    ggplot2::annotate('rect', xmin = p_int[1], xmax = p_int[2],
                      ymin = 0, ymax = nrow(CI_matrix), alpha=.2, fill='dodgerblue3') +
    ggplot2::geom_vline(xintercept = hyp, linetype = 'dashed') +
    ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                         "Agnostic" = "goldenrod",
                                         "Reject" = "darkred")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 35)))+
    ggplot2::coord_cartesian(ylim=c(0 + 0.5, nrow(CI_matrix) - 0.5), clip = "off") +
    ggplot2::xlim(arr_inf, arr_sup) +
    ggplot2::annotate("text", x = hyp - CI_sd,
                      y = -1.5, label = "Favors Control",
                      colour = "orangered")+
    ggplot2::annotate("segment", x = hyp - CI_sd/4, xend = arr_inf,
                      y = -1, yend = -1, colour = "orangered", size = 1, alpha=0.6,
                      arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
    ggplot2::annotate("text", x = hyp + CI_sd,
                      y = -1.5, label = "Favors Treatment",
                      colour = "slateblue4")+
    ggplot2::annotate("segment", x = hyp + CI_sd/4, xend = arr_sup,
                      y = -1, yend = -1, colour = "slateblue4", size = 1, alpha=0.6,
                      arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
    ggplot2::labs(title = "NNT-based REACT Forestplot",
                  x = "Mean Difference", y = "Study", color = "Decision")

  methods::show(p)

  invisible(p)
}
