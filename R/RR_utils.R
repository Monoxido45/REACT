#' @title REACT Relative Risk forestplot.
#' @description REACT adapted forestplot for agnostic and pragmatic meta analysis of relative risk.
#' @param CI_matrix Matrix of lower and upper bounds of confidence intervals related to each covariate
#' with columns named CIlower and CIupper.
#' @param point_estim Logical or vector of point estimations related to each covariate. If TRUE,
#' it uses the mean of the confidence intervals as the point estimate. Default is FALSE.
#' @param covar_names Vector of covariate names. Default is NULL.
#' @param hyp precise hypothesis of interest - default is 1.
#' @param NNT Number Necessary to Treat used to build pragmatic hypothesis. Default is NULL.
#' @param CER Control Error Rate used to build pragmatic hypothesis. Default is 1.
#' @param RRR Relative Risk Reduction used to build pragmatic hypothesis. Default is NULL.
#' @return REACT forestplot ggplot object.
#' @export
REACT_RRplot <- function(CI_matrix,
                         NNT = NULL,
                         CER = 1,
                         RRR = NULL,
                         hyp = 1,
                         point_estim = FALSE,
                         covar_names = NULL){
  if(is.null(NNT) & is.null(RRR)) stop("Please specify either the NNT or the RRR.")
  if(!is.null(NNT)) RRR <- 1/(NNT*CER)
  tol <- abs(hyp - abs(1 - RRR))
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

  if(is.null(covar_names)){
    covar_names <- as.factor(1:nrow(CI_matrix))
  }
  if(is.logical(point_estim)){
    if(!point_estim){
      point_estim <- NA
    }else{
      point_estim <- rowSums(CI_matrix)/2
    }
  }
  RR_data <- data.frame(CIlower = CI_matrix[,1],
                          CIupper = CI_matrix[,2],
                          points = point_estim,
                          color = cols,
                          covarlab = factor(covar_names, levels = rev(sort(unique(covar_names)))))

  # arrows limits
  arr_inf <- min(c(CI_matrix[,1], p_int[1])) - 0.05
  arr_sup <- max(c(CI_matrix[, 2], p_int[2])) + 0.05

  # adding annotations
  p <- ggplot2::ggplot(data = RR_data, ggplot2::aes(y = covarlab, xmin = CIlower, xmax = CIupper, col = color)) +
    ggplot2::geom_errorbarh(ggplot2::aes(height=.2)) +
    {if(!all(is.na(RR_data$points))) ggplot2::geom_point(ggplot2::aes(x = points))} +
    ggplot2::annotate('rect', xmin = p_int[1], xmax = p_int[2],
                      ymin = 0, ymax = nrow(CI_matrix)+.5, alpha=.2, fill='dodgerblue3') +
    ggplot2::geom_vline(xintercept = hyp, linetype = 'dashed') +
    ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                         "Agnostic" = "goldenrod",
                                         "Reject" = "darkred")) +
    # ggplot2::annotate("text", x = p_int[1] - 0.1,
    #                   y = -.1, label = "Favors Control",
    #                   colour = "orangered")+
    # ggplot2::annotate("segment", x = hyp - 0.05, xend = arr_inf,
    #                   y = 0, yend = 0, colour = "orangered", size = 1, alpha=0.6,
    #                   arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
    # ggplot2::annotate("text", x = p_int[2] + 0.1,
    #                   y = -.1, label = "Favors Treatment",
    #                   colour = "slateblue4")+
    # ggplot2::annotate("segment", x = hyp + 0.05, xend = arr_sup,
    #                   y = 0, yend = 0, colour = "slateblue4", size = 1, alpha=0.6,
    #                   arrow=arrow(length = ggplot2::unit(0.1, "inches")))+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 35)))+
    ggplot2::coord_cartesian(ylim=c(0 + 0.5, nrow(CI_matrix)), clip = "off") +
    ggplot2::labs(title = "RRR-based REACT Procedure",
                  x = "Relative Risk (RR)", y = "Covariate", color = "Decision")

  methods::show(p)

  invisible(p)
}
