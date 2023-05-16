#' @title REACT forestplot.
#' @description REACT adapted forestplot for agnostic and pragmatic meta analysis data.
#' @param CI_matrix Matrix of lower and upper bounds of confidence intervals related to each study
#' with columns named CIlower and CIupper.
#' @param point_estim Logical or vector of point estimations related each study parameter. If TRUE,
#' it uses the mean of the confidence intervals as the point estimate. Default is FALSE.
#' @param study_names Vector of study names. Default is NULL.
#' @param hip precise hypothesis of interest - default is 0.
#' @param NNT Number Necessary to Treat used to build pragmatic hypothesis.
#' @export
REACT_forestplot <- function(CI_matrix,
                             NNT,
                             hip = 0,
                             point_estim = FALSE,
                             study_names = NULL){
  epsilon <- 1/NNT
  idxs <- matrix(findInterval(CI_matrix,
                              c(hip - epsilon,hip + epsilon)), ncol = 2)

  cols <- factor(ifelse(idxs[,1] == 1 & idxs[,2] == 1, 0,
                        ifelse((idxs[,1] == 0 & idxs[,2] == 0) |
                                 (idxs[,1] == 2 & idxs[,2] == 2), 1,
                               1/2)),
                 levels = c(0, 1/2, 1),
                 labels = c("accept", "agnostic", "reject"))

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
  p <- meta_data %>%
    ggplot2::ggplot(ggplot2::aes(y = studlab, xmin = CIlower, xmax = CIupper, col = color)) +
    ggplot2::geom_errorbarh(height=.2, linewdith = 1) +
    {if(!all(is.na(meta_data$points))) ggplot2::geom_point(ggplot2::aes(x = points))} +
    ggplot2::annotate('rect', xmin = hip - epsilon, xmax = hip + epsilon,
                      ymin = 0, ymax = nrow(CI_matrix), alpha=.2, fill='dodgerblue3') +
    ggplot2::geom_vline(xintercept = hip, linetype = 'dashed') +
    ggplot2::scale_color_manual(labels = c("Accept", "Agnostic", "Reject"),
                                values=c("darkgreen", "goldenrod", "darkred")) +
    labs(title = "NNT-based REACT Forestplot", x = "Mean Difference", y = "Study", color = "Decision")

  methods::show(p)

  invisible(p)
}
