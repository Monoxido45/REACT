#' @title REACT for multiple parameters comparisons with logical consistency
#' @description REACT plotting for multiple parameters comparisons. It uses asymptotic properties from
#' MLE to plot pairwise confidence ellipsis and conclude if we reject, accept or remain agnostic about
#' differences $|\theta_i - \theta_j|$ between each pair of parameter $\theta_i$ and $\theta_j$.
#' @param alpha alpha level (default is 0.05)
#' @param par vector of point estimates of each parameter
#' @param f_matrix fisher information matrix associated to the vector of parameters
#' @param nrow Number of rows for multiple plotting. If NA, we choose based on vector length
#' (default is NA)
#' @param ncol Number of columns for multiple plotting. If NA, we choose based on vector length
#' (default is NA)
#' @param verbose set whether text output should be generated (verbose = TRUE) or not (verbose = FALSE)
#' @export

m_comparisons <- function(alpha = 0.05, nrow = NA, ncol = NA, par, f_matrix){
  # variance covariance matrix
  f_var_cov <- solve(f_matrix)
  # counting
  c <- 1
  q <- qchisq(1 - alpha, length(par))

  data_list = list()
  for(i in 1:(nrow(f_var_cov) - 1)){
    for(j in (i + 1):nrow(f_var_cov)){
      new_cov <- f_var_cov[c(i, j), c(i, j)]
      new_par <- par[c(i, j)]

      eig <- eigen(new_cov)
      # maxvalue lambda
      idx <- which.max(eig$values)
      # minor and major axis
      a <- sqrt(eig$values[idx]*q)
      b <- sqrt(eig$values[-idx]*q)

      larg_eigevec <- eig$vectors[,idx]
      angle = atan2(larg_eigevec[2], larg_eigevec[1])

      points <- DescTools::DrawEllipse(x = new_par[1], y = new_par[2],
                                       radius.x = a, radius.y = b,
                                       rot = angle, plot = FALSE, nv = 1000)

      # points data
      tex_label = paste0("$\\theta_", i,  "$ and $\\theta_", j, "$")
      data_list[[c]] <- data.frame(x = points$x,
                                   y = points$y,
                                   label = tex_label)

      c = c + 1
    }
  }

  # concatenating data frame
  all_data <- do.call("rbind", data_list)

  x_min <- min(all_data$x) - delta
  x_max <- max(all_data$x) + delta

  # adding y to each x between x_min and x_max
  y_1 <- seq(x_min - delta, x_max - delta, length.out = 50)
  y_2 <- seq(x_min + delta, x_max + delta, length.out = 50)
  x <- seq(x_min, x_max, length.out = 50)
  prag_data <- data.frame(x = x,
                          y_1 = y_1,
                          y_2 = y_2)

  appender <- function(string){
   latex2exp::TeX(string)}


  p <- ggplot2::ggplot(ggplot2::aes(x = x, y = y), data = all_data) +
    ggplot2::geom_polygon(colour = "black", fill=NA)+
    ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_1, ymax = y_2),
                data = prag_data, inherit.aes = FALSE,
                fill = "dodgerblue3", alpha = 0.2) +
    ggplot2::facet_wrap(~label, scales = "free",
                        labeller = ggplot2::as_labeller(appender,
                                        default = ggplot2::label_parsed))

  show(p)
  invisible(p)
}
