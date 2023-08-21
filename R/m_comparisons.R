#' @title REACT for multiple parameters comparisons with logical consistency
#' @description REACT plotting for multiple parameters comparisons. It uses asymptotic properties from
#' MLE to plot pairwise confidence ellipsis and conclude if we reject, accept or remain agnostic about
#' each differences \eqn{|\theta_i - \theta_j|} between each pair of parameter \eqn{\theta_i} and \eqn{\theta_j}.
#' Returns a plot list with all plottings of each comparisson.
#' @param alpha alpha level (default is 0.05)
#' @param par vector of point estimates of each parameter
#' @param tol tolerance to build pragmatics.
#' @param f_matrix fisher information matrix associated to the vector of parameters
#' @param invert Boolean indicating whether fisher matrix needs to be inverted or no.
#' Default is FALSE
#' @param nrow Number of rows for multiple plotting. If NA, we choose based on vector length
#' (default is NA)
#' @param ncol Number of columns for multiple plotting. If NA, we choose based on vector length
#' (default is NA)
#' @return A list with each parameter comparisons ggplot-based plots
#' @export
#'

m_comparisons <- function(alpha = 0.05,
                          nrow = NA,
                          ncol = NA,
                          invert = FALSE,
                          tol, par,
                          f_matrix){
  # variance covariance matrix
  if(invert){
  f_var_cov <- solve(f_matrix)
  }else{
    f_var_cov <- f_matrix
  }
  # counting
  c <- 1
  q <- stats::qchisq(1 - alpha, length(par))
  plot_list <- list()

  for(i in 1:(nrow(var_cov) - 1)){
    for(j in (i + 1):nrow(var_cov)){
      new_cov <- f_var_cov[c(i, j), c(i, j)]
      new_par <- par[c(i, j)]
      q <- qchisq(1 - alpha, nrow(var_cov))

      eig <- Matrix::Schur(new_cov)
      vecs <- eig$Q
      eigs <- eig$EValues

      main_eigenvec <- vecs[,1]

      # minor and major axis
      a <- sqrt(eigs[1]*q)
      b <- sqrt(eigs[2]*q)

      angle = atan2(main_eigenvec[2], main_eigenvec[1])

      # diferrences between two means to increase xlab and ylab

      points <- DescTools::DrawEllipse(x = new_par[1], y = new_par[2],
                                       radius.x = a, radius.y = b,
                                       rot = angle, plot = FALSE, nv = 1000)

      # using the points to determine whether the ellipse is intercepted by one of the lines
      # above line
      check_ellipse <- function(c1, c2, a, b, theta, tol){
        A <- ((cos(theta)^2/a^2) + (sin(theta)^2/b^2))
        B <- (2*cos(theta)*sin(theta)*((1/a^2) - (1/b^2)))
        C <- ((cos(theta)^2/b^2) + (sin(theta)^2/a^2))

        # checking if it has solutions for each pragmatic bound
        C_1 <- A + B + C
        C_2 <- ((tol - c1 - c2)*B) - (2*A*c1) + (C*((2*tol) - (2 * c2)))
        C_3 <- ((A * (c1^2)) + (((-c1*tol) + (c1*c2))*B) +
                  (C * (tol - c2)^2) - 1)

        return((C_2^2) - (4*C_1*C_3))
      }
      abv_roots <- check_ellipse(new_par[1], new_par[2], a, b, angle, tol)
      blw_roots <- check_ellipse(new_par[1], new_par[2], a, b, angle, -tol)

      res <- ifelse(any(blw_roots > 0, abv_roots > 0), 1/2, 0)

      if(res == 0){
        # checking if center of ellipse is inside pragmatic region
        inside <- ((new_par[1] - tol) <= new_par[2] & (new_par[1] + tol) >= new_par[2])

        # if is not inside, and is tangent, then agnostic
        if(!inside & any(blw_roots == 0, abv_roots == 0)){
          res <- 1/2
          # if is inside, then accept
        }else if(inside){
          res <- 0
          # if not inside, then reject
        }else{
          res <- 1
        }
      }

      cols <- factor(res,
                     levels = c(0, 1/2, 1),
                     labels = c("Accept", "Agnostic", "Reject"))

      # plotting based on res
      # points data
      tex_title = latex2exp::TeX(paste0("$\\theta_", i,  "$ and $\\theta_", j, "$"))
      tex_xlab = latex2exp::TeX(paste0("$\\theta_", i, "$"))
      tex_ylab = latex2exp::TeX(paste0("$\\theta_", j, "$"))

      data_used <- data.frame(x = points$x,
                              y = points$y,
                              col = cols)

      # building data to construct pragmatic region
      x_min <- min(data_used$x) - tol
      x_max <- max(data_used$x) + tol

      # adding y to each x between x_min and x_max
      y_1 <- seq(x_min - tol, x_max - tol, length.out = 300)
      y_2 <- seq(x_min + tol, x_max + tol, length.out = 300)
      x <- seq(x_min, x_max, length.out = 300)
      prag_data <- data.frame(x = x,
                              y_1 = y_1,
                              y_2 = y_2)

      plot_list[[c]] <- ggplot2::ggplot(
        ggplot2::aes(x = x, y = y, fill = col, colour = col), data = data_used) +
        ggplot2::geom_polygon(alpha = 0.4)+
        ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                             "Agnostic" = "goldenrod",
                                             "Reject" = "darkred"),
                                    drop = FALSE) +
        ggplot2::scale_fill_manual(values=c("Accept" = "darkgreen",
                                            "Agnostic" = "goldenrod",
                                            "Reject" = "darkred"),
                                   drop = FALSE)+
        ggplot2::geom_ribbon(ggplot2::aes(x = x, ymin = y_1, ymax = y_2),
                    data = prag_data, inherit.aes = FALSE,
                    fill = "dodgerblue3", alpha = 0.2)+
        ggplot2::theme_bw()+
        ggplot2::labs(x = tex_xlab,
             y = tex_ylab,
             title = tex_title,
             colour = "Decision",
             fill = "Decision")

      c <- c + 1
    }
  }

  p <- ggpubr::ggarrange(plotlist = plot_list, nrow = nrow, ncol = ncol, common.legend = TRUE)
  show(p)

  invisible(plot_list)
}
