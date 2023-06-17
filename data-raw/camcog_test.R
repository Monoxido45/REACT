library(dplyr)
library(ggplot2)
library(ggforce)
library(latex2exp)

# importing data
camcog <- read.csv("data-raw/CAMCOG.csv")
camcog |> count(Diagnostico)
my_labs <- (camcog |> count(Diagnostico))$Diagnostico

# different counts for each group of diagnostic
# computing means for each group
stats <- camcog |>
  group_by(Diagnostico) |>
  summarise(mu = mean(CAMCOG, na.rm = TRUE),
            sigma = var(CAMCOG, na.rm = TRUE)/n())

var_cov <- diag(stats$sigma, nrow = 3, ncol = 3)
mu <- stats$mu
plot_list <- list()
c <- 1
alpha <- 0.05
delta <- 15


# Testing one ellipsis ----------------------------------------------------
# testing ggplot2
new_cov <- var_cov[c(1, 2), c(1, 2)]
new_mu <- mu[c(1, 2)]
q <- qchisq(1 - alpha, nrow(var_cov))

eig <- Matrix::Schur(new_cov)
eigs = eig$EValues
# maxvalue lambda
idx <- which.max(eig$EValues)
vecs <- eig$Q
main_eigenvec <- eig$Q[,1]

# minor and major axis
a <- sqrt(eigs[1]*q)
b <- sqrt(eigs[2]*q)

angle = atan2(main_eigenvec[2], main_eigenvec[1])

# diferrences between two means to increase xlab and ylab
mu_dif <- abs(new_mu[1] - new_mu[2])

points <- DescTools::DrawEllipse(x = new_mu[1], y = new_mu[2],
                                 radius.x = a, radius.y = b,
                       rot = angle , plot = FALSE, nv = 4000)


uniroot(check_ellipse,
        c(min(points$x) - 10, max(points$x) + 10),
        extendInt = "yes",
        delta = delta, c1 = new_mu[1], c2 = new_mu[2], a = a, b = b,
        tol = 0.0001)



print(check_ellipse(x_cand, x_cand - delta,
                    c1 = new_mu[1], c2 = new_mu[2],
                    a = a, b = b))

print(check_ellipse(x_cand, x_cand + delta,
                    c1 = new_mu[1], c2 = new_mu[2],
                    a = a, b = b))

# points data
points_data <- data.frame(x = points$x,
                          y = points$y)

ggplot(aes(x = x, y = y), data = points_data) +
  geom_polygon(colour="black", fill=NA)+
  geom_abline(intercept = delta, linetype = "dashed")+
  geom_abline(intercept = -delta, linetype = "dashed")+
  coord_cartesian(xlim = c(new_mu[1] - mu_dif - delta, new_mu[1] + mu_dif + delta),
                  ylim = c(new_mu[2] - mu_dif - delta, new_mu[2] + mu_dif + delta))



# Multiple comparisons ----------------------------------------------------
c = 1
plot_list = list()
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    new_cov <- var_cov[c(i, j), c(i, j)]
    new_mu <- mu[c(i, j)]
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

    points <- DescTools::DrawEllipse(x = new_mu[1], y = new_mu[2],
                                     radius.x = a, radius.y = b,
                                     rot = angle, plot = FALSE, nv = 1000)

    # using the points to determine whether the ellipse is intercepted by one of the lines
    # above line
    check_ellipse <- function(c1, c2, a, b, theta, delta){
      A <- ((cos(theta)^2/a^2) + (sin(theta)^2/b^2))
      B <- (2*cos(theta)*sin(theta)*((1/a^2) - (1/b^2)))
      C <- ((cos(theta)^2/b^2) + (sin(theta)^2/a^2))

      # checking if it has solutions for each pragmatic bound
      C_1 <- A + B + C
      C_2 <- ((delta - c1 - c2)*B) - (2*A*c1) + (C*((2*delta) - (2 * c2)))
      C_3 <- ((A * (c1^2)) + (((-c1*delta) + (c1*c2))*B) +
                (C * (delta - c2)^2) - 1)

      return((C_2^2) - (4*C_1*C_3))
    }
    abv_roots <- check_ellipse(new_mu[1], new_mu[2], a, b, angle, delta)
    blw_roots <- check_ellipse(new_mu[1], new_mu[2], a, b, angle, -delta)

    res <- ifelse(any(blw_roots > 0, abv_roots > 0), 1/2, 0)

    if(res == 0){
      # checking if center of ellipse is inside pragmatic region
      inside <- ((new_mu[1] - delta) <= new_mu[2] & (new_mu[1] + delta) >= new_mu[2])

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
    tex_title = TeX(paste0("$\\mu_", i,  "$ and $\\mu_", j, "$"))
    tex_xlab = TeX(paste0("$\\mu_", i, "$"))
    tex_ylab = TeX(paste0("$\\mu_", j, "$"))

    data_used <- data.frame(x = points$x,
                            y = points$y,
                            col = cols)

    # building data to construct pragmatic region
    x_min <- min(data_used$x) - delta
    x_max <- max(data_used$x) + delta

    # adding y to each x between x_min and x_max
    y_1 <- seq(x_min - delta, x_max - delta, length.out = 300)
    y_2 <- seq(x_min + delta, x_max + delta, length.out = 300)
    x <- seq(x_min, x_max, length.out = 300)
    prag_data <- data.frame(x = x,
                            y_1 = y_1,
                            y_2 = y_2)

    plot_list[[c]] <- ggplot(aes(x = x, y = y, fill = col, colour = col), data = data_used) +
      geom_polygon(alpha = 0.4)+
      ggplot2::scale_color_manual(values=c("Accept" = "darkgreen",
                                           "Agnostic" = "goldenrod",
                                           "Reject" = "darkred"),
                                  drop = FALSE) +
      ggplot2::scale_fill_manual(values=c("Accept" = "darkgreen",
                                           "Agnostic" = "goldenrod",
                                           "Reject" = "darkred"),
                                 drop = FALSE)+
      geom_ribbon(aes(x = x, ymin = y_1, ymax = y_2),
                  data = prag_data, inherit.aes = FALSE,
                  fill = "dodgerblue3", alpha = 0.2)+
      theme_bw()+
      labs(x = tex_xlab,
           y = tex_ylab,
           title = tex_title,
           colour = "Decision",
           fill = "Decision")

    c = c + 1
  }
}

p <- ggpubr::ggarrange(plotlist = plot_list, nrow = 1, ncol = 3, common.legend = TRUE)



