library(dplyr)
library(ggplot2)
library(ggforce)

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

# testing ggplot2
new_cov <- var_cov[c(1, 2), c(1, 2)]
new_mu <- mu[c(1, 2)]
q <- qchisq(1 - alpha, nrow(var_cov))

eig <- eigen(new_cov)
# maxvalue lambda
idx <- which.max(eig$values)
# minor and major axis
a <- sqrt(eig$values[idx]*q)
b <- sqrt(eig$values[-idx]*q)

larg_eigevec <- eig$vectors[,idx]
angle = atan2(larg_eigevec[2], larg_eigevec[1])

# diferrences between two means to increase xlab and ylab
mu_dif <- abs(new_mu[1] - new_mu[2])

points <- DescTools::DrawEllipse(x = new_mu[1], y = new_mu[2], radius.x = a, radius.y = b,
                       rot = angle, plot = FALSE, nv = 4000)

# points data
points_data <- data.frame(x = points$x,
                          y = points$y)

ggplot(aes(x = x, y = y), data = points_data) +
  geom_polygon(colour="black", fill=NA)+
  geom_abline(intercept = delta, linetype = "dashed")+
  geom_abline(intercept = -delta, linetype = "dashed")+
  coord_cartesian(xlim = c(new_mu[1] - mu_dif - delta, new_mu[1] + mu_dif + delta),
                  ylim = c(new_mu[2] - mu_dif - delta, new_mu[2] + mu_dif + delta))



c = 1
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    new_cov <- var_cov[c(i, j), c(i, j)]
    new_mu <- mu[c(i, j)]
    q <- qchisq(1 - alpha, nrow(var_cov))

    eig <- eigen(new_cov)
    # maxvalue lambda
    idx <- which.max(eig$values)
    # minor and major axis
    a <- sqrt(eig$values[idx]*q)
    b <- sqrt(eig$values[-idx]*q)

    larg_eigevec <- eig$vectors[,idx]
    angle = atan2(larg_eigevec[2], larg_eigevec[1])

    # diferrences between two means to increase xlab and ylab
    mu_dif <- abs(new_mu[1] - new_mu[2])

    points <- DescTools::DrawEllipse(x = new_mu[1], y = new_mu[2], radius.x = a, radius.y = b,
                                     rot = angle, plot = FALSE, nv = 4000)

    # points data
    points_data <- data.frame(x = points$x,
                              y = points$y)

    plot_list[[c]] <- ggplot(aes(x = x, y = y), data = points_data) +
      geom_polygon(colour="black", fill=NA)+
      geom_abline(intercept = delta, linetype = "dashed")+
      geom_abline(intercept = -delta, linetype = "dashed")+
      coord_cartesian(xlim = c(new_mu[1] - mu_dif - delta, new_mu[1] + mu_dif + delta),
                      ylim = c(new_mu[2] - mu_dif - delta, new_mu[2] + mu_dif + delta))+
      labs(x = bquote(mu[.(my_labs[i])]),
           y = bquote(mu[.(my_labs[j])]))

    c = c + 1
  }
}

ggpubr::ggarrange(plotlist = plot_list, nrow = 1, ncol = 3)


