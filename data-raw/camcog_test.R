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



# Multiple comparisons ----------------------------------------------------
c = 1
data_list = list()
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
  TeX(string)}


ggplot(aes(x = x, y = y), data = all_data) +
  geom_polygon(colour = "black", fill=NA)+
  geom_ribbon(aes(x = x, ymin = y_1, ymax = y_2),
              data = prag_data, inherit.aes = FALSE,
              fill = "dodgerblue3", alpha = 0.2) +
  facet_wrap(~label, scales = "free", labeller = as_labeller(appender,
                                                             default = label_parsed))


ggpubr::ggarrange(plotlist = plot_list, nrow = 1, ncol = 3)


