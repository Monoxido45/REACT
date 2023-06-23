library(dplyr)
library(ggplot2)
library(ggforce)
library(latex2exp)

# importing data (CG - Control, DA - Alzheimer, CCL - MCI?)
camcog <- read.csv("data-raw/CAMCOG.csv")
camcog |> count(Diagnostico)
my_labs <- (camcog |> count(Diagnostico))$Diagnostico

# different counts for each group of diagnostic
# computing means for each group
stats <- camcog |>
  group_by(Diagnostico) |>
  summarise(mu = mean(CAMCOG, na.rm = TRUE),
            sigma = var(CAMCOG, na.rm = TRUE)/n(),
            n = n())

# computing pooled variance
pooled_var <- stats |>
  mutate(new_sigma = sigma*(n - 1)/(sum(n) - 3)) |>
  pull(new_sigma) |> sum()

var_cov <- diag(stats$sigma, nrow = 3, ncol = 3) #1 - CCL, 2 - DA, 3 - GC
pooled_var_cov <- diag(pooled_var, nrow = 3, ncol = 3)
mu <- stats$mu
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
names_list <- c("MCI", "AD", "CG")
# supposing heterocedastic behavior:
p <- REACT::m_comparisons(alpha = 0.05, nrow = 1, ncol = 3,
                     delta = delta, par = mu, f_matrix = var_cov)

# modifying labs
c <- 1
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    tex_xlab = TeX(paste0("$\\mu_{", names_list[i], "}$"))
    tex_ylab = TeX(paste0("$\\mu_{", names_list[j], "}$"))

    p[[c]] <- p[[c]] + labs(x = tex_xlab,
                  y = tex_ylab,
                  title = paste0(names_list[i], " vs ", names_list[j]),
                  colour = "Decision",
                  fill = "Decision")
    c <- c + 1
  }}

ggpubr::ggarrange(plotlist = p, nrow = 1, ncol = 3, common.legend = TRUE)

# supposing homocedastic behavior:
p <- REACT::m_comparisons(alpha = 0.05, nrow = 1, ncol = 3,
                          delta = delta, par = mu, f_matrix = pooled_var_cov)

# modifying labs
c <- 1
for(i in 1:(nrow(var_cov) - 1)){
  for(j in (i + 1):nrow(var_cov)){
    tex_xlab = TeX(paste0("$\\mu_{", names_list[i], "}$"))
    tex_ylab = TeX(paste0("$\\mu_{", names_list[j], "}$"))

    p[[c]] <- p[[c]] + labs(x = tex_xlab,
                            y = tex_ylab,
                            title = paste0(names_list[i], " vs ", names_list[j]),
                            colour = "Decision",
                            fill = "Decision")
    c <- c + 1
  }}

ggpubr::ggarrange(plotlist = p, nrow = 1, ncol = 3, common.legend = TRUE)



