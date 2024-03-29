
<!-- README.md is generated from README.Rmd. Please edit that file -->

# REACT: Region of Equivalence Agnostic Confidence-based Test

<!-- badges: start -->
<!-- badges: end -->

This package offers a user-friendly interface to REACT, a framework that
combines equivalence testing (evaluating both practical significance and
statistical significance) and three-way testing (allowing for a
hypothesis to be accepted, rejected or for the user to remain
undecided/agnostic). Generally speaking, REACT is performed in three
steps:

1.  **Use the information of what are negligible differences (for
    example, results that are different only due to measurement errors)
    to establish the null hypothesis of actual practical interest,
    called the region of equivalence.**
2.  **Build a multivariate confidence set of all parameters to be
    tested. If a specific parameter will never be tested, it does not
    have to be included in the confidence set.**
3.  **Compare the null hypothesis with the confidence set, which is
    simply checking if the confidence set is either: fully contained by
    the null hypothesis (accept $H_0$), completely outside of it (reject
    $H_0$) or somewhere in between (remain agnostic).**

<!-- Continue writing more later -->

REACT ideas, properties and applications are more detailed in the paper:

[Izbicki R.,Cabezas L. M. C., Colugnatti F. A. B., Lassance R. F. L., de
Souza A. A. L., Stern R. B. (2023). Rethinking Hypothesis Tests. arXiv
preprint arXiv:2308.09112](https://arxiv.org/abs/2308.09112).

## Installation

You can install the development version of REACT from
[GitHub](https://github.com/Monoxido45/REACT) with:

``` r
# install.packages("devtools")
devtools::install_github("Monoxido45/REACT")
```

## Base test example

Any extended simple hypothesis can be tested through REACT by providing
a confidence interval, a chosen tolerance and the original simple
hypothesis in the *base_test* function. For example, if the extended
hypothesis is of the form $H_0: |\mu_1 - \mu_2| \le tol$,

``` r
library(REACT)
## REACT t-test

set.seed(125)
obs1 <- rnorm(n = 30, mean = 1)
obs2 <- rnorm(n = 30, mean = 1.1)

# building confidence set
ci <- t.test(obs1, obs2, var.equal=TRUE, conf.level = 0.95)$conf.int
# tolerance
tol <- 1.5

# performing base test and getting the output
test <- base_test(ci, tol = 1.5, hyp = 0, verbose = TRUE)
#> REACT results:
#> Pragmatic lower bound:  -1.50
#> Pragmatic upper bound:  1.50
#> Confidence interval:
#> lower bound: -0.55
#> upper bound: 0.598
#> REACT conclusion:
#> Based on the provided confidence interval we accept the null hypothesis.
```

We can also plot the CI compared to the region of equivalence:

``` r
plot(test)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## Pairwise comparisons example

One can also perform pairwise comparisons of multiple parameters while
maintaining logical coherence by using *m_comparisons*. To do this, we
assume that the estimators were obtained through MLE and use Fisher’s
information as follows:

``` r
# vector of point estimations
par <- c(1.1, 3.5, 1.5)

# fisher matrix
var_cov <- diag(c(0.05, 0.05, 0.1), nrow = 3, ncol = 3)

# tolerance
tol <- 1.5

REACT::m_comparisons(alpha = 0.05, nrow = 1, ncol = 3,
                     tol = tol, par = par, f_matrix = var_cov)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
