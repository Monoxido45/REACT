
<!-- README.md is generated from README.Rmd. Please edit that file -->

# REACT: Region of Equivalence Agnostic Confidence-based Test

<!-- badges: start -->
<!-- badges: end -->

The goal of REACT is to perform equivalence three-way hypothesis testing
in a user-friendly manner. In its simplest form REACT consists in:

<!-- enrich this steps more later -->

1.  **Establishing the region of equivalence.**
2.  **Build a Confidence Set.**
3.  **Test null hypothesis using confidence set following a three-way
    rule.**

<!-- Continue writing more later -->

## Installation

You can install the development version of REACT from
[GitHub](https://github.com/Monoxido45/REACT) with:

``` r
# install.packages("devtools")
devtools::install_github("Monoxido45/REACT")
```

## Base test example

One can conduct any REACT simple hypothesis testing by providing a
confidence interval, a chosen tolerance and the original simple
hypothesis in the *base_test* function:

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

We can have three kinds of conclusions.
