# some tests using simple test function
# CI falling inside region
CI <- c(-0.5, 0.5)
REACT::base_test(CI, tol = 0.3)

CI <- c(0, 0.25)
REACT::base_test(CI, tol = 0.5)

CI <- c(-0.25, 0.25)
REACT::base_test(CI, tol = 0.5)

CI <- c(0.5, 0.75)
REACT::base_test(CI, tol = 0.3)

