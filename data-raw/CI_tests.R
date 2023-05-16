# some tests using simple test function
# CI falling inside region
library(REACT)
CI <- c(-0.5, 0.5)
obj <- REACT::base_test(CI, tol = 0.3)
plot(obj)

CI <- c(0, 0.25)
obj <- REACT::base_test(CI, tol = 0.5)
plot(obj)

CI <- c(-0.25, 0.25)
obj <- REACT::base_test(CI, tol = 0.5)
plot(obj)

CI <- c(0.5, 0.75)
obj <- REACT::base_test(CI, tol = 0.3)
plot(obj)
