#===== libraries =====#
library(hierbasis2)
library(glmnet)
#=====================#



#===== parameters =====#
set.seed(420)
n <- 1000
p <- 5
theta <- 0 # covariance parameter

#===== generate data =====#
Z <- matrix(rnorm(n * p), ncol = p)
SIGMA <- outer(1:p, 1:p, function(i, j) theta^(abs(i - j)))
C <- chol(SIGMA)
X <- Z %*% C

round(cov(X), 4)









