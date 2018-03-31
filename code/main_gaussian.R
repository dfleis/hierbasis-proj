#==============================#
#
#
#
#==============================#

#=====================#
#===== LIBRARIES =====#
#=====================#
library(HierBasis)
library(glmnet)

#=====================#
#===== FUNCTIONS =====#
#=====================#

#=========================#
#===== GENERATE DATA =====#
#=========================#
set.seed(124)
n <- 1e2
p <- 1e2
PCT_SPARSE <- 0.75
sigma <- 1 # response noise

# DESIGN MATRIX (add a column of 1s for the intercept)
X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
# LINEAR COEFFICIENTS
beta.init <- sample((-1)^(1:p) * (1:p))
beta.0 <- sample(c(0, 1), size = p, replace = T, prob = c(PCT_SPARSE, 1 - PCT_SPARSE))
beta <- beta.init * beta.0
# NOISE (response noise)
eps <- rnorm(n, 0, sigma)

# RESPONSE
ytrue <- X %*% beta
y <- ytrue + eps







