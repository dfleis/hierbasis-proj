#==================================================#
#
# Literally just replicate the examples listed in 
# the HierBasis readme
#
# Univariate case
#
#==================================================#

#===== load libraries =====#
#devtools::install_github("asadharis/HierBasis")
library(HierBasis)

#===== parameters =====#
set.seed(124)
n <- 250
nl <- 20
sigma <- 0.35
beta <- c(0, -4, -1, 3, 9, -6, -1, 7, 9)

K <- length(beta) - 1

#===== generate data =====#
x <- sort(runif(n, -1, 1))
X <- sapply(0:K, function(d) x^d)
eps <- rnorm(n, 0, sigma)

y <- X %*% beta + eps

#===== fit model =====#
pt <- proc.time()
fit <- HierBasis(x, y, nlam = nl)
proc.time() - pt


#===== figures =====#
xseq <- seq(-1, 1, length.out = 1e3)
Xseq <- sapply(0:K, function(d) xseq^d)
ytrue <- Xseq %*% beta # noiseless mean function

# visualize fit for k^th value of lambda
k <- 9
plot(NA, xlim = range(xseq), ylim = range(y), 
     xlab = 'x', ylab = 'y')
abline(h = 0, col = 'gray60', lty = 'longdash')
abline(v = 0, col = 'gray60', lty = 'longdash')
lines(ytrue ~ xseq, lwd = 1.5, col = 'gray20')
points(y ~ x, type = 'p', pch = 21, bg = 'white', cex = 0.55)
lines(fit$fitted.values[,k] ~ x, lwd = 1.5, col = 'red')



LAM <- matrix(rep(fit$lambdas, each = n), nrow = n)
plot(NA, xlim = range(fit$lambdas), ylim = range(fit$beta), log = 'x',
     xlab = expression(lambda), ylab = expression(hat(beta)))
abline(h = 0, col = 'gray60')
matlines(t(LAM), t(fit$beta))


plot(fit$active ~ fit$lambdas, pch = 21, bg = 'white',
     ylab = "Nonzero Coefficients", 
     xlab = expression(lambda))





