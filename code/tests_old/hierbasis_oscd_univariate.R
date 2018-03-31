#========================================#
#
# Computing the one-step coordinate
# descent algorithm outlined in
# Haris et al. 2018
#
#========================================#
#===== libraries/external tools =====#
library(CVXR)
source("~/projects/matrix-plot-tools/code/image.mat.R")

#===== functions =====#
# mean response function
ftrue <- function(x) {
  x #sin(2 * pi * x)
}

redblue_grad <- colorRampPalette(
  c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25)))

#===== parameters =====#
set.seed(124)
n <- 300
K <- 25; m <- 3
sigma <- 0.1
lams <- 10^seq(-6, 3, length.out = 50)

#===== generate data =====#
x <- matrix(sort(runif(n, -1, 1)), nrow = n)
eps <- matrix(rnorm(n, sd = sigma), nrow = n)
y <- ftrue(x) + eps

plot(y ~ x, pch = 21, bg = 'white')

#===== DO WORK =====#
PSI <- sapply(1:K, function(k) x^k)



s <- svd(PSI)
U <- sqrt(n) * s$u
V <- tcrossprod(diag(s$d), s$v)/sqrt(n)

#### Do one step coordinate descent
Uty <- crossprod(U, y)
w <- (1:K)^m - (1:K - 1)^m

beta_hats <- vector(mode = 'list', length = length(lams))
pt <- proc.time()
for (j in 1:length(beta_hats)) {
  # initialize
  betas <- sapply(1:K, function(i) Uty/n)
  
  for (k in K:1) {
    a <- max((1 - w[k] * lams[j]/sqrt(sum((betas[k:K,k])^2))), 0)
    betas[k:K, k - 1] <- a * betas[k:K,k]
  }
  
  beta_hats[[j]] <- betas[,1]
}
proc.time() - pt

fhat <- function(z, nbasis = K, lam_idx) {
  PSIz <- sapply(1:nbasis, function(k) z^k)
  PSIz %*% beta_hats[[lam_idx]]
}

#===== figures =====#
xs <- seq(min(x), max(x), length.out = 5e2)

# image.mat(PSI, legend = T)
# image.mat(cor(PSI), legend = T)
# image.mat(abs(cor(PSI)), legend = T)

# plot fits
lam_grad <- redblue_grad(length(lams))
plot(NA, xlim = range(x), ylim = range(y),
     xlab = 'x', ylab = 'y')
lines(ftrue(xs) ~ xs, lwd = 2, col = 'darkblue')
points(y ~ x, pch = 21, cex = 0.75, bg = 'white')
for (j in 1:length(lams)) {
  lines(fhat(xs, lam_idx = j) ~ xs, lwd = 1.5, col = lam_grad[j])
}

# plot coefficient paths as a function of lambda
beta_mat <- matrix(unlist(beta_hats), nrow = K)
plot(NA, xlim = range(lams), ylim = range(beta_hats),
     xlab = expression(lambda), ylab = expression(hat(beta)),
     log = 'x')
for (k in 1:nrow(beta_mat)) {
  lines(beta_mat[k,] ~ lams, lwd = 1.5, col = 'darkblue')
}

image.mat(beta_mat,
          xlab = expression(lambda),
          ylab = expression(abs(hat(beta))))








