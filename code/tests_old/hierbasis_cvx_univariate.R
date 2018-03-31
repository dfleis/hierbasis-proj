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
  x^2 #sin(2 * pi * x)
}

# reformulated heirbasis penalty (see equation 14)
heirbasis_penalty14 <- function(betatilde, nbasis, lambda, m = 2) {
  w <- (1:nbasis)^m - ((1:nbasis) - 1)^m
  
  lnorms <- lapply(1:nbasis, function(k) w[k] * norm2(betatilde[k:nbasis]))
  lambda * sum(Reduce("+", lnorms))
}

# reformulated heirbasis penalty (see equation 15)
heirbasis_penalty15 <- function(beta, nbasis, lambda, m = 2) {
  # SOMETHING WRONG HERE... check relationship between beta and betatilde
  w <- (1:nbasis)^m - ((1:nbasis) - 1)^m
  
  lnorms <- lapply(1:nbasis, function(k) w[k] * norm2(beta[k:nbasis]))
  lambda * sum(Reduce("+", lnorms))
}

redblue_grad <- colorRampPalette(
  c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25)))

#===== parameters =====#
set.seed(124)
n <- 100
K <- 10
sigma <- 0.1
lams <- 10^seq(-8, 0, length.out = 8)

#===== generate data =====#
x <- matrix(sort(runif(n, -1, 1)), nrow = n)
eps <- matrix(rnorm(n, sd = sigma), nrow = n)
y <- ftrue(x) + eps

plot(y ~ x)

#===== DO WORK =====#
PSI <- sapply(1:K, function(k) x^k)

s <- svd(PSI)
U <- sqrt(n) * s$u
V <- tcrossprod(diag(s$d), s$v)/sqrt(n)

### Try CVXR on the reformulated problem
beta <- Variable(K)
Uty <- crossprod(U, y)
loss <- 0.5 * norm2(Uty - beta)^2

beta_hats <- vector(mode = 'list', length = length(lams))
for (i in 1:length(lams)) {
  
  pt <- proc.time()
  obj <- loss + heirbasis_penalty15(beta, nbasis = K, lambda = lams[i], m = 2)
  prob <- Problem(Minimize(obj))
  res <- solve(prob)
  print(c(i, (proc.time() - pt)[3]))
  
  beta_hats[[i]] <- res$getValue(beta)
}

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
plot(NA, xlim = range(x), ylim = c(-5, 2),
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
for (k in 1:ncol(beta_mat)) {
  lines(beta_mat[k,] ~ lams, lwd = 1.5, col = 'darkblue')
}

image.mat(abs(beta_mat),
          xlab = expression(lambda),
          ylab = expression(abs(hat(beta))))








