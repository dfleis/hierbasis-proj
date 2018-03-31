#==============================#
#
# Fitting additive models.
# (solving equations (3) or (4) 
# in Haris et al. 2018)
#
#==============================#
#===== libraries =====#
library(Matrix)
library(CVXR)
library(ggplot2)
source("~/projects/parula/code/parula.R")

#===== functions =====#
f0 <- function(z) rep(0, length(z))
f1 <- function(z) 2 * sin(2 * z)
f2 <- function(z) z
f3 <- function(z) 0.1 * z^3
ftrue <- function(z) z^2

heirbasis_penalty5 <- function(PSI, beta, m = 2) {
  n <- nrow(PSI); K <- ncol(PSI)
  w <- (1:K)^m - ((1:K) - 1)^m
  
  OMEGA <- lapply(1:K, function(k) {
    w[k] * norm2(PSI[,k:K] %*% beta[k:K])
    #w[k] * p_norm(PSI[,k:K] %*% beta[k:K], 2)
  })
  
  sum(Reduce("+", OMEGA))/sqrt(n)
}
heirbasis_penalty14 <- function(betatilde, nbasis, m = 2) {
  w <- (1:nbasis)^m - ((1:nbasis) - 1)^m
  
  lnorm <- lapply(1:K, function(k) {
    w[k] * norm2(betatilde[k:K])
  })
  
  sum(Reduce("+", lnorm))
}
heirbasis_penalty15 <- function(beta, nbasis, m = 2) {
  w <- (1:nbasis)^m - ((1:nbasis) - 1)^m
  
  lnorm <- lapply(1:K, function(k) {
    w[k] * norm2(beta[k:K])
  })
  
  sum(Reduce("+", lnorm))
}


#===== parameters =====#
set.seed(124)
n <- 40 # observations
K <- n
sigma <- 0.25 # noise dispersion
lams <- 10^seq(-3, 1, length.out = 10)

#===== generate data =====#
x <- matrix(sort(runif(n, -1, 1)), nrow = n)
eps <- matrix(rnorm(n, sd = sigma), nrow = n)
y <- ftrue(x) + eps

plot(y ~ x, pch = 21, bg = 'white')

#===== fit =====#
ybar <- mean(y)
yc <- y - ybar

# generate bases of order n
PSI <- sapply(1:K, function(d) x^d)
psibar <- apply(PSI, 2, mean)
PSIc <- as(scale(PSI, scale = FALSE), Class = "dgCMatrix")

beta <- Variable(K)  

myqr <- Matrix::qr(PSIc)
U <- Matrix::qr.Q(myqr) * sqrt(n)
V <- Matrix::qrR(myqr) / sqrt(n)

Utyn <- crossprod(U, yc)/n

loss <- 1/(2 * n) * sum((y - PSI %*% beta)^2)
loss2 <- 1/2 * sum((Utyn - beta)^2)
b0 <- vector(mode = 'numeric', length = length(lams))
beta_hats <- vector(mode = 'list', length = length(lams))
for (i in 1:length(lams)) {
  pt <- proc.time()
  #obj <- loss + lams[i] * heirbasis_penalty5(PSI, beta, m = 2)
  #obj <- loss2 + lams[i] * heirbasis_penalty14(betatilde, nbasis = K, m = 2)
  obj <- loss2 + lams[i] * heirbasis_penalty15(beta, nbasis = K, m = 2)
  prob <- Problem(Minimize(obj))
  res <- solve(prob)
  print(c(i, (proc.time() - pt)[3]))
  
  betatmp <- res$getValue(beta)
  beta_hats[[i]] <- backsolve(V, betatmp)
  b0 <- ybar - psibar %*% beta_hats[[i]]
}



yhats <- lapply(beta_hats, function(bh) PSI %*% bh)
beta_hats_mat <- matrix(as.numeric(unlist(beta_hats)), 
                        nrow = n, ncol = length(lams))

fhat <- function(z, nbasis, lam_idx) {
  PSI <- sapply(1:nbasis, function(d) z^d)
  PSI %*% beta_hats[[lam_idx]]
}

#===== figures =====#
par(mfrow = c(1,1), pty = 'm')
xs <- seq(min(x), max(x), length.out = 1e3)

# par(mfrow = c(1,1), pty = 'm')
# beta_hats_flip <- t(apply(beta_hats_mat, 2, rev))
# image(lams, 1:n, beta_hats_flip, yaxt = 'n',
#       xlab = expression(lambda), ylab = "Basis Order",
#       col = parula(1e2))
# axis(2, at = 1:n, labels = n:1)

lam_ind <- 1
plot(y ~ x, pch = 21, cex = 0.75, bg = 'white',
     main = substitute(paste(lambda, " = ", lam),
                       list(lam = lams[lam_ind])))
lines(ftrue(x) ~ x, lwd = 2, col = 'darkblue')
lines(fhat(xs, n, lam_ind) ~ xs, col = 'red', lwd = 1.5)


#### PLOT AGAINST LAMBDA ####
plot(NA, xlim = range(lams), ylim = range(beta_hats_mat),
     log = 'x',
     #main = substitute(paste("j = ", j_id), list(j_id = j)),
     ylab = "", #ylab = expression(hat(beta)),
     xlab = expression(lambda))
abline(h = 0, lwd = 1.5, col = 'gray70', lty = 'dashed')
for (i in 1:nrow(beta_hats_mat)) {
  lines(beta_hats_mat[i,] ~ lams, lwd = 1.5, col = 'darkblue')
}



