#===== libraries =====#
library(HierBasis)
source("~/projects/matrix-plot-tools/code/image.mat.R")

#===== functions =====#
rmvnorm <- function(n, p = 1, mu = 0, SIGMA = diag(p)) {
  # Generate n random p-dimensional normally 
  # distributed deviates with mean vector mu 
  # and covariance matrix SIGMA.
  #
  # To do: 
  # Include some PSD tolerances and sanitization 
  # on SIGMA.
  # Confirm whether sweep() is the quickest way 
  # of adding mu to our deviates.
  
  if (length(mu) != 1 & length(mu) != p) {
    stop("Error in rmvnorm: Mean mu must be either a scalar
         or a p-vector.")
  }
  if (nrow(SIGMA) != ncol(SIGMA)) {
    stop("Error in rmvnorm: Covariance matrix SIGMA must be
         a square p * p matrix.")
  }
  if (p != nrow(SIGMA)) {
    stop("Error in rmvnorm: Dimension p does not agree with
         covariance matrix SIGMA.")
  } 
  if (length(mu) != 1 & nrow(SIGMA) != length(mu)) {
    stop("Error in rmvnorm: Dimension of mean vector mu does
         not agree with covariance matrix SIGMA.")
  }
  
  if (p == 1) {
    # unvariate case
    rnorm(n, mu, sqrt(SIGMA))
  } else {
    # multivariate case
    C <- chol(SIGMA)
    Z <- matrix(rnorm(n * p), nrow = n)
    sweep(Z %*% C, 2, mu, "+")
  }
}

cv_ahb_lam <- function(mod, K_folds = 5) {
  # Perform K-fold cross-validation on an additive 
  # hierbasis object over tuning parameters lambda.
  #
  # To do: 
  # Consider other loss functions.
  # Check K_folds is an integer within 2 and N
  
  # extract relevant components
  X <- mod$x; y <- mod$y; lams <- mod$lam
  n <- nrow(X); p <- ncol(X)
  
  # assign observations to cross-validation
  # fold 1, ..., K_folds
  folds <- sample(cut(1:n, breaks = K_folds, labels = F))
  
  # get training subset indices
  train_idx <- lapply(1:K_folds, function(i) !(folds %in% i))
  
  cv_folds_err <- sapply(train_idx, function(trn) {
    tst <- !trn # test subset indices
    
    Xtrn <- X[trn,]; ytrn <- y[trn]
    Xtst <- X[tst,]; ytst <- y[tst]
    
    # compute training model and training error
    mod_cv <- AdditiveHierBasis(x = Xtrn, y = ytrn, nbasis = mod$nbasis,
                                nlam = nlam, max.lambda = max(lams), 
                                lam.min.ratio = min(lams)/max(lams), 
                                m.const = mod$m.const)
    trnerr <- colSums((mod_cv$fitted.values - ytrn)^2)/sum(trn)
    
    # fit test data with training data estimates
    yhattst <- predict(mod_cv, new.x = Xtst)
    
    # compute test error across each lambda
    # (squared error loss)
    tsterr <- colSums((yhattst - ytst)^2)/sum(tst)
    
    list("train" = trnerr, "test" = tsterr)
  })
  cv_folds_trn_err <- do.call("cbind", cv_folds_err["train",])
  cv_folds_tst_err <- do.call("cbind", cv_folds_err["test",])
 
  trn_err <- rowSums(cv_folds_trn_err)/K_folds
  tst_err <- rowSums(cv_folds_tst_err)/K_folds
  
  min_err_idx <- which(tst_err == min(tst_err))
  best_lam <- lams[min_err_idx]
  
  out <- vector(mode = 'list')
  out$train_err    <- trn_err
  out$test_err     <- tst_err
  out$lams         <- lams
  out$best_lam_idx <- min_err_idx
  out$best_lam     <- best_lam
  out
}

fadd <- function(Z, funs) {
  # Generate some additive response given 
  # covariates Z and component functions funs.
  n <- nrow(Z); p <- ncol(Z)
  
  # If there are more functions than covariates,
  # consider on the first p functions
  if (length(funs) > p) 
    funs <- funs[1:p]
  
  yj <- matrix(0, nrow = n, ncol = p)
  
  for (j in 1:length(funs)) {
    yj[,j] <- funs[[j]](Z[,j])
  } 
  
  list("y" = rowSums(yj), "yj" = yj, "funs" = funs)
} 

f1 <- function(z) 2 * sin(2 * z)
f2 <- function(z) z
f3 <- function(z) 0.1 * z^3

#===== parameters =====#
set.seed(124)
n <- 1e3
p <- 2
sigma <- 0.3 # noise dispersion
theta <- 0 # covariate correlation term
funs <- list(f1, f2, f3)
K <- 10
cv_folds <- 10
max_lam <- 0.1
min_lam <- 1e-4
nlam <- 100

#===== generate data =====#
THETA <- outer(1:p, 1:p, FUN = function(i, j) (-1)^(i - j) * theta^abs(i - j))
X <- rmvnorm(n, p, SIGMA = THETA)
eps <- matrix(rnorm(n, sd = sigma), nrow = n)

ytrue_obj <- fadd(X, funs)
y <- ytrue_obj$y + eps

#===== fit model =====#
# fit entire model
pt <- proc.time()
mod <- HierBasis::AdditiveHierBasis(X, y, nbasis = K, nlam = nlam, 
                                    max.lambda = max_lam, 
                                    lam.min.ratio = min_lam/max_lam)
proc.time() - pt

err <- colSums((mod$fitted.values - y)^2)/n

# cross validate across lambda
pt <- proc.time()
cv <- cv_ahb_lam(mod, K_folds = cv_folds)
proc.time() - pt

plot(NA, xlab = expression(lambda), ylab = "Error", log = 'xy',
     xlim = range(cv$lams), 
     ylim = range(cv$train_err))
     #ylim = range(c(cv$test_err, cv$train_err)))
lines(cv$train_err ~ cv$lams, lwd = 2, col = 'darkblue')
lines(cv$test_err ~ cv$lams, lwd = 2, col = 'red')

### extract model components ###
lams <- mod$lam
best_lam <- cv$best_lam
best_lam_idx <- cv$best_lam_idx
betahat <- mod$beta
betahat_arr <- array(as.numeric(betahat), dim = c(K, p, nlam))
yhat <- predict(mod)

best_betahat <- betahat_arr[,,best_lam_idx]
best_yhat <- yhat[,best_lam_idx]



#===== figures =====#

j <- 1
fj <- ytrue_obj$funs[[j]]
xseq <- seq(min(X[,j]), max(X[,j]), length.out = 1e3)
yjseq <- fj(xseq)

plot(NA, xlim = range(X[,j]), ylim = range(y),
     xlab = substitute(paste(X[varidx]),
                       list(varidx = j)), 
     ylab = 'y')
points(y ~ X[,j], pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.35))
lines(yjseq ~ xseq, lwd = 2, col = 'darkblue')
points(best_yhat ~ X[,j], pch = 19, cex = 0.5, col = rgb(1, 0, 0, 0.35))


PSIseq <- sapply(1:mod$nbasis, function(k) xseq^k)
yjseqhat <- PSIseq %*% best_betahat[,j]

plot(NA, xlim = range(xseq), ylim = range(yjseq),
     xlab = substitute(paste(X[varidx]),
                       list(varidx = j)), 
     ylab = 'y')
lines(yjseq ~ xseq, lwd = 2, col = 'darkblue')
lines(yjseqhat ~ xseq, lwd = 2, col = 'red')





# pairs(X, pch = 21, cex = 0.5, bg = 'white') 
# image.mat(THETA,
#           xlab = expression(X[i]),
#           ylab = expression(X[j]),
#           main = "Covariate Correlation Matrix")





































