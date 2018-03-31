#===== libraries =====#
library(HierBasis)
library(glmnet)
library(fields)
source("~/projects/parula/code/parula.R")

#===== functions =====#
cv_ahb_lam <- function(mod, alpha = 0.5, K_folds = 5) {
  # Perform K-fold cross-validation on an additive 
  # hierbasis object over tuning parameters lambda.
  #
  # To do: 
  # Add alternative loss functions for non regression tasks.
  # Check K_folds is an valid integer.
  
  # extract relevant components
  X <- mod$x; y <- mod$y; lams <- mod$lam; nlam <- length(lams)
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
                                m.const = mod$m.const,
                                alpha = alpha)
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
unscale <- function(x, scale_obj) {
  xout <- x * attr(scale_obj, "scaled:scale") + 
    attr(scale_obj, "scaled:center")
  attributes(xout) <- NULL
  xout
}
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
fadd <- function(Z, funs) {
  # Generate some additive response given 
  # covariates Z and component functions funs.
  n <- nrow(Z); p <- ncol(Z)
  
  yj <- matrix(0, nrow = n, ncol = p)
  
  for (j in 1:length(funs)) {
    yj[,j] <- funs[[j]](Z[,j])
  } 
  
  list("y" = rowSums(yj), "yj" = yj, "funs" = funs)
} 

f0 <- function(z) rep(0, length(z))
f1 <- function(z) 5 * sin(2 * z)
f2 <- function(z) -3 * z
f3 <- function(z) 4 * z^3
#f4 <- function(z) exp(z)/(1 + exp(z))

#===== parameters =====#
set.seed(124)
n <- 100
p <- 1e3
sparse_pct <- 0.9 # proportion of zero predictors
sigma <- 0 # noise dispersion
theta <- 0 # control the correlation between predictors
funs_unique <- c(f0, f1, f2, f3)

act_varid <- sample(1:p, size = floor(p * (1 - sparse_pct)), replace = F)
zero_varid <- (1:p)[!((1:p) %in% act_varid)]
act_funid <- sample(2:(length(funs_unique)), size = length(act_varid), replace = T)
zero_funid <- sample(c(1), size = p - length(act_varid), replace = T)

varid <- c(act_varid, zero_varid)
funid <- c(act_funid, zero_funid)[order(varid)]
funs <- funs_unique[funid]

#===== generate data ====#
THETA <- outer(1:p, 1:p, function(i, j) (-1)^(i - j) * theta^abs(i - j))
X <- rmvnorm(n, p, SIGMA = THETA)
eps <- rnorm(n, 0, sigma)

yadd <- fadd(X, funs)
ytrue <- yadd$y
y <- ytrue + eps

trnidx <- sample(1:n, size = floor(n/2), replace = F)
tstidx <- (1:n)[!((1:n) %in% trnidx)]

Xtrn <- X[trnidx,]
Xtst <- X[tstidx,]

Xstrn <- scale(Xtrn)
xbar <- attr(Xstrn, "scaled:center")
sdx <- attr(Xstrn, "scaled:scale")
Xstst <- sapply(1:ncol(X), function(i) (Xtst[,i] - xbar[i])/sdx[i])

ytrn <- y[trnidx]
ytst <- y[tstidx]

ystrn <- scale(ytrn)

#===== fit models =====#
m_lm <- lm(ystrn ~ Xstrn)
yshat_lm <- m_lm$fitted.values
yhat_lm <- unscale(yshat_lm, ystrn)
rmse_lm <- sqrt(mean((yhat_lm - ytst)^2))

cv_glmnet <- cv.glmnet(Xstrn, ystrn, nfolds = 5, alpha = 1)
glmnet_lamidx <- which(cv_glmnet$lambda == cv_glmnet$lambda.min)
yshat_glmnet <- predict(cv_glmnet, newx = Xstst, s = cv_glmnet$lambda.min)
yhat_glmnet <- unscale(yshat_glmnet, ystrn)
rmse_glmnet <- sqrt(mean((yhat_glmnet - y)^2))
plot(cv_glmnet)
rmse_glmnet

alpha <- 0.5
pt <- proc.time()
m_ahb <- AdditiveHierBasis(Xstrn, ystrn, nbasis = 6,
                           alpha = alpha, m.const = 3)
proc.time() - pt
pt <- proc.time()
cv_ahb <- cv_ahb_lam(m_ahb, alpha = alpha, K_folds = 5)
proc.time() - pt
plot(NA, log = 'x', xlab = expression(lambda),
     ylab = "Error",
     xlim = range(cv_ahb$lams), 
     ylim = range(cv_ahb$test_err, cv_ahb$train_err))
lines(cv_ahb$train_err ~ cv_ahb$lams, lwd = 1.5, col = 'darkblue')
lines(cv_ahb$test_err ~ cv_ahb$lams, lwd = 1.5, col = 'red')

yshat_ahb <- predict(m_ahb, new.x = Xstst)[,cv_ahb$best_lam_idx]
yhat_ahb <- unscale(yshat_ahb, ystrn)
rmse_ahb <- sqrt(mean((yhat_ahb - y)^2))


rmse_lm
rmse_glmnet
rmse_ahb

m_glmnet_beta <- cv_glmnet$glmnet.fit$beta[,glmnet_lamidx]
m_ahb_beta <- matrix(m_ahb$beta[,cv_ahb$best_lam_idx], nrow = m_ahb$nbasis)
sum(m_glmnet_beta == 0)/p
sum(apply(m_ahb_beta, 2, function(x) sum(abs(x)) == 0))/p








