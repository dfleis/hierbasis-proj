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
logit <- function(p) log(p/(1 - p))
expit <- function(z) 1/(1 + exp(-z))
rmse <- function(z, zhat) {
  # root-mean-squared error loss
  n <- length(z)
  sqrt(sum((z - zhat)^2)/n)
}
se <- function(z) {
  # standard error
  sd(z)/sqrt(length(z))
}
cv.AdditiveHierBasis <- function(mod, 
                                 nfolds   = 10, 
                                 alpha    = NULL, 
                                 max.iter = 100,
                                 tol      = 1e-04) {
  # Perform n-fold cross-validation on an AdditiveHierBasis 
  # object.
  # Only defined for Gaussian family AdditiveHierBasis models.
  #
  # Note: Since alpha, max.iter, and tol are not stored in
  # the AdditiveHierBasis object it is critical to set the
  # parameters to the same parameters that were used for the
  # input.
  
  # extract relevant components
  X    <- mod$x
  y    <- mod$y
  lams <- mod$lam
  nlam <- length(lams)
  n    <- nrow(X)
  p    <- ncol(X)
  type <- mod$type
  
  # assign observations to cross-validation
  # fold 1, ..., nfolds
  folds <- sample(cut(1:n, breaks = nfolds, labels = F))
  
  # get training subset indices
  train.idx <- lapply(1:nfolds, function(i) !(folds %in% i))
  
  cv.err <- sapply(train.idx, function(trn) {
    tst <- !trn # test subset indices
    
    X.trn <- X[trn,]; y.trn <- y[trn]
    X.tst <- X[tst,]; y.tst <- y[tst]
    
    # compute training model and training error
    mod.cv <- AdditiveHierBasis(x = X.trn, y = y.trn, 
                                nbasis        = mod$nbasis,
                                nlam          = nlam, 
                                max.lambda    = max(lams), 
                                lam.min.ratio = min(lams)/max(lams), 
                                m.const       = mod$m.const,
                                type          = type,
                                tol           = tol,
                                max.iter      = max.iter,
                                alpha         = alpha)
    yhat.trn <- mod.cv$fitted.values
    # fit test data using the trained model
    yhat.tst <- predict(mod.cv, new.x = X.tst)
    
    # compute errors
    if (type == "binomial") {
      warning("Warning: 'type' not yet defined for 'binomial' models.")
    } else {
      trn.err <- apply(yhat.trn, 2, function(yhat) rmse(y.trn, yhat))
      tst.err <- apply(yhat.tst, 2, function(yhat) rmse(y.tst, yhat))
    }
    
    list("train" = trn.err, "test" = tst.err)
  })
  cv.trn.err <- do.call("cbind", cv.err["train",])
  cv.tst.err <- do.call("cbind", cv.err["test",])
  
  # average over all folds 
  # entries correspond to an error for a unique lambda
  trn.err <- rowSums(cv.trn.err)/nfolds
  tst.err <- rowSums(cv.tst.err)/nfolds
  
  tst.err.se <- apply(cv.tst.err, 1, se)
  tst.err.hi <- tst.err + tst.err.se
  tst.err.lo <- tst.err - tst.err.se
  
  min.err.idx <- which(tst.err == min(tst.err))
  min.lam <- lams[min.err.idx]
  
  out <- vector(mode = 'list')
  out$train.err      <- trn.err
  out$test.err       <- tst.err
  out$test.err.se    <- tst.err.se
  out$test.err.hi    <- tst.err.hi
  out$test.err.lo    <- tst.err.lo
  out$lams           <- lams
  out$min.lambda.idx <- min.err.idx
  out$min.lambda     <- min.lam
  out
}



#=========================#
#===== GENERATE DATA =====#
#=========================#
set.seed(124)
n <- 50
p <- 50
PCT_SPARSE <- 0.75 # proportion of zero coefficients
sigma <- 1 # response noise

# DESIGN MATRIX (exclude a column of 1s for the intercept)
X <- matrix(rnorm(n * (p - 1)), nrow = n)
# LINEAR COEFFICIENTS
beta.init <- sample((-1)^(1:p) * (1:p))/p
beta.0 <- sample(c(0, 1), size = p, replace = T, prob = c(PCT_SPARSE, 1 - PCT_SPARSE))
beta <- beta.init * beta.0
# NOISE (response noise)
eps <- rnorm(n, 0, sigma)

# RESPONSE
ytrue <- cbind(1, X) %*% beta
y <- ytrue + eps


#======================#
#===== FIT MODELS =====#
#======================#
nfolds <- 5
pct.train <- 0.5
pct.valid <- 1 - pct.train
alpha.glmnet <- seq(0, 1, length.out = 50)
# remove 0 from alpha.ahb since It caused problems with the plotting later on
alpha.ahb    <- seq(0, 1, length.out = 25); alpha.ahb <- alpha.ahb[alpha.ahb != 0]
nbasis.ahb   <- 10
m.const.ahb  <- 3

#==== split data into training and validation sets ====#
# We we use the training data set to do our cross validation over
# and the validation set to compute the FINAL error estimates.
# Note: This training set is not the same as the training sets
# we produce for cross validation. Our cross validation training
# sets are subsets of this primary training set, and the test sets
# are taken to be the other portion of these subsets.
# After cross-validation is completed with the training and 
# test sets we do a final fit of the validation data set.
# This fit is what we compare our model errors to.

train.idx <- sample(1:n, size = pct.train * n, replace = F)
valid.idx <- (1:n)[!((1:n) %in% train.idx)]

X.train <- X[train.idx,]
X.valid <- X[valid.idx,]
y.train <- y[train.idx]
y.valid <- y[valid.idx]

#==== standardize ====#
# do we want to do this?
# glmnet does it automatically and I don't think we need to
# do it for HierBasis (I don't remember thinking it was a huge
# deal for Hierbasis...)
#X.train.s <- scale(X.train)
#X.valid.s <- copyscale(X.valid, X.train.s)

#==== linear regression (baseline) ====#
# Note: The OLS estimate will not make sense for 
# p > n (infinite solutions)

# fit model
m.OLS <- lm(y.train ~ X.train)
# fit validation responses
beta.OLS <- coef(m.OLS)
yhat.valid.OLS <- cbind(1, X.valid) %*% beta.OLS

#==== glmnet (alpha = 0 => ridge, alpha = 1 => lasso) ====#

yhat.valid.glmnet <- sapply(alpha.glmnet, function(alpha) {
  cv.glmnet <- cv.glmnet(X.train, y.train,
                         alpha       = alpha,
                         nfolds      = nfolds,
                         family      = "gaussian",
                         standardize = T,
                         intercept   = T)
  predict(cv.glmnet, X.valid, s = "lambda.1se")
})

#==== AdditiveHierBasis =====#

yhat.valid.ahb <- matrix(NA, nrow = nrow(X.valid), ncol = length(alpha.ahb)) 
for (i in 1:ncol(yhat.valid.ahb)) {
  alpha <- alpha.ahb[i]
  
  m.ahb <- AdditiveHierBasis(X.train, y.train,
                             nbasis        = nbasis.ahb,
                             alpha         = alpha, 
                             m.const       = m.const.ahb, 
                             type          = "gaussian")
  pt <- proc.time()
  cv.ahb <- cv.AdditiveHierBasis(m.ahb, 
                                 nfolds = nfolds, 
                                 alpha  = alpha)
  tm <- proc.time() - pt
  print(c(i, tm))
  
  yhat.valid.ahb[,i] <- predict(m.ahb, new.x = X.valid)[,cv.ahb$min.lambda.idx[1]]
}


#=====================================#
#===== COMPARE VALIDATION ERRORS =====#
#=====================================#
rmse.OLS <- rmse(y.valid, yhat.valid.OLS)
rmse.glmnet <- apply(yhat.valid.glmnet, 2, function(yhat) rmse(y.valid, yhat))
rmse.ahb <- apply(yhat.valid.ahb, 2, function(yhat) rmse(y.valid, yhat))

plot(rmse.glmnet ~ alpha.glmnet, type = 'o', col = 'green2',
     pch = 19, cex = 0.5, lty = 'dotdash',
     ylim = range(rmse.glmnet, rmse.ahb))
points(rmse.ahb ~ alpha.ahb, type = 'o', col = 'darkblue',
       pch = 21, bg = 'white', cex = 0.75)
abline(h = rmse.OLS, col = 'red', lty = 'dashed')




