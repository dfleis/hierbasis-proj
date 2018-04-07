set.seed(124)

#===== libraries =====#
library(HierBasis)
library(glmnet)
library(scatterplot3d)
library(fields)
source("~/projects/parula/code/parula.R")

#===== load data =====#
# response: `Class` (1 for fraud, 0 otherwise)
dat <- read.csv("data/creditcard.csv", nrow = 1e3)

#===== functions =====#
unscale <- function(x, scale_obj) {
  xout <- x * attr(scale_obj, "scaled:scale") + 
    attr(scale_obj, "scaled:center")
  attributes(xout) <- NULL
  xout
}

agg <- function(x, y) {
  d <- aggregate(x ~ y, FUN = function(z) {
    c(
      "n"      = length(z),
      "min"    = min(z),
      "q1"     = unname(quantile(z, probs = 0.25)),
      "mean"   = mean(z), 
      "med"    = unname(quantile(z, probs = 0.5)), 
      "q3"     = unname(quantile(z, probs = 0.75)),
      "max"    = max(z),
      "sd"     = sd(z),
      "se"     = sd(z)/sqrt(length(z)),
      "iqr"    = unname(diff(quantile(z, probs = c(0.25, 0.75)))),
      "lo.sd"  = mean(z) - sd(z),
      "hi.sd"  = mean(z) + sd(z),
      "lo.se"  = mean(z) - sd(z)/sqrt(length(z)),
      "hi.se"  = mean(z) + sd(z)/sqrt(length(z))
    )
  })
  d$x <- data.frame(d$x)
  d$x$lo.iqr <- d$x$med - 1.5 * d$x$iqr
  d$x$hi.iqr <- d$x$med + 1.5 * d$x$iqr
  list("d" = d, "X" = x)
}

plot.agg <- function(agg.list, ...) {
  agg.obj <- agg.list$d
  
  # quick hack to solve plotting range issue
  plotlim <- quantile(agg.list$X, probs = c(0.05, 0.95))
  
  plot(NA, xlim = range(agg.obj$y), 
       ylim = range(agg.obj$x$lo.iqr, agg.obj$x$hi.iqr, plotlim), ...)
  segments(x0 = agg.obj$y, x1 = agg.obj$y, lwd = 2,
           y0 = agg.obj$x$lo.sd,
           y1 = agg.obj$x$hi.sd)
  segments(x0 = agg.obj$y - 0.1, x1 = agg.obj$y + 0.1, lwd = 1.5,
           y0 = agg.obj$x$lo.sd, 
           y1 = agg.obj$x$lo.sd)
  segments(x0 = agg.obj$y - 0.1, x1 = agg.obj$y + 0.1, lwd = 1.5,
           y0 = agg.obj$x$hi.sd, 
           y1 = agg.obj$x$hi.sd)
  points(agg.obj$y, agg.obj$x$mean, cex = 0.75, pch = 19, col = 'red2')
  points(agg.obj$y, agg.obj$x$med, cex = 0.75, pch = 19, col = 'blue2')
  points(agg.obj$y, agg.obj$x$lo.iqr, cex = 0.75, pch = 21, bg = 'white')
  points(agg.obj$y, agg.obj$x$hi.iqr, cex = 0.75, pch = 21, bg = 'white')
}


cv_ahb_lam <- function(mod, K_folds = 5) {
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

#===== clean data ====#
dat <- dat[dat$Amount != 0,] # remove 0 amounts for log transforms

# define splits into training and testing subsets 
tstidx <- sample(x = 1:nrow(dat), size = nrow(dat)/5, replace = F)
trnidx <- (1:nrow(dat))[!(1:nrow(dat) %in% tstidx)]

X <- as.matrix(dat[,2:29])
y <- log(dat$Amount)

Xtrn <- X[trnidx,]
Xtst <- X[tstidx,]
ytrn <- y[trnidx]
ytst <- y[tstidx]

# standardize
# note: we do not standardize the test response
# since this will (in principle) be left unknown
# to us in practice when generating model
# (i.e. we are using ytst for model comparison)
Xstrn <- scale(Xtrn)
xbar <- attr(Xstrn, "scaled:center")
sdx <- attr(Xstrn, "scaled:scale")
Xstst <- sapply(1:ncol(X), function(i) (Xtst[,i] - xbar[i])/sdx[i])
ystrn <- scale(ytrn)

#===== set parameters =====#
n <- nrow(dat)
ntrn <- length(trnidx)
ntst <- length(tstidx)
nbasis <- 10
K_folds <- 5
alphas <- seq(0, 1, length.out = 10)
ncvsims <- 5

#==== fit models =====#
m_lm <- lm(ystrn ~ -1 + Xstrn)
yshat_lm <- Xstst %*% coef(m_lm)
yhat_lm <- unscale(yshat_lm, ystrn)
rmse_lm <- sqrt(mean((yhat_lm - ytst)^2))

rmse_glmnet_mat <- matrix(NA, nrow = ncvsims, ncol = length(alphas))
for (j in 1:length(alphas)) {
  a <- alphas[j]
  
  rmse_glmnet_mat[,j] <- replicate(ncvsims, {
    cv <- cv.glmnet(x = Xstrn, y = ystrn, nfolds = K_folds, alpha = a)
    yshat_glmnet <- predict(cv, newx = Xstst, s = cv$lambda.1se)
    yhat_glmnet <- unscale(yshat_glmnet, scale_obj = ystrn)
    sqrt(mean((yhat_glmnet - ytst)^2))
  })
}
rmse_glmnet <- apply(rmse_glmnet_mat, 2, function(x) {
  c(se.lo = mean(x) - sd(x)/sqrt(length(x)),
    mean  = mean(x), 
    se.hi = mean(x) + sd(x)/sqrt(length(x)))
})

plot(NA, xlim = range(alphas), ylim = range(rmse_glmnet),
     xlab = expression(alpha), ylab = "RMSE")
lines(rmse_glmnet["se.lo",] ~ alphas, lwd = 1.5, col = 'gray50')
lines(rmse_glmnet["mean",] ~ alphas, lwd = 1.5, col = 'black')
points(rmse_glmnet["mean",] ~ alphas, pch = 19, cex = 0.5, col = 'red')
lines(rmse_glmnet["se.hi",] ~ alphas, lwd = 1.5, col = 'gray50')


cv <- cv.glmnet(x = Xstrn, y = ystrn, nfolds = K_folds, alpha = 1)
plot(cv)

Xstrn2 <- Xstrn; dimnames(Xstrn2) <- NULL
ystrn2 <- ystrn; dimnames(ystrn2) <- NULL
pt <- proc.time()
m_ahb <- AdditiveHierBasis(x = Xstrn2, y = ystrn2, nbasis = 3,
                           alpha = 0.5, m.const = 3)
proc.time() - pt
pt <- proc.time()
cv_ahb <- cv_ahb_lam(m_ahb, K_folds = K_folds)
proc.time() - pt

plot(NA, log = 'x', xlab = expression(lambda),
     ylab = "Error",
     xlim = range(cv_ahb$lams), 
     ylim = range(cv_ahb$test_err, cv_ahb$train_err))
lines(cv_ahb$train_err ~ cv_ahb$lams, lwd = 1.5, col = 'darkblue')
lines(cv_ahb$test_err ~ cv_ahb$lams, lwd = 1.5, col = 'red')

cv_ahb_lam_idx <- cv_ahb$best_lam_idx
yshat_ahb <- predict(m_ahb, new.x = Xstst)[,cv_ahb_lam_idx]
yhat_ahb <- unscale(yshat_ahb, ystrn)
rmse_ahb <- sqrt(mean((yhat_ahb - ytst)^2))

rmse_lm
min(rmse_glmnet["mean",])
rmse_ahb

m_ahb_beta <- matrix(m_ahb$beta[,cv_ahb_lam_idx], ncol = ncol(X))
m_ahb_act_vars <- sum(apply(m_ahb_beta, 2, function(mb) sum(abs(mb)) != 0))
m_ahb_act_vars

#===== figures =====#
# check correlation structure
varsubset <- !(colnames(dat) %in% c("Time", "Amount", "Class"))
dat2 <- dat[,varsubset]

C <- cor(dat2)
image(1:nrow(C), 1:ncol(C), C, col = parula(64))

CU <- C[upper.tri(C)]
mx <- max(CU)
mn <- min(CU)
mxidx <- which(C == mx, arr.ind = T)
mnidx <- which(C == mn, arr.ind = T)
C[mxidx]; C[mnidx]

plot(sort(CU), pch = 21, bg = 'white', cex = 0.5,
     ylab = "corr", ylim = c(-1,1))

hist(CU, breaks = 1e2)









