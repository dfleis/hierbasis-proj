set.seed(124)

#===== libraries =====#
library(HierBasis)
library(glmnet)
library(ade4)

#===== load data =====#
dat <- read.csv("data/train.csv")

#===== functions =====#
cv_ahb_lam <- function(mod, K_folds = 5) {
  # Perform K-fold cross-validation on an additive 
  # hierbasis object over tuning parameters lambda.
  #
  # To do: 
  # Consider other loss functions.
  # Check K_folds is an integer within 2 and N
  
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


unscale <- function(x, scale_obj) {
  xout <- x * attr(scale_obj, "scaled:scale") + attr(scale_obj, "scaled:center")
  attributes(xout) <- NULL
  xout
}
dummy <- function(df) {
  ISFACT <- sapply(df, is.factor)
  FACTS <- acm.disjonctif(df[, ISFACT, drop = FALSE])
  NONFACTS <- df[, !ISFACT,drop = FALSE]
  data.frame(NONFACTS, FACTS)
}

#===== organize data =====#
excludeNA <- c("Alley", "LotFrontage", "FireplaceQu", "PoolQC", 
               "Fence", "MiscFeature")
exclude <- c("Id", "TotalBsmtSF", "GrLivArea", 
             "Electrical", "Exterior2nd", excludeNA) 
dat2 <- dat[,!(colnames(dat) %in% exclude)]

cc <- complete.cases(dat2)
dat3 <- dat2[cc,]

prednames <- colnames(dat2)[!(colnames(dat2) %in% "SalePrice")]
form <- paste0("SalePrice ~ ", paste(prednames, collapse = " + "))
mtmp <- lm(form, data = dat3)

X <- model.matrix(mtmp)
y <- log(dat[cc,"SalePrice"])

tstidx <- sample(1:nrow(X), size = floor(nrow(X)/4), replace = F)
trnidx <- (1:nrow(X))[!(1:nrow(X) %in% tstidx)]

Xs <- X#scale(X)
Xstrn <- Xs[trnidx,]
Xstst <- Xs[tstidx,]

ytrn <- y[trnidx]
ytst <- y[tstidx]
ystrn <- ytrn#scale(ytrn)

#===== fit models =====#
nfolds <- 10
nbasis <- 10
alphas <- seq(0, 1, length.out = 50)

m_lm <- lm(ystrn ~ -1 + Xstrn)
yshat_lm <- Xstst %*% coef(m_lm)
yhat_lm <- yshat_lm# unscale(yshat_lm, scale_obj = ystrn)
rmse_lm <- sqrt(mean((yhat_lm - ytst)^2))
rmse_lm

library(MASS)
m <- lm.ridge(ystrn ~ -1 + Xstrn, lambda = 1)



pt <- proc.time()
cv_glmnet <- lapply(alphas, function(a) {
  cv <- cv.glmnet(Xstrn, ystrn, nfolds = nfolds, alpha = a)
  
  yshat_glmnet_min <- predict(cv, newx = Xstst, s = cv$lambda.min)
  yshat_glmnet_1se <- predict(cv, newx = Xstst, s = cv$lambda.1se)
  
  yhat_glmnet_min <- unscale(yshat_glmnet_min, scale_obj = ystrn)
  yhat_glmnet_1se <- unscale(yshat_glmnet_1se, scale_obj = ystrn)
  
  rmse_glmnet_min <- sqrt(mean((yhat_glmnet_min - ytst)^2))
  rmse_glmnet_1se <- sqrt(mean((yhat_glmnet_1se - ytst)^2))
  
  out <- list()
  out$cv <- cv
  out$yhat_min <- yhat_glmnet_min
  out$yhat_1se <- yhat_glmnet_1se
  out$rmse_min <- rmse_glmnet_min
  out$rmse_1se <- rmse_glmnet_1se
  out
})
proc.time() - pt

rmse_min_glmnet <- sapply(cv_glmnet, function(cv) cv$rmse_min)
rmse_1se_glmnet <- sapply(cv_glmnet, function(cv) cv$rmse_1se)

cv_lasso <- cv.glmnet(Xstrn, ystrn, nfolds = nfolds, alpha = 1)
plot(cv_lasso)
lmin_idx <- which(cv_lasso$lambda == cv_lasso$lambda.min)
l1se_idx <- which(cv_lasso$lambda == cv_lasso$lambda.1se)
b_lasso <- cv_lasso$glmnet.fit$beta[,lmin_idx]
var_lasso <- which(b_lasso != 0)
Xstrn_lasso <- Xstrn[,var_lasso] 
Xstst_lasso <- Xstst[,var_lasso]
m_lm_lasso <- lm(ystrn ~ Xstrn_lasso)
yshat_lasso <- cbind(1, Xstst_lasso) %*% coef(m_lm_lasso)
yhat_lasso <- unscale(yshat_lasso, scale_obj = ystrn)
rmse_lm_lasso <- sqrt(mean((yhat_lasso - ytst)^2))

min(rmse_min_glmnet)
min(rmse_1se_glmnet)
rmse_lm
rmse_lm_lasso

pt <- proc.time()
m_ahb <- AdditiveHierBasis(x = Xstrn, y = ystrn, nbasis = nbasis, 
                         nlam = length(cv_lasso$lambda),
                         max.lambda = 1,
                         lam.min.ratio = 0.001,
                         alpha = 0.5, m.const = 35)
proc.time() - pt
pt <- proc.time()
cv_ahb <- cv_ahb_lam(m_ahb, K_folds = nfolds)
proc.time() - pt

yshat_ahb <- predict(m_ahb, new.x = Xstst)[,cv_ahb$best_lam_idx]
yhat_ahb <- unscale(yshat_ahb, ystrn)
sqrt(mean((yhat_ahb - ytst)^2))

plot(cv_ahb$test_err ~ cv_ahb$lams, pch = 19, cex = 0.5, 
     log = 'xy', col = 'red')

