#===== libraries =====#
library(HierBasis)
library(glmnet)
library(fields)

#===== functions =====#
logit <- function(z) log(z/(1 - z))
expit <- function(z) exp(z)/(1 + exp(z))

classif_err <- function(y, yhat) {
  #mean(u != v)
  mean(y != yhat)
  #mean(log(1 + exp(-y * yhat)))
}
l2loss <- function(z) {
  sum(z^2)
}
l2nloss <- function(z) {
  mean(z^2)
}

cv_ahb_lams <- function(mod, K_folds = 5, alpha, tol, max.iter) {
  # Perform K-fold cross-validation on an additive 
  # hierbasis object over tuning parameters lambda.
  #
  # To do: 
  # Check K_folds is an valid integer.
  
  # extract relevant components
  X    <- mod$x
  y    <- mod$y
  lams <- mod$lam
  nlam <- length(lams)
  n    <- nrow(X)
  p    <- ncol(X)
  type <- mod$type
  
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
    mod_cv <- AdditiveHierBasis(x = Xtrn, y = ytrn, 
                                nbasis        = mod$nbasis,
                                nlam          = nlam, 
                                max.lambda    = max(lams), 
                                lam.min.ratio = min(lams)/max(lams), 
                                m.const       = mod$m.const,
                                type          = type,
                                tol           = tol,
                                max.iter      = max.iter,
                                alpha         = alpha)
    # fit test data with training data estimates
    yhattst0 <- predict(mod_cv, new.x = Xtst)
    
    # compute errors
    if (type == "binomial") {
      yhattrn <- ifelse(mod_cv$fitted.values > 0.75, 1, 0)
      yhattst <- ifelse(yhattst0 > 0.75, 1, 0)
      trnerr <- apply(yhattrn, 2, function(yh) classif_err(ytrn, yh))
      tsterr <- apply(yhattst, 2, function(yh) classif_err(ytst, yh))
    } else {
      warning("test")
      #trnerr <- apply(yhattrn, 2, function(yh) l2nloss(yh - ytrn))
      #tsterr <- apply(yhattst, 2, function(yh) l2nloss(yh - ytst))
    }
    
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
copyscale <- function(x, scale_obj) {
  m <- attr(scale_obj, "scaled:center")
  s <- attr(scale_obj, "scaled:scale")
  
  sapply(1:ncol(x), function(i) (x[,i] - m[i])/s[i])
}


#===== load & restructure data =====#
load(file = "data/colitis.rda")

dat <- data.frame(y = colitis$y, X = t(colitis$x))
colnames(dat) <- c("y", paste0("X", 1:nrow(colitis$x)))

# subset only colitis or crohn's patients
dat2 <- dat[dat$y != 1,] 

#===== exploratory work =====#
X <- as.matrix(dat2[,2:ncol(dat2)])
y <- dat2$y - 2 # binary {0, 1}

xbars <- colMeans(X)
xvars <- apply(X, 2, var)

hist(log10(xbars), breaks = 1e2)
hist(log10(xvars), breaks = 1e2)
plot(sort(log10(xvars)), cex = 0.5, pch = 21, bg = 'white')

xvarid <- which(xvars >= rev(sort(xvars))[100])
X <- X[,xvarid]

#===== parameters =====#
n <- length(y)
p <- ncol(X)
ntrain <- 55
nsplits <- 30
K_folds <- 5
alpha_glmnet <- 1 # 1 = lasso, 0 = ridge
alpha_ahb <- 0.5 # in [0, 1] controlling the balance between the sparsity penalty and the hierarchical penalty

splits <- lapply(1:nsplits, function(i) {
  trnid <- sample(1:n, size = ntrain, replace = F)
  tstid <- (1:n)[!((1:n) %in% trnid)]
  
  train <- list(X = X[trnid,], y = y[trnid])
  test  <- list(X = X[tstid,], y = y[tstid])
  
  list(train = train, test = test)
})

i <- 1
Xtrn <- splits[[i]]$train$X
ytrn <- splits[[i]]$train$y
Xtst <- splits[[i]]$test$X
ytst <- splits[[i]]$test$y

Xstrn <- scale(Xtrn)
Xstst <- copyscale(Xtst, Xstrn)

pt <- proc.time()
cv_glmnet <- cv.glmnet(x = Xstrn, y = ytrn, 
                       family = "binomial", 
                       alpha = alpha_glmnet, 
                       nfolds = K_folds)
proc.time() - pt

best_lam_glmnet <- cv_glmnet$lambda.min
lamid_glmnet <- which(cv_glmnet$lambda == best_lam_glmnet)

yhat0_glmnet <- predict(cv_glmnet, newx = Xstst, s = best_lam_glmnet)
yhat_glmnet <- ifelse(yhat0_glmnet > 0.75, 1, 0)

plot(ytst - yhat_glmnet)
plot(ytst)
points(yhat_glmnet, col = 'red')
plot(cv_glmnet)

cv_glmnet$beta



alpha_ahb <- 0.5
tol <- 2.5e-4
max.iter <- 1e3
ytrn[ytrn == 0] <- 0.5

pt <- proc.time()
m_ahb <- AdditiveHierBasis(x = Xstrn, y = ytrn,
                           type = "binomial",
                           nbasis = 6,
                           max.lambda = 0.2,
                           lam.min.ratio = 1e-1,
                           m.const = 3,
                           tol = tol, 
                           max.iter = max.iter,
                           alpha = alpha_ahb)
proc.time() - pt
m_ahb$fitted.values
err <- apply(ifelse(m_ahb$fitted.values < 0.75, 1, 0), 2, function(yh) classif_err(ytrn, yh))
plot(err ~ m_ahb$lam, log = 'x')

colgrad <- colorRampPalette(c("red", "blue"))(ncol(m_ahb$fitted.values))
plot(NA, xlim = range(m_ahb$lam), ylim = c(0, 1), log = 'x')
for (i in 1:nrow(m_ahb$fitted.values)) {
  lines(m_ahb$fitted.values[i,] ~ m_ahb$lam, lwd = 1.5, col = colgrad[i])
}

nz_ahb_basis <- apply(m_ahb$beta, 2, function(b) {
  apply(matrix(b, nrow = m_ahb$nbasis), 2, function(x) sum(x != 0))
})
nz_ahb <- apply(nz_ahb_basis, 2, function(x) sum(x != 0))
image.plot(1:ncol(X), -log10(m_ahb$lam), nz_ahb_basis)
plot(nz_ahb ~ m_ahb$lam, log = 'x', type = 'l', lwd = 2, col = 'darkblue')

image.plot(as.matrix(abs(m_ahb$beta)))


pt <- proc.time()
cv_ahb <- cv_ahb_lams(mod = m_ahb, K_folds = K_folds, 
                      tol = tol, max.iter = max.iter,
                      alpha = alpha_ahb)
proc.time() - pt

lamid_ahb <- cv_ahb$best_lam_idx[3]
cv_ahb$best_lam_idx

yhat0_ahb <- predict(m_ahb, new.x = Xstst)[,lamid_ahb]
yhat0_ahb
yhat_ahb <- ifelse(yhat0_ahb > 0.75, 1, 0)

plot(cv_ahb$train_err ~ cv_ahb$lams, 
     ylim = range(cv_ahb$train_err, cv_ahb$test_err),
     log = 'x',
     type = 'l', lwd = 2, col = 'darkblue')
lines(cv_ahb$test_err ~ cv_ahb$lams, lwd = 2, col = 'red')

# glmnet results
1 - cv_glmnet$nzero[lamid_glmnet]/p # sparsity
mean(ytst != yhat_glmnet)
# ahb results
1 - nz_ahb[lamid_ahb]/p # sparsity
mean(ytst != yhat_ahb)


plot(1 - nz_ahb/p ~ m_ahb$lam, type = 'l', log = 'x', lwd = 2)
plot(1 - cv_glmnet$nzero/p ~ cv_glmnet$lambda, type = 'l', log = 'x', lwd = 2)

K <- m_ahb$nbasis
X <- m_ahb$x
y <- m_ahb$y
beta0 <- m_ahb$intercept
betahats <- m_ahb$beta

PSIs <- lapply(1:ncol(X), function(j) sapply(1:K, function(k) X[,j]^k))
betahats_list <- lapply(1:ncol(betahats), function(i) {
  matrix(betahats[,i], nrow = K)
})
etas <- sapply(1:ncol(betahats), function(lamidx) {
  PSIbetahat <- sapply(1:ncol(X), function(j) {
    PSIs[[j]] %*% betahats_list[[lamidx]][,j]
  })
  eta <- beta0[lamidx] + rowSums(PSIbetahat)
})
yhats <- apply(etas, 2, function(e) 1/(1 + exp(-e)))

etas_glmnet <- Xstrn %*% as.matrix(cv_glmnet$glmnet.fit$beta)

bahb <- as.matrix(betahats)
bglmn <- as.matrix(cv_glmnet$glmnet.fit$beta)

image.plot(bahb)
image.plot(bglmn)

image.plot(as.matrix(m_ahb$beta))
image.plot(as.matrix(cv_glmnet$glmnet.fit$beta))
image.plot(yhats)
range(betahats)
range(cv_glmnet$glmnet.fit$beta)



image.plot(etas)
range(etas_glmnet)



image.plot(1:nrow(X), 1:ncol(betahats), ifelse(yhats - 0.25 > 0.5, 1, 0))
image.plot(x = 1:ncol(X), 
           z = as.matrix(cv_glmnet$glmnet.fit$beta),
           col = colorRampPalette(c("red", "white", "blue"))(3),
           nlevel = 3,
           breaks = c(-5, -0.001, 0.001, 5))






