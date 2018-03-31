set.seed(124)

#===== libraries =====#
library(HierBasis)
library(glmnet)
library(scatterplot3d)
source("~/projects/parula/code/parula.R")

#===== load data =====#
# response: wine quality (`quality`)
dat <- read.csv("data/winequality-red.csv")
# predictor names
prednames <- names(dat)[!(names(dat) %in% "quality")]

# define splits into training and testing subsets 
tstidx <- sample(x = 1:nrow(dat), size = nrow(dat)/5, replace = F)
trnidx <- (1:nrow(dat))[!(1:nrow(dat) %in% tstidx)]

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

#===== structure data =====#
X <- as.matrix(dat[,prednames])
y <- dat$quality
ytrn <- y[trnidx]
ytst <- y[tstidx]

# standardize
# note: we do not standardize the test response
# since this will (in principle) be left unknown
# to us in practice when generating model
# (i.e. we are using ytst for model comparison)
Xs <- scale(X)
Xstrn <- Xs[trnidx,]
Xstst <- Xs[tstidx,]
ystrn <- scale(ytrn)

#===== set parameters =====#
n <- nrow(dat)
ntrn <- length(trnidx)
ntst <- length(tstidx)
nbasis <- 10

#===== fit models =====#
m_lm <- lm(ystrn ~ Xstrn)
yshat_lm <- cbind(1, Xstst) %*% coef(m_lm)
yhat_lm <- unscale(yshat_lm, scale_obj = ystrn)
sqrt(mean((yhat_lm - ytst)^2))

plot(ytst)
points(yhat_lm, col = 'red')
m_ahb <- AdditiveHierBasis(X, y, nbasis = 10)
head(dat)

p <- 8
agg1 <- agg(dat[,p], dat[,"quality"])
plot.agg(agg1, main = colnames(dat)[p],
         xlab = "Quality Score", ylab = "Level")


p1 <- 2; p2 <- 3
ncols <- length(unique(dat$quality))
cols <- parula(ncols)[as.numeric(cut(dat$quality, breaks = ncols))]
scatterplot3d(x = dat[,p1], 
              y = dat[,p2],
              z = dat$quality, 
              xlab = colnames(dat)[p1],
              ylab = colnames(dat)[p2],
              zlab = "Quality Level",
              pch = 21, cex.symbols = 0.5,
              angle = 40, color = cols)
