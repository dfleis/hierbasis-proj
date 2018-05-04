#===== libraries =====#
#install_github("dfleis/hierbasis2")
library(hierbasis2)
library(glmnet)
#library(viridis)
#=====================#

#===== parameters =====#
set.seed(680)
nsims <- 25
n <- 500 # number of observations
p <- 9 # number of predictors (excluding intercept)
SPARSE_PCT <- 0.5 # proportion of sparse predictors
# which predictors are sparse?
SPARSE_IDX <- sample(2:(p + 1), size = floor(SPARSE_PCT * p)) 
SNR <- 3 # signal-to-noise ratio (controls noise dispersions)
beta <- rnorm(p + 1, 0, 10)
beta[SPARSE_IDX] <- 0

lambdas <- exp(seq(1.25, -3.25, length.out = 50))


#===== generate data =====#
pt <- proc.time()
sim <- replicate(nsims, {
  X <- matrix(rnorm(n * p), ncol = p) # iid normal deviates
  # compute noiseless response
  ytrue <- cbind(1, X) %*% beta
  # compute noise dispersion satisfying the SNR ratio
  sigma2 <- sum(ytrue^2)/(n - 1) * 1/SNR 
  eps <- rnorm(n, 0, sqrt(sigma2))
  
  # compute perturbed response
  y <- ytrue + eps
  
  ### split data into training and validation sets ###
  X_train <- X[1:(n/2),]
  X_valid <- X[(n/2 + 1):n,]
  y_train <- y[1:(n/2)]
  y_valid <- y[(n/2 + 1):n]
  
  #===== fit models =====#
  # lasso
  mod.glmnet.cv <- cv.glmnet(x = X_train, y = y_train, alpha = 1)
  # additive hierbasis
  mod.ahb.cv <- cv.additivehierbasis(X = X_train, y = y_train, lambdas = lambdas)
  
  # predict
  yhat.glmnet.min <- predict(mod.glmnet.cv, X_valid, s = mod.glmnet.cv$lambda.1se)
  yhat.glmnet.1se <- predict(mod.glmnet.cv, X_valid, s = mod.glmnet.cv$lambda.min)
  yhat.ahb.min <- predict(mod.ahb.cv, X_valid, lam.idx = mod.ahb.cv$lambda.min.idx)
  yhat.ahb.1se <- predict(mod.ahb.cv, X_valid, lam.idx = mod.ahb.cv$lambda.1se.idx)
  
  # compute errors
  err.glmnet.min <- yhat.glmnet.min - y_valid
  err.glmnet.1se <- yhat.glmnet.1se - y_valid
  err.ahb.min <- yhat.ahb.min - y_valid
  err.ahb.1se <- yhat.ahb.1se - y_valid
  
  list("glmnet.min" = err.glmnet.min,
       "glmnet.1se" = err.glmnet.1se,
       "ahb.min" = err.ahb.min,
       "ahb.1se" = err.ahb.1se)
})
proc.time() - pt

glmnet.min <- sapply(sim["glmnet.min",], function(x) mean(x^2))
glmnet.1se <- sapply(sim["glmnet.1se",], function(x) mean(x^2))
ahb.min <- sapply(sim["ahb.min",], function(x) mean(x^2))
ahb.1se <- sapply(sim["ahb.1se",], function(x) mean(x^2))

d.glmnet.min <- density(glmnet.min)
d.glmnet.1se <- density(glmnet.1se)
d.ahb.min <- density(ahb.min)
d.ahb.1se <- density(ahb.1se)

plot(d.glmnet.min, lwd = 2, lty = 'dashed', col = 'gray60',
     xlim = quantile(c(d.glmnet.min$x, d.glmnet.1se$x,
                  d.ahb.min$x, d.ahb.1se$x), probs = c(0.1, 0.7)),
     ylim = range(d.glmnet.min$y, d.glmnet.1se$y,
                  d.ahb.min$y, d.ahb.1se$y),
     xlab = "Mean-Square Error", main = "MSE Density")
lines(d.glmnet.1se, lwd = 2, lty = 'dotted', col = 'gray60')
lines(d.ahb.min, lwd = 2, lty = 'dashed')
lines(d.ahb.1se, lwd = 2, lty = 'dotted')
legend("topright", legend = c("gnet.min", "gnet.1se", 
                              "HB.min", "HB.1se"),
       lty = c("dashed", "dotted", "dashed", "dotted"),
       col = c("gray60", "gray60", "black", "black"),
       bty = "n", cex = 0.75, lwd = 2)

