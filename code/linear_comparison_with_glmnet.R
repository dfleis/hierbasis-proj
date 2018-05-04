#===== libraries =====#
#install_github("dfleis/hierbasis2")
library(hierbasis2)
library(glmnet)
#library(viridis)
#=====================#

#===== parameters =====#
set.seed(680)

n <- 2500 # number of observations
p <- 9 # number of predictors (excluding intercept)
SNR <- 3 # signal-to-noise ratio (controls noise dispersions)

SPARSE_PCT <- 0.5 # proportion of sparse predictors
SPARSE_IDX <- sample(2:(p + 1), size = floor(SPARSE_PCT * p)) 

beta <- rnorm(p + 1, 0, 10)
beta[SPARSE_IDX] <- 0

lams.ahb <- exp(seq(1, -3, length.out = 50))

#===== generate data =====#
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
pt <- proc.time()
mod.glmnet.cv <- cv.glmnet(x = X_train, y = y_train, alpha = 1)
proc.time() - pt

# additive hierbasis
pt <- proc.time()
mod.ahb.cv <- cv.additivehierbasis(X = X_train, y = y_train)
proc.time() - pt

plot(mod.glmnet.cv)
plot(mod.ahb.cv)



### predict validation sets ###
yhat.glmnet.min <- predict(mod.glmnet.cv, X_valid, s = mod.glmnet.cv$lambda.1se)
yhat.glmnet.1se <- predict(mod.glmnet.cv, X_valid, s = mod.glmnet.cv$lambda.min)
yhat.ahb.min <- predict(mod.ahb.cv, X_valid, lam.idx = mod.ahb.cv$lambda.min.idx)
yhat.ahb.1se <- predict(mod.ahb.cv, X_valid, lam.idx = mod.ahb.cv$lambda.1se.idx)

mse.glmnet.min <- mean((yhat.glmnet.min - y_valid)^2)
mse.glmnet.1se <- mean((yhat.glmnet.1se - y_valid)^2)
mse.ahb.min <- mean((yhat.ahb.min - y_valid)^2)
mse.ahb.1se <- mean((yhat.ahb.1se - y_valid)^2)

mse.glmnet.min
mse.glmnet.1se
mse.ahb.min
mse.ahb.1se

plot(mod.glmnet.cv)
plot(mod.ahb.cv)

coef(mod.glmnet.cv, s = mod.glmnet.cv$lambda.min)
coef(mod.glmnet.cv, s = mod.glmnet.cv$lambda.1se)
plot(coef(mod.ahb.cv, lam.idx = mod.ahb.cv$lambda.min.idx))
coef(mod.ahb.cv, lam.idx = mod.ahb.cv$lambda.1se.idx)
plot


plot(coef(mod.glmnet.cv$glmnet.fit))
plot(mod.ahb.cv$model.fit, pred.idx = 1)
plot(mod.ahb.cv$model.fit, plot.stat = 'active', legend = T)
plot(mod.ahb.cv$model.fit, pred.idx = 1, plot.type = 'image', plot.stat = 'coef', ylab = " ")


beta.glmnet <- as.matrix(coef(mod.glmnet.cv$glmnet.fit))
beta.glmnet <- beta.glmnet[,2:ncol(beta.glmnet)]
plot(beta.glmnet[1,] ~ mod.glmnet.cv$lambda, type = 'l', ylim = range(beta.glmnet))
for(j in 2:p) {
  lines(beta.glmnet[j,] ~ mod.glmnet.cv$lambda)
}
#===== figures =====#
SIGMA.flip <- apply(t(SIGMA), 2, rev)
ax <- pretty(1:p)
rwb_grad <- colorRampPalette(c("blue", "white", "red"))(31)
image.plot(x = 1:p, y = 1:p, z = SIGMA.flip, 
           xlab = "Covariate", 
           ylab = "Covariate",
           zlim = c(-1, 1),
           col = rwb_grad, 
           axes = F,
           main = "Covariance Matrix")
axis(1, at = ax, labels = ax)
axis(2, at = ax, labels = ax)

