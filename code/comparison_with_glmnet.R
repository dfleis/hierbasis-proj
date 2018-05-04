#===== libraries =====#
#install_github("dfleis/hierbasis2")

library(hierbasis2)
library(glmnet)
#=====================#

#===== functions =====#
f1 <- function(x) x
f2 <- function(x) (2 * x - 1)^2
f3 <- function(x) 2 * sin(2 * pi * x)/(2 - sin(2 * pi * x))
f4 <- function(x) 0.1 * sin(2 * pi * x) + 0.2 * cos(2 * pi * x) + 
  0.3 * sin(2 * pi * x)^2 + 
  0.4 * cos(2 * pi * x)^3 + 
  0.5 * sin(2 * pi * x)^3

#===== data parameters =====#
set.seed(680)
n <- 5000 # number of observations
p <- 9 # number of predictors (excluding intercept)
SNR <- 10 # signal to noise ratio


#===== generate data =====#
X <- matrix(rnorm(n * p), ncol = p) # iid normal deviates
y1 <- 5 * f1(X[, 1])
y2 <- 3 * f2(X[, 2])
y3 <- 4 * f3(X[, 3])
y4 <- 6 * f4(X[, 4])
ytrue <- 2 + y1 + y2 + y3 + y4

# compute noise dispersion satisfying the SNR ratio
sigma2 <- sum(ytrue^2)/(n - 1) * 1/SNR 
eps <- 0#rnorm(n, 0, sqrt(sigma2))

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
mod.ahb.cv <- cv.additivehierbasis(X = X_train, y = y_train, 
                                   lambdas = exp(seq(1.25, -3, length.out = 50)))
proc.time() - pt


b.ahb <- coef(mod.ahb.cv, lam.idx = mod.ahb.cv$lambda.1se.idx)

n_valid <- nrow(X_valid)
K <- mod.ahb.cv$model.fit$nbasis
# each slice of the array has the orthogonal design
# corresponding to each covariate
PSI.array <- array(NA, dim = c(n_valid, K, p))

# generate and center basis exansion PSI over every predictor
# j = 1, ..., p of order 1 to nbasis
for (j in 1:p) {
  # POLYNOMIAL basis expansion
  PSI.array[,, j] <- outer(X_valid[,j], 1:K, "^")
}

j <- 3
plot(y3[(n/2+1):n] ~ X_valid[,j])
yhat <- PSI.array[,,j] %*% b.ahb$X[,j]
points(yhat ~ X_valid[,j], col = 'red')




plot(mod.glmnet.cv)
plot(mod.ahb.cv)

lammin.glmnet.idx <- which(mod.glmnet.cv$lambda == mod.glmnet.cv$lambda.min)
lam1se.glmnet.idx <- which(mod.glmnet.cv$lambda == mod.glmnet.cv$lambda.1se)
lammin.ahb.idx <- mod.ahb.cv$lambda.min.idx
lam1se.ahb.idx <- mod.ahb.cv$lambda.1se.idx

coefmin.glmnet <- coef(mod.glmnet.cv, mod.glmnet.cv$lambda.min)
coef1se.glmnet <- coef(mod.glmnet.cv, mod.glmnet.cv$lambda.1se)
coefmin.ahb <- coef(mod.ahb.cv, mod.ahb.cv$lambda.min)
coef1se.ahb <- coef(mod.ahb.cv, mod.ahb.cv$lambda.1se)

plot(mod.ahb.cv$model.fit,
     plot.stat = "active", plot.type = "image", legend = T)

#===== figures =====#

plot(y1 ~ X[,1], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y2 ~ X[,2], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y3 ~ X[,3], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y4 ~ X[,4], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))

plot(y ~ X[,1], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y ~ X[,2], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y ~ X[,3], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))
plot(y ~ X[,4], pch = 19, cex = 0.75,
     col = rgb(0, 0, 0, 0.5))







