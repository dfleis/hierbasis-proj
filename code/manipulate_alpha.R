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
nsims <- 5
alpha <- seq(0, 1, length.out = 21); alpha <- alpha[alpha != 0]
n <- 1500 # number of observations
p <- 9 # number of predictors (excluding intercept)
SNR <- 4 # signal to noise ratio

is_nonzero <- c(rep(T, 4), rep(F, p - 4))

#===== SIMULATE =====#
pt <- proc.time()
sim <- replicate(nsims, {
  #===== generate data =====#
  X <- matrix(rnorm(n * p), ncol = p) # iid normal deviates
  y1 <- 5 * f1(X[, 1])
  y2 <- 3 * f2(X[, 2])
  y3 <- 4 * f3(X[, 3])
  y4 <- 6 * f4(X[, 4])
  ytrue <- 2 + y1 + y2 + y3 + y4
  # compute noise dispersion satisfying the SNR ratio
  sigma2 <- sum(ytrue^2)/(n - 1) * 1/SNR 
  eps <- rnorm(n, 0, sqrt(sigma2))
  # compute perturbed response
  y <- ytrue + eps
  
  ### split data into training and validation sets
  X_train <- X[1:(n/2),]; X_valid <- X[(n/2 + 1):n,]
  y_train <- y[1:(n/2)]; y_valid <- y[(n/2 + 1):n]
  
  #===== fit models =====#
  mods <- vector(mode = "list", length = length(alpha))
  select.tab <- matrix(nrow = p, ncol = length(alpha))
  err <- vector(mode = 'numeric', length = length(alpha))
  
  for (i in 1:length(alpha)) {
    mods[[i]] <- cv.additivehierbasis(X = X_train, y = y_train, alpha = alpha[i])
    
    beta.1se <- coef(mods[[i]], lam.idx = mods[[i]]$lambda.1se.idx)
    yhat.1se <- predict(mods[[i]]$model.fit, new.X = X_valid, 
                        lam.idx = mods[[i]]$lambda.1se.idx)
    beta.select <- apply(beta.1se$X, 2, function(b) sum(abs(b)) != 0)
    
    err[i] <- mean((yhat.1se - y_valid)^2)
    select.tab[,i] <- beta.select
  }
  
  list("select" = select.tab, "err" = err)
})
proc.time() - pt


p_00i <- Reduce("+", lapply(sim["select",], function(ss) apply(ss, 2, function(s) (s == F) & (is_nonzero == F))))/nsims
p_11i <- Reduce("+", lapply(sim["select",], function(ss) apply(ss, 2, function(s) (s == T) & (is_nonzero == T))))/nsims
p_01i <- Reduce("+", lapply(sim["select",], function(ss) apply(ss, 2, function(s) (s == F) & (is_nonzero == T))))/nsims
p_10i <- Reduce("+", lapply(sim["select",], function(ss) apply(ss, 2, function(s) (s == T) & (is_nonzero == F))))/nsims
p_1i <- p_11i + p_00i
p_0i <- p_01i + p_10i

p_00 <- colMeans(p_00i)
p_11 <- colMeans(p_11i)
p_01 <- colMeans(p_01i) # type 2
p_10 <- colMeans(p_10i) # type 1
p_1 <- colMeans(p_1i) # total correct classif
p_0 <- colMeans(p_0i) # total misclassif

v_00 <- apply(p_00i * (1 - p_00i), 2, sum)/p^2
v_11 <- apply(p_11i * (1 - p_11i), 2, sum)/p^2
v_01 <- apply(p_01i * (1 - p_01i), 2, sum)/p^2
v_10 <- apply(p_10i * (1 - p_10i), 2, sum)/p^2
v_1 <- apply(p_1i * (1 - p_1i), 2, sum)/p^2
v_0 <- apply(p_0i * (1 - p_0i), 2, sum)/p^2

par(mfrow = c(2, 2))
plot(NA, main = "Correct Variable\nSelection/Rejection",
     ylab = "Rate",
     xlab = expression(alpha),
     xlim = range(alpha),
     ylim = range(p_1 + sqrt(v_1), p_1 - sqrt(v_1)))
segments(x0 = alpha, 
         x1 = alpha,
         y0 = p_1 - sqrt(v_1),
         y1 = p_1 + sqrt(v_1))
points(p_1 ~ alpha, pch = 19, col = 'red', cex = 0.75)

plot(NA, main = "Incorrect Variable\nSelection/Rejection",
     ylab = "Rate",
     xlab = expression(alpha),
     xlim = range(alpha),
     ylim = range(p_0 + sqrt(v_0), p_0 - sqrt(v_0)))
segments(x0 = alpha, 
         x1 = alpha,
         y0 = p_0 - sqrt(v_0),
         y1 = p_0 + sqrt(v_0))
points(p_0 ~ alpha, pch = 19, col = 'red', cex = 0.75)

plot(NA, main = "Type I Error Rates\nVariable Selection/Rejection",
     ylab = "Rate",
     xlab = expression(alpha),
     xlim = range(alpha),
     ylim = range(p_10 + sqrt(v_10), p_10 - sqrt(v_10)))
segments(x0 = alpha, 
         x1 = alpha,
         y0 = p_10 - sqrt(v_10),
         y1 = p_10 + sqrt(v_10))
points(p_10 ~ alpha, pch = 19, col = 'red', cex = 0.75)

plot(NA, main = "Type II Error Rates\nVariable Selection/Rejection",
     ylab = "Rate",
     xlab = expression(alpha),
     xlim = range(alpha),
     ylim = range(p_01 + sqrt(v_01), p_01 - sqrt(v_01)))
segments(x0 = alpha, 
         x1 = alpha,
         y0 = p_01 - sqrt(v_01),
         y1 = p_01 + sqrt(v_01))
points(p_01 ~ alpha, pch = 19, col = 'red', cex = 0.75)



