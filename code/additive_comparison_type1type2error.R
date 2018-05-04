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
nsims <- 10

n <- 50 # number of observations
p <- 100 # number of predictors (excluding intercept)
SNR <- 3 # signal to noise ratio

is_nonzero <- c(rep(T, 4), rep(F, p - 4))
#===== generate data =====#
pt <- proc.time()
sim <- replicate(nsims, {
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
  
  #===== fit models =====#
  # lasso
  mod.glmnet.cv <- cv.glmnet(x = X, y = y, alpha = 1)
  # additive hierbasis
  mod.ahb.cv <- cv.additivehierbasis(X = X, y = y, tol = 1e-03)
  
  sparse_test.glmnet <- as.numeric(coef(mod.glmnet.cv)[2:(p + 1)]) != 0
  sparse_test.ahb <- apply(coef(mod.ahb.cv)$X, 2, function(bcol) sum(abs(bcol)) != 0)
  cbind(sparse_test.glmnet, sparse_test.ahb)
})
proc.time() - pt

# glmnet table for probability of type 1 and type 2 errors
glmnet.tabs <- lapply(1:nsims, function(i) 
  table(is_nonzero, factor(sim[,1,i], levels = c(F, T))))
glmnet.arr <- simplify2array(glmnet.tabs)
glmnet.prob <- apply(glmnet.arr, c(1, 2), sum)/(nsims * p)

# ahb table for probability of type 1 and type 2 errors
ahb.tabs <- lapply(1:nsims, function(i) 
  table(is_nonzero, factor(sim[,2,i], levels = c(F, T))))
ahb.arr <- simplify2array(ahb.tabs)
ahb.prob <- apply(ahb.arr, c(1, 2), sum)/(nsims * p)

glmnet.correct <- sum(diag(glmnet.prob))
glmnet.type1 <- glmnet.prob[1,2] # false pos
glmnet.type2 <- glmnet.prob[2,1] # false neg
ahb.correct <- sum(diag(ahb.prob))
ahb.type1 <- ahb.prob[1,2] # false pos
ahb.type2 <- ahb.prob[2,1] # false neg

prob.correct <- c(glmnet.correct, ahb.correct)
names(prob.correct) <- c("glmnet", "HierBasis")
prob.errs <- cbind(c(glmnet.type1, glmnet.type2), 
                   c(ahb.type1, ahb.type2))
rownames(prob.errs) <- c("Type 1", "Type 2")
colnames(prob.errs) <- c("glmnet", "HierBasis")

round(rbind(prob.correct, prob.errs), 4)


