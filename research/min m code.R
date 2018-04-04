library(devtools)
#library(microbenchmark)

library(glmnet)
library(HierBasis)
#install_github("dfleis/hierbasis2")
library(hierbasis2)

#=====================#
#===== functions =====#
#=====================#
f1 <- function(z) 2 * z
f2 <- function(z) -1 * z^3
f3 <- function(z) 5 * sin(0.25 * pi * z)
f4 <- function(z) -2 * cos(1 * pi * z)

rmse <- function(y, yhat) {
  n <- length(y)
  sqrt(sum((y - yhat)^2)/n)
}

#=========================#
#===== generate data =====#
#=========================#
set.seed(124)
n <- 250
sigma <- 20

x <- sort(runif(n, -3, 3)) # sort just to make plotting easier later
eps <- rnorm(n, 0, sigma)

ytrue <- f1(x) + f2(x) + f3(x) + f4(x)
y <- ytrue + eps

plot(ytrue ~ x, ylim = range(y), pch = 19, cex = 0.25,
     col = rgb(0, 0, 0, 1), type = 'o')
points(y ~ x, pch = 19, col = rgb(0, 0, 0.75, 0.5), cex = 0.5)

#=================================#
#===== plots with varying m  =====#
#=================================#
c <-1 
m.vec <- seq(0,10, 0.1)
store <- matrix(0,nrow=length(m.vec), ncol=2)

for (l in m.vec){
  mod.hb2 <- hierbasis2::hierbasis(x = x, y = y, nbasis = 10, m.const=l)
  yhat.hb2 <- predict(mod.hb2)
  store[c,1] <- l
  store[c,2] <- rmse(ytrue, yhat.hb2)
  c <- c+1
}
plot(store, xlab="m values", ylab="RMSE", pch = 19, col = rgb(0, 0, 0, 0.5),
     cex = 0.5)
lines(store)

#====================================#
#===== min m for varying sigma  =====#
#====================================#
set.seed(124)
d <- 1
sig.vec <- 1:20
out <- matrix(0,nrow=length(sig.vec), ncol=2)

for (s in sig.vec){
c <- 1 
m.vec <- seq(0,10, 0.1)
store <- matrix(0,nrow=length(m.vec), ncol=2)
eps <- rnorm(n, 0, s)
y <- ytrue + eps

  for (l in m.vec){
    mod.hb2 <- hierbasis2::hierbasis(x = x, y = y, nbasis = 10, m.const=l)
    yhat.hb2 <- predict(mod.hb2)
    store[c,1] <- l
    store[c,2] <- rmse(ytrue, yhat.hb2)
    c <- c+1
  }
out[d,1] <- s
out[d,2] <- store[which.min(store[,2]),1]

d <- d+1
}
out
plot(out)
