
#===== functions =====#
logit <- function(z) log(z/(1 - z))
expit <- function(z) exp(z)/(1 + exp(z))

classloss <- function(y, yhat, type = NULL, avg = T) {
  loss <- vector(mode = 'numeric', length = length(y))
  
  if (type == "01" || is.null(type)) {
    loss <- yhat != y
  } else if (type == "logistic") {
    loss <- 1/log(2) * log(1 + exp(-y * yhat))
  } else if (type == "hinge") {
    loss <- pmax(0, 1 - y * yhat)
  } else if (type == "entropy") {
    s <- (1 + y)/2 
    loss <- -s * log(yhat) + (1 - s) * log(1 - yhat)
  } else if (type == "square") {
    loss <- (1 + y * yhat)^2
  } else {
    stop("Error in classloss: No loss function name provided.")
  }
  
  if (avg) {
    mean(loss)
  } else {
    loss
  }
}

#===== parameters =====#
n <- 500
ntrn <- 400
beta <- c(0, -2)

#===== generate data =====#
p <- length(beta)
X <- matrix(rnorm(n * (p - 1)), nrow = n)
eta <- cbind(1, X) %*% beta
pi <- expit(eta)

y <- rbinom(n = n, size = 1, prob = pi)

#===== models =====#
m_glm <- glm(y ~ X, family = binomial)
betahat <- coef(m_glm)
etahat <- cbind(1, X) %*% betahat
pihat <- expit(etahat)
pihat


type <- "entropy"
classloss(ytst, yhat_01, type = type)
classloss(ytst, yhat_11, type = type)
classloss(ytst, yhat_12, type = type)



#===== figures =====#
xp <- X[,1]
plot(y ~ xp, pch = 19, cex = 1, col = rgb(0, 0, 0, 0.5))
lines(yhat_01[order(xp)] ~ sort(xp), col = rgb(1, 0, 0, 0.75), 
      lwd = 5)
lines(yhat_11[order(xp)] ~ sort(xp), col = rgb(0, 1, 0, 0.75), 
      lwd = 5, lty = 'dashed')
lines(yhat_12[order(xp)] ~ sort(xp), col = rgb(0, 0, 1, 0.75), 
      lwd = 5, lty = 'dotdash')

  
  





