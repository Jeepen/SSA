# ---------------------------------------------------
#
# Author: Jeppe
# Date: 2023-01-24
#  
# Script Description: Simulation playground and demography
#
# ---------------------------------------------------


# Parameters --------------------------------------------------------------
library(MortalityTables)
nsim <- 100
n <- 10000
set.seed(24012023)
delta <- Inf

# All cause mortality -----------------------------------------------------
mortalityTables.load("Austria_Annuities_EROMF")
deathhazards <- mort.AT.observed.unisex@deathProbs[,ncol(mort.AT.observed.unisex@deathProbs)]
# deathprobs <- deathSurv <- rep(1, 101)
# deathprobs[1] <- deathhazards[1]
# deathSurv[2] <- 1 - deathprobs[1]
# for(i in 2:101){
#   deathprobs[i] <- deathhazards[i] * deathSurv[i-1]
#   deathSurv[i] <- deathSurv[i-1] - deathprobs[i]
# }
# plot(0:100, deathprobs, type = "l", xlab = "", ylab = "", las = 1)
# Z <- sample(0:100, 1e3, prob = deathprobs, replace = TRUE)
a <- exp(coef(lm(log(deathhazards) ~ I(0:100)))[1])
b <- coef(lm(log(deathhazards) ~ I(0:100)))[2]
dat <- data.frame(trt = rep(0, n))
sim <- simsurv(dist = "gompertz", lambdas = a, gammas = b, betas = c(trt = 0), x = dat)
Z <- sim$eventtime

# Depression --------------------------------------------------------------
# probDepressed <- c(1.1, 5.1, 5.0, 4.7, 4, 2)
# probDepressed <- probDepressed / sum(probDepressed)
# y <- sample(1:6, 1e3, replace = TRUE, prob = probDepressed)
# y[y == 1] <- sample(20:29, sum(y == 1), replace = TRUE)
# y[y == 2] <- sample(30:39, sum(y == 2), replace = TRUE)
# y[y == 3] <- sample(40:49, sum(y == 3), replace = TRUE)
# y[y == 4] <- sample(50:59, sum(y == 4), replace = TRUE)
# y[y == 5] <- sample(60:69, sum(y == 5), replace = TRUE)
# y[y == 6] <- sample(70:99, sum(y == 6), replace = TRUE)
rc <- rn <- ra <- deltas <- nsub <- numeric(nsim)
for(i in 1:nsim){
  cat(i, "\n")
  sim <- simsurv(dist = "gompertz", lambdas = a, gammas = b, betas = c(trt = 0), x = dat)
  Z <- sim$eventtime
  Y <- runif(n, 15, 100)
  X <- runif(n, 40, 100)
  d <- data.frame(x = X[pmax(X,Y) < Z & abs(X-Y) < delta], 
                  y = Y[pmax(X,Y) < Z & abs(X-Y) < delta])
  nsub[i] <- nrow(d)
  rc[i] <- with(d, sum(x < y) / sum(y < x))
  # delta <- with(dsub, max(abs(y-x)))
  helpers1 <- helpers2 <- numeric(100)
  for(j in 1:100){
    dsample <- d
    dsample[,2] <- sample(d[,2])
    if(nsubs != Inf){
      dsample <- dsample[abs(dsample$y-dsample$x) < delta,]
    }
    helpers1[j] <- with(dsample, sum(y>x))
    helpers2[j] <- with(dsample, sum(x>y))
  }
  rn[i] <- sum(helpers1) / sum(helpers2)
  ra[i] <- rc[i] / rn[i]
  cat("crude ", rc[i], "\n")
  cat("null ", rn[i], "\n")
  cat("adjusted ", ra[i], "\n")
}
summary(rc)
summary(rn)
summary(ra)
# summary(deltas)
summary(nsub)

# condMean <- numeric(101)
# for(i in 1:101) condMean[i] <- sum(deathSurv[i:101]) / deathSurv[i]
# plot(0:100, condMean, type = "l", xlab = "", ylab = "", las = 1)
# condMean[30]
# plot(0:100, condMean + (0:100), type = "l", xlab = "", ylab = "", las = 1)
