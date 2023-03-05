library(survival)
library(MortalityTables)
library(foreach)
library(doParallel)
registerDoParallel(makeCluster(3))
library(doRNG)

rgompertz <- function(n, a = 0.0001755355, b = 0.0649845, x = Inf, HR = 0){
  u <- z <- runif(n)
  limit <- Inf
  if(HR == 0){
    z <- rep(Inf, n)
  }
  else if(all(x < Inf) & HR != 1){
    limit <- 1 - exp(-a/b * (exp(b*x)))
    z[u < limit] <- log(1 - b / a * log(1-u[u < limit])) / b
    z[u > limit] <- log((-b/a*log(1-u[u>limit]) + 1 - exp(b*x[u>limit]) + HR * exp(b*x[u>limit])) / HR) / b
  } 
  else{
    z <-  log(1 - b / a * log(1-u)) / b
  }
  z
}

simsubs <- function(n, l1 = .01, l2 = .01, HR = 1, frailtyVar = 0, delta = Inf, 
                    a = 0.0001755355, b = 0.0649845, HR2 = 1){
  if(frailtyVar > 0){
    frail <- rgamma(n, shape = 1/frailtyVar, scale = frailtyVar)            
    l1frail <- l1 * frail
    l2frail <- l2 * frail
  }
  else{
    l1frail <- rep(l1, n)
    l2frail <- rep(l2, n)
  }
  X <- rexp(n, rate = l1frail)
  Y <- numeric(n)
  u <- runif(n)
  cond <- u < (1-exp(-l2frail*X)) 
  Y[cond] <- -log(1-u[cond]) / l2frail[cond]
  Y[!cond] <- (l2frail[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2frail[!cond] * HR)
  Z <- rgompertz(n, a = a, b = b, x = X, HR = HR2)
  data.frame(x = X[pmax(X,Y) < Z],  
             y = Y[pmax(X,Y) < Z])
}

simfunction <- function(n = 1e5, nsim = 100, l1 = .01, l2 = .01, HR = 1, frailtyVar = 0, delta = Inf,
                        a = 0.0001755355, b = 0.0649845, HR2 = 1){
    rc <- rn <- ra <- rcox <- nsub <- deltas <- p <- numeric(nsim)
    out <- foreach(i = 1:nsim, .combine = "rbind", .export = c("rgompertz", "simsubs")) %dorng% {
        if(frailtyVar > 0){
            frail <- rgamma(n, shape = 1/frailtyVar, scale = frailtyVar)
            l1frail <- l1 * frail
            l2frail <- l2 * frail
        }
        else{
            l1frail <- rep(l1, n)
            l2frail <- rep(l2, n)
        }
        X <- rexp(n, rate = l1frail)
        Y <- numeric(n)
        u <- runif(n)
        cond <- u < (1-exp(-l2frail*X))
        Y[cond] <- -log(1-u[cond]) / l2frail[cond]
        Y[!cond] <- (l2frail[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2frail[!cond] * HR)
        Z <- rgompertz(n, a = a, b = b, x = X, HR = HR2)
        d <- data.frame(x = X[pmax(X,Y) < Z],
                        y = Y[pmax(X,Y) < Z])
      # ntotal <- 0
      # ncounter <- 0
      # d <- NULL
      # while(ncounter < n){
      #   d <- rbind(d, simsubs(n = n, l1 = l1, l2 = l2, HR = HR, frailtyVar = frailtyVar, 
      #                         delta = delta, a = a, b = b, HR2 = HR2))
      #   ntotal <- ntotal + n
      #   ncounter <- with(d, sum(abs(y-x) < delta))
      # }
      # d <- with(d, d[1:match(n, cumsum(abs(x-y) < delta)),])
      # cat("Right?", with(d, sum(abs(x-y) < delta)), "\n")
      dsub <- subset(d, abs(x-y) < delta)
      nsub <- nrow(dsub)
      rc <- with(dsub, sum(x < y) / sum(y < x))

      helpers1 <- helpers2 <- numeric(100)
      for(j in 1:100){
        dsample <- d
        dsample[,2] <- sample(d[,2])
        dsample <- dsample[abs(dsample$y-dsample$x) < delta,]
        helpers1[j] <- with(dsample, sum(y>x))
        helpers2[j] <- with(dsample, sum(x>y))
      }
      rn <- sum(helpers1) / sum(helpers2)
      ra <- rc / rn
      # dCox <- data.frame(time1 = rep(0, sum(Y<X)), time2 = Y[Y < X], status = 1, treat = 0)
      # dCox <- rbind(dCox, data.frame(time1 = rep(0, sum(Y>X)), time2 = X[Y>X], status = 0, treat = 0))
      # dCox <- rbind(dCox, data.frame(time1 = X[Y>X], time2 = Y[Y>X], status = 1, treat = 1))
      # rcox <- coef(coxph(Surv(time1, time2, status) ~ treat, data = dCox, timefix = FALSE))
      # c(rc, rn, ra, rcox, nsub)
        c(rc, rn, ra, nsub)
    }
    out
}

