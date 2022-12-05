# Parameters --------------------------------------------------------------
library(survival)
l1 <- 5
l2 <- 10
frailtyVar <- 0
nsim <- 1000
set.seed(05122022)


# Sim ---------------------------------------------------------------------
rc <- rn <- ra <- rcox <- nsub <- numeric(nsim)
for (i in 1:nsim) {
  cat(i, "\n")
  if (frailtyVar > 0) {
    frail <- rgamma(n, shape = 1 / frailtyVar, scale = frailtyVar)
    l1frail <- l1 * frail
    l2frail <- l2 * frail
  }
  else{
    l1frail <- l1
    l2frail <- l2
  }
  X <- rexp(n, rate = l1frail)
  Y <- rexp(n, rate = l2frail)
  d <- data.frame(x = X, y = Y)
  dCox <- data.frame(time1 = rep(0, sum(Y<X)), time2 = Y[Y < X], status = 1, treat = 0)
  dCox <- rbind(dCox, data.frame(time1 = rep(0, sum(Y>X)), time2 = X[Y>X], status = 0, treat = 0))
  dCox <- rbind(dCox, data.frame(time1 = X[Y>X], time2 = Y[Y>X], status = 1, treat = 1))
  rcox[i] <- coef(coxph(Surv(time1, time2, status) ~ treat, data = dCox, timefix = FALSE))
  rc[i] <- with(d, sum(x < y) / sum(y < x))
  d <- d[1:min(1000, length(X)),]
  a <- sum(sapply(1:nrow(d), function(w) sum(d$y > d$x[w]))) / (min(1000, length(X)) ^ 2)
  rn[i] <- a / (1 - a)
  ra[i] <- rc[i] / rn[i]
}
mean(rc)
mean(rn)
mean(ra)
exp(mean(rcox))

