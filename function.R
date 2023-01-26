library(survival)
library(MortalityTables)

rgompertz <- function(n, alpha = exp(-10), beta = .085){
  u <- runif(n)
  log(1 - beta / alpha * log(1-u)) / beta
}


simfunction <- function(n = 1e5, nsim = 100, l1 = .072, l2 = .072, t1 = 0, t2 = Inf, HR = 1, frailtyVar = 0, delta = Inf, comprisk = FALSE){
    rc <- rn <- ra <- rcox <- nsub <- deltas <- p <- numeric(nsim)
    for(i in 1:nsim){
        # cat(i, "\n")
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
        Z <- rep(Inf, n)
        if(comprisk){
          Z <- rgompertz(n)
          # dat <- data.frame(trt = rep(0,n))
          # sim <- simsurv(dist = "gompertz", lambdas = exp(-10), gammas = 0.085, betas = c(trt = 0), x = dat)
          # Z <- sim$eventtime
        } 
        d <- data.frame(x = X[pmin(X,Y)>t1 & pmax(X,Y)<pmin(t2,Z)], 
                        y = Y[pmin(X,Y)>t1 & pmax(X,Y)<pmin(t2,Z)])
        dsub <- subset(d, abs(x-y) < delta)
        # if(t1 == 0 & t2 == Inf){
          # d <- data.frame(x = X[Y < (2*X)], y = Y[Y < (2*X)])                    # Virker ikke pga. tidstrend
        # }
        # if(nsubs == Inf){
          # dsub <- d
        # }
        # else{
          # dsub <- d[order(abs(d$y-d$x)),][1:nsubs,]
        # }
        # delta <- with(dsub, max(abs(y-x)))
        # deltas[i] <- delta
        nsub[i] <- nrow(dsub)
        rc[i] <- with(dsub, sum(x < y) / sum(y < x))
        # test <- with(d, (sum(x < y) - sum(y < x))^2 / nrow(d))
        # p[i] <- 1 - pchisq(test, df = 1)
        # if(adjust == "good"){
          # da <- data.frame(x = X, y = Y)[1:min(nsub[i], 1000),]
        # }
        # else{
          # da <- d[1:min(nsub[i], 1000),]
        # }
        
        helpers1 <- helpers2 <- numeric(100)
        for(j in 1:100){
          dsample <- d
          dsample[,2] <- sample(d[,2])
          # if(nsubs != Inf){
            dsample <- dsample[abs(dsample$y-dsample$x) < delta,]
          # }
          helpers1[j] <- with(dsample, sum(y>x))
          helpers2[j] <- with(dsample, sum(x>y))
        }
        rn[i] <- sum(helpers1) / sum(helpers2)
        ra[i] <- rc[i] / rn[i]
        dCox <- data.frame(time1 = rep(0, sum(Y<X)), time2 = Y[Y < X], status = 1, treat = 0)
        dCox <- rbind(dCox, data.frame(time1 = rep(0, sum(Y>X)), time2 = X[Y>X], status = 0, treat = 0))
        dCox <- rbind(dCox, data.frame(time1 = X[Y>X], time2 = Y[Y>X], status = 1, treat = 1))
        rcox[i] <- coef(coxph(Surv(time1, time2, status) ~ treat, data = dCox, timefix = FALSE))
        # cat("crude ", rc[i], "\n")
        # cat("null ", rn[i], "\n")
        # cat("adjusted ", ra[i], "\n")
        # cat("cox ", rcox[i], "\n")
        # cat("nsub ", nsub[i], "\n")
    }
    # data.frame(rc = rc, rn = rn, ra = ra, rcox = rcox, delta = deltas, p = p)
    # data.frame(rc = rc, rn = rn, ra = ra, rcox = rcox, delta = deltas)
    data.frame(rc = rc, rn = rn, ra = ra, rcox = rcox, nsub = nsub)
}

