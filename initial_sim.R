rm(list=ls())
source("~/Dropbox/phd/symmetry/function.R")

## Simulations
badwindow <- simfunction(window = "bad") # Normal, bad window
write.csv(badwindow, file = "~/Dropbox/phd/symmetry/results/badwindow.csv")
goodwindow <- simfunction() # Normal, good window
write.csv(goodwindow, file = "~/Dropbox/phd/symmetry/results/goodwindow.csv")
goodwindowFrailty <- simfunction(frailty = TRUE) # Normal, good window
write.csv(goodwindowFrailty, file = "~/Dropbox/phd/symmetry/results/goodwindowFrailty.csv")
timeTrend <- simfunction(l2 = 5, HR = 1)
write.csv(timeTrend, file = "~/Dropbox/phd/symmetry/results/timeTrend.csv")

## Artikel
HR <- 2
u <- runif(n)
X <- rexp(n, rate = l1)
Y <- numeric(n)
cond <- u < (1-exp(-l2*X)) 
Y[cond] <- -log(1-u[cond]) / l2
Y[!cond] <- (l2 * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2 * HR)

## Normal, dårligt vindue
l1 <- 10
l2 <- 10
t1 <- 0
t2 <- .001
tt <- quantile(pmax(X,Y), c(.001, .01, .1, .3, .5, .75, .9))
rc <- rn <- ra <- nsub <- matrix(NA, nrow = nsim, ncol = length(tt))
for(i in 1:nsim){
    cat(i, "\n")
    u <- runif(n)
    X <- rexp(n, rate = l1)
    Y <- numeric(n)
    cond <- u < (1-exp(-l2*X)) 
    Y[cond] <- -log(1-u[cond]) / l2
    Y[!cond] <- (l2 * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2 * HR)
    for(j in 1:length(tt)){
        d <- data.frame(x = X, y = Y)
        d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < tt[j], ]
        nsub[i,j] <- nrow(d)
        rc[i,j] <- with(d, sum(x < y) / sum(y < x))
        d <- d[1:min(1000, nsub[i,j]),]
        a <- sum(sapply(1:nrow(d), function(i) sum(d$y > d$x[i]))) / (nsub[i,j]^2)
        rn[i,j] <- a / (1 - a)
        ra[i,j] <- rc[i,j] / rn[i,j]
    }
    cat("crude ", rc[i,], "\n")
    cat("null ", rn[i,], "\n")
    cat("adjusted ", ra[i,], "\n")
}
cat("rc", apply(rc, 2, mean), "\n")
cat("rn", apply(rn, 2, mean), "\n")
cat("ra", apply(ra, 2, mean), "\n")
cat("nsub", apply(nsub, 2, mean), "\n")
print(xtable(matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)), include.rownames=FALSE)

mat <- matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)
rownames(mat) <- c("rc", "rn", "ra", "nsub")
colnames(mat) <- c(.001, .01, .1, .3, .5, .75, .9)
write.csv(mat, file = "~/Dropbox/phd/symmetry/results/badwindow.csv")

## Normal, godt vindue
l1 <- 10
l2 <- 10
t1 <- 0
t2 <- .001
tt <- quantile(pmax(X,Y), c(.001, .01, .1, .3, .5, .75, .9))
rc <- rn <- ra <- nsub <- matrix(NA, nrow = nsim, ncol = length(tt))
for(i in 1:nsim){
    cat(i, "\n")
    u <- runif(n)
    X <- rexp(n, rate = l1)
    Y <- numeric(n)
    cond <- u < (1-exp(-l2*X)) 
    Y[cond] <- -log(1-u[cond]) / l2
    Y[!cond] <- (l2 * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2 * HR)
    for(j in 1:length(tt)){
        d <- data.frame(x = X, y = Y)
        d <- d[pmin(d$x, d$y) > (pmax(d$x,d$y) - tt[j]), ]
        nsub[i,j] <- nrow(d)
        rc[i,j] <- with(d, sum(x < y) / sum(y < x))
        d <- d[1:min(1000, nsub[i,j]),]
        a <- sum(sapply(1:nrow(d), function(i) sum(d$y > d$x[i]))) / (nsub[i,j]^2)
        rn[i,j] <- a / (1 - a)
        ra[i,j] <- rc[i,j] / rn[i,j]
    }
    cat("crude ", rc[i,], "\n")
    cat("null ", rn[i,], "\n")
    cat("adjusted ", ra[i,], "\n")
}
cat("rc", apply(rc, 2, mean), "\n")
cat("rn", apply(rn, 2, mean), "\n")
cat("ra", apply(ra, 2, mean), "\n")
cat("nsub", apply(nsub, 2, mean), "\n")
print(xtable(matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)), include.rownames=FALSE)

mat <- matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)
rownames(mat) <- c("rc", "rn", "ra", "nsub")
colnames(mat) <- c(.001, .01, .1, .3, .5, .75, .9)
write.csv(mat, file = "~/Dropbox/phd/symmetry/results/goodwindow.csv")



## Frailty, godt vindue
HR <- 2
u <- runif(n)
frailty <- rgamma(n, shape = 1, scale = 1)
l1 <- frailty * 10
l2 <- frailty * 10
X <- rexp(n, rate = l1)
Y <- numeric(n)
cond <- u < (1-exp(-l2*X)) 
Y[cond] <- -log(1-u[cond]) / l2[cond]
Y[!cond] <- (l2[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2[!cond] * HR)

l1 <- 10
l2 <- 10
t1 <- 0
t2 <- .001
tt <- quantile(pmax(X,Y), c(.001, .01, .1, .3, .5, .75, .9))
rc <- rn <- ra <- nsub <- matrix(NA, nrow = nsim, ncol = length(tt))
for(i in 1:nsim){
    cat(i, "\n")
    HR <- 2
    u <- runif(n)
    frailty <- rgamma(n, shape = 1, scale = 1)
    l1 <- frailty * 10
    l2 <- frailty * 10
    X <- rexp(n, rate = l1)
    Y <- numeric(n)
    cond <- u < (1-exp(-l2*X)) 
    Y[cond] <- -log(1-u[cond]) / l2[cond]
    Y[!cond] <- (l2[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2[!cond] * HR)
    for(j in 1:length(tt)){
        d <- data.frame(x = X, y = Y)
        d <- d[pmin(d$x, d$y) > (pmax(d$x,d$y) - tt[j]), ]
        nsub[i,j] <- nrow(d)
        rc[i,j] <- with(d, sum(x < y) / sum(y < x))
        d <- d[1:min(1000, nsub[i,j]),]
        a <- sum(sapply(1:nrow(d), function(i) sum(d$y > d$x[i]))) / (nsub[i,j]^2)
        rn[i,j] <- a / (1 - a)
        ra[i,j] <- rc[i,j] / rn[i,j]
    }
    cat("crude ", rc[i,], "\n")
    cat("null ", rn[i,], "\n")
    cat("adjusted ", ra[i,], "\n")
}
cat("rc", apply(rc, 2, mean), "\n")
cat("rn", apply(rn, 2, mean), "\n")
cat("ra", apply(ra, 2, mean), "\n")
cat("nsub", apply(nsub, 2, mean), "\n")
print(xtable(matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)), include.rownames=FALSE)

mat <- matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)
rownames(mat) <- c("rc", "rn", "ra", "nsub")
colnames(mat) <- c(.001, .01, .1, .3, .5, .75, .9)
write.csv(mat, file = "~/Dropbox/phd/symmetry/results/goodwindowFrailty.csv")


HR <- 2
u <- runif(n)
frailty <- rgamma(n, shape = 1, scale = 1)
l1 <- frailty * 10
l2 <- frailty * 10
X <- rexp(n, rate = l1)
Y <- numeric(n)
cond <- u < (1-exp(-l2*X)) 
Y[cond] <- -log(1-u[cond]) / l2[cond]
Y[!cond] <- (l2[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2[!cond] * HR)


tt <- quantile(pmax(X,Y), c(.001, .01, .1, .3, .5, .75, .9))
rc <- rn <- ra <- matrix(NA, nrow = nsim, ncol = length(tt))
for(i in 1:nsim){
    cat(i, "\n")
    for(j in 1:length(tt)){  
        d <- data.frame(x = X, y = Y)
        d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < tt[j], ]
        rc[i,j] <- with(d, sum(x < y) / sum(y < x))
        d <- d[1:min(10000, nrow(d)),]
        a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
        rn[i,j] <- a / (1 - a)
        ra[i,j] <- rc[i,j] / rn[i,j]
    }
    cat("crude ", rc[i,], "\n")
    cat("null ", rn[i,], "\n")
    cat("adjusted ", ra[i,], "\n")
}
rc
rn
ra

## Cox model
d <- do.call("rbind", lapply(1:1e4, function(i){
    if(X[i] > Y[i]) data.frame(id = i, time1 = 0, time2 = Y[i], status = 1, treat = 0)
    else data.frame(id = i, time1 = c(0, X[i]), time2 = c(X[i], Y[i]), status = c(0,1), treat = c(0,1))
}))
summary(coxph(Surv(time1, time2, status) ~ treat, data = d))

## for(i in 1:6){
##     cat(i, "\n")
##     d <- data.frame(x = X, y = Y)
##     d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < tt[i], ]
##     d <- d[1:min(10000, nrow(d)),]
##     ## a <- nullSR(d$x, d$y)
##     a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
##     rn[i] <- a / (1 - a)
##     ra[i] <- rc[i] / rn[i]
##     cat("null ", rn[i], "\n")
##     cat("adjusted ", ra[i], "\n")
## }
rc
rn
ra

d <- data.frame(x = X, y = Y)
d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < t2, ]
(rc <- with(d, sum(x < y) / sum(y < x)))



a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
rn <- a / (1 - a)
rc / rn
exp(confint(glm((x<y) ~ offset(rep(log(rn), nrow(d))), family = "binomial", data = d)))

integrate(function(t) dexp(t, rate = l2 * HR), t1, t2)$value / integrate(function(t) dexp(t, rate = l2), t1, t2)$value



## Independent
X <- rexp(n, rate = 10)
Y <- rexp(n, rate = 5)
d <- data.frame(x = X, y = Y)
d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < t2, ]
(rc <- with(d, sum(x < y) / sum(y < x)))
sum(X<Y) / sum(Y<X)


## Old
Z1 <- rexp(n, l2*HR)
Z0 <- rexp(n, l2)
sum(Z1 < .005) / sum(Z0 < .005)




rc <- sum(X<Y) / sum(Y>X)
a <- sum(sapply(1:n, function(i) sum(Y>X[i]))) / (n^2)
rn <- a / (1 - a)
rc / rn
coef(glm((X<Y) ~ offset(rep(log(rn), n)), family = "binomial"))










