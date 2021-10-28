rm(list=ls())
set.seed(26102021)
n <- 1e4
l1 <- 5
l2 <- 10
rho <- .5
t1 <- .1
t2 <- .5

u1 <- rnorm(n)
u2 <- u1 * rho + sqrt(1-rho^2) * rnorm(n)
x <- qexp(pnorm(u1), rate = l1)
y <- qexp(pnorm(u2), rate = l2)
d <- data.frame(x, y)
d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < t2, ]
rc <- with(d, sum(x < y) / sum(y < x))

a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
rn <- a / (1 - a)
rc / rn
summary(glm((x<y) ~ offset(rep(log(rn), nrow(d))), family = "binomial", data = d))


integrate(function(t) pexp(t, rate = l1) * dexp(t, rate = l2), lower = t1, upper = t2)$value  /
integrate(function(t) pexp(t, rate = l2) * dexp(t, rate = l1), lower = t1, upper = t2)$value

(integrate(function(t) pexp(t, rate = l2, lower.tail=FALSE) * dexp(t, rate = l1), lower = t1, upper = t2)$value -
                 pexp(t2, rate = l2, lower.tail=FALSE) * (pexp(t2, rate = l1) - pexp(t1, rate = l1))) /
    (integrate(function(t) pexp(t, rate = l2) * dexp(t, rate = l1), lower = t1, upper = t2)$value -
       pexp(t1, rate = l1) * (pexp(t2, rate = l1) - pexp(t1, rate = l1)))



## Artikel
HR <- 1
u <- runif(n)
X <- rexp(n, rate = l1)
Y <- numeric(n)
cond <- u < (1-exp(-l2*X)) 
Y[cond] <- -log(1-u[cond]) / l2
Y[!cond] <- (l2 * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2 * HR)

d <- data.frame(x = X, y = Y)
d <- d[d$x < tau & d$y < tau, ]
rc <- with(d, sum(x < y) / sum(y < x))

a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
rn <- a / (1 - a)
rc / rn
summary(glm((x<y) ~ offset(rep(log(rn), nrow(d))), family = "binomial", data = d))


rc <- sum(X<Y) / sum(Y>X)
a <- sum(sapply(1:n, function(i) sum(Y>X[i]))) / (n^2)
rn <- a / (1 - a)
rc / rn
summary(glm((X<Y) ~ offset(rep(log(rn), n)), family = "binomial"))










