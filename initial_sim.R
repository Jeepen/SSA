rm(list=ls())
set.seed(26102021)
n <- 1e5
l1 <- 5
l2 <- 10
rho <- 0
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

integrate(function(t) (pexp(t2, rate = l2) - pexp(t, rate = l2)) * dexp(t, rate = l1), lower = t1, upper = t2)$value /
integrate(function(t) (pexp(t, rate = l2) - pexp(t1, rate = l2)) * dexp(t, rate = l1), lower = t1, upper = t2)$value



## Artikel
HR <- 10
u <- runif(n)
X <- rexp(n, rate = l1*re)
Y <- numeric(n)
cond <- u < (1-exp(-l2*X)) 
Y[cond] <- -log(1-u[cond]) / (l2)
Y[!cond] <- (l2 * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2 * HR)

d <- data.frame(x = X, y = Y)
d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < t2, ]
(rc <- with(d, sum(x < y) / sum(y < x)))

a <- sum(sapply(1:nrow(d), function(i) with(d, sum(y>x[i])))) / (nrow(d)^2)
rn <- a / (1 - a)
rc / rn
summary(glm((x<y) ~ offset(rep(log(rn), nrow(d))), family = "binomial", data = d))




Z1 <- rexp(n, l2*HR)
Z0 <- rexp(n, l2)
sum(Z1 < .2) / sum(Z0 < .2)




rc <- sum(X<Y) / sum(Y>X)
a <- sum(sapply(1:n, function(i) sum(Y>X[i]))) / (n^2)
rn <- a / (1 - a)
rc / rn
coef(glm((X<Y) ~ offset(rep(log(rn), n)), family = "binomial"))










