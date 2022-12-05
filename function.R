simfunction <- function(n = 1e6, nsim = 100, l1 = 10, l2 = 10, t1 = 0, frailty = FALSE, HR = 2, seed = 25032022,
                        window = "good", frailtyVar = 1){
    set.seed(seed)
    l1 <- rep(l1, n)
    l2 <- rep(l2, n)
    if(frailty){
        frail <- rgamma(n, shape = 1/frailtyVar, scale = frailtyVar)            
        l1frail <- l1 * frail
        l2frail <- l2 * frail
    }
    else{
        l1frail <- l1
        l2frail <- l2
    }
    u <- runif(n)
    X <- rexp(n, rate = l1frail)
    Y <- numeric(n)
    cond <- u < (1-exp(-l2frail*X)) 
    Y[cond] <- -log(1-u[cond]) / l2frail[cond]
    Y[!cond] <- (l2frail[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2frail[!cond] * HR)
    # maxXY <- pmax(X,Y)
    # X <- X[maxXY < 1]
    # Y <- Y[maxXY < 1]
    tt <- quantile(pmax(X,Y), c(.0001, .001, .01, .1, .3, .5, .75, .9))
    rc <- rn <- ra <- rcox <- nsub <- matrix(NA, nrow = nsim, ncol = length(tt))
    for(i in 1:nsim){
        cat(i, "\n")
        u <- runif(n)
        if(frailty){
            frail <- rgamma(n, shape = 1/frailtyVar, scale = frailtyVar)            
            l1frail <- l1 * frail
            l2frail <- l2 * frail
        }
        else{
            l1frail <- l1
            l2frail <- l2
        }
        X <- rexp(n, rate = l1frail)
        Y <- numeric(n)
        cond <- u < (1-exp(-l2frail*X)) 
        Y[cond] <- -log(1-u[cond]) / l2frail[cond]
        Y[!cond] <- (l2frail[!cond] * X[!cond] * (HR-1) - log(1-u[!cond])) / (l2frail[!cond] * HR)
        # maxXY <- pmax(X,Y)
        # X <- X[maxXY < 1]
        # Y <- Y[maxXY < 1]
        
        # d <- data.frame(tim1 = rep(0, sum(cond)), time2 = Y[cond], status = 1)
        # d <- rbind(d, time0 = rep(0, sum(!cond)), time2 = X[!cond], status = 0)
        
        for(j in 1:length(tt)){
            d <- data.frame(x = X, y = Y)
            ifelse(window == "good", d <- d[(pmax(d$x,d$y) - pmin(d$x, d$y)) < tt[j], ],
                   d <- d[pmin(d$x, d$y) > t1 & pmax(d$x, d$y) < tt[j], ])
            nsub[i,j] <- nrow(d)
            rc[i,j] <- with(d, sum(x < y) / sum(y < x))
            d <- d[1:min(1000, nsub[i,j]),]
            a <- sum(sapply(1:nrow(d), function(w) sum(d$y > d$x[w]))) / (min(1000, nsub[i,j])^2)
            rn[i,j] <- a / (1 - a)
            ra[i,j] <- rc[i,j] / rn[i,j]
        }
        cat("crude ", rc[i,], "\n")
        cat("null ", rn[i,], "\n")
        cat("adjusted ", ra[i,], "\n")
        cat("nsub ", nsub[i,], "\n")
    }
    mat <- matrix(c(apply(rc, 2, mean), apply(rn, 2, mean),apply(ra, 2, mean),apply(nsub, 2, mean)), nrow = 4, byrow = TRUE)
    rownames(mat) <- c("rc", "rn", "ra", "nsub")
    colnames(mat) <- c(.0001, .001, .01, .1, .3, .5, .75, .9)
    mat
}

