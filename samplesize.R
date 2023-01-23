n <- 1e3
nsim <- 1e4
odds <- numeric(nsim)
for(i in 1:nsim){
  y <- rbinom(n, 1, .5)
  odds[i] <- mean(y) / (1-mean(y))
}
sd(log(odds))
2 / sqrt(n)


power <- function(effect = NULL, SE = NULL, power = NULL){
  if(is.null(power)){
    return(pnorm(-qnorm(.975), mean = -abs(effect) / SE) * 2 )
  }
  if(is.null(SE)){
    return(uniroot(function(SE) pnorm(-qnorm(.975), mean = -abs(effect) / SE) * 2 - power, interval = c(0,100))$root)
  }
}

samplesize <- function(effect, power){
  4 / power(effect = effect, power = power)^2
}

samplesize(effect = .1, power = .9)

power(effect = .1, SE = 2 / sqrt(1200))
