rc <- rn <- ra <- numeric(1e4)
for(i in 1:1e4){
  x <- rexp(1e4, rate = 2)
  y <- rexp(1e4, rate = 1)
  z <- rexp(1e4, rate = 1)
  d <- data.frame(x = x, y = y)[pmax(x,y) < z,]
  dsub <- d
  rc[i] <- with(dsub, sum(x < y) / sum(y < x))
  
  helpers1 <- helpers2 <- numeric(100)
  for(j in 1:100){
    dsample <- d
    dsample[,2] <- sample(d[,2])
    helpers1[j] <- with(dsample, sum(y>x))
    helpers2[j] <- with(dsample, sum(x>y))
  }
  rn[i] <- sum(helpers1) / sum(helpers2)
  ra[i] <- rc[i] / rn[i]
}
summary(ra)