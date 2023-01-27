rm(list=ls())
library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))
source("function.R")
source("demography.R")
set.seed(25012023)

# Parameters --------------------------------------------------------------
deltas <- c(Inf, 10, 1.5)
lambda1 <- c(.01, .02)
frailtyVar <- c(0,1)
HR <- c(1,2)                                                                    # Between X and Y
HR2 <- c(0, 1, 5)                                                               # Between X and Z, note that 0 is a hack - it results in Z=Inf
totalscenarios <- length(deltas) * length(lambda1) * length(frailtyVar) * length(HR) * length(HR2)

# All simulations ---------------------------------------------------------
results <- list()
i <- 1
starttime <- Sys.time()
for(hr in HR){
  for(hr2 in HR2){
    for(lambda in lambda1){
      for(vars in frailtyVar){
        for(delta in deltas){
          cat(i, "\n")
          results[[i]] <- list(name = data.frame(delta = delta, HR2 = hr2, lambda = lambda, var = vars, HR = hr), 
                               ests = simfunction(n = 1e4, nsim = 1000, l1 = lambda, HR = hr, frailtyVar = vars, delta = delta, HR2 = hr2))
          cat("Expected time left:", (totalscenarios-i) * difftime(Sys.time(), starttime, units = "mins") / i, "\n")
          cat("Total time so far:", difftime(Sys.time(), starttime, units = "mins"), "\n")
          cat("Expected total time:", difftime(Sys.time(), starttime, units = "mins") + 
                (totalscenarios-i) * difftime(Sys.time(), starttime, units = "mins") / i, "\n")
          cat("Results", apply(results[[i]]$ests, 2, mean), "\n")
          i <- i + 1
        }
      }
    }
  }
}
endtime <- Sys.time()
difftime(endtime, starttime)
saveRDS(results, "allresults2.rds")

# Save results ------------------------------------------------------------
smallresults <- meanresults <- data.frame()
for(i in 1:length(results)){
  smallresults <- rbind(smallresults, results[[i]]$name)
  meanresults <- rbind(meanresults, apply(results[[i]]$ests, 2, mean))
}
res <- cbind(smallresults, meanresults)
names(res)[6:9] <- c("rc", "rn", "ra", "nsub")
saveRDS(res, "meanresults2.rds")
