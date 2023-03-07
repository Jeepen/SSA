rm(list=ls())
library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))
source("function.R")
source("demography.R")
set.seed(25012023)


# Simulations -------------------------------------------------------------
# No competing risks, no unmeasured confounding ---------------------------
simple_simple <- simfunction(n = 1e4, nsim = 1000, HR2 = 0)
round(apply(simple_simple, 2, mean), 2)
simple_effect <- simfunction(n = 1e4, nsim = 1000, HR2 = 0, HR = 10)
round(apply(simple_effect, 2, mean), 2)
simple_timetrend <- simfunction(n = 1e4, nsim = 1000, HR2 = 0, l1 = .05)
round(apply(simple_timetrend, 2, mean), 2)

# Competing risks, no unmeasured confounding ------------------------------
competingrisk <- simfunction(n = 1e4, nsim = 1000, HR2 = 1, l1 = .05)
round(apply(competingrisk, 2, mean), 2)
competingrisk_dependent <- simfunction(n = 1e4, nsim = 1000, HR2 = 10, l1 = .01)
round(apply(competingrisk_dependent, 2, mean), 2)

# No competing risks, unmeasured confounding ------------------------------
unmeasured_confounding <- simfunction(n = 1e4, nsim = 1000, HR2 = 0, l1 = .05, frailtyVar = 1)
round(apply(unmeasured_confounding, 2, mean), 2)



effect_competingrisk <- simfunction(n = 1e4, nsim = 1000, HR = 2)
round(apply(effect_competingrisk, 2, mean), 2)
noeffect_competingrisk_timetrend <- simfunction(n = 1e4, nsim = 1000, l1 = .04)
round(apply(noeffect_competingrisk_timetrend, 2, mean), 2)
noeffect_dependentcompetingrisk_timetrend <- simfunction(n = 1e4, nsim = 1000, l1 = .04, HR2 = .2)
round(apply(noeffect_dependentcompetingrisk_timetrend, 2, mean), 2)
noeffect_competingrisk_timetrend_confounding <- simfunction(n = 1e4, nsim = 1000, l1 = .04, frailtyVar = 1)
round(apply(noeffect_competingrisk_timetrend_confounding, 2, mean), 2)
noeffect_competingrisk_timetrend_confounding_smallwindow <- 
  simfunction(n = 1e4, nsim = 1000, l1 = .04, frailtyVar = 1, delta = .5)
round(apply(noeffect_competingrisk_timetrend_confounding_smallwindow, 2, mean), 2)
effect_competingrisk_timetrend_confounding_smallwindow <- 
  simfunction(n = 1e4, nsim = 1000, l1 = .04, frailtyVar = 1, delta = .5, HR = 2)
round(apply(effect_competingrisk_timetrend_confounding_smallwindow, 2, mean), 2)

effect_timetrend_confounding <- 
  simfunction(n = 1e4, nsim = 1000, l1 = .02, frailtyVar = 1, HR2 = 0)
round(apply(effect_timetrend_confounding, 2, mean), 2)


effect_timetrend_confounding_smallwindow <- 
  simfunction(n = 1e4, nsim = 1000, l1 = .04, frailtyVar = 1, delta = .5, HR2 = 0)
round(apply(effect_timetrend_confounding_smallwindow, 2, mean), 2)



effect_smallwindow <- simfunction(n = 1e4, nsim = 1000, HR = 2, delta = 1.5)
apply(effect_smallwindow, 2, mean)


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
                               ests = simfunction(n = 1e2, nsim = 1000, l1 = lambda, HR = hr, frailtyVar = vars, delta = delta, HR2 = hr2))
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
saveRDS(results, "allresults3.rds")

# Save results ------------------------------------------------------------
smallresults <- meanresults <- data.frame()
for(i in 1:length(results)){
  smallresults <- rbind(smallresults, results[[i]]$name)
  meanresults <- rbind(meanresults, apply(results[[i]]$ests, 2, mean))
}
res <- cbind(smallresults, meanresults)
names(res)[6:10] <- c("rc", "rn", "ra", "nsub", "totalsubs")
saveRDS(res, "meanresults3.rds")
