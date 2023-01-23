rm(list=ls())
library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))
source("function.R")

# No effect of treatment, no frailty, no timetrend, and infinite delta
simNoEffNoFrailtyNoTTInfDelta <- simfunction(nsim = 1000)
round(apply(simNoEffNoFrailtyNoTTInfDelta, 2, mean), 2)                         # Everything works
round(exp(mean(simNoEffNoFrailtyNoTTInfDelta[,4])), 2)

# Effect of treatment, no frailty, no timetrend, and infinite delta
simEffNoFrailtyNoTTInfDelta <- simfunction(HR = 2, nsim = 100)
round(apply(simEffNoFrailtyNoTTInfDelta, 2, mean), 2)                           # crude SR = 1 due to infinite delta. Adjusted SR is weird.
round(exp(mean(simEffNoFrailtyNoTTInfDelta[,4])), 2)

# Effect of treatment, no frailty, no timetrend, small delta
simEffNoFrailtyNoTTSmallDelta <- simfunction(HR = 2, nsubs = 1200, nsim = 1000)
round(apply(simEffNoFrailtyNoTTSmallDelta, 2, mean), 2)                         # Everything works
round(exp(mean(simEffNoFrailtyNoTTSmallDelta[,4])), 2)

# No effect of treatment, no frailty, timetrend, and infinite delta
simNoEffNoFrailtyTTInfDelta <- simfunction(l1 = 2, nsim = 1000)
round(apply(simNoEffNoFrailtyTTInfDelta, 2, mean), 2)                           # Bias due to time-trend. Adjusted makes sense here as only case
round(exp(mean(simNoEffNoFrailtyTTInfDelta[,4])), 2)

# No effect of treatment, no frailty, timetrend and small delta
simNoEffNoFrailtyTTSmallDelta <- simfunction(l1 = 2, nsubs = 1200, nsim = 1000)
round(apply(simNoEffNoFrailtyTTSmallDelta, 2, mean), 2)                         # Works due to small delta, adjustment doesn't matter
round(exp(mean(simNoEffNoFrailtyTTSmallDelta[,4])), 2)

# No effect of treatment, frailty, timetrend, and infinite delta
simNoEffFrailtyTTInfDelta <- simfunction(HR = 1, l1 = 2, frailtyVar = 1, nsim = 1000)
round(apply(simNoEffFrailtyTTInfDelta, 2, mean), 2)                             # Works but now Cox is broken due to unmeasured confounding
round(exp(mean(simNoEffFrailtyTTInfDelta[,4])), 2)


# No effect of treatment, frailty, timetrend, and small delta
simNoEffFrailtyTTSmallDelta <- simfunction(HR = 1, l1 = 2, nsubs = 1200, frailtyVar = 1, nsim = 1000)
round(apply(simNoEffFrailtyTTSmallDelta, 2, mean), 2)                           # Works but now Cox is broken due to unmeasured confounding
round(exp(mean(simNoEffFrailtyTTSmallDelta[,4])), 2)


# Effect of treatment, frailty, timetrend, and small delta
simEffFrailtyTTSmallDelta <- simfunction(HR = 2, l1 = 2, nsubs = 1200, frailtyVar = 1, nsim = 1000)
round(apply(simEffFrailtyTTSmallDelta, 2, mean), 2)                             # Works but now Cox is broken due to unmeasured confounding
round(exp(mean(simEffFrailtyTTSmallDelta[,4])), 2)

# Competing risk
simComp <- simfunction(nsim = 1000, comprisk = TRUE, l1 = 4)
round(apply(simComp, 2, mean), 2)                                               # Works but now Cox is broken due to unmeasured confounding
round(exp(mean(simComp[,4])), 2)


# 
# # Effect of treatment, no frailty, timetrend, small delta
# simEffNoFrailtyTTSmallDelta <- simfunction(HR = 2, l2 = 5, delta = 0.005)
# apply(simEffNoFrailtyTTSmallDelta, 2, mean)                                     # More or less the same conclusion as above

