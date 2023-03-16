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
simple_effect <- simfunction(n = 1e4, nsim = 1000, HR2 = 0, HR = 2)
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

# Short window ------------------------------------------------------------
unmeasured_confounding_shortwindow <- simfunction(n = 1e4, nsim = 1000, HR2 = 1, frailtyVar = 1, HR = 2, delta = 1)
round(apply(unmeasured_confounding_shortwindow, 2, mean), 2)



