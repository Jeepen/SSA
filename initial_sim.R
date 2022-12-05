rm(list=ls())
library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))
source("function.R")

## Simulations
badwindow <- simfunction(window = "bad") # Normal, bad window
write.csv(badwindow, file = "~/Dropbox/phd/symmetry/results/badwindow.csv")
goodwindow <- simfunction() # Normal, good window
write.csv(goodwindow, file = "~/Dropbox/phd/symmetry/results/goodwindow.csv")
goodwindowFrailty <- simfunction(frailty = TRUE) # Normal, good window
write.csv(goodwindowFrailty, file = "~/Dropbox/phd/symmetry/results/goodwindowFrailty.csv")
timeTrend <- simfunction(l2 = 5, HR = 1)
write.csv(timeTrend, file = "~/Dropbox/phd/symmetry/results/timeTrend.csv")
