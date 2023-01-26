# ---------------------------------------------------
#
# Author: Jeppe
# Date: 2023-01-26
#  
# Script Description: Look at results from simulation
#
# ---------------------------------------------------

library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))

d <- readRDS("meanresults.rds")
subset(d, var == 0 & comprisk == FALSE & lambda == .072)
subset(d, HR == 2 & delta == 1)
