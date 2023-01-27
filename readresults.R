# ---------------------------------------------------
#
# Author: Jeppe
# Date: 2023-01-26
#  
# Script Description: Look at results from simulation
#
# ---------------------------------------------------
rm(list = ls())
library(rstudioapi)
setwd(dirname(getSourceEditorContext()$path))

d <- readRDS("meanresults2.rds")
subset(d, var == 0 & HR2 == 5 & HR == 1)
subset(d, HR == 2 & delta == 1)
subset(d, HR == 1 & lambda == .036 & var == 0 & delta == Inf)
subset(d, HR2 == 1 & HR == 1 & lambda == .036 & var == 0)
subset(d, lambda == 0.036)
