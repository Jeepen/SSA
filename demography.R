# ---------------------------------------------------
#
# Author: Jeppe
# Date: 2023-01-24
#  
# Script Description: Simulation playground and demography
#
# ---------------------------------------------------


# Parameters --------------------------------------------------------------
library(MortalityTables)

# All cause mortality -----------------------------------------------------
mortalityTables.load("Austria_Annuities_EROMF")
deathhazards <- EROF.G1950.female@deathProbs
a <- exp(coef(lm(log(deathhazards) ~ I(0:100)))[1])
b <- coef(lm(log(deathhazards) ~ I(0:100)))[2]
