rm(list = ls())

# install.packages("mvtnorm")
# install.packages("nlme")

library(mvtnorm)
library(nlme)
PATH <- getwd()

######## ######## ######## ######## ######## ########
# Simulation I
## SIM1:
n <- 25
Rep <- 1
Repp <- 195
seednum <- n * 1000 + 1000000 + Rep
source(paste(PATH, "/Code/simSEM25.R", sep = ""))

## SIM2:
n <- 50
Rep <- 1
Repp <- 131
seednum <- n * 1000 + 1000000 + Rep
source(paste(PATH, "/Code/simSEM25.R", sep = ""))

## SIM3:
n <- 100
Rep <- 1
Repp <- 113
seednum <- n * 1000 + 1000000 + Rep
source(paste(PATH, "/Code/simSEM25.R", sep = ""))

## SIM4:
n <- 200
Rep <- 1
Repp <- 101
seednum <- n * 1000 + 1000000 + Rep
source(paste(PATH, "/Code/simSEM25.R", sep = ""))

## SIM5:
n <- 400
Rep <- 1
Repp <- 102
seednum <- n * 1000 + 1000000 + Rep
source(paste(PATH, "/Code/simSEM25.R", sep = ""))

######## ######## ######## ######## ######## ########
# Simulation II
## SIM1:
n <- 25
Rep <- 1
Repp <- 181
seednum <- n * 1000 + 2000000 + Rep
source(paste(PATH, "/Code/simSEM50.R", sep = ""))

## SIM2:
n <- 50
Rep <- 1
Repp <- 136
seednum <- n * 1000 + 2000000 + Rep
source(paste(PATH, "/Code/simSEM50.R", sep = ""))

## SIM3:
n <- 100
Rep <- 1
Repp <- 116
seednum <- n * 1000 + 2000000 + Rep
source(paste(PATH, "/Code/simSEM50.R", sep = ""))

## SIM4:
n <- 200
Rep <- 1
Repp <- 103
seednum <- n * 1000 + 2000000 + Rep
source(paste(PATH, "/Code/simSEM50.R", sep = ""))

## SIM5:
n <- 400
Rep <- 1
Repp <- 101
seednum <- n * 1000 + 2000000 + Rep
source(paste(PATH, "/Code/simSEM50.R", sep = ""))

######## ######## ######## ######## ######## ########
# Simulation III
## SIM1:
n <- 25
Rep <- 1
Repp <- 157
seednum <- n * 1000 + 4000000 + Rep
source(paste(PATH, "/Code/simSEM75.R", sep = ""))

## SIM2:
n <- 50
Rep <- 1
Repp <- 152
seednum <- n * 1000 + 4000000 + Rep
source(paste(PATH, "/Code/simSEM75.R", sep = ""))

## SIM3:
n <- 100
Rep <- 1
Repp <- 119
seednum <- n * 1000 + 4000000 + Rep
source(paste(PATH, "/Code/simSEM75.R", sep = ""))

## SIM4:
n <- 200
Rep <- 1
Repp <- 105
seednum <- n * 1000 + 4000000 + Rep
source(paste(PATH, "/Code/simSEM75.R", sep = ""))

## SIM5:
n <- 400
Rep <- 1
Repp <- 105
seednum <- n * 1000 + 4000000 + Rep
source(paste(PATH, "/Code/simSEM75.R", sep = ""))