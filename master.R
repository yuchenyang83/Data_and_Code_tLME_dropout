rm(list=ls())
library(ggplot2)
library(cowplot)
library(rlang)
library(ggtext)
library(glue)
library(gridExtra)
library(grid)
library(pROC)
library(mvtnorm)
library(nlme)
PATH = getwd()

# Re-produce Figure 1
source(paste(PATH, '/Code/Fig1.R', sep=''))

# Re-produce Figure 2
source(paste(PATH, '/Code/Fig2.R', sep=''))

# Re-produce Figure 3
source(paste(PATH, '/Code/Fig3.R', sep=''))

# Re-produce Figure 4
source(paste(PATH, '/Code/Fig4.R', sep=''))

# Re-produce Figure 5
source(paste(PATH, '/Code/Fig5.R', sep=''))

# Re-produce Figure 6
source(paste(PATH, '/Code/Fig6.R', sep=''))



# Re-produce Table 1
source(paste(PATH, '/Code/Tab1.R', sep=''))

# Re-produce Table 2
source(paste(PATH, '/Code/Tab2.R', sep=''))

# Re-produce Table 3
source(paste(PATH, '/Code/Tab3.R', sep=''))

# Re-produce Table 4
source(paste(PATH, '/Code/Tab4.R', sep=''))

# Re-produce Table 5
source(paste(PATH, '/Code/Tab5.R', sep=''))

