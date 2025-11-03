# data to set up model simulations -------------------------------------------------------------
# load packages
library(deSolve)

setwd("E:/fdu/PhD project/Infectious disease/spatial/part2/codes/Mobility and Dengue")
# load climate data
load("data/climateData.RData")

# set immigration and emmigration rate
#ie <- 0.01
ie <- 0

# set up list of sites
sites <- c("Xishuangbanna")

# set human population numbers for each site 
population <- c(1333000)

# set birth and death rates
BRs <- c(6.89) # birth rates per 1000 persons
DRs <- c(7.27) # death rates per 1000 persons

# model timestep
timestep = 1/12

