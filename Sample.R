
rm(list=ls())
gc()

######## Loading Packages #######################
library("mvtnorm")
library("lattice")
library("igraph")
library("pracma")
library("Rcpp")
library("RcppArmadillo")

######## Call Functions #########################
# Need to set the working directory
source("Functions_BSS.R")
#cpp code for block fused lasso and block lasso
sourceCpp("Functions_BSS.cpp")

########## General Settings ######################
T <- (10^3*4); #number of time points
p <- 10; # number of time series components
brk <- c(floor(T/3),floor(2*T/3),T+1); # true break points with T+1 as the last element
m0 <- length(brk) -1; # number of break points
q.t <- 1; # the true AR order
m <- m0+1 #number of segments
sp_density <- c(0.05, 0.05, 0.05) #sparsity level (5%)

#####################################################
########## Data Generation       ####################
#####################################################
#randomly generate transition matrix and timeseries data
try <- simu_var(method = 'sparse', nobs = T, k=p, sp_pattern= 'random', arlags=seq(1,q.t,1),
                brk=brk, sp_density = sp_density , seed = 1 )
print(plot.matrix(do.call("cbind",try$model_param), m ))
data <- try$series
data <- as.matrix(data)

######################################################
######## block segmentation scheme (BSS)    ##########
######################################################
#run the bss method
ptm <- proc.time()
temp <- bss(data)
proc.time() - ptm

#display the estimated break points
print("Estimated break points:")
print(temp$final.selected.points)
#display the true break points
print("True break points:")
print(brk[-length(brk)])

