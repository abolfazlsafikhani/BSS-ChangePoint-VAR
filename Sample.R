
rm(list=ls())
gc()

######## Loading Packages #######################
library("mvtnorm")
library("lattice")
library("Rcpp")
library("RcppArmadillo")

######## Call Functions #########################
source("Functions_BSS.R")
#cpp code for block fused lasso and block lasso
sourceCpp("Functions_BSS.cpp")

########## General Settings ######################
T <- (10^4); #number of time points
p <- 10; # number of time series components
brk <- c(floor(T/3),floor(2*T/3),T+1); # true break points with T+1 as the last element
m0 <- length(brk) -1; # number of break points
q.t <- 1; # the true AR order

########## Phi Gneration #########################
m <- m0+1
phi.full <- matrix(0, p, p*q.t*m)
aa <- 0.8
set.seed(123456)
for(mm in 1 : m){
  phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)] <- 0
  for (j in 1:(p-1)){
    bool_1 <- sample(0:2,1,prob = c(0.2,0.6,0.2))
    x_shift = sample(0:4,2)
    if (bool_1 > 0 &&  (j +  x_shift[1:bool_1] <= p) ){
      phi.full[j,((mm-1)*q.t*p+j +  x_shift[1:bool_1])] <- -aa
    }
  }
  if(mm %% 2 == 0 ){
    phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)]  <- -phi.full[,((mm-1)*q.t*p+1):((mm)*q.t*p)]
  }
}  

#display the true AR coefficients matrice
print(plot.matrix((phi.full), m))

#####################################################
########## Data Generation       ####################
#####################################################
set.seed(123456)
e.sigma <- as.matrix(1*diag(p));
try=var.sim.break(T, arlags=q.t, malags=NULL, phi=phi.full, sigma=e.sigma, brk = brk)
data <- try$series
data <- as.matrix(data)

######################################################
######## block segmentation scheme (BSS)    ##########
######################################################
#run the bss method
ptm <- proc.time()
temp <- bss(data, refine = FALSE)
proc.time() - ptm

#display the estimated break points
print("Estimated break points:")
print(temp$final.selected.points)
#display the true break points
print("True break points:")
print(brk[-length(brk)])

#run the bss method
ptm <- proc.time()
temp.2 <- bss(data, refine = TRUE)
proc.time() - ptm

#display the estimated break points
print("Estimated break points:")
print(temp.2$final.selected.points)
#display the true break points
print("True break points:")
print(brk[-length(brk)])

