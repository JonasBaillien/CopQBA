# CopQBA

# NOTE
## Before using anything from Part 2 or Part 3, always run Part 1 as it contains the functions used in the other parts

### Required packages ###
#########################

library(QBAsyDist)
library(plot3D)
library(rgl)
library(misc3d)
library(lattice)
library(cowplot)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(snpar)
library(goftest)
library(copula)
library(nloptr)

### Contents ###
################

# Part 1: functions
## derivatives of the log-copula density of Clayton, Gumbel, AMH, Frank, Joe, Gaussian and Student's t-copulas
### - first and second order derivative to the copula parameter
### - first derivative to the argments u_j   
### - mixed derivate to u_j and the copula parameter                                           

## function to simulate from a given copula with specified QBA-marginals 
## functions for fitting the four univariate quantile based distributions
## function to fit a given copula structure with specified quantile-based margins using IFM
## function to calculate the asymptotic variance-covariance matrix of the univariate QBA-normal, QBA-Laplace, QBA-Student's t and QBA-logistic distribution
## function to calculate the empirical covariance matrix of the parameters estimated by means of IFM of five Archimedean copulas 
## functions to calculate the empirical covariance matrix of the parameters estimated by means of IFM of the Gaussian and Student's t-copula 
## function for determining the best fitting QBA-marginal to fit a marginal of the data
## function to determine the best fitting copula based on AIC using the empirical pseudo-observation 
## construction of the contribution of the copula with given marginals to the CIC

# More information on input/output can be found in comment in the code itself


# Part 2: Simulations
## Everything related to the simulation study, including figures


# Part 3: Data Examples
## Everything related to the data examples, including figures
