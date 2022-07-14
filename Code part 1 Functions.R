########################################################################
########################################################################
### Code accompagnying the paper:                                    ###
###     "Inference for copulas with two-piece margins"               ###
### Authors: Baillien Jonas, Gijbels irène and Verhasselt Anneleen   ###
########################################################################
########################################################################

#########################
### Part 1: Functions ###
#########################

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




## functions for the ratios of the symmetric reference densities of the margins
## with respect to its first and second derivatives.                           
#################################################################################

# f'(x)/f(x)
#############

ratiofirst=function(x,basefunc,df=NULL){
  if(basefunc=="t" & is.null(df)){
    stop('no degrees of freedom provided')
  }
  out=switch(basefunc,'normal'=-x,
             't'=-(df+1)*x/(df+x^2),
             'logistic'=-(1-exp(-x))/(1+exp(-x)),
             'laplace'=-sign(x)
  )
  return(out)
}


# f''(x)/f(x) 
##############

ratiosecond=function(x,basefunc,df=NULL){
  if(basefunc=="t" & is.null(df)){
    stop('no degrees of freedom provided')
  }
  out=switch(basefunc,'normal'=x^2-1,
             't'=((df+1)^2*x^2-(df-x^2)*(df+1))/(df+x^2)^2,
             'logistic'=((1-exp(-x))/exp(-x)-2)*exp(-x)/(1+exp(-x))^2,
             'laplace'=1
  )
  return(out)
}












## derivatives of the log-copula density of Archimedean  
## copulas included: first and second order derivative   
## to the copula parameter, first derivative to the     
## argments u_j and mixed derivate to u_j and the copula 
## parameter.                                            
## Copulas included are: Clayton, Gumbel, AMH, Frank and 
## Joe                                                  
#########################################################
#########################################################


### clayton copula ###
######################


derivsClayton=function(theta,U){
  
  # function to calculate the derivatives of log-density of the Clayton copula
  # theta is the copula parameter theta_C, numeric
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### k in the sum
  k = 0:(d-1)
  
  
  ### t(u) and its derivatives
  
  # t(u)
  tu = apply(U^(-theta)-1,1,sum)
  
  # first derivative of t(u) w.r.t. the copula parameter
  dt_tu = apply(-U^(-theta)*log(U),1,sum)
  
  # second derivative of t(u) w.r.t. the copula parameter
  dtt_tu = apply(U^(-theta)*log(U)^2,1,sum)
  
  # derivative of t(u) w.r.t. u_j (in column j)
  du_tu = -theta*U^(-theta-1)
  
  # mixed derivative of t(u) (w.r.t. u_j in column j)
  dut_tu  =U^(-theta-1)*(theta*log(U)-1)
  
  
  
  ### first derivative to the copula parameter
  partialcop = sum(k/(theta*k+1)) - apply(log(U),1,sum) + 1/theta^2*log(1+tu) -  
    (d+1/theta)*dt_tu/(1+tu)
  
  
  
  ### second derivative to the copula parameter
  partialcopcop = -sum(k^2/(theta*k+1)^2) - 2/theta^3*log(1+tu) + 2/theta^2*dt_tu/(1+tu) - 
    (d+1/theta)*(dtt_tu*(1+tu)-dt_tu^2)/(1+tu)^2
  
  
  ### first derivative to u_j (in column j)
  partialu = -(1+theta)/U - (d+1/theta)*du_tu/(1+tu)
  
  
  ### mixed derivative to u_j (in column j) and the copula parameter
  partialcopu = -1/U + 1/theta^2*du_tu/(1+tu) - 
    (d+1/theta)*(dut_tu*(1+tu)-dt_tu*du_tu)/(1+tu)^2
  
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}



### AMH copula ###
##################


derivsAMH=function(theta,U){
  
  # function to calculate the derivatives of log-density of the AMH copula
  # theta is the copula parameter theta_C, numeric
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### h_theta^A(u)
  hu = theta*apply(U/(1-theta*(1-U)),1,prod)
  
  
  ### polylogarithm functions of order -d, -d-1 and -d-2
  LiD0 = copula::polylog(z = hu,s = -d,method = "default")
  LiD1 = copula::polylog(z = hu,s = -d-1,method = "default")
  LiD2 = copula::polylog(z = hu,s = -d-2,method = "default")
  
  
  ### terms that are frequent in appearance
  
  # sum over j of (1-u_j)/(1-theta*(1-u_j))
  fa1 = apply((1-U)/(1-theta*(1-U)),1,sum)
  
  # sum over j of [(1-u_j)/(1-theta*(1-u_j))]^2
  fa2 = apply(((1-U)/(1-theta*(1-U)))^2,1,sum)
  
  
  ### first derivative w.r.t. the copula parameter
  partialcop = -(d+1)/(1-theta) - 1/theta + fa1 + LiD1*(1/theta + fa1)/LiD0
  
  
  
  ### second derivative w.r.t. the copula parameter
  partialcopcop = -(d+1)/(1-theta)^2 + 1/theta^2 + fa2 + 
    ( LiD0*( LiD2*(1/theta + fa1)^2 + LiD1*(-1/theta^2 + fa2) ) -
        (1/theta + fa1)^2*LiD1^2)/LiD0^2
  
  
  ### first derivative w.r.t. u_j (in column j)
  partialu = -theta/(1-theta*(1-U)) - 1/U + LiD1/LiD0*(1-theta)/(U*(1 - theta*(1 - U)))
  
  
  ### mixed derivative w.r.t. u_j (in column j) and the copula parameter
  partialcopu = -1/(1 - theta*(1-U))^2 + 1/LiD0^2*( -LiD0*LiD1/(1 - theta*(1-U))^2 + 
                                                      LiD0*LiD2*(1 - theta)/(U*(1 - theta*(1-U)))*(1/theta + fa1) - 
                                                      LiD1^2*(1/theta + fa1)*(1 - theta)/(U*(1 - theta*(1-U))) )
  
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}



### Frank copula ###
####################


derivsFrank=function(theta,U){
  
  # function to calculate the derivatives of log-density of the Frank copula
  # theta is the copula parameter theta_C, numeric
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### h_theta^F(u)
  hu = (1 - exp(-theta))^(1 - d)*apply(1-exp(-theta*U),1,prod)
  
  
  ### polylogarithm functions of order -d, -d-1 and -d-2
  LiD0 = copula::polylog(z = hu,s = -(d-1),method = "default")
  LiD1 = copula::polylog(z = hu,s = -d,method = "default")
  LiD2 = copula::polylog(z = hu,s = -(d+1),method = "default")
  
  
  ### terms that are frequent in appearance
  
  # theta*exp(-theta*u_j)/(1-exp(-theta*u_j))
  fa1 = theta*exp(-theta*U)/(1-exp(-theta*U))
  
  # sum over j of u_j*exp(-theta*u_j)/(1-exp(-theta*u_j))
  fa2 = apply(U*exp(-theta*U)/(1-exp(-theta*U)),1,sum)
  
  # sum over j of u_j^2*exp(-theta*u_j)/(1-exp(-theta*u_j))^2
  fa3 = apply(U^2*exp(-theta*U)/(1-exp(-theta*U))^2,1,sum)
  
  
  
  ### first derivative w.r.t. the copula parameter
  partialcop = (d-1)/theta + LiD1*( (1-d)*exp(-theta)/(1-exp(-theta)) + fa2)/LiD0 - 
    apply(U,1,sum) - fa2
  
  
  
  ### second derivative w.r.t. the copula parameter
  partialcopcop = -(d-1)/theta^2 + fa3 + 1/LiD0^2*( (LiD2*((1-d)*exp(-theta)/(1-exp(-theta)) + fa2)^2 -
                                                       LiD1*((1-d)*exp(-theta)/(1-exp(-theta))^2 + fa3))*LiD0 - 
                                                      (LiD1*((1-d)*exp(-theta)/(1-exp(-theta)) + fa2))^2 )
  
  
  ### first derivative w.r.t. u_j (in column j)
  partialu = LiD1*fa1/LiD0 - theta - fa1
  
  
  ### mixed derivative w.r.t. u_j (in column j) and the copula parameter
  partialcopu = -1 - exp(-theta*U)/(1-exp(-theta*U))^2*(1-exp(-theta*U)-theta*U) + 
    1/LiD0^2*( (exp(-theta*U)/(1-exp(-theta*U))^2*(1-exp(-theta*U)-theta*U)*LiD1 + 
                  ((1-d)*exp(-theta)/(1-exp(-theta)) + fa2)*LiD2*fa1)*LiD0 - 
                 LiD1^2*fa1*(((1-d)*exp(-theta)/(1-exp(-theta)) + fa2)) )
  
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}


### Joe copula ###
##################

derivsJoe=function(theta,U){
  
  # function to calculate the derivatives of log-density of the Joe copula
  # theta is the copula parameter theta_C, numeric
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### h_theta^J(u) and its derivatives
  
  # h_theta^J(u)
  hu = apply(1-(1-U)^theta,1,prod)
  
  # first derivative to the copula parameter
  dt_hu = -hu*apply((1-U)^theta*log(1-U)/(1-(1-U)^theta),1,sum) 
  
  # second derivative w.r.t. the copula parameter
  dtt_hu = hu*( apply((1-U)^theta*log(1-U)/(1-(1-U)^theta),1,sum)^2 - 
                  apply((1-U)^theta*log(1-U)^2/(1-(1-U)^theta)^2,1,sum) )
  
  # derivative w.r.t. u_j (in the j-th column)
  du_hu = hu*theta*(1-U)^(theta-1)/(1-(1-U)^theta) 
  
  # mixed derivative w.r.t. the copula parameter and u_j (in the j-th column)
  dut_hu = hu*( (1-U)^(theta-1)*(1-(1-U)^theta+theta*log(1-U))/(1-(1-U)^theta)^2 - 
                  theta*(1-U)^(theta-1)/(1-(1-U)^theta)*apply((1-U)^theta*log(1-U)/(1-(1-U)^theta),1,sum) )
  
  
  
  ### P_(d,1/theta)^J function and its derivatives
  
  # k in the sum of p
  k = 0:(d-1)
  
  # a_(d,k)^J(1/theta)
  a = sapply(k+1,Stirling2,n=d)*gamma(k+1-1/theta)/gamma(1-1/theta)
  
  # derivative of a_(d,k)^J(1/theta) w.r.t. the copula parameter
  dt_a = sapply(k+1,Stirling2,n=d)*gamma(k+1-1/theta)*( digamma(k+1-1/theta) - 
                                                          digamma(1-1/theta))/(theta^2*gamma(1-1/theta) )
  
  # second derivative of a_(d,k)^J(1/theta) w.r.t. the copula parameter
  dtt_a = -sapply(k+1,Stirling2,n=d)*gamma(k+1-1/theta)/(theta^4*gamma(1-1/theta))*( -digamma(k+1-1/theta)^2 + 
                                                                                       2*theta*digamma(k+1-1/theta) + 2*digamma(1-1/theta)*digamma(k+1-1/theta) - trigamma(k+1-1/theta) - 
                                                                                       digamma(1-1/theta)^2 - 2*theta*digamma(1-1/theta) + trigamma(1-1/theta) )
  
  # x = h_theta^J/(1-h_theta^J)
  x = hu/(1-hu)
  
  # derivative of x w.r.t. the copula parameter
  dt_x = dt_hu/(1-hu)^2
  
  # second derivative of x w.r.t. the copula parameter
  dtt_x = (dtt_hu*(1-hu) + 2*dt_hu^2)/(1-hu)^3
  
  # derivative of x w.r.t. u_j (in the j-th column)
  du_x = du_hu/(1-hu)^2
  
  # mixed derivative of x w.r.t. the copula parameter and u_j (in the j-th column)
  dut_x = (dut_hu*(1-hu) + 2*du_hu*dt_hu)/(1-hu)^3
  
  # P_(d,1/theta)^J(x)  
  P = apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d),MARGIN = 2,STATS = k,FUN = '^'),2,a,'*'),1,sum)
  
  # first derivative of P w.r.t. the copula parameter
  dt_P = dt_x*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*a)[-1],'*'),1,sum) + 
    apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d),MARGIN = 2,STATS = k,FUN = '^'),2,dt_a,'*'),1,sum)
  
  # second derivative of P w.r.t. the copula parameter
  dtt_P = 1*(d>2)*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-2),MARGIN = 2,STATS = k[-c(1,2)]-2,FUN = '^'),2,(k*(k-1)*a)[-c(1,2)],'*'),1,sum)*dt_x^2 +
    dtt_x*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*a)[-1],'*'),1,sum) + 
    apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d),MARGIN = 2,STATS = k,FUN = '^'),2,dtt_a,'*'),1,sum) + 
    2*dt_x*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*dt_a)[-1],'*'),1,sum)
  
  # derivative w.r.t. u_j of P (in the j-th column)
  du_P = apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*a)[-1],'*'),1,sum)*du_x
  
  # mixed derivative w.r.t. the copula parameter and u_j of P (in the j-th column)
  dut_P = du_x*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*dt_a)[-1],'*'),1,sum) + 
    1*(d>2)*apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-2),MARGIN = 2,STATS = k[-c(1,2)]-2,FUN = '^'),2,(k*(k-1)*a)[-c(1,2)],'*'),1,sum)*dt_x*du_x +
    apply(sweep(sweep(x = matrix(x,nrow=n,ncol=d-1),MARGIN = 2,STATS = k[-1]-1,FUN = '^'),2,(k*a)[-1],'*'),1,sum)*dut_x
  
  
  
  
  ### first derivative w.r.t. the copula parameter
  partialcop = (d-1)/theta + apply(log(1-U),1,sum) - 1/theta^2*log(1-hu) +
    (1-1/theta)*dt_hu/(1-hu) + dt_P/P
  
  ### second derivative w.r.t. the copula parameter 
  partialcopcop = -(d-1)/theta^2 + 2/theta^3*log(1-hu) + 2/theta^2*dt_hu/(1-hu) +
    (1-1/theta)*(dtt_hu*(1-hu) + dt_hu^2)/(1-hu)^2 + 
    (dtt_P*P - dt_P^2)/P^2
  
  ### first derivative w.r.t. u_j (in column j)
  partialu = -(theta-1)/(1-U) + (1-1/theta)*du_hu/(1-hu) + du_P/P
  
  
  ### mixed derivative w.r.t. u_j (in column j) and the copula parameter
  partialcopu = -1/(1-U) + 1/theta^2*du_hu/(1-hu) + (1-1/theta)*(dut_hu*(1-hu)+du_hu*dt_hu)/(1-hu)^2 + 
    (dut_P*P-du_P*dt_P)/P^2
  
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}




### Gumbel copula ###
#####################

derivsGumbel=function(theta,U){
  
  # function to calculate the derivatives of log-density of the Gumbel copula
  # theta is the copula parameter theta_C, numeric
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### t_theta^G(u) and its derivatives
  
  # t_theta^G(u)
  tu = apply((-log(U))^theta,1,sum)
  
  # first derivative to the copula parameter
  dt_tu = apply((-log(U))^theta*log(-log(U)),1,sum)
  
  # second derivative to the copula parameter
  dtt_tu = apply((-log(U))^theta*log(-log(U))^2,1,sum)
  
  # derivative to u_j (in the j-th column)
  du_tu = -theta*(-log(U))^(theta-1)*1/U
  
  # mixed derivative to the copula parameter and u_j (in the j-th column)
  dut_tu = -(-log(U))^(theta-1)/U*(theta*log(-log(U))+1)
  
  
  
  ### P_(d,1/theta)^G function and its derivatives
  
  # k in the sum of p
  k = 1:d
  
  # a_(d,k)^J(1/theta) 
  a = (-1)^(d-k)
  for(j in 1:d){
    a[j]=a[j]*sum((1/theta)^(k[j]:d)*sapply(k[j]:d,Stirling1,n=d)*sapply(k[j]:d,Stirling2,k=k[j]))
  }
  
  # derivative of a_(d,k)^J(1/theta) to the copula parameter
  dt_a = -(-1)^(d-k)
  for(j in 1:d){
    dt_a[j] = dt_a[j]*sum((k[j]:d)*(1/theta)^((k[j]:d)+1)*sapply(k[j]:d,Stirling1,n=d)*sapply(k[j]:d,Stirling2,k=k[j]))
  }
  
  # second derivative of a_(d,k)^J(1/theta) to the copula parameter
  dtt_a = (-1)^(d-k)
  for(j in 1:d){
    dtt_a[j] = dtt_a[j]*sum((k[j]:d)*((k[j]:d)+1)*(1/theta)^((k[j]:d)+2)*sapply(k[j]:d,Stirling1,n=d)*sapply(k[j]:d,Stirling2,k=k[j]))
  }
  
  
  # P_(d,1/theta)^J(x) 
  P = apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,a,'*'),1,sum)
  
  # first derivative of P to the copula parameter
  dt_P = apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,dt_a,'*'),1,sum) +
    dt_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*a,'*'),1,sum) - 
    log(tu)*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,k/theta^2*a,'*'),1,sum)
  
  # second derivative of P to the copula parameter --> niet correct
  dtt_P = apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,dtt_a,'*'),1,sum) +
    2*dt_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*dt_a,'*'),1,sum) -
    2*log(tu)*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,k/theta^2*dt_a,'*'),1,sum)  -
    2*dt_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta^2*a,'*'),1,sum) +
    dt_tu^2*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-2,FUN = '^'),2,k/theta*(k/theta-1)*a,'*'),1,sum) + 
    dtt_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*a,'*'),1,sum) -
    2*dt_tu*log(tu)*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k^2/theta^3*a,'*'),1,sum) +
    log(tu)^2*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,k^2/theta^4*a,'*'),1,sum) + 
    2*log(tu)*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta,FUN = '^'),2,k/theta^3*a,'*'),1,sum)  
  
  
  # derivative to u_j of P (in the j-th column)
  du_P = du_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*a,'*'),1,sum)
  
  
  # mixed derivative to the copula parameter and u_j of P (in the j-th column) --> niet correct
  dut_P = -du_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta^2*a,'*'),1,sum) + 
    du_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*dt_a,'*'),1,sum) +
    du_tu*dt_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-2,FUN = '^'),2,k/theta*(k/theta-1)*a,'*'),1,sum) -
    log(tu)*du_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k^2/theta^3*a,'*'),1,sum) + 
    dut_tu*apply(sweep(sweep(x = matrix(tu,nrow=n,ncol=d),MARGIN = 2,STATS = k/theta-1,FUN = '^'),2,k/theta*a,'*'),1,sum)
  
  
  
  
  ### first derivative to the copula parameter
  partialcop = d/theta + apply(log(-log(U)),1,sum) - tu^(1/theta-1)*(theta*dt_tu - tu*log(tu))/theta^2 -
    d*dt_tu/tu + dt_P/P
  
  ### second derivative to the copula parameter 
  partialcopcop = -d/theta^2 + 2*tu^(1/theta-1)*(theta*dt_tu - tu*log(tu))/theta^3 - 
    tu^(1/theta-1)/theta^2*(theta*dt_tu - tu*log(tu))*((1/theta-1)/tu*dt_tu - log(tu)/theta^2) - 
    tu^(1/theta-1)/theta^2*(theta*dtt_tu - dt_tu*log(tu)) - d*(dtt_tu*tu-dt_tu^2)/tu^2 + 
    (dtt_P*P - dt_P^2)/P^2
  
  ### first derivative to u_j (in column j)
  partialu = -1/theta*tu^(1/theta-1)*du_tu + (theta-1)/(U*log(U)) - d*du_tu/tu - 1/U + du_P/P
  
  
  ### mixed derivative to u_j (in column j) and the copula parameter --> niet correct
  partialcopu = 1/(U*log(U)) - d/tu^2*(dut_tu*tu - dt_tu*du_tu) + (dut_P*P - dt_P*du_P)/P^2 - 
    1/theta^2*( (1/theta-1)*tu^(1/theta-2)*du_tu*(theta*dt_tu - tu*log(tu)) + 
                  tu^(1/theta-1)*(theta*dut_tu-du_tu*(log(tu)+1)) )
  
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}










## derivatives of the log-copula density of Gaussian    
## copula included: first and second order derivative    
## to the copula parameter, first derivative to the     
## argments u_j and mixed derivate to u_j and the copula 
## parameter.                                            
#########################################################
#########################################################

derivsgaussian=function(theta,U){
  
  # function to calculate the derivatives of log-density of the Gaussian copula
  # theta are the copula parameters Sigma, numeric vector (S21,S31,...,S32,...,...) of length d(d-1)/2
  # U is a matrix of size nxd containing the pseudo observations
  # n=number of observations; d=dimension
  
  n = nrow(U)
  d = ncol(U)
  
  
  ### creation of Sigma, the correlation matrix
  Sigma=diag(0.5,nrow=d)
  Sigma[lower.tri(Sigma,diag=F)]=theta
  Sigma=Sigma+t(Sigma)
  
  SigmaInv=solve(Sigma)
  
  ### determinant of Sigma
  dets0out=det(Sigma)
  
  
  ### matrix to hold the determinants of matrices with 1 row and 1 column missing
  dets1out=matrix(0,nrow=d,ncol=d)
  if(d==2){
    for(j in 1:d){
      for(i in 1:j){
        
        dets1out[j,i]=dets1out[i,j]=Sigma[-j,-i]
        
      }
    }
  } else {
    for(j in 1:d){
      for(i in 1:j){
        
        dets1out[j,i]=dets1out[i,j]=det(Sigma[-j,-i])
        
      }
    }
  }
  
  
  ### array to hold the determinants of matrices with 2 rows and 2 columns missing
  # dets2out[j,i,k,l] is the determinant with the j- and k-th row and i- and l-th column missing
  
  if(d==2){
    dets2out=array(1,dim = rep(d,4))
  } else if(d==3){
    dets2out=array(1,dim = rep(d,4))
    for(j in 2:d){
      for(i in 2:d){
        for(k in 1:(j-1)){
          for(l in 1:(i-1)){
            
            dets2out[j,i,k,l]=Sigma[-c(j,k),-c(i,l)]
            dets2out[k,i,j,l]=dets2out[j,i,k,l]
            dets2out[j,l,k,i]=dets2out[j,i,k,l]
            dets2out[k,l,j,i]=dets2out[j,i,k,l]
            
          }
        }
      }
    }
  } else {
    dets2out=array(1,dim = rep(d,4))
    for(j in 2:d){
      for(i in 2:d){
        for(k in 1:(j-1)){
          for(l in 1:(i-1)){
            
            dets2out[j,i,k,l]=det(Sigma[-c(j,k),-c(i,l)])
            dets2out[k,i,j,l]=dets2out[j,i,k,l]
            dets2out[j,l,k,i]=dets2out[j,i,k,l]
            dets2out[k,l,j,i]=dets2out[j,i,k,l]
            
          }
        }
      }
    }
  }
  
  ### array to hold the determinants of matrices with 3 rows and 3 columns missing
  # dets3out[j,i,k,l,r,s] is the determinant with the j-, k-, r-th row and i-,l-,s-th column missing
  if(d==2){
    dets3out=array(0,dim=rep(d,6))
  } else if(d==3){
    dets3out=array(1,dim=rep(d,6))
  } else if(d==4){
    dets3out=array(1,dim=rep(d,6))
    for(j in 3:d){
      for(i in 3:d){
        for(k in 2:(j-1)){
          for(l in 2:(i-1)){
            for(r in 1:(k-1)){
              for(s in 1:(l-1)){
                
                posrow=gtools::permutations(n = 3,r = 3,v = c(j,k,r))
                poscol=gtools::permutations(n = 3,r = 3,v = c(i,l,s))
                pos=matrix(0,nrow=36,ncol=6)
                for(q in 1:6){
                  pos[(1:6)+(q-1)*6,c(1,3,5)]=posrow
                  pos[(1:6)+(q-1)*6,c(2,4,6)]=sweep(x=pos[(1:6)+(q-1)*6,c(2,4,6)],MARGIN = 2,STATS =-poscol[q,])
                }
                dets3out[pos]=Sigma[-c(j,k,r),-c(i,l,s)]
                
              }
            }
          }
        }
      }
    }
  } else {
    dets3out=array(1,dim=rep(d,6))
    for(j in 3:d){
      for(i in 3:d){
        for(k in 2:(j-1)){
          for(l in 2:(i-1)){
            for(r in 1:(k-1)){
              for(s in 1:(l-1)){
                
                posrow=gtools::permutations(n = 3,r = 3,v = c(j,k,r))
                poscol=gtools::permutations(n = 3,r = 3,v = c(i,l,s))
                pos=matrix(0,nrow=36,ncol=6)
                for(q in 1:6){
                  pos[(1:6)+(q-1)*6,c(1,3,5)]=posrow
                  pos[(1:6)+(q-1)*6,c(2,4,6)]=sweep(x=pos[(1:6)+(q-1)*6,c(2,4,6)],MARGIN = 2,STATS =-poscol[q,])
                }
                dets3out[pos]=det(Sigma[-c(j,k,r),-c(i,l,s)])
                
              }
            }
          }
        }
      }
    }
  }
  
  
  ### construction of the elements B_{j,i}^{k,l} and {r,s}^D_{j,i}^{k,l}
  B=array(NA,dim = rep(d,4))
  D=array(NA,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 1:d){
        for(l in 1:d){
          
          # B_{j,i}^{k,l}
          B[j,i,k,l] = (1-1*(j!=k)*1*(i!=l))*1*(i!=j)*(-1)^(k+l+1*(k>i)+1*(l>j))*dets2out[j,i,l,k] +
            (1-1*(j!=l)*1*(i!=k))*1*(i!=j)*(-1)^(k+l+1*(k>j)+1*(l>i))*dets2out[j,i,k,l] +
            1*(j!=k)*1*(j!=l)*1*(i!=k)*1*(i!=l)*((-1)^(k+l+1*(k>i)+1*(l>j))*dets2out[j,i,l,k]
                                                 +(-1)^(k+l+1*(k>j)+1*(l>i))*dets2out[j,i,k,l] )
          
          
          
          # {r,s}^D_{j,i}^{k,l}
          for(r in 1:d){
            for(s in 1:d){
              
              D[j,i,k,l,r,s]= 1*(r!=j)*1*(r!=i)*1*(r!=k)*1*(r!=l)*1*(s!=j)*1*(s!=i)*1*(s!=k)*1*(s!=l)*(
                (-1)^(r+s+1*(s>i)+1*(s>l)+1*(r>j)+1*(r>k))*dets3out[j,i,k,l,r,s] +
                  (-1)^(r+s+1*(r>i)+1*(r>l)+1*(s>j)+1*(s>k))*dets3out[j,i,k,l,s,r] ) +
                (1-1*(r!=i)*1*(r!=l)*1*(s!=j)*1*(s!=k))*1*(s!=i)*1*(s!=l)*1*(r!=j)*1*(r!=k)*(
                  (-1)^(r+s+1*(s>i)+1*(s>l)+1*(r>j)+1*(r>k))*dets3out[j,i,k,l,r,s] ) +
                (1-1*(r!=j)*1*(r!=k)*1*(s!=i)*1*(s!=l))*1*(s!=j)*1*(s!=k)*1*(r!=i)*1*(r!=l)*(
                  (-1)^(r+s+1*(r>i)+1*(r>l)+1*(s>j)+1*(s>k))*dets3out[j,i,k,l,s,r] )
            }
          }
          
        }
      }
    }
  }
  
  
  ### second derivatives w.r.t. Sigma_{k,l} and Sigma_{r,s} of det(Sigma_{-j,-i})
  B2=array(NA,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          for(r in 2:d){
            for(s in 1:(r-1)){
              
              B2[j,i,k,l,r,s] = (1-1*(j!=k)*1*(i!=l))*1*(i!=j)*(-1)^(k+l+1*(k>i)+1*(l>j))*D[j,i,l,k,r,s] +
                (1-1*(j!=l)*1*(i!=k))*1*(i!=j)*(-1)^(k+l+1*(k>j)+1*(l>i))*D[j,i,k,l,r,s] +
                1*(j!=k)*1*(j!=l)*1*(i!=k)*1*(i!=l)*((-1)^(k+l+1*(k>i)+1*(l>j))*D[j,i,l,k,r,s]
                                                     +(-1)^(k+l+1*(k>j)+1*(l>i))*D[j,i,k,l,r,s] )
              
              
              B2[j,i,l,k,r,s]=B2[j,i,k,l,r,s]
              B2[j,i,k,l,s,r]=B2[j,i,k,l,r,s]
              B2[j,i,l,k,s,r]=B2[j,i,k,l,r,s]
              
            }
          }
        }
      }
    }
  }
  
  ### first derivatives w.r.t. Sigma_{k,l} of Sigma^{-1}_{i,j}
  D1SigmaInv=array(0,dim = rep(d,4))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          
          D1SigmaInv[i,j,k,l] = (-1)^(i+j)*( dets0out*B[j,i,k,l] 
                                             - 2*dets1out[j,i]*(-1)^(k+l)*dets1out[k,l] )
          D1SigmaInv[i,j,l,k]=D1SigmaInv[i,j,k,l]
          
        }
      }
    }
  }
  D1SigmaInv=D1SigmaInv/dets0out^2
  
  
  ### second derivatives w.r.t. Sigma_{k,l} and Sigma_{r,s} of Sigma^{-1}_{i,j}
  D2SigmaInv=array(0,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          for(r in 2:d){
            for(s in 1:(r-1)){
              
              D2SigmaInv[j,i,k,l,r,s] = (-1)^(i+j)*(
                B2[j,i,k,l,r,s]*dets0out^2 -
                  2*B[j,i,k,l]*(-1)^(r+s)*dets1out[r,s]*dets0out -
                  2*B[j,i,r,s]*(-1)^(k+l)*dets1out[k,l]*dets0out -
                  2*dets1out[j,i]*(-1)^(k+l)*B[k,l,r,s]*dets0out +
                  8*dets1out[j,i]*(-1)^(k+l+r+s)*dets1out[k,l]*dets1out[r,s]
              )
              
              D2SigmaInv[j,i,l,k,r,s]=D2SigmaInv[j,i,k,l,r,s]
              D2SigmaInv[j,i,k,l,s,r]=D2SigmaInv[j,i,k,l,r,s]
              D2SigmaInv[j,i,l,k,s,r]=D2SigmaInv[j,i,k,l,r,s]
              
            }
          }
        }
      }
    }
  }
  D2SigmaInv=D2SigmaInv/dets0out^3
  
  
  ### terms involving U
  
  # Gaussian quantile function evaluated at U
  qU=qnorm(U)
  
  # derivative of the Gaussian quantile function evaluated at U
  dqU=1/dnorm(qnorm(U))
  
  
  
  ### indices of the element of Sigma
  indices=matrix(NA,nrow=d*(d-1)/2,ncol=2)
  count=1
  for(l in 1:(d-1)){
    for(k in (l+1):d){
      
      indices[count,]=c(k,l)
      
      count=count+1
    }
  }
  
  
  ### first derivative w.r.t. the copula parameters
  # rows for observations
  # columns for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcop = matrix(NA,nrow=n,ncol=d*(d-1)/2)
  
  for(count1 in 1:(d*(d-1)/2)){    
    
    k=indices[count1,1]
    l=indices[count1,2]
    
    # term with the sums over i and j
    temp=matrix(NA,nrow=n,ncol=1)
    for(q in 1:n){
      temp[q]=qU[q,]%*%D1SigmaInv[,,k,l]%*%qU[q,]
    }
    
    # the derivative in the n observations
    partialcop[,count1]=-(-1)^(k+l)*dets1out[k,l]/dets0out - 1/2*temp
    
  }
  
  
  
  ### second derivative w.r.t. the copula parameter
  # first dimension is for observations
  # second and third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcopcop = array(NA,dim = c(n,d*(d-1)/2,d*(d-1)/2))
  
  for(count1 in 1:(d*(d-1)/2)){
    for(count2 in count1:(d*(d-1)/2)){
      
      k=indices[count1,1]
      l=indices[count1,2]
      r=indices[count2,1]
      s=indices[count2,2]
      
      # term involving the sum over i and j
      temp=matrix(NA,nrow=n,ncol=1)
      for(q in 1:n){
        temp[q]=qU[q,]%*%D2SigmaInv[,,k,l,r,s]%*%qU[q,]
      }
      
      # second derivative in the observations with respect to Sigma_{k,l} and Sigma_{r,s}
      partialcopcop[,count1,count2]= -(-1)^(k+l)*B[k,l,r,s]/dets0out +
        2*(-1)^(k+l+r+s)*dets1out[k,l]*dets1out[r,s]/dets0out^2 -
        1/2*temp 
      partialcopcop[,count2,count1]=partialcopcop[,count1,count2]
      
    }
  }
  
  
  ### first derivative w.r.t. u_j (in column j)
  # rows for observations
  # column for u_1,...,u_d
  partialu = matrix(NA,nrow=n,ncol=d)
  for(q in 1:d){
    
    partialu[,q]=-dqU[,q]*(qU%*%SigmaInv[,q]-qU[,q])
    
  }
  
  
  ### mixed derivative w.r.t. u_j (in column j) and the copula parameter
  # first dimension is for the observations
  # second dimension is for u_1,...,u_d
  # third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcopu = array(NA,dim=c(n,d,d*(d-1)/2))
  for(count1 in 1:(d*(d-1)/2)){
    for(q in 1:d){
      
      k=indices[count1,1]
      l=indices[count1,2]
      
      partialcopu[,q,count1]=-dqU[,q]*(qU%*%(D1SigmaInv[,q,k,l]))
      
    }
  }
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}










## derivatives of the log-copula density of Student's t  
## copula included: first and second order derivative    
## to the copula parameter, first derivative to the      
## arguments u_j and mixed derivate to u_j and the      
## copula parameter.                                     
#########################################################
#########################################################

derivst=function(theta,U,nu){
  
  # function to calculate the derivatives of log-density of the Gaussian copula
  # theta are the copula parameters Sigma, numeric vector (S21,S31,...,S32,...,...,nu) of length d(d-1)/2 + 1 
  # U is a matrix of size nxd containing the pseudo observations in the rows
  # n=number of observations; d=dimension; nc=number of copula parameters
  
  n = nrow(U)
  d = ncol(U)
  nc=d*(d-1)/2+1
  
  
  ### creation of Sigma, the correlation matrix
  Sigma=diag(0.5,nrow=d)
  Sigma[lower.tri(Sigma,diag=F)]=theta[-nc]
  Sigma=Sigma+t(Sigma)
  
  SigmaInv=solve(Sigma)
  
  
  
  
  ### calculation of terms found in the derivatives ###
  #####################################################
  
  ### determinant of Sigma
  dets0out=det(Sigma)
  
  
  ### matrix to hold the determinants of matrices with 1 row and 1 column missing
  dets1out=matrix(0,nrow=d,ncol=d)
  if(d==2){
    for(j in 1:d){
      for(i in 1:j){
        
        dets1out[j,i]=dets1out[i,j]=Sigma[-j,-i]
        
      }
    }
  } else {
    for(j in 1:d){
      for(i in 1:j){
        
        dets1out[j,i]=dets1out[i,j]=det(Sigma[-j,-i])
        
      }
    }
  }
  
  
  ### array to hold the determinants of matrices with 2 rows and 2 columns missing
  # dets2out[j,i,k,l] is the determinant with the j- and k-th row and i- and l-th column missing
  
  if(d==2){
    dets2out=array(1,dim = rep(d,4))
  } else if(d==3){
    dets2out=array(1,dim = rep(d,4))
    for(j in 2:d){
      for(i in 2:d){
        for(k in 1:(j-1)){
          for(l in 1:(i-1)){
            
            dets2out[j,i,k,l]=Sigma[-c(j,k),-c(i,l)]
            dets2out[k,i,j,l]=dets2out[j,i,k,l]
            dets2out[j,l,k,i]=dets2out[j,i,k,l]
            dets2out[k,l,j,i]=dets2out[j,i,k,l]
            
          }
        }
      }
    }
  } else {
    dets2out=array(1,dim = rep(d,4))
    for(j in 2:d){
      for(i in 2:d){
        for(k in 1:(j-1)){
          for(l in 1:(i-1)){
            
            dets2out[j,i,k,l]=det(Sigma[-c(j,k),-c(i,l)])
            dets2out[k,i,j,l]=dets2out[j,i,k,l]
            dets2out[j,l,k,i]=dets2out[j,i,k,l]
            dets2out[k,l,j,i]=dets2out[j,i,k,l]
            
          }
        }
      }
    }
  }
  
  ### array to hold the determinants of matrices with 3 rows and 3 columns missing
  # dets3out[j,i,k,l,r,s] is the determinant with the j-, k-, r-th row and i-,l-,s-th column missing
  if(d==2){
    dets3out=array(0,dim=rep(d,6))
  } else if(d==3){
    dets3out=array(1,dim=rep(d,6))
  } else if(d==4){
    dets3out=array(1,dim=rep(d,6))
    for(j in 3:d){
      for(i in 3:d){
        for(k in 2:(j-1)){
          for(l in 2:(i-1)){
            for(r in 1:(k-1)){
              for(s in 1:(l-1)){
                
                posrow=gtools::permutations(n = 3,r = 3,v = c(j,k,r))
                poscol=gtools::permutations(n = 3,r = 3,v = c(i,l,s))
                pos=matrix(0,nrow=36,ncol=6)
                for(q in 1:6){
                  pos[(1:6)+(q-1)*6,c(1,3,5)]=posrow
                  pos[(1:6)+(q-1)*6,c(2,4,6)]=sweep(x=pos[(1:6)+(q-1)*6,c(2,4,6)],MARGIN = 2,STATS =-poscol[q,])
                }
                dets3out[pos]=Sigma[-c(j,k,r),-c(i,l,s)]
                
              }
            }
          }
        }
      }
    }
  } else {
    dets3out=array(1,dim=rep(d,6))
    for(j in 3:d){
      for(i in 3:d){
        for(k in 2:(j-1)){
          for(l in 2:(i-1)){
            for(r in 1:(k-1)){
              for(s in 1:(l-1)){
                
                posrow=gtools::permutations(n = 3,r = 3,v = c(j,k,r))
                poscol=gtools::permutations(n = 3,r = 3,v = c(i,l,s))
                pos=matrix(0,nrow=36,ncol=6)
                for(q in 1:6){
                  pos[(1:6)+(q-1)*6,c(1,3,5)]=posrow
                  pos[(1:6)+(q-1)*6,c(2,4,6)]=sweep(x=pos[(1:6)+(q-1)*6,c(2,4,6)],MARGIN = 2,STATS =-poscol[q,])
                }
                dets3out[pos]=det(Sigma[-c(j,k,r),-c(i,l,s)])
                
              }
            }
          }
        }
      }
    }
  }
  
  
  ### construction of the elements B_{j,i}^{k,l} and {r,s}^D_{j,i}^{k,l}
  B=array(NA,dim = rep(d,4))
  D=array(NA,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 1:d){
        for(l in 1:d){
          
          # B_{j,i}^{k,l}
          B[j,i,k,l] = (1-1*(j!=k)*1*(i!=l))*1*(i!=j)*(-1)^(k+l+1*(k>i)+1*(l>j))*dets2out[j,i,l,k] +
            (1-1*(j!=l)*1*(i!=k))*1*(i!=j)*(-1)^(k+l+1*(k>j)+1*(l>i))*dets2out[j,i,k,l] +
            1*(j!=k)*1*(j!=l)*1*(i!=k)*1*(i!=l)*((-1)^(k+l+1*(k>i)+1*(l>j))*dets2out[j,i,l,k]
                                                 +(-1)^(k+l+1*(k>j)+1*(l>i))*dets2out[j,i,k,l] )
          
          
          
          # {r,s}^D_{j,i}^{k,l}
          for(r in 1:d){
            for(s in 1:d){
              
              D[j,i,k,l,r,s]= 1*(r!=j)*1*(r!=i)*1*(r!=k)*1*(r!=l)*1*(s!=j)*1*(s!=i)*1*(s!=k)*1*(s!=l)*(
                (-1)^(r+s+1*(s>i)+1*(s>l)+1*(r>j)+1*(r>k))*dets3out[j,i,k,l,r,s] +
                  (-1)^(r+s+1*(r>i)+1*(r>l)+1*(s>j)+1*(s>k))*dets3out[j,i,k,l,s,r] ) +
                (1-1*(r!=i)*1*(r!=l)*1*(s!=j)*1*(s!=k))*1*(s!=i)*1*(s!=l)*1*(r!=j)*1*(r!=k)*(
                  (-1)^(r+s+1*(s>i)+1*(s>l)+1*(r>j)+1*(r>k))*dets3out[j,i,k,l,r,s] ) +
                (1-1*(r!=j)*1*(r!=k)*1*(s!=i)*1*(s!=l))*1*(s!=j)*1*(s!=k)*1*(r!=i)*1*(r!=l)*(
                  (-1)^(r+s+1*(r>i)+1*(r>l)+1*(s>j)+1*(s>k))*dets3out[j,i,k,l,s,r] )
            }
          }
          
        }
      }
    }
  }
  
  
  ### second derivatives w.r.t. Sigma_{k,l} and Sigma_{r,s} of det(Sigma_{-j,-i})
  B2=array(NA,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          for(r in 2:d){
            for(s in 1:(r-1)){
              
              B2[j,i,k,l,r,s] = (1-1*(j!=k)*1*(i!=l))*1*(i!=j)*(-1)^(k+l+1*(k>i)+1*(l>j))*D[j,i,l,k,r,s] +
                (1-1*(j!=l)*1*(i!=k))*1*(i!=j)*(-1)^(k+l+1*(k>j)+1*(l>i))*D[j,i,k,l,r,s] +
                1*(j!=k)*1*(j!=l)*1*(i!=k)*1*(i!=l)*((-1)^(k+l+1*(k>i)+1*(l>j))*D[j,i,l,k,r,s]
                                                     +(-1)^(k+l+1*(k>j)+1*(l>i))*D[j,i,k,l,r,s] )
              
              
              B2[j,i,l,k,r,s]=B2[j,i,k,l,r,s]
              B2[j,i,k,l,s,r]=B2[j,i,k,l,r,s]
              B2[j,i,l,k,s,r]=B2[j,i,k,l,r,s]
              
            }
          }
        }
      }
    }
  }
  
  ### first derivatives w.r.t. Sigma_{k,l} of Sigma^{-1}_{i,j}
  D1SigmaInv=array(0,dim = rep(d,4))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          
          D1SigmaInv[i,j,k,l] = (-1)^(i+j)*( dets0out*B[j,i,k,l] 
                                             - 2*dets1out[j,i]*(-1)^(k+l)*dets1out[k,l] )
          D1SigmaInv[i,j,l,k]=D1SigmaInv[i,j,k,l]
          
        }
      }
    }
  }
  D1SigmaInv=D1SigmaInv/dets0out^2
  
  
  ### second derivatives w.r.t. Sigma_{k,l} and Sigma_{r,s} of Sigma^{-1}_{i,j}
  D2SigmaInv=array(0,dim = rep(d,6))
  for(j in 1:d){
    for(i in 1:d){
      for(k in 2:d){
        for(l in 1:(k-1)){
          for(r in 2:d){
            for(s in 1:(r-1)){
              
              D2SigmaInv[j,i,k,l,r,s] = (-1)^(i+j)*(
                B2[j,i,k,l,r,s]*dets0out^2 -
                  2*B[j,i,k,l]*(-1)^(r+s)*dets1out[r,s]*dets0out -
                  2*B[j,i,r,s]*(-1)^(k+l)*dets1out[k,l]*dets0out -
                  2*dets1out[j,i]*(-1)^(k+l)*B[k,l,r,s]*dets0out +
                  8*dets1out[j,i]*(-1)^(k+l+r+s)*dets1out[k,l]*dets1out[r,s]
              )
              
              D2SigmaInv[j,i,l,k,r,s]=D2SigmaInv[j,i,k,l,r,s]
              D2SigmaInv[j,i,k,l,s,r]=D2SigmaInv[j,i,k,l,r,s]
              D2SigmaInv[j,i,l,k,s,r]=D2SigmaInv[j,i,k,l,r,s]
              
            }
          }
        }
      }
    }
  }
  D2SigmaInv=D2SigmaInv/dets0out^3
  
  
  
  # Student's t quantile function evaluated at U
  S=qt(p = U,df = nu)
  
  
  # Student's t density evaluated at S
  s=dt(x = S,df = nu)
  
  # derivative of S w.r.t. the d.f.
  derivqt=function(u,nu){ # numerical differentiation of first order of the quantile function in a single element u
    funcD1=function(x,u){ # function which accepts as first argument the d.f. and returns the quantile function
      return(qt(p = u,df = x))
    }
    return(pracma::fderiv(f = funcD1,x = nu,n = 1,method = "central",u=u))
  }
  dnuS=apply(U,c(1,2),derivqt,nu=nu)
  
  
  # second order derivative of S w.r.t. the d.f.
  deriv2qt=function(u,nu){ # numerical differentiation of second order of the quantile function in a single element u
    funcD2=function(x,u){ # function which accepts as first argument the d.f. and returns the quantile function
      return(qt(p = u,df = x))
    }
    return(pracma::fderiv(f = funcD2,x = nu,n = 2,method = "central",u=u))
  }
  d2nuS=apply(U,c(1,2),deriv2qt,nu=nu)
  
  
  # derivative of S w.r.t. its argument 
  duS=1/s
  
  # derivative of S w.r.t. its argument and the d.f.
  derivduS=function(u,nu){
    funcduS=function(nu,u){ 
      return(1/dt(qt(p = u,df = nu),df = nu))
    }
    return(pracma::fderiv(f = funcduS,x = nu,n = 1,method = "central",u=u))
  }
  dunuS=apply(U,c(1,2),derivduS,nu=nu)
  
  
  ### indices of the element of Sigma
  indices=matrix(NA,nrow=d*(d-1)/2,ncol=2)
  count=1
  for(l in 1:(d-1)){
    for(k in (l+1):d){
      
      indices[count,]=c(k,l)
      
      count=count+1
    }
  }
  
  N=1+diag(S%*%SigmaInv%*%t(S))/nu
  
  
  
  ### derivatives of the log-copula density ###
  #############################################
  
  ### first derivative to elements of Sigma
  # rows for observations
  # columns for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1} (and one for nu, see below)
  partialcop = matrix(NA,nrow=n,ncol=nc)
  
  for(count1 in 1:(d*(d-1)/2)){    
    
    k=indices[count1,1]
    l=indices[count1,2]
    
    # term with the sums over i and j
    temp=matrix(NA,nrow=n,ncol=1)
    for(q in 1:n){
      temp[q]=S[q,]%*%D1SigmaInv[,,k,l]%*%S[q,]
    }
    
    # the derivative in the n observations
    partialcop[,count1]=-(-1)^(k+l)*dets1out[k,l]/dets0out - (nu+d)/(2*nu)*temp/N
    
  }
  
  
  
  ### first derivative to nu
  partialcop[,nc]=1/2*digamma((nu+d)/2) + (d-1)/2*digamma(nu/2) - d/2*digamma((nu+1)/2) - 1/2*log(N) - 
    (nu+d)/(2*nu^2)*(2*nu*diag(S%*%SigmaInv%*%t(dnuS)) - diag(S%*%SigmaInv%*%t(S)))/N +
    1/2*apply(log(1+S^2/nu),1,sum) + (nu+1)/(2*nu^2)*apply((2*nu*S*dnuS-S^2)/(1+S^2/nu),1,sum)
  
  
  ### second derivative to the copula parameter
  # first dimension is for observations
  # second and third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1} and the d.f.
  partialcopcop = array(NA,dim = c(n,nc,nc))
  
  for(count1 in 1:(d*(d-1)/2)){
    for(count2 in count1:(d*(d-1)/2)){
      
      k=indices[count1,1]
      l=indices[count1,2]
      r=indices[count2,1]
      s=indices[count2,2]
      
      # term involving the sum over i and j
      temp1=matrix(NA,nrow=n,ncol=1)
      for(q in 1:n){
        temp1[q]=S[q,]%*%D2SigmaInv[,,k,l,r,s]%*%S[q,]
      }
      
      temp2=matrix(NA,nrow=n,ncol=1)
      for(q in 1:n){
        temp2[q]=S[q,]%*%D1SigmaInv[,,k,l]%*%S[q,]*S[q,]%*%D1SigmaInv[,,r,s]%*%S[q,]
      }
      
      # second derivative in the observations with respect to Sigma_{k,l} and Sigma_{r,s}
      partialcopcop[,count1,count2]= -(-1)^(k+l)*B[k,l,r,s]/dets0out +
        2*(-1)^(k+l+r+s)*dets1out[k,l]*dets1out[r,s]/dets0out^2 -
        (nu+d)/(2*nu)*temp1/N + (nu+d)/(2*nu^2)*temp2/N^2 
      partialcopcop[,count2,count1]=partialcopcop[,count1,count2]
      
    }
  }
  
  
  ### second derivative to nu
  partialcopcop[,nc,nc]=1/4*trigamma((nu+d)/2) + (d-1)/4*trigamma(nu/2) - d/4*trigamma((nu+1)/2) - 1/(2*nu^2)*( 2*nu*diag(S%*%SigmaInv%*%t(dnuS)) - diag(S%*%SigmaInv%*%t(S)) )/N + 
    (nu+2*d)/(2*nu^3)*( 2*nu*diag(S%*%SigmaInv%*%t(dnuS)) - diag(S%*%SigmaInv%*%t(S)) )/N - 
    (nu+d)/nu*(diag(dnuS%*%SigmaInv%*%t(dnuS)) + diag(S%*%SigmaInv%*%t(d2nuS)))/N +
    (nu+d)/(2*nu^4)*(2*nu*diag(S%*%SigmaInv%*%t(dnuS)) - diag(S%*%SigmaInv%*%t(S)))^2/N^2 -
    1/(nu^3)*apply((2*nu*S*dnuS-S^2)/(1+S^2/nu),1,sum) + 
    (nu+1)/nu*apply((dnuS^2+S*d2nuS)/(1+S^2/nu),1,sum) -
    (nu+1)/(2*nu^4)*apply(((2*nu*S*dnuS-S^2)/(1+S^2/nu))^2,1,sum)
  
  
  ### mixed second order derivatives w.r.t Sigma_{k,l} and the d.f.
  for(count1 in 1:(d*(d-1)/2)){    
    
    k=indices[count1,1]
    l=indices[count1,2]
    
    # term with the sums over i and j
    temp1=matrix(NA,nrow=n,ncol=1)
    temp2=matrix(NA,nrow=n,ncol=1)
    for(q in 1:n){
      temp1[q]=S[q,]%*%D1SigmaInv[,,k,l]%*%S[q,]
      temp2[q]=S[q,]%*%D1SigmaInv[,,k,l]%*%dnuS[q,]
    }
    
    # the derivative in the n observations
    partialcopcop[,count1,nc]=-1/(2*nu)*temp1/N - (nu+d)/(2*nu^2)*(2*nu*temp2-temp1)/N +
      (nu+d)/(2*nu^3)*temp1*(2*nu*diag(S%*%SigmaInv%*%t(dnuS))-diag(S%*%SigmaInv%*%t(S)))/N^2
    partialcopcop[,nc,count1]=partialcopcop[,count1,nc]
    
  }
  
  ### first derivative to u_j (in column j)
  # rows for observations
  # column for u_1,...,u_d
  partialu = -(nu+d)/nu*(S%*%SigmaInv*duS)/N + (nu+1)/nu*S*duS/(1+S^2/nu)
  
  
  
  ### mixed derivative to u_j (in column j) and the copula parameter
  # first dimension is for the observations
  # second dimension is for u_1,...,u_d
  # third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1} and nu
  partialcopu = array(NA,dim=c(n,d,nc))
  
  for(count1 in 1:(d*(d-1)/2)){
    
    k=indices[count1,1]
    l=indices[count1,2]
    
    # term with the sums over i and j
    temp=matrix(NA,nrow=n,ncol=1)
    for(q in 1:n){
      temp[q]=S[q,]%*%D1SigmaInv[,,k,l]%*%S[q,]
    }
    
    partialcopu[,,count1]=-(nu+d)/nu*S%*%D1SigmaInv[,,k,l]*duS/N +(nu+d)/nu^2*sweep(S%*%SigmaInv*duS,1,temp,"*")/N^2
    
  }
  
  partialcopu[,,nc]=d/nu^2*S%*%SigmaInv*duS/N - (nu+d)/nu*(dnuS%*%SigmaInv*duS + S%*%SigmaInv*dunuS)/N + 
    (nu+d)/nu^3*(S%*%SigmaInv*duS)*(2*nu*diag(S%*%SigmaInv%*%t(dnuS))-diag(S%*%SigmaInv%*%t(S)))/N^2 - 
    1/nu^2*S*duS/(1+S^2/nu) + (nu+1)/nu*(dnuS*duS+S*dunuS)/(1+S^2/nu) -
    (nu+1)/nu^3*(2*nu*S*dnuS-S^2)*S*duS/(1+S^2/nu)^2
  
  
  return(list("partialcop"=partialcop,
              "partialcopcop"=partialcopcop,
              "partialu"=partialu,
              "partialcopu"=partialcopu)
  )
  
}










## function to simulate from a given copula with     
## specified quantile-based margins                  
#####################################################

samplecopula=function(n,alpha,mu,phi,df=NULL,margfuncs,cop,seed=NULL){
  # a function to simulate n data points from a specified copula structure cop with
  # margins margfuncs. The margins are assumed to be from the quantile based family  
  # with skewing parameters alpha, location parameters mu, scale parameters phi and
  # possibly df degrees of freedom when a student-t distribution is involved
  
  
  # n: numeric
  # alpha: numeric vector of length d with elements in [0,1]
  # mu: numeric vector of length d
  # phi: numeric vector of length d with positive elements
  # df: numeric vector with elements larger than 2
  # margfuncs: character vector of length d. Options for elements are "normal", "laplace", "logistic" or "t"
  # cop: object of class "Copula" with dimensions d
  
  
  # dimensionality of the problem
  d=cop@dimension
  if(length(alpha)!=d | length(mu)!=d | length(phi)!=d | length(margfuncs)!=d){
    
    stop('incorrect length of parameter vector')
    
  }
  
  # check if margfuncss are ok
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(margfuncs,possibilities)
  if(sum(test)>0){
    
    stop("incorrect margin function detected")
    
  }
  
  # check if correct amount of different df are given
  if(is.element("t",margfuncs)){
    df=na.omit(df)
    if(is.null(df) | length(df)!=sum(margfuncs=="t")){
      
      stop("no/incorrect number of degrees of freedom for student-t distribution provided")
      
    } else { 
      
      # check if df are all larger than 2
      if(sum(df<2)>0){
        stop('df out of bounds')
      }
      
      # create helpvector with degrees of freedom
      index=which(margfuncs=="t") 
      tpars=as.matrix(rep(NA,d),nrow=1)
      tpars[index]=df
      df=tpars
      
    }
  }
  
  # check elements of phi and alpha
  if(sum(alpha<0|alpha>1)>0){
    
    stop('alpha out of bounds')
    
  }
  if(sum(phi<=0)>0){
    
    stop('phi out of bounds')
    
  }
  
  # set seed if supplied
  if(!is.null(seed)){
    
    set.seed(seed)
    
  }
  
  # generate a uniform sample from the copula structure
  U=rCopula(n=n,copula=cop)
  
  # transform U to data scale using the inverse of the marginals
  X=matrix(NA,nrow=n,ncol=d)
  for(i in 1:d){
    if(margfuncs[i]=="normal"){
      X[,i]=qAND(beta = U[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else if(margfuncs[i]=="logistic"){
      X[,i]=qALoD(beta = U[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else if(margfuncs[i]=="laplace"){
      X[,i]=qALaD(beta = U[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    } else {
      X[,i]=qATD(beta = U[,i],mu=mu[i],phi=phi[i],alpha=alpha[i],nu=df[i])
    }
  }
  
  # returns the uniform numbers from the copula U and the generated sample X
  return(list("U"=U,"X"=X))
}











## functions for fitting the four univariate quantile based distributions 
##########################################################################

# normal
#########

# log-likelihood function in the correct form
dQBN=function(pars,x){
  # pars contains in order: alpha, mu and phi
  # x is the data
  
  # extracting parameters
  alpha=pars[1]
  mu=pars[2]
  phi=pars[3]
  
  # calculating likelihood in each datapoint
  f=1*(x<=mu)*(alpha*(1-alpha)*sqrt(2/(pi*phi^2))*exp( -1/2*(1-alpha)^2*((mu-x)/phi)^2))+
    1*(x>mu)*(alpha*(1-alpha)*sqrt(2/(pi*phi^2))*exp( -1/2*alpha^2*((x-mu)/phi)^2)) 
  
  # returning minus log-likelihood
  return(-sum(log(f)))
}

# function for fitting
fitAND=function(data,start=NULL,nstart=10,seed=NULL){
  # function for fitting a quantile-based normal distribution to the data
  # using maximum likelihood
  
  # data is a numeric vector containing the data
  # start is an optional set of starting values for the optimizer
  # nstart is the number of different random starting values for the parameter fits
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  lowerbounds=c(0,-Inf,0)
  upperbounds=c(1,Inf,Inf)
  
  
  if(is.null(start)){
    
    # set seed if supplied
    if(!is.null(seed)){
      
      set.seed(seed)
      
    } else {
      
      seed=sample(1:10^8,1)
      set.seed(seed)
      
    }
    
    # generate starting values
    startalpha=runif(nstart)
    startmu=runif(nstart,min=min(data),max=max(data))
    startphi=runif(nstart,min=0,max=sd(data))
    
    # combine starting values in matrix
    x0=cbind(startalpha,startmu,startphi)
    
    # holding vectors for parameter estimates and log-likelihood
    parsfit=matrix(NA,nrow=nstart,ncol=3)
    loglfit=rep(NA,nstart)
    
    # main loop for parameter estimation using the bobyqa function from the nloptr package
    for(i in 1:nstart){
      
      # in case of error for certain starting values, they are suppressed
      try({
        # set seed for consistency of results
        seed=seed+i
        
        # minimization of minus the log likelihood
        output=bobyqa(x0=x0[i,],fn=dQBN,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
        i=i+1
        
        # optimal parameters and minus log-likelihood
        parsfit[i,]=output$par
        loglfit[i]=output$value
      },silent=T)
    }
    
    # returning best fit
    indmin=which.min(loglfit)
    bestpars=parsfit[indmin,]
    
    return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
    
  } else {
    
    # set seed for consistency of results
    seed=seed+i
    
    # minimization of minus log-likelihood
    output=bobyqa(x0=start,fn=dQBN,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=5000,xtol_rel=10^-5))
    
    # returning parameters and log-likelihood
    return(list("alpha"=output$par[1],"mu"=output$par[2],"phi"=output$par[3],"LogLikelihood"=-output$value))
    
  }
}

# Laplace 
##########

# log-likelihood function in the correct form
dQBLa=function(pars,x){
  # pars contains in order: alpha, mu and phi
  # x is the data
  
  # extracting parameters
  alpha=pars[1]
  mu=pars[2]
  phi=pars[3]
  
  # calculating likelihood in each datapoint
  f=1*(x<=mu)*(alpha*(1-alpha)/phi*exp( -(1-alpha)*(mu-x)/phi))+
    1*(x>mu)*(alpha*(1-alpha)/phi*exp( -alpha*(x-mu)/phi))  
  
  # returning minus log-likelihood
  return(-sum(log(f)))
}

# function for fitting
fitALaD=function(data,start=NULL,nstart=10,seed=NULL){
  # function for fitting a quantile-based normal distribution to the data
  # using maximum likelihood
  
  # data is a numeric vector containing the data
  # start is an optional set of starting values for the optimizer
  # nstart is the number of different random starting values for the parameter fits
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  lowerbounds=c(0,-Inf,0)
  upperbounds=c(1,Inf,Inf)
  
  
  if(is.null(start)){
    
    # set seed if supplied
    if(!is.null(seed)){
      
      set.seed(seed)
      
    } else {
      
      seed=sample(1:10^8,1)
      set.seed(seed)
      
    }
    
    # generate starting values
    startalpha=runif(nstart)
    startmu=runif(nstart,min=min(data),max=max(data))
    startphi=runif(nstart,min=0,max=sd(data))
    
    # combine starting values in matrix
    x0=cbind(startalpha,startmu,startphi)
    
    # holding vectors for parameter estimates and log-likelihood
    parsfit=matrix(NA,nrow=nstart,ncol=3)
    loglfit=rep(NA,nstart)
    
    # main loop for parameter estimation using the bobyqa function from the nloptr package
    for(i in 1:nstart){
      
      # in case of error for certain starting values, they are suppressed
      try({
        # set seed for consistency of results
        seed=seed+i
        
        # minimization of minus the log likelihood
        output=bobyqa(x0=x0[i,],fn=dQBLa,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
        
        # optimal parameters and minus log-likelihood
        parsfit[i,]=output$par
        loglfit[i]=output$value
      },silent=T)
    }
    
    # returning best fit
    indmin=which.min(loglfit)
    bestpars=parsfit[indmin,]
    
    return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
    
  } else {
    
    # set seed for consistency of results
    seed=seed+i
    
    # minimization of minus log-likelihood
    output=bobyqa(x0=start,fn=dQBLa,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
    
    # returning parameters and log-likelihood
    return(list("alpha"=output$par[1],"mu"=output$par[2],"phi"=output$par[3],"LogLikelihood"=-output$value))
    
  }
}

# Logistic
###########

# log-likelihood function in the correct form
dQBLo=function(pars,x){
  # pars contains in order: alpha, mu and phi
  # x is the data
  
  # extracting parameters
  alpha=pars[1]
  mu=pars[2]
  phi=pars[3]
  
  # calculating likelihood in each datapoint
  f=1*(x<=mu)*(2*alpha*(1-alpha)/phi*exp(-(1-alpha)*(mu-x)/phi)/(1+exp(-(1-alpha)*(mu-x)/phi))^2)+
    1*(x>mu)*(2*alpha*(1-alpha)/phi*exp(-alpha*(x-mu)/phi)/(1+exp(-alpha*(x-mu)/phi))^2)  
  
  # returning minus log-likelihood
  return(-sum(log(f)))
}

# function for fitting
fitALoD=function(data,start=NULL,nstart=10,seed=NULL){
  # function for fitting a quantile-based normal distribution to the data
  # using maximum likelihood
  
  # data is a numeric vector containing the data
  # start is an optional set of starting values for the optimizer
  # nstart is the number of different random starting values for the parameter fits
  
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi)
  lowerbounds=c(0,-Inf,0)
  upperbounds=c(1,Inf,Inf)
  
  
  if(is.null(start)){
    
    # set seed if supplied
    if(!is.null(seed)){
      
      set.seed(seed)
      
    } else {
      
      seed=sample(1:10^8,1)
      set.seed(seed)
      
    }
    
    # generate starting values
    startalpha=runif(nstart)
    startmu=runif(nstart,min=min(data),max=max(data))
    startphi=runif(nstart,min=0,max=sd(data))
    
    # combine starting values in matrix
    x0=cbind(startalpha,startmu,startphi)
    
    # holding vectors for parameter estimates and log-likelihood
    parsfit=matrix(NA,nrow=nstart,ncol=3)
    loglfit=rep(NA,nstart)
    
    # main loop for parameter estimation using the bobyqa function from the nloptr package
    for(i in 1:nstart){
      
      # in case of error for certain starting values, they are suppressed
      try({
        # set seed for consistency of results
        seed=seed+i
        
        # minimization of minus the log likelihood
        output=bobyqa(x0=x0[i,],fn=dQBLo,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
        
        # optimal parameters and minus log-likelihood
        parsfit[i,]=output$par
        loglfit[i]=output$value
      },silent=T)
    }
    
    # returning best fit
    indmin=which.min(loglfit)
    bestpars=parsfit[indmin,]
    
    return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"LogLikelihood"=-loglfit[indmin]))
    
  } else {
    
    # set seed for consistency of results
    seed=seed+i
    
    # minimization of minus log-likelihood
    output=bobyqa(x0=start,fn=dQBLo,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
    
    # returning parameters and log-likelihood
    return(list("alpha"=output$par[1],"mu"=output$par[2],"phi"=output$par[3],"LogLikelihood"=-output$value))
    
  }
}

# Student's t 
##############

# log-likelihood function in the correct form
dQBT=function(pars,x){
  # pars contains in order: alpha, mu, phi and degrees of freedom
  # x is the data
  
  # extracting parameters
  alpha=pars[1]
  mu=pars[2]
  phi=pars[3]
  nu=pars[4]
  
  # calculating likelihood in each datapoint
  f=1*(x<=mu)*(2*alpha*(1-alpha)/(phi*sqrt(nu)*beta(1/2,nu/2))*(1+(1-alpha)^2/nu*((mu-x)/phi)^2)^(-(nu+1)/2))+
    1*(x>mu)*(2*alpha*(1-alpha)/(phi*sqrt(nu)*beta(1/2,nu/2))*(1+alpha^2/nu*((x-mu)/phi)^2)^(-(nu+1)/2))  
  
  # returning minus log-likelihood
  return(-sum(log(f)))
}

# function for fitting
fitATD=function(data,start=NULL,nstart=10,seed=NULL){
  # function for fitting a quantile-based normal distribution to the data
  # using maximum likelihood
  
  # data is a numeric vector containing the data
  # start is an optional set of starting values for the optimizer
  # nstart is the number of different random starting values for the parameter fits
  
  # upper and lower bounds for parameters in optimization (alpha,mu,phi,df)
  lowerbounds=c(0,-Inf,0,2)
  upperbounds=c(1,Inf,Inf,2000)
  
  
  if(is.null(start)){
    
    # set seed if supplied
    if(!is.null(seed)){
      
      set.seed(seed)
      
    } else {
      
      seed=sample(1:10^8,1)
      set.seed(seed)
      
    }
    
    # generate starting values
    startalpha=runif(nstart)
    startmu=runif(nstart,min=min(data),max=max(data))
    startphi=runif(nstart,min=0,max=sd(data))
    startnu=runif(nstart,min=2,max=200)
    
    # combine starting values in matrix
    x0=cbind(startalpha,startmu,startphi,startnu)
    
    # holding vectors for parameter estimates and log-likelihood
    parsfit=matrix(NA,nrow=nstart,ncol=4)
    loglfit=rep(NA,nstart)
    
    # main loop for parameter estimation using the bobyqa function from the nloptr package
    for(i in 1:nstart){
      
      # in case of error for certain starting values, they are suppressed
      try({
        # set seed for consistency of results
        seed=seed+i
        
        # minimization of minus the log likelihood
        output=bobyqa(x0=x0[i,],fn=dQBT,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
        
        # optimal parameters and minus log-likelihood
        parsfit[i,]=output$par
        loglfit[i]=output$value
      },silent=T)
    }
    
    # returning best fit
    indmin=which.min(loglfit)
    bestpars=parsfit[indmin,]
    
    return(list("alpha"=bestpars[1],"mu"=bestpars[2],"phi"=bestpars[3],"nu"=bestpars[4],"LogLikelihood"=-loglfit[indmin]))
    
  } else {
    
    # set seed for consistency of results
    seed=seed+i
    
    # minimization of minus log-likelihood
    output=bobyqa(x0=start,fn=dQBT,lower=lowerbounds,upper=upperbounds,nl.info=F,x=data,control = list(maxeval=50000,xtol_rel=10^-5))
    
    # returning parameters and log-likelihood
    return(list("alpha"=output$par[1],"mu"=output$par[2],"phi"=output$par[3],"nu"=output$par[4],"LogLikelihood"=-output$value))
    
  }
}












## function to fit a given copula structure        
## with specified quantile-based margins using IFM 
###################################################

fitcopula=function(data,margfuncs,cop,nstart=10,seed=NULL,method="L-BFGS-B"){
  # function to fit a copula object to a dataset with quantile based margins
  # the used technique is the "inference functions for margins" (IFM) technique 
  # as this is similar to full MLE in terms of performance but faster and 
  # easier to implement. It is assumed that the margins are known and correct, 
  # otherwise theoretical results do not hold anymore. First the margins are fit
  # and these fitted parameters are used in the fitting of the copula-object
  
  # data: numeric nxd matrix with 1 data point per row
  # margfuncs: character vector of length d. options are "normal", "laplace","t", "logistic"
  # cop: object of class "Copula" with dimensions d
  # method: optimization method for "optim" in copula fitting
  
  
  # dimensions of the data
  n=length(data[,1])
  d=length(data[1,])
  
  
  # check if dimensions are ok
  if(length(margfuncs)!=d){
    
    stop('incorrect number of margins given')
    
  }
  if(cop@dimension!=d){
    
    stop('incorrect copula dimension')
    
  }
  
  
  # check if margins are of correct type
  possibilities=c("t","normal","logistic","laplace")
  test=!is.element(margfuncs,possibilities)
  if(sum(test)>0){
    
    stop("incorrect margin function detected")
    
  }
  
  
  ### step 1: fitting the margins
  ###############################
  
  # required files
  source(file="~/PhD/code/copulas/univariatefit.R")
  
  # matrix to store uniform transformed data
  U=matrix(NA,nrow=n,ncol=d)
  
  # vectors to store parameters
  alpha=rep(NA,d)
  mu=rep(NA,d)
  phi=rep(NA,d)
  df=rep(NA,d)
  
  # vector to store the log-likelihood of the fitted margins
  logl=rep(NA,d)
  
  # fitting the margins
  for(i in 1:d){
    
    if(margfuncs[i]=="normal"){
      
      fit=fitAND(data=data[,i],nstart=nstart,seed=seed)
      
      U[,i]=pAND(q=data[,i],mu=fit$mu,phi=fit$phi,alpha=fit$alpha)
      
      alpha[i]=fit$alpha
      mu[i]=fit$mu
      phi[i]=fit$phi
      
      logl[i]=fit$LogLikelihood
      
    } else if(margfuncs[i]=="logistic"){
      
      fit=fitALoD(data=data[,i],nstart=nstart,seed=seed)
      
      U[,i]=pALoD(q=data[,i],mu=fit$mu,phi=fit$phi,alpha=fit$alpha)
      
      alpha[i]=fit$alpha
      mu[i]=fit$mu
      phi[i]=fit$phi
      
      logl[i]=fit$LogLikelihood
      
    } else if(margfuncs[i]=="laplace"){
      
      fit=fitALaD(data=data[,i],nstart=nstart,seed=seed)
      
      U[,i]=pALaD(q=data[,i],mu=fit$mu,phi=fit$phi,alpha=fit$alpha)
      
      alpha[i]=fit$alpha
      mu[i]=fit$mu
      phi[i]=fit$phi 
      
      logl[i]=fit$LogLikelihood
      
    } else {
      
      fit=fitATD(data=data[,i],nstart=nstart,seed=seed)
      
      U[,i]=pATD(q=data[,i],mu=fit$mu,phi=fit$phi,alpha=fit$alpha,nu = fit$nu)
      
      alpha[i]=fit$alpha
      mu[i]=fit$mu
      phi[i]=fit$phi 
      df[i]=fit$nu
      
      logl[i]=fit$LogLikelihood
      
    }
    
  }
  
  ### step 2: fitting the copula to the uniform transformed data
  ##############################################################
  
  # fittin the copula object cop
  copulafit=fitCopula(copula = cop,data = U,method="ml",optim.method=method) # for frank replace to "ml", otherwise error
  
  # parameter estimates for the copula parameter and log-likelihood
  coppar=copulafit@estimate
  coplogl=copulafit@loglik
  
  
  
  
  # return uniformly transformed data, margin parameters, copula parameters and total log-likelihood
  loglikelihood=sum(logl,coplogl)
  return(list("U"=U,
              "alpha"=alpha,
              "mu"=mu,
              "phi"=phi,
              "df"=df,
              "copulapar"=coppar,
              "logl"=loglikelihood))
  
}













## Function to calculate the asymptotic variance-covariance matrix of the univariate  
## QBA-normal, QBA-Laplace, QBA-Student's t and QBA-logistic distribution             
#####################################################################################

asymptoticvariance=function(alpha,mu,phi,nu,margfunc){
  # alpha:    the skewing parameter in (0,1)
  # mu:       the location parameter in (-inf,inf)
  # phi:      the scaling parameter in (0,inf)
  # nu:       is the degrees of freedom for the QBA-Student's t in (0,inf)
  # margfunc: is a string containing either "normal", "laplace", logistic" or "t" representing
  #           the appropriate member of the QBA-family
  
  # returns the variance-covariance matrix of the parameters concerning the specified margin 
  # (either a 3x3 matrix or a 4x4 matrix in case of margfunc=="t")
  
  
  # order of output QBAsyDist is: mu; phi; alpha; nu
  
  if(margfunc=="t"){
    
    FI=matrix( c((alpha*(1-alpha)*(nu+1))/(phi^2*(nu+3)), 0, -(8*gamma(nu/2+3/2))/(phi*sqrt(nu*pi)*(nu+3)*gamma(nu/2)), 0,
                 0, (2*nu)/(phi^2*(nu+3)), -(2*(1-2*alpha)*nu)/(phi*alpha*(1-alpha)*(nu+3)), -2/(phi*(nu+1)*(nu+3)),
                 -(8*gamma(nu/2+3/2))/(phi*sqrt(nu*pi)*(nu+3)*gamma(nu/2)), -(2*(1-2*alpha)*nu)/(phi*alpha*(1-alpha)*(nu+3)), (5*alpha^2*nu-3*alpha^2-5*alpha*nu+3*alpha+2*nu)/(alpha^2*(1-alpha)^2*(nu+3)), (2*(1-2*alpha))/(alpha*(1-alpha)*(nu+1)*(nu+3)),
                 0, -2/(phi*(nu+1)*(nu+3)), (2*(1-2*alpha))/(alpha*(1-alpha)*(nu+1)*(nu+3)), -1/2*(1/2*trigamma((nu+1)/2)-1/2*trigamma(nu/2))-(nu+5)/(2*nu*(nu+1)*(nu+3))),
               nrow=4,ncol=4,byrow=T)
    
  } else if(margfunc=="normal"){
    
    FI=matrix( c( (2*alpha*(1-alpha)*1/2)/(phi^2), 0, -2*sqrt(2/pi)/phi,
                  0, 1/(phi^2)*(2*3/2-1), -(1-2*alpha)*(2*3/2-1)/(alpha*(1-alpha)*phi),
                  -2*sqrt(2/pi)/phi, -(1-2*alpha)*(2*3/2-1)/(alpha*(1-alpha)*phi), ((alpha^3+(1-alpha)^3)*2*3/2-(1-2*alpha)^2)/(alpha^2*(1-alpha)^2))
               ,nrow=3,ncol=3,byrow=T)
    
    
  } else if(margfunc=="logistic"){
    
    FI=matrix( c( (2*alpha*(1-alpha)*1/6)/(phi^2), 0, -2*(1/6+log(2)/3)/phi,
                  0, 1/(phi^2)*(2*(2/3+pi^2/18)-1), -(1-2*alpha)*(2*(2/3+pi^2/18)-1)/(alpha*(1-alpha)*phi),
                  -2*(1/6+log(2)/3)/phi, -(1-2*alpha)*(2*(2/3+pi^2/18)-1)/(alpha*(1-alpha)*phi), ((alpha^3+(1-alpha)^3)*2*(2/3+pi^2/18)-(1-2*alpha)^2)/(alpha^2*(1-alpha)^2))
               ,nrow=3,ncol=3,byrow=T)
    
  } else if(margfunc=="laplace"){
    
    FI=matrix( c( (2*alpha*(1-alpha)*1/2)/(phi^2), 0, -2*1/2/phi,
                  0, 1/(phi^2)*(2*1-1), -(1-2*alpha)*(2*1-1)/(alpha*(1-alpha)*phi),
                  -2*1/2/phi, -(1-2*alpha)*(2*1-1)/(alpha*(1-alpha)*phi), ((alpha^3+(1-alpha)^3)*2*1-(1-2*alpha)^2)/(alpha^2*(1-alpha)^2))
               ,nrow=3,ncol=3,byrow=T)
    
  } else {
    
    stop("incorrect margin")
    
  }
  
  var=diag(solve(FI))
  return(var)
  
  
}




















## Function to calculate the empirical covariance matrix 
## of the parameters estimated by means of IFM of five   
## Archimedean copulas                                   
#########################################################

FI_IFM=function(Z,alpha,mu,phi,df=NULL,coppars,margfuncs,copula){
  
  
  # in this, the empirical FI matrix of IFM using QBA margins is calculated.
  # This is the Godambe information matrix where expectations are replaced
  # by their empirical counterpart (averages of sums of derivatives)
  
  
  # Z is a nxd matrix containing the observations in the rows
  # alpha, mu, phi and optional df are the quantile-based asymmetric margin parameters,
  #                                numerical vectors of length d
  # margfuncs is a character vector containing the type of quantile-based margin
  # options for elements are "normal", "logistic", "laplace" or "t"
  # coppars is a numeric vector containing the copula parameters
  # copula is the family of archimedean copulas. Options are: "clayton", "AMH"
  #      "gumbel", "joe", "frank
  
  
  
  ### extras: dimensions and pseudo-observations ###
  ##################################################
  
  # lengths of vectors
  n=length(Z[,1])
  d=length(Z[1,])
  indt=which(margfuncs=='t')
  lt=length(indt)
  
  # degrees of freedom
  df=na.omit(df)
  if(lt==length(df)){
    df2=rep(NA,d)
    df2[indt]=df
    df=df2
  } else {
    stop("incorrect number of degrees of freedom provided")
  }
  
  # calculation of pseudo-observations
  U=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    U[,j]=switch(margfuncs[j],'normal'=QBAsyDist::pAND(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'logistic'=QBAsyDist::pALoD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'t'=QBAsyDist::pATD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                 ,'laplace'=QBAsyDist::pALaD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
    )
  }
  
  
  
  ### calculation of the require derivatives w.r.t. the copula parameter and margin parameters ###
  ################################################################################################
  
  ### first derivatives of the margins
  
  partialalpha=matrix(NA,nrow=n,ncol=d)
  partialmu=matrix(NA,nrow=n,ncol=d)
  partialphi=matrix(NA,nrow=n,ncol=d)
  partialdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivatives of the log-margin with respect to alpha, mu and phi
    partialalpha[,j]=(1-2*alpha[j])/(alpha[j]*(1-alpha[j])) +
      arg*(1*(Z[,j]<=mu[j])*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
             1*(Z[,j]>mu[j])*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))
    partialmu[,j]=1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    partialphi[,j]= -1/phi[j] + 1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*(-alpha[j])/phi[j]*arg*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    
    # when basefunc is student-t
    if(margfuncs[j]=="t"){
      partialdf[,j] = -1/(2*df[j]) + 1/2*digamma((df[j]+1)/2) - 1/2*digamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*(-1/2*log(1+(1-alpha[j])^2*arg^2/df[j])+(df[j]+1)/2*((1-alpha[j])^2*arg^2/df[j]^2)/(1+(1-alpha[j])^2*arg^2/df[j])) +
        1*(Z[,j]>mu[j])*(-1/2*log(1+alpha[j]^2*arg^2/df[j])+(df[j]+1)/2*(alpha[j]^2*arg^2/df[j]^2)/(1+alpha[j]^2*arg^2/df[j]))
    }
  }
  
  ### second derivatives of the margins
  partialalphaalpha=matrix(NA,nrow=n,ncol=d)
  partialalphamu=matrix(NA,nrow=n,ncol=d)
  partialalphaphi=matrix(NA,nrow=n,ncol=d)
  partialalphadf=matrix(NA,nrow=n,ncol=d)
  partialmumu=matrix(NA,nrow=n,ncol=d)
  partialmuphi=matrix(NA,nrow=n,ncol=d)
  partialmudf=matrix(NA,nrow=n,ncol=d)
  partialphiphi=matrix(NA,nrow=n,ncol=d)
  partialphidf=matrix(NA,nrow=n,ncol=d)
  partialdfdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # second order derivatives of the log-margin with respect to alpha, mu and phi
    
    # alpha-alpha
    partialalphaalpha[,j]=(-2*alpha[j]^2+2*alpha[j]-1)/(alpha[j]^2*(1-alpha[j])^2)+
      1*(Z[,j]<=mu[j])*arg^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*arg^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                               (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialalphaalpha[,j]=(2*((alpha[j])^3+(1-alpha[j])^3)-(1-2*alpha[j])^2)/(alpha[j]*(1-alpha[j]))^2
    }
    
    # alpha-mu
    partialalphamu[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                           arg*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                           (1-alpha[j])/phi[j]*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # alpha-phi
    partialalphaphi[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                            arg^2*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                            (1-alpha[j])/phi[j]*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg^2*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # mu-mu
    partialmumu[,j]=1*(Z[,j]<=mu[j])*((1-alpha[j])/phi[j])^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j])^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                                             (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialmumu[,j]=-alpha[j]*(1-alpha[j])/phi[j]^2
    }
    
    
    # mu-phi
    partialmuphi[,j]=1*(Z[,j]<=mu[j])*(-(1-alpha[j])/phi[j]^2*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                         arg*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                         (1-alpha[j])^2/phi[j]^2*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j]^2*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # phi-phi
    partialphiphi[,j]=1/phi[j]^2 + 
      1*(Z[,j]<=mu[j])*(-2*(1-alpha[j])/phi[j]^2*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                          arg^2*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                          (1-alpha[j])^2/phi[j]^2*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(2*alpha[j]/phi[j]^2*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg^2*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # related to degrees of freedom (df) when basefunc is student-t 
    if(margfuncs[j]=='t'){
      
      L=1+(1-alpha[j])^2/df[j]*arg^2
      R=1+alpha[j]^2/df[j]*arg^2
      
      # alpha-df
      partialalphadf[,j]=1*(Z[,j]<=mu[j])*( (1-alpha[j])/df[j]^2*arg^2/L + (1-alpha[j])^3*(df[j]+1)/df[j]^3*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( alpha[j]/df[j]^2*arg^2/R - alpha[j]^3*(df[j]+1)/df[j]^3*arg^4/R^2 )
      
      # mu-df
      partialmudf[,j]=1*(Z[,j]<=mu[j])*( (-1/phi[j]+(df[j]+1)/(phi[j]*df[j]))*(1-alpha[j])^2/df[j]*arg/L + (1-alpha[j])^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/L^2 ) +
        1*(Z[,j]>mu[j])*( (1/phi[j]-(df[j]+1)/(phi[j]*df[j]))*alpha[j]^2/df[j]*arg/R + alpha[j]^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/R^2 )
      
      # phi-df
      partialphidf[,j]=1*(Z[,j]<=mu[j])*( -(1-alpha[j])^2/(df[j]^2*phi[j])*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( -alpha[j]^2/(df[j]^2*phi[j])*arg^2/R + alpha[j]^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/R^2 )
      
      # df-df
      partialdfdf[,j]=1/(2*df[j]^2) + 1/4*trigamma((df[j]+1)/2) - 1/4*trigamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*( ( (1-alpha[j])^2/df[j]^2 - (1-alpha[j])^2*(df[j]+1)/df[j]^3 )*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(2*df[j]^4)*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( ( alpha[j]^2/df[j]^2 - alpha[j]^2*(df[j]+1)/df[j]^3 )*arg^2/R + alpha[j]^4*(df[j]+1)/(2*df[j]^4)*arg^4/R^2 )
      
    }
    
  }
  
  ### derivatives of the log-copula density 
  
  # additional function for the student-t case as it requires the calculation of an integral
  derivIBeta=function(Z,nu){
    
    # function to numerically approximate the derivate of the incomplete beta function
    # approximation is in the integral term, other terms are calculated exactly
    
    # number of observations
    n=length(Z)
    
    # function passed down to integrate
    intfun=function(s,nu){
      return(log(s)*s^(nu/2-1)*(1-s)^(-1/2))
    }
    
    # exact part of the derivative
    gz=nu/(nu+Z)
    out=(1-gz)^(-1/2)*gz^(nu/2-1)*Z/(nu+Z)^2
    
    # calculation of the integral part per observation
    for(j in 1:n){
      x=1/2*integrate(f = intfun,lower = 0,upper = gz[j],nu=nu,subdivisions = 1000)$value
      out[j]=out[j]+x
    }
    
    return(out)
  }
  
  
  
  # derivatives of pseudo-obvervations U w.r.t the margin parameters
  Ualpha=matrix(NA,nrow=n,ncol=d)
  Umu=matrix(NA,nrow=n,ncol=d)
  Uphi=matrix(NA,nrow=n,ncol=d)
  Udf=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivative of u_j w.r.t eta_j
    temp=switch(margfuncs[j],
                'normal'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pnorm(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dnorm(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*pnorm(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dnorm(x = alpha[j]*arg))
                  , -QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'logistic'=cbind(
                  1*(Z[,j]<=mu[j])*(2*plogis(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dlogis(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*plogis(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dlogis(x = alpha[j]*arg))
                  , -QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'t'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pt(q = (1-alpha[j])*arg,df = df[j])-2*alpha[j]*arg*dt(x = (1-alpha[j])*arg,df = df[j])) +
                    1*(Z[,j]>mu[j])*(2-2*pt(q = alpha[j]*arg,df = df[j])+2*(1-alpha[j])*arg*dt(x = alpha[j]*arg,df = df[j]))
                  , -QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , -arg*QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , 1*(Z[,j]<=mu[j])*alpha[j]*( beta(df[j]/2,1/2)*derivIBeta(Z=((1-alpha[j])*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+((1-alpha[j])*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2 -
                    1*(Z[,j]>mu[j])*(1-alpha[j])*( beta(df[j]/2,1/2)*derivIBeta(Z=(alpha[j]*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+(alpha[j]*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2
                )
                ,'laplace'=cbind(
                  1*(Z[,j]<=mu[j])*(2*LaplacesDemon::plaplace(q = (1-alpha[j])*arg)-2*alpha[j]*arg*LaplacesDemon::dlaplace(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*LaplacesDemon::plaplace(q = alpha[j]*arg)+2*(1-alpha[j])*arg*LaplacesDemon::dlaplace(x = alpha[j]*arg))
                  , -QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
    )
    
    Ualpha[,j]=temp[,1]
    Umu[,j]=temp[,2]
    Uphi[,j]=temp[,3]
    if(margfuncs[j]=='t'){
      Udf[,j]=temp[,4]
    }
    
  }
  
  derivatives=switch(copula,"clayton"=derivsClayton(theta = coppars,U = U)
                     ,"joe"=derivsJoe(theta = coppars,U = U)
                     ,"frank"=derivsFrank(theta = coppars,U = U)
                     ,"gumbel"=derivsGumbel(theta = coppars,U = U)
                     ,"AMH"=derivsAMH(theta = coppars,U = U)
  )
  
  # first derivative of the log-copula density w.r.t. the copula parameter
  partialcop=derivatives$partialcop
  
  # second derivative of the log-copula density w.r.t. the copula parameter
  partialcopcop=derivatives$partialcopcop
  
  # mixed derivative of the log-copula density w.r.t. copula parameter and u_j (in the j-th column)
  partialcopu=derivatives$partialcopu
  
  # observations which do not give NA in the derivatives
  ii=which(!is.na(partialcop))
  
  
  # second derivative of the copula with respect to the copula parameters and margin parameters
  partialcopeta=matrix(NA,nrow=length(ii),ncol=3*d+lt)
  counter=1
  for(j in 1:d){
    
    cc=2
    
    # this term is a common denominator in all partial derivatives, 
    # the derivate of the log-copula w.r.t u_j
    repeating=partialcopu[ii,j]
    
    
    # copula parameter-alpha
    partialcopalpha=repeating*Ualpha[ii,j]
    
    # copula parameter-mu
    partialcopmu=repeating*Umu[ii,j]
    
    # copula parameter-phi
    partialcopphi=repeating*Uphi[ii,j]
    
    # combining the previous three
    temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
    
    # when basefunc is student-t
    if(margfuncs[j]=='t'){
      cc=3
      
      # copula parameter-df
      partialcopdf=repeating*Udf[ii,j]
      temp=cbind(temp,partialcopdf)
    }
    
    
    
    partialcopeta[,counter:(counter+cc)]=temp
    counter=counter+cc+1
  }
  
  
  
  
  
  ##############################
  ### calculation of the FIM ###
  ##############################
  
  # number of margin parameters
  nm=3*d+sum(margfuncs=='t')
  
  
  
  ### construction of K_eta ###
  #############################
  
  # block by block corresponding to the j-th margin
  temp=c()
  for(j in 1:d){
    if(margfuncs[j]=='t'){
      temp=cbind(temp,partialalpha[ii,j],partialmu[,j],partialphi[ii,j],partialdf[ii,j])
    } else {
      temp=cbind(temp,partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
    }
  }
  K_eta=1/(length(ii))*t(temp)%*%temp
  
  # block by block corresponding to the j-th margin
  # K_eta=matrix(0,nrow = nm,ncol = nm)
  # counter=1
  # for(j in 1:d){
  # 
  #   cc=2
  #   if(margfuncs[j]=='t'){
  #     cc=3
  #   }
  #   temp2=t(temp[,counter:(counter+cc)])%*%temp[,counter:(counter+cc)]
  # 
  #   K_eta[counter:(counter+cc),counter:(counter+cc)]=temp2
  #   counter=counter+cc+1
  # }
  # K_eta=1/(length(ii))*K_eta
  # 
  
  ### construction of K_eta_thetaC ###
  ####################################
  
  # counter for creation of the full matrix
  counter=1
  
  # full matrix
  K_eta_thetaC=matrix(0,nrow=nm,ncol=1)
  
  # block by block corresponding to the j-th margin
  for(j in 1:d){
    if(margfuncs[j]=='t'){
      cc=3
      temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
    } else {
      cc=2
      temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
    }
    K_eta_thetaC[counter:(counter+cc),1]=1/(length(ii))*t(temp)%*%partialcop[ii]
    counter=counter+cc+1
  }
  
  
  
  ### construction of K_thetaC ###
  ################################
  K_thetaC=1/(length(ii))*t(partialcop[ii])%*%partialcop[ii]
  
  
  
  ### construction of K ###
  #########################
  K=matrix(NA,nrow=nm+1,ncol=nm+1)
  K[1:nm,1:nm]=K_eta
  K[nm+1,1:nm]=t(K_eta_thetaC) # when the model is correctly specified, asymptotically this is 0
  K[1:nm,nm+1]=K_eta_thetaC    # when the model is correctly specified, asymptotically this is 0
  K[nm+1,nm+1]=K_thetaC
  
  
  
  ### construction of I_eta ###
  #############################
  
  # counter for creation of the full matrix
  counter=1
  
  # full matrix
  I_eta=matrix(0,nrow=nm,ncol=nm)
  
  # minus average of second order derivatives
  paa=-1/(length(ii))*apply(partialalphaalpha[ii,],2,sum)
  pam=-1/(length(ii))*apply(partialalphamu[ii,],2,sum)
  pap=-1/(length(ii))*apply(partialalphaphi[ii,],2,sum)
  pad=-1/(length(ii))*apply(partialalphadf[ii,],2,sum)
  pmm=-1/(length(ii))*apply(partialmumu[ii,],2,sum)
  pmp=-1/(length(ii))*apply(partialmuphi[ii,],2,sum)
  pmd=-1/(length(ii))*apply(partialmudf[ii,],2,sum)
  ppp=-1/(length(ii))*apply(partialphiphi[ii,],2,sum)
  ppd=-1/(length(ii))*apply(partialphidf[ii,],2,sum)
  pdd=-1/(length(ii))*apply(partialdfdf[ii,],2,sum)
  
  
  # block by block corresponding to the j-th margin
  for(j in 1:d){
    
    cc=2
    temp=matrix(NA,nrow=3,ncol=3)
    
    if(margfuncs[j]=='t'){
      cc=3
      temp=matrix(0,nrow=4,ncol=4)
      temp[1,4]=temp[4,1]=pad[j]
      temp[2,4]=temp[4,2]=pmd[j]
      temp[3,4]=temp[4,3]=ppd[j]
      temp[4,4]=pdd[j]
    } 
    
    temp[1,1]=paa[j]
    temp[1,2]=temp[2,1]=pam[j]
    temp[1,3]=temp[3,1]=pap[j]
    temp[2,2]=pmm[j]
    temp[2,3]=temp[3,2]=pmp[j]
    temp[3,3]=ppp[j]
    
    I_eta[counter:(counter+cc),counter:(counter+cc)]=temp
    counter=counter+cc+1
  }
  
  
  
  ### construction of I_eta_thetaC ###
  ####################################
  
  # full matrix
  I_eta_thetaC=-1/(length(ii))*apply(partialcopeta,2,sum)
  
  
  
  ### construction of I_thetaC ###
  ################################
  I_thetaC=-1/(length(ii))*sum(partialcopcop[ii])
  
  
  
  ### construction of G ###
  #########################
  G=matrix(0,nrow=nm+1,ncol=nm+1)
  
  G[1:nm,1:nm]=I_eta
  G[nm+1,1:nm]=t(I_eta_thetaC)
  G[nm+1,nm+1]=I_thetaC
  
  
  FIM=t(G)%*%solve(K)%*%G
  
  return(list('FIM'=FIM,'G'=G,'K'=K,'not used'=(1:n)[-ii]))
  
}








## Function to calculate the empirical covariance matrix 
## of the parameters estimated by means of IFM of the   
## Gaussian copula                                       
#########################################################

FI_IFM_Gaussian=function(Z,alpha,mu,phi,df=NULL,coppars,margfuncs){
  
  
  # in this, the empirical FI matrix of IFM using QBA margins is calculated.
  # This is the Godambe information matrix where expectations are replaced
  # by their empirical counterpart (averages of sums of derivatives)
  
  
  # Z is a nxd matrix containing the observations in the rows
  # alpha, mu, phi and optional df are the quantile-based asymmetric margin parameters,
  #                                numerical vectors of length d
  # margfuncs is a character vector containing the type of quantile-based margin
  # options are "normal", "logistic", "laplace" or "t"
  # coppars is a numeric vector containing the copula parameters
  
  
  
  ### extras: dimensions and pseudo-observations ###
  ##################################################
  
  # lengths of vectors
  n=length(Z[,1])
  d=length(Z[1,])
  indt=which(margfuncs=='t')
  lt=length(indt)
  
  # degrees of freedom
  df=na.omit(df)
  if(lt==length(df)){
    df2=rep(NA,d)
    df2[indt]=df
    df=df2
  } else {
    stop("incorrect number of degrees of freedom provided")
  }
  
  # calculation of pseudo-observations
  U=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    U[,j]=switch(margfuncs[j],'normal'=QBAsyDist::pAND(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'logistic'=QBAsyDist::pALoD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'t'=QBAsyDist::pATD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                 ,'laplace'=QBAsyDist::pALaD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
    )
  }
  
  
  
  ### calculation of the require derivatives w.r.t. the copula parameter and margin parameters ###
  ################################################################################################
  
  ### first derivatives of the margins
  
  partialalpha=matrix(NA,nrow=n,ncol=d)
  partialmu=matrix(NA,nrow=n,ncol=d)
  partialphi=matrix(NA,nrow=n,ncol=d)
  partialdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivatives of the log-margin with respect to alpha, mu and phi
    partialalpha[,j]=(1-2*alpha[j])/(alpha[j]*(1-alpha[j])) +
      arg*(1*(Z[,j]<=mu[j])*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
             1*(Z[,j]>mu[j])*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))
    partialmu[,j]=1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    partialphi[,j]= -1/phi[j] + 1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*arg*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    
    # when basefunc is student-t
    if(margfuncs[j]=="t"){
      partialdf[,j] = -1/(2*df[j]) + 1/2*digamma((df[j]+1)/2) - 1/2*digamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*(-1/2*log(1+(1-alpha[j])^2*arg^2/df[j])+(df[j]+1)/2*((1-alpha[j])^2*arg^2/df[j]^2)/(1+(1-alpha[j])^2*arg^2/df[j])) +
        1*(Z[,j]>mu[j])*(-1/2*log(1+alpha[j]^2*arg^2/df[j])+(df[j]+1)/2*(alpha[j]^2*arg^2/df[j]^2)/(1+alpha[j]^2*arg^2/df[j]))
    }
  }
  
  ### second derivatives of the margins
  partialalphaalpha=matrix(NA,nrow=n,ncol=d)
  partialalphamu=matrix(NA,nrow=n,ncol=d)
  partialalphaphi=matrix(NA,nrow=n,ncol=d)
  partialalphadf=matrix(NA,nrow=n,ncol=d)
  partialmumu=matrix(NA,nrow=n,ncol=d)
  partialmuphi=matrix(NA,nrow=n,ncol=d)
  partialmudf=matrix(NA,nrow=n,ncol=d)
  partialphiphi=matrix(NA,nrow=n,ncol=d)
  partialphidf=matrix(NA,nrow=n,ncol=d)
  partialdfdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # second order derivatives of the log-margin with respect to alpha, mu and phi
    
    # alpha-alpha
    partialalphaalpha[,j]=(-2*alpha[j]^2+2*alpha[j]-1)/(alpha[j]^2*(1-alpha[j])^2)+
      1*(Z[,j]<=mu[j])*arg^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*arg^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                               (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialalphaalpha[,j]=(2*((alpha[j])^3+(1-alpha[j])^3)-(1-2*alpha[j])^2)/(alpha[j]*(1-alpha[j]))^2
    }
    
    # alpha-mu
    partialalphamu[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                           arg*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                           (1-alpha[j])/phi[j]*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # alpha-phi
    partialalphaphi[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                            arg^2*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                            (1-alpha[j])/phi[j]*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg^2*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # mu-mu
    partialmumu[,j]=1*(Z[,j]<=mu[j])*((1-alpha[j])/phi[j])^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j])^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                                             (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialmumu[,j]=-alpha[j]*(1-alpha[j])/phi[j]^2
    }
    
    
    # mu-phi
    partialmuphi[,j]=1*(Z[,j]<=mu[j])*(-(1-alpha[j])/phi[j]^2*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                         arg*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                         (1-alpha[j])^2/phi[j]^2*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j]^2*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # phi-phi
    partialphiphi[,j]=1/phi[j]^2 + 
      1*(Z[,j]<=mu[j])*(-2*(1-alpha[j])/phi[j]^2*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                          arg^2*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                          (1-alpha[j])^2/phi[j]^2*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(2*alpha[j]/phi[j]^2*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg^2*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # related to degrees of freedom (df) when basefunc is student-t 
    if(margfuncs[j]=='t'){
      
      L=1+(1-alpha[j])^2/df[j]*arg^2
      R=1+alpha[j]^2/df[j]*arg^2
      
      # alpha-df
      partialalphadf[,j]=1*(Z[,j]<=mu[j])*( (1-alpha[j])/df[j]^2*arg^2/L + (1-alpha[j])^3*(df[j]+1)/df[j]^3*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( alpha[j]/df[j]^2*arg^2/R - alpha[j]^3*(df[j]+1)/df[j]^3*arg^4/R^2 )
      
      # mu-df
      partialmudf[,j]=1*(Z[,j]<=mu[j])*( (-1/phi[j]+(df[j]+1)/(phi[j]*df[j]))*(1-alpha[j])^2/df[j]*arg/L + (1-alpha[j])^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/L^2 ) +
        1*(Z[,j]>mu[j])*( (1/phi[j]-(df[j]+1)/(phi[j]*df[j]))*alpha[j]^2/df[j]*arg/R + alpha[j]^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/R^2 )
      
      # phi-df
      partialphidf[,j]=1*(Z[,j]<=mu[j])*( -(1-alpha[j])^2/(df[j]^2*phi[j])*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( -alpha[j]^2/(df[j]^2*phi[j])*arg^2/R + alpha[j]^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/R^2 )
      
      # df-df
      partialdfdf[,j]=1/(2*df[j]^2) + 1/4*trigamma((df[j]+1)/2) - 1/4*trigamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*( ( (1-alpha[j])^2/df[j]^2 - (1-alpha[j])^2*(df[j]+1)/df[j]^3 )*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(2*df[j]^4)*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( ( alpha[j]^2/df[j]^2 - alpha[j]^2*(df[j]+1)/df[j]^3 )*arg^2/R + alpha[j]^4*(df[j]+1)/(2*df[j]^4)*arg^4/R^2 )
      
    }
    
  }
  
  ### derivatives of the log-copula density 
  
  # additional function for the student-t case as it requires the calculation of an integral
  derivIBeta=function(Z,nu){
    
    # function to numerically approximate the derivate of the incomplete beta function
    # approximation is in the integral term, other terms are calculated exactly
    
    # number of observations
    n=length(Z)
    
    # function passed down to integrate
    intfun=function(s,nu){
      return(log(s)*s^(nu/2-1)*(1-s)^(-1/2))
    }
    
    # exact part of the derivative
    gz=nu/(nu+Z)
    out=(1-gz)^(-1/2)*gz^(nu/2-1)*Z/(nu+Z)^2
    
    # calculation of the integral part per observation
    for(j in 1:n){
      x=1/2*integrate(f = intfun,lower = 0,upper = gz[j],nu=nu,subdivisions = 1000)$value
      out[j]=out[j]+x
    }
    
    return(out)
  }
  
  
  
  # derivatives of pseudo-obvervations U w.r.t the margin parameters
  Ualpha=matrix(NA,nrow=n,ncol=d)
  Umu=matrix(NA,nrow=n,ncol=d)
  Uphi=matrix(NA,nrow=n,ncol=d)
  Udf=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivative of u_j w.r.t eta_j
    temp=switch(margfuncs[j],
                'normal'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pnorm(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dnorm(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*pnorm(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dnorm(x = alpha[j]*arg))
                  , -QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'logistic'=cbind(
                  1*(Z[,j]<=mu[j])*(2*plogis(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dlogis(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*plogis(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dlogis(x = alpha[j]*arg))
                  , -QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'t'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pt(q = (1-alpha[j])*arg,df = df[j])-2*alpha[j]*arg*dt(x = (1-alpha[j])*arg,df = df[j])) +
                    1*(Z[,j]>mu[j])*(2-2*pt(q = alpha[j]*arg,df = df[j])+2*(1-alpha[j])*arg*dt(x = alpha[j]*arg,df = df[j]))
                  , -QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , -arg*QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , 1*(Z[,j]<=mu[j])*alpha[j]*( beta(df[j]/2,1/2)*derivIBeta(Z=((1-alpha[j])*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+((1-alpha[j])*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2 -
                    1*(Z[,j]>mu[j])*(1-alpha[j])*( beta(df[j]/2,1/2)*derivIBeta(Z=(alpha[j]*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+(alpha[j]*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2
                )
                ,'laplace'=cbind(
                  1*(Z[,j]<=mu[j])*(2*LaplacesDemon::plaplace(q = (1-alpha[j])*arg)-2*alpha[j]*arg*LaplacesDemon::dlaplace(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*LaplacesDemon::plaplace(q = alpha[j]*arg)+2*(1-alpha[j])*arg*LaplacesDemon::dlaplace(x = alpha[j]*arg))
                  , -QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
    )
    
    Ualpha[,j]=temp[,1]
    Umu[,j]=temp[,2]
    Uphi[,j]=temp[,3]
    if(margfuncs[j]=='t'){
      Udf[,j]=temp[,4]
    }
    
  }
  
  derivatives=derivsgaussian(theta = coppars,U = U)
  
  # first derivative of the log-copula density w.r.t. the copula parameter
  # rows for observations
  # columns for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}  
  partialcop=derivatives$partialcop
  
  # second derivative of the log-copula density w.r.t. the copula parameter
  # first dimension is for observations
  # second and third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcopcop=derivatives$partialcopcop
  
  # mixed derivative of the log-copula density w.r.t. copula parameter and u_j
  # first dimension is for the observations
  # second dimension is for u_1,...,u_d
  # third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcopu=derivatives$partialcopu
  
  # observations which do not give NA in the derivatives
  ii=which(!is.na(apply(partialcop,1,sum)))
  
  
  # second derivative of the copula with respect to the copula parameters and margin parameters
  # first dimension for observations
  # second dimension for alpha_1,..,kappa_1,alpha_2,...,kappa_d
  # third dimension for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
  partialcopeta=array(NA,dim =c(length(ii),ncol=3*d+lt,d*(d-1)/2))
  for(k in 1:(d*(d-1)/2)){
    counter=1
    for(j in 1:d){
      cc=2
      
      # this term is a common denominator in all partial derivatives, 
      # the derivate of the log-copula w.r.t u_j
      repeating=partialcopu[ii,j,k]
      
      
      # copula parameter-alpha
      partialcopalpha=repeating*Ualpha[ii,j]
      
      # copula parameter-mu
      partialcopmu=repeating*Umu[ii,j]
      
      # copula parameter-phi
      partialcopphi=repeating*Uphi[ii,j]
      
      # combining the previous three
      temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
      
      # when basefunc is student-t
      if(margfuncs[j]=='t'){
        cc=3
        
        # copula parameter-df
        partialcopdf=repeating*Udf[ii,j]
        temp=cbind(temp,partialcopdf)
      }
      
      
      
      partialcopeta[,counter:(counter+cc),k]=temp
      counter=counter+cc+1
    }
    
  }
  
  
  
  
  
  ##############################
  ### calculation of the FIM ###
  ##############################
  
  # number of margin parameters
  nm=3*d+sum(margfuncs=='t')
  
  
  
  ### construction of K_eta ###
  #############################
  
  # block by block corresponding to the j-th margin
  temp=c()
  for(j in 1:d){
    if(margfuncs[j]=='t'){
      temp=cbind(temp,partialalpha[ii,j],partialmu[,j],partialphi[ii,j],partialdf[ii,j])
    } else {
      temp=cbind(temp,partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
    }
  }
  K_eta=1/(length(ii))*t(temp)%*%temp
  
  
  
  ### construction of K_eta_thetaC ###
  ####################################
  
  
  # full matrix
  K_eta_thetaC=matrix(0,nrow=nm,ncol=d*(d-1)/2)
  
  for(k in 1:(d*(d-1)/2)){
    
    # counter for creation of the full matrix
    counter=1
    
    # block by block corresponding to the j-th margin
    for(j in 1:d){
      if(margfuncs[j]=='t'){
        cc=3
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
      } else {
        cc=2
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
      }
      K_eta_thetaC[counter:(counter+cc),k]=1/(length(ii))*t(temp)%*%partialcop[ii,k]
      counter=counter+cc+1
    }
    
  }
  
  
  
  
  ### construction of K_thetaC ###
  ################################
  K_thetaC=1/(length(ii))*t(partialcop[ii,])%*%partialcop[ii,]
  
  
  
  ### construction of K ###
  #########################
  K=matrix(NA,nrow=nm+d*(d-1)/2,ncol=nm+d*(d-1)/2)
  K[1:nm,1:nm]=K_eta
  K[(nm+1):(nm+d*(d-1)/2),1:nm]=t(K_eta_thetaC) 
  K[1:nm,(nm+1):(nm+d*(d-1)/2)]=K_eta_thetaC    
  K[(nm+1):(nm+d*(d-1)/2),(nm+1):(nm+d*(d-1)/2)]=K_thetaC
  
  
  
  ### construction of I_eta ###
  #############################
  
  # counter for creation of the full matrix
  counter=1
  
  # full matrix
  I_eta=matrix(0,nrow=nm,ncol=nm)
  
  # minus average of second order derivatives
  paa=-1/(length(ii))*apply(partialalphaalpha[ii,],2,sum)
  pam=-1/(length(ii))*apply(partialalphamu[ii,],2,sum)
  pap=-1/(length(ii))*apply(partialalphaphi[ii,],2,sum)
  pad=-1/(length(ii))*apply(partialalphadf[ii,],2,sum)
  pmm=-1/(length(ii))*apply(partialmumu[ii,],2,sum)
  pmp=-1/(length(ii))*apply(partialmuphi[ii,],2,sum)
  pmd=-1/(length(ii))*apply(partialmudf[ii,],2,sum)
  ppp=-1/(length(ii))*apply(partialphiphi[ii,],2,sum)
  ppd=-1/(length(ii))*apply(partialphidf[ii,],2,sum)
  pdd=-1/(length(ii))*apply(partialdfdf[ii,],2,sum)
  
  
  # block by block corresponding to the j-th margin
  for(j in 1:d){
    
    cc=2
    temp=matrix(NA,nrow=3,ncol=3)
    
    if(margfuncs[j]=='t'){
      cc=3
      temp=matrix(0,nrow=4,ncol=4)
      temp[1,4]=temp[4,1]=pad[j]
      temp[2,4]=temp[4,2]=pmd[j]
      temp[3,4]=temp[4,3]=ppd[j]
      temp[4,4]=pdd[j]
    } 
    
    temp[1,1]=paa[j]
    temp[1,2]=temp[2,1]=pam[j]
    temp[1,3]=temp[3,1]=pap[j]
    temp[2,2]=pmm[j]
    temp[2,3]=temp[3,2]=pmp[j]
    temp[3,3]=ppp[j]
    
    I_eta[counter:(counter+cc),counter:(counter+cc)]=temp
    counter=counter+cc+1
  }
  
  
  
  ### construction of I_eta_thetaC ###
  ####################################
  if(d==2){
    I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2),sum)
  } else {
    I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2,3),sum)
  }
  
  
  
  
  ### construction of I_thetaC ###
  ################################
  if(d==2){
    I_thetaC=-1/(length(ii))*sum(partialcopcop[ii,,])
  } else {
    I_thetaC=-1/(length(ii))*apply(partialcopcop[ii,,],c(2,3),sum)
  }
  
  
  
  
  ### construction of G ###
  #########################
  G=matrix(0,nrow=nm+d*(d-1)/2,ncol=nm+d*(d-1)/2)
  
  G[1:nm,1:nm]=I_eta
  G[(nm+1):(nm+d*(d-1)/2),1:nm]=t(I_eta_thetaC)
  G[(nm+1):(nm+d*(d-1)/2),(nm+1):(nm+d*(d-1)/2)]=I_thetaC
  
  
  FIM=t(G)%*%solve(K)%*%G
  
  return(list('FIM'=FIM,'G'=G,'K'=K,'not used'=(1:n)[-ii]))
  
}


### Function to calculate the empirical covariance matrix 
### of the parameters estimated by means of IFM of the    
### Gaussian copula                                       
#########################################################

FI_IFM_Student=function(Z,alpha,mu,phi,df=NULL,coppars,margfuncs){
  
  
  # in this, the empirical FI matrix of IFM using QBA margins is calculated.
  # This is the Godambe information matrix where expectations are replaced
  # by their empirical counterpart (averages of sums of derivatives)
  
  
  # Z is a nxd matrix containing the observations in the rows
  # alpha, mu, phi and optional df are the quantile-based asymmetric margin parameters,
  #                                numerical vectors of length d
  # margfuncs is a character vector containing the type of quantile-based margin
  # options are "normal", "logistic", "laplace" or "t"
  # coppars is a numeric vector containing the copula parameters + degrees of freedom
  # of the copula
  
  
  
  ### extras: dimensions and pseudo-observations ###
  ##################################################
  
  # lengths of vectors
  n=length(Z[,1])
  d=length(Z[1,])
  indt=which(margfuncs=='t')
  lt=length(indt)
  
  # degrees of freedom
  df=na.omit(df)
  if(lt==length(df)){
    df2=rep(NA,d)
    df2[indt]=df
    df=df2
  } else {
    stop("incorrect number of degrees of freedom provided")
  }
  
  # degrees of freedom of copula
  nuc=coppars[d*(d-1)/2+1]
  
  # calculation of pseudo-observations
  U=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    U[,j]=switch(margfuncs[j],'normal'=QBAsyDist::pAND(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'logistic'=QBAsyDist::pALoD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'t'=QBAsyDist::pATD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                 ,'laplace'=QBAsyDist::pALaD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
    )
  }
  
  
  
  ### calculation of the require derivatives w.r.t. the copula parameter and margin parameters ###
  ################################################################################################
  
  ### first derivatives of the margins
  
  partialalpha=matrix(NA,nrow=n,ncol=d)
  partialmu=matrix(NA,nrow=n,ncol=d)
  partialphi=matrix(NA,nrow=n,ncol=d)
  partialdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivatives of the log-margin with respect to alpha, mu and phi
    partialalpha[,j]=(1-2*alpha[j])/(alpha[j]*(1-alpha[j])) +
      arg*(1*(Z[,j]<=mu[j])*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
             1*(Z[,j]>mu[j])*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))
    partialmu[,j]=1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    partialphi[,j]= -1/phi[j] + 1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*arg*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    
    # when basefunc is student-t
    if(margfuncs[j]=="t"){
      partialdf[,j] = -1/(2*df[j]) + 1/2*digamma((df[j]+1)/2) - 1/2*digamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*(-1/2*log(1+(1-alpha[j])^2*arg^2/df[j])+(df[j]+1)/2*((1-alpha[j])^2*arg^2/df[j]^2)/(1+(1-alpha[j])^2*arg^2/df[j])) +
        1*(Z[,j]>mu[j])*(-1/2*log(1+alpha[j]^2*arg^2/df[j])+(df[j]+1)/2*(alpha[j]^2*arg^2/df[j]^2)/(1+alpha[j]^2*arg^2/df[j]))
    }
  }
  
  ### second derivatives of the margins
  partialalphaalpha=matrix(NA,nrow=n,ncol=d)
  partialalphamu=matrix(NA,nrow=n,ncol=d)
  partialalphaphi=matrix(NA,nrow=n,ncol=d)
  partialalphadf=matrix(NA,nrow=n,ncol=d)
  partialmumu=matrix(NA,nrow=n,ncol=d)
  partialmuphi=matrix(NA,nrow=n,ncol=d)
  partialmudf=matrix(NA,nrow=n,ncol=d)
  partialphiphi=matrix(NA,nrow=n,ncol=d)
  partialphidf=matrix(NA,nrow=n,ncol=d)
  partialdfdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # second order derivatives of the log-margin with respect to alpha, mu and phi
    
    # alpha-alpha
    partialalphaalpha[,j]=(-2*alpha[j]^2+2*alpha[j]-1)/(alpha[j]^2*(1-alpha[j])^2)+
      1*(Z[,j]<=mu[j])*arg^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*arg^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                               (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialalphaalpha[,j]=(2*((alpha[j])^3+(1-alpha[j])^3)-(1-2*alpha[j])^2)/(alpha[j]*(1-alpha[j]))^2
    }
    
    # alpha-mu
    partialalphamu[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                           arg*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                           (1-alpha[j])/phi[j]*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # alpha-phi
    partialalphaphi[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                            arg^2*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                            (1-alpha[j])/phi[j]*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg^2*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # mu-mu
    partialmumu[,j]=1*(Z[,j]<=mu[j])*((1-alpha[j])/phi[j])^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j])^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                                             (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialmumu[,j]=-alpha[j]*(1-alpha[j])/phi[j]^2
    }
    
    
    # mu-phi
    partialmuphi[,j]=1*(Z[,j]<=mu[j])*(-(1-alpha[j])/phi[j]^2*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                         arg*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                         (1-alpha[j])^2/phi[j]^2*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j]^2*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # phi-phi
    partialphiphi[,j]=1/phi[j]^2 + 
      1*(Z[,j]<=mu[j])*(-2*(1-alpha[j])/phi[j]^2*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                          arg^2*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                          (1-alpha[j])^2/phi[j]^2*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(2*alpha[j]/phi[j]^2*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg^2*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # related to degrees of freedom (df) when basefunc is student-t 
    if(margfuncs[j]=='t'){
      
      L=1+(1-alpha[j])^2/df[j]*arg^2
      R=1+alpha[j]^2/df[j]*arg^2
      
      # alpha-df
      partialalphadf[,j]=1*(Z[,j]<=mu[j])*( (1-alpha[j])/df[j]^2*arg^2/L + (1-alpha[j])^3*(df[j]+1)/df[j]^3*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( alpha[j]/df[j]^2*arg^2/R - alpha[j]^3*(df[j]+1)/df[j]^3*arg^4/R^2 )
      
      # mu-df
      partialmudf[,j]=1*(Z[,j]<=mu[j])*( (-1/phi[j]+(df[j]+1)/(phi[j]*df[j]))*(1-alpha[j])^2/df[j]*arg/L + (1-alpha[j])^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/L^2 ) +
        1*(Z[,j]>mu[j])*( (1/phi[j]-(df[j]+1)/(phi[j]*df[j]))*alpha[j]^2/df[j]*arg/R + alpha[j]^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/R^2 )
      
      # phi-df
      partialphidf[,j]=1*(Z[,j]<=mu[j])*( -(1-alpha[j])^2/(df[j]^2*phi[j])*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( -alpha[j]^2/(df[j]^2*phi[j])*arg^2/R + alpha[j]^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/R^2 )
      
      # df-df
      partialdfdf[,j]=1/(2*df[j]^2) + 1/4*trigamma((df[j]+1)/2) - 1/4*trigamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*( ( (1-alpha[j])^2/df[j]^2 - (1-alpha[j])^2*(df[j]+1)/df[j]^3 )*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(2*df[j]^4)*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( ( alpha[j]^2/df[j]^2 - alpha[j]^2*(df[j]+1)/df[j]^3 )*arg^2/R + alpha[j]^4*(df[j]+1)/(2*df[j]^4)*arg^4/R^2 )
      
    }
    
  }
  
  ### derivatives of the log-copula density 
  
  # additional function for the student-t case as it requires the calculation of an integral
  derivIBeta=function(Z,nu){
    
    # function to numerically approximate the derivate of the incomplete beta function
    # approximation is in the integral term, other terms are calculated exactly
    
    # number of observations
    n=length(Z)
    
    # function passed down to integrate
    intfun=function(s,nu){
      return(log(s)*s^(nu/2-1)*(1-s)^(-1/2))
    }
    
    # exact part of the derivative
    gz=nu/(nu+Z)
    out=(1-gz)^(-1/2)*gz^(nu/2-1)*Z/(nu+Z)^2
    
    # calculation of the integral part per observation
    for(j in 1:n){
      x=1/2*integrate(f = intfun,lower = 0,upper = gz[j],nu=nu,subdivisions = 1000)$value
      out[j]=out[j]+x
    }
    
    return(out)
  }
  
  
  
  # derivatives of pseudo-obvervations U w.r.t the margin parameters
  Ualpha=matrix(NA,nrow=n,ncol=d)
  Umu=matrix(NA,nrow=n,ncol=d)
  Uphi=matrix(NA,nrow=n,ncol=d)
  Udf=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivative of u_j w.r.t eta_j
    temp=switch(margfuncs[j],
                'normal'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pnorm(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dnorm(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*pnorm(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dnorm(x = alpha[j]*arg))
                  , -QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'logistic'=cbind(
                  1*(Z[,j]<=mu[j])*(2*plogis(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dlogis(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*plogis(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dlogis(x = alpha[j]*arg))
                  , -QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'t'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pt(q = (1-alpha[j])*arg,df = df[j])-2*alpha[j]*arg*dt(x = (1-alpha[j])*arg,df = df[j])) +
                    1*(Z[,j]>mu[j])*(2-2*pt(q = alpha[j]*arg,df = df[j])+2*(1-alpha[j])*arg*dt(x = alpha[j]*arg,df = df[j]))
                  , -QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , -arg*QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , 1*(Z[,j]<=mu[j])*alpha[j]*( beta(df[j]/2,1/2)*derivIBeta(Z=((1-alpha[j])*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+((1-alpha[j])*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2 -
                    1*(Z[,j]>mu[j])*(1-alpha[j])*( beta(df[j]/2,1/2)*derivIBeta(Z=(alpha[j]*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+(alpha[j]*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2
                )
                ,'laplace'=cbind(
                  1*(Z[,j]<=mu[j])*(2*LaplacesDemon::plaplace(q = (1-alpha[j])*arg)-2*alpha[j]*arg*LaplacesDemon::dlaplace(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*LaplacesDemon::plaplace(q = alpha[j]*arg)+2*(1-alpha[j])*arg*LaplacesDemon::dlaplace(x = alpha[j]*arg))
                  , -QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
    )
    
    Ualpha[,j]=temp[,1]
    Umu[,j]=temp[,2]
    Uphi[,j]=temp[,3]
    if(margfuncs[j]=='t'){
      Udf[,j]=temp[,4]
    }
    
  }
  
  derivatives=derivst(theta = coppars,U = U,nu=nuc)
  
  # first derivative of the log-copula density w.r.t. the copula parameter
  # rows for observations
  # columns for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1},nu  
  partialcop=derivatives$partialcop
  
  # second derivative of the log-copula density w.r.t. the copula parameter
  # first dimension is for observations
  # second and third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}, nu
  partialcopcop=derivatives$partialcopcop
  
  # mixed derivative of the log-copula density w.r.t. copula parameter and u_j
  # first dimension is for the observations
  # second dimension is for u_1,...,u_d
  # third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}, nu
  partialcopu=derivatives$partialcopu
  
  # observations which do not give NA in the derivatives
  ii=which(!is.na(apply(partialcop,1,sum)))
  
  
  # second derivative of the copula with respect to the copula parameters and margin parameters
  # first dimension for observations
  # second dimension for alpha_1,..,kappa_1,alpha_2,...,kappa_d
  # third dimension for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1},nu
  partialcopeta=array(NA,dim =c(length(ii),ncol=3*d+lt,d*(d-1)/2+1))
  for(k in 1:(d*(d-1)/2+1)){
    counter=1
    for(j in 1:d){
      cc=2
      
      # this term is a common denominator in all partial derivatives, 
      # the derivate of the log-copula w.r.t u_j
      repeating=partialcopu[ii,j,k]
      
      
      # copula parameter-alpha
      partialcopalpha=repeating*Ualpha[ii,j]
      
      # copula parameter-mu
      partialcopmu=repeating*Umu[ii,j]
      
      # copula parameter-phi
      partialcopphi=repeating*Uphi[ii,j]
      
      # combining the previous three
      temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
      
      # when basefunc is student-t
      if(margfuncs[j]=='t'){
        cc=3
        
        # copula parameter-df
        partialcopdf=repeating*Udf[ii,j]
        temp=cbind(temp,partialcopdf)
      }
      
      
      
      partialcopeta[,counter:(counter+cc),k]=temp
      counter=counter+cc+1
    }
    
  }
  
  
  
  
  
  ##############################
  ### calculation of the FIM ###
  ##############################
  
  # number of margin parameters
  nm=3*d+sum(margfuncs=='t')
  
  
  
  ### construction of K_eta ###
  #############################
  
  # block by block corresponding to the j-th margin
  temp=c()
  for(j in 1:d){
    if(margfuncs[j]=='t'){
      temp=cbind(temp,partialalpha[ii,j],partialmu[,j],partialphi[ii,j],partialdf[ii,j])
    } else {
      temp=cbind(temp,partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
    }
  }
  K_eta=1/(length(ii))*t(temp)%*%temp
  
  
  
  ### construction of K_eta_thetaC ###
  ####################################
  
  
  # full matrix
  K_eta_thetaC=matrix(0,nrow=nm,ncol=d*(d-1)/2+1)
  
  for(k in 1:(d*(d-1)/2+1)){
    
    # counter for creation of the full matrix
    counter=1
    
    # block by block corresponding to the j-th margin
    for(j in 1:d){
      if(margfuncs[j]=='t'){
        cc=3
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
      } else {
        cc=2
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
      }
      K_eta_thetaC[counter:(counter+cc),k]=1/(length(ii))*t(temp)%*%partialcop[ii,k]
      counter=counter+cc+1
    }
    
  }
  
  
  
  
  ### construction of K_thetaC ###
  ################################
  K_thetaC=1/(length(ii))*t(partialcop[ii,])%*%partialcop[ii,]
  
  
  
  ### construction of K ###
  #########################
  K=matrix(NA,nrow=nm+d*(d-1)/2+1,ncol=nm+d*(d-1)/2+1)
  K[1:nm,1:nm]=K_eta
  K[(nm+1):(nm+d*(d-1)/2+1),1:nm]=t(K_eta_thetaC) 
  K[1:nm,(nm+1):(nm+d*(d-1)/2+1)]=K_eta_thetaC    
  K[(nm+1):(nm+d*(d-1)/2+1),(nm+1):(nm+d*(d-1)/2+1)]=K_thetaC
  
  
  
  ### construction of I_eta ###
  #############################
  
  # counter for creation of the full matrix
  counter=1
  
  # full matrix
  I_eta=matrix(0,nrow=nm,ncol=nm)
  
  # minus average of second order derivatives
  paa=-1/(length(ii))*apply(partialalphaalpha[ii,],2,sum)
  pam=-1/(length(ii))*apply(partialalphamu[ii,],2,sum)
  pap=-1/(length(ii))*apply(partialalphaphi[ii,],2,sum)
  pad=-1/(length(ii))*apply(partialalphadf[ii,],2,sum)
  pmm=-1/(length(ii))*apply(partialmumu[ii,],2,sum)
  pmp=-1/(length(ii))*apply(partialmuphi[ii,],2,sum)
  pmd=-1/(length(ii))*apply(partialmudf[ii,],2,sum)
  ppp=-1/(length(ii))*apply(partialphiphi[ii,],2,sum)
  ppd=-1/(length(ii))*apply(partialphidf[ii,],2,sum)
  pdd=-1/(length(ii))*apply(partialdfdf[ii,],2,sum)
  
  
  # block by block corresponding to the j-th margin
  for(j in 1:d){
    
    cc=2
    temp=matrix(NA,nrow=3,ncol=3)
    
    if(margfuncs[j]=='t'){
      cc=3
      temp=matrix(0,nrow=4,ncol=4)
      temp[1,4]=temp[4,1]=pad[j]
      temp[2,4]=temp[4,2]=pmd[j]
      temp[3,4]=temp[4,3]=ppd[j]
      temp[4,4]=pdd[j]
    } 
    
    temp[1,1]=paa[j]
    temp[1,2]=temp[2,1]=pam[j]
    temp[1,3]=temp[3,1]=pap[j]
    temp[2,2]=pmm[j]
    temp[2,3]=temp[3,2]=pmp[j]
    temp[3,3]=ppp[j]
    
    I_eta[counter:(counter+cc),counter:(counter+cc)]=temp
    counter=counter+cc+1
  }
  
  
  
  ### construction of I_eta_thetaC ###
  ####################################
  if(d==2){
    I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2),sum)
  } else {
    I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2,3),sum)
  }
  
  
  
  
  ### construction of I_thetaC ###
  ################################
  if(d==2){
    I_thetaC=-1/(length(ii))*sum(partialcopcop[ii,,])
  } else {
    I_thetaC=-1/(length(ii))*apply(partialcopcop[ii,,],c(2,3),sum)
  }
  
  
  
  
  ### construction of G ###
  #########################
  G=matrix(0,nrow=nm+d*(d-1)/2+1,ncol=nm+d*(d-1)/2+1)
  
  G[1:nm,1:nm]=I_eta
  G[(nm+1):(nm+d*(d-1)/2+1),1:nm]=t(I_eta_thetaC)
  G[(nm+1):(nm+d*(d-1)/2+1),(nm+1):(nm+d*(d-1)/2+1)]=I_thetaC
  
  
  FIM=t(G)%*%solve(K)%*%G
  
  return(list('FIM'=FIM,'G'=G,'K'=K,'not used'=(1:n)[-ii]))
  
}







## function for determining the best fitting quantile-
## based margin to fit a margin of the data            
#######################################################

bestmargin=function(data,crit="AIC",nstart=NULL,all=F,seed=NULL){
  # accepts a numeric vector as input which forms a single margin of the data on which 
  # the model is fit. The best suited margin is determined from a list of 4 and based
  # on either the correlation of a uniform QQ-plot (crit="QQ") of the fitted distribution
  # or AIC/BIC (crit="AIC"/"BIC") or the best one according to the p-value of a Kolmogorov-
  # Smirnov test. If all=T the criterion for all possibilities is returned so the user can 
  # choose on his/her preference, otherwise only the best one toghether with the uniform 
  # transformation of the data, the fitted parameters and the log-likelihood is returned
  
  
  # length of data vector
  n=length(data)
  
  # number of starting values
  if(is.null(nstart)){
    nstart=sqrt(n)
  }
  
  
  # choices for margin functions
  mf=c("normal","laplace","logistic","t")
  
  # fitting the the possible margins
  fitN=fitAND(data,nstart = nstart,seed = seed)
  fitLa=fitALaD(data,nstart = nstart,seed = seed)
  fitLo=fitALoD(data,nstart= nstart,seed = seed)
  fitT=fitATD(data,nstart = nstart,seed = seed)
  
  # uniform tranform (CDF values) of the data
  U=matrix(NA,nrow=n,ncol=4)
  U[,1]=QBAsyDist::pAND(q=data,mu=fitN$mu,phi=fitN$phi,alpha=fitN$alpha)
  U[,2]=QBAsyDist::pALaD(q=data,mu=fitLa$mu,phi=fitLa$phi,alpha=fitLa$alpha)
  U[,3]=QBAsyDist::pALoD(q=data,mu=fitLo$mu,phi=fitLo$phi,alpha=fitLo$alpha)
  U[,4]=QBAsyDist::pATD(q=data,mu=fitT$mu,phi=fitT$phi,alpha=fitT$alpha,nu = fitT$nu)
  
  # fitted parameters
  pars=matrix(NA,nrow=4,ncol=4)
  rownames(pars)=c('alpha','mu','phi','nu')
  pars[1:3,1]=c(fitN$alpha,fitN$mu,fitN$phi)
  pars[1:3,2]=c(fitLa$alpha,fitLa$mu,fitLa$phi)
  pars[1:3,3]=c(fitLo$alpha,fitLo$mu,fitLo$phi)
  pars[1:4,4]=c(fitT$alpha,fitT$mu,fitT$phi,fitT$nu)
  
  # log-likelihood
  logl=rep(NA,4)
  logl[1]=fitN$LogLikelihood
  logl[2]=fitLa$LogLikelihood
  logl[3]=fitLo$LogLikelihood
  logl[4]=fitT$LogLikelihood
  
  if(crit=="AIC"){
    
    AICfit=-2*c(fitN$LogLikelihood,fitLa$LogLikelihood,fitLo$LogLikelihood,fitT$LogLikelihood)+2*c(3,3,3,4)
    if(all==T){
      
      return(cbind("margin"=mf,"AIC"=AICfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.min(AICfit)],
                  "U"=U[,which.min(AICfit)],
                  "parameters"=pars[,which.min(AICfit)],
                  "logl"=logl[which.min(AICfit)]))
      
    }
    
  } else if(crit=="BIC"){
    
    BICfit=-2*c(fitN$LogLikelihood,fitLa$LogLikelihood,fitLo$LogLikelihood,fitT$LogLikelihood)+log(n)*c(3,3,3,4)
    
    if(all==T){
      
      return(cbind("margin"=mf,"BIC"=BICfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.min(BICfit)],
                  "U"=U[,which.min(BICfit)],
                  "parameters"=pars[,which.min(BICfit)],
                  "logl"=logl[which.min(BICfit)]))
      
    }
    
    
  } else if(crit=="QQ"){
    
    corfit=c(
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,1])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,2])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,3])),
      cor(seq(1/(n+1),n/(n+1),length.out=n),sort(U[,4]))
    )
    
    if(all==T){
      
      return(cbind("margin"=mf,"Correlation QQ-plot"=corfit,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.max(corfit)],
                  "U"=U[,which.max(corfit)],
                  "parameters"=pars[,which.max(corfit)],
                  "logl"=logl[which.max(corfit)]))
      
    }
    
  } else if(crit=="KS"){
    
    suppressWarnings({
      pvalue=c(ks.test(data,pAND,mu=fitN$mu,phi=fitN$phi,alpha=fitN$alpha)$p.value,
               ks.test(data,pALaD,mu=fitLa$mu,phi=fitLa$phi,alpha=fitLa$alpha)$p.value,
               ks.test(data,pALoD,mu=fitLo$mu,phi=fitLo$phi,alpha=fitLo$alpha)$p.value,
               ks.test(data,pATD,mu=fitT$mu,phi=fitT$phi,alpha=fitT$alpha,nu = fitT$nu)$p.value
      )})
    
    if(all==T){
      
      return(cbind("margin"=mf,"KS-test p-value"=pvalue,"parameters"=t(pars)))
      
    } else {
      
      return(list("margin"=mf[which.max(pvalue)],
                  "U"=U[,which.max(pvalue)],
                  "parameters"=pars[,which.max(pvalue)],
                  "logl"=logl[which.max(pvalue)]))
      
    }
  } else {
    
    stop("invalid criterion")
    
  }
  
}  






## function to determine the best fitting copula based on 
## AIC using the empirical pseudo-observation             
##########################################################

bestcopulafit=function(U,all=F,crit="AIC"){
  # accepts a nxd matrix U with the pseudo-observations from the data and determines the best
  # copula from a list of possibilities based on AIC or BIC. 
  # if all=T, both criteria for all copulas are returned, otherwise, only the best copula object
  # together with its parameters is returned
  
  # possibilities are: Gumbel, Clayton, Frank, Joe, AMH, Gaussian and Student's t
  
  
  # dimensions
  d=length(U[1,])
  n=length(U[,1])
  
  # names and copulas
  names=c("Gumbel","Clayton","Frank","Joe","normal","t")
  cops=list(gumbelCopula(dim=d),
            claytonCopula(dim=d),
            frankCopula(dim=d),
            joeCopula(dim=d),
            normalCopula(dim=d,dispstr = "un"),
            tCopula(dim=d,dispstr = "un")
  )
  
  
  # fitting all possibilites and storing the AIC, BIC, parameters and log-likelihood
  AICfit=rep(NA,d)
  BICfit=rep(NA,d)
  pars=list()
  logl=rep(NA,d)
  for(i in 1:6){
    
    try({
      fit=fitCopula(copula = cops[[i]],data=U)
      AICfit[i]=-2*fit@loglik+2*length(fit@estimate)
      BICfit[i]=-2*fit@loglik+log(n)*length(fit@estimate)
      pars[[i]]=fit@estimate
      logl[i]=fit@loglik
    },silent=T)
    
    
  }
  
  if(all==T){
    
    names(pars)=names
    # returns copula names with corresponding AIC and BIC
    return(list("info"=cbind("name"=names,"logl"=logl,"AIC"=AICfit,"BIC"=BICfit),"parameters"=pars))
    
  } else {
    
    # returns the best copula object
    if(crit=="AIC"){
      
      return(list("copula"=names[which.min(AICfit)],
                  "parameters"=pars[[which.min(AICfit)]],
                  "logl"=logl[which.min(AICfit)]))
      
    } else if(crit=="BIC"){
      
      return(list("copula"=names[which.min(BICfit)],
                  "parameters"=pars[[which.min(BICfit)]],
                  "logl"=logl[which.min(BICfit)]))
      
    }
    
    
  }
  
}








## Construction of the contribution of the copula with 
## given margins to the CIC 
########################################################

CICpenalty=function(Z,margfuncs,alpha,mu,phi,df,copula,theta,nu){
  
  ######################################################################################
  ### function to calculate the penalty constant for the copula contribution to CIC  ###
  ### Z = nxd matrix containing raw data                                             ###
  ### margfuncs = d-vector containing fitted margin functions belonging              ###
  ###             to the QBA-family ("normal", "Laplace", "t" or "logistic")         ###
  ### alpha, mu and phi = d-vectors containing fitted skewing, location and          ###
  ###                     scale parameters                                           ###
  ### df = vector containing degrees of freedom of QBA-Student's t distributions     ###
  ### copula = the class of copula ("normal", "t", "gumbel", "clayton",              ###
  ###          "frank" or "joe")                                                     ###
  ### theta = vector containing the fitted copula parameters                         ###
  ### nu = degrees of freedom of Student's t-copula                                  ###
  ######################################################################################
  
  ### extras: dimensions and pseudo-observations ###
  ##################################################
  
  # lengths of vectors
  n=length(Z[,1])
  d=length(Z[1,])
  indt=which(margfuncs=='t')
  lt=length(indt)
  
  # degrees of freedom
  df=na.omit(df)
  if(lt==length(df)){
    df2=rep(NA,d)
    df2[indt]=df
    df=df2
  } else {
    stop("incorrect number of degrees of freedom provided")
  }
  
  # calculation of pseudo-observations
  U=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    U[,j]=switch(margfuncs[j],'normal'=QBAsyDist::pAND(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'logistic'=QBAsyDist::pALoD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                 ,'t'=QBAsyDist::pATD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                 ,'laplace'=QBAsyDist::pALaD(q = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
    )
  }
  
  
  
  ### calculation of the require derivatives w.r.t. the copula parameter and margin parameters ###
  ################################################################################################
  
  ### first derivatives of the log-margins w.r.t. the margin parameters ###
  
  partialalpha=matrix(NA,nrow=n,ncol=d)
  partialmu=matrix(NA,nrow=n,ncol=d)
  partialphi=matrix(NA,nrow=n,ncol=d)
  partialdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # derivatives of the log-margin with respect to alpha, mu and phi
    partialalpha[,j]=(1-2*alpha[j])/(alpha[j]*(1-alpha[j])) +
      arg*(1*(Z[,j]<=mu[j])*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
             1*(Z[,j]>mu[j])*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))
    
    partialmu[,j]=1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    
    partialphi[,j]= -1/phi[j] + 1*(Z[,j]<=mu[j])*(1-alpha[j])/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
      1*(Z[,j]>mu[j])*-alpha[j]/phi[j]*arg*ratiofirst(x =  alpha[j]*arg,basefunc = margfuncs[j],df = df[j])
    
    # when basefunc is student-t
    if(margfuncs[j]=="t"){
      partialdf[,j] = -1/(2*df[j]) + 1/2*digamma((df[j]+1)/2) - 1/2*digamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*(-1/2*log(1+(1-alpha[j])^2*arg^2/df[j])+(df[j]+1)/2*((1-alpha[j])^2*arg^2/df[j]^2)/(1+(1-alpha[j])^2*arg^2/df[j])) +
        1*(Z[,j]>mu[j])*(-1/2*log(1+alpha[j]^2*arg^2/df[j])+(df[j]+1)/2*(alpha[j]^2*arg^2/df[j]^2)/(1+alpha[j]^2*arg^2/df[j]))
    }
  }
  
  ### second derivatives of the log-margins w.r.t. the margin parameters
  partialalphaalpha=matrix(NA,nrow=n,ncol=d)
  partialalphamu=matrix(NA,nrow=n,ncol=d)
  partialalphaphi=matrix(NA,nrow=n,ncol=d)
  partialalphadf=matrix(NA,nrow=n,ncol=d)
  partialmumu=matrix(NA,nrow=n,ncol=d)
  partialmuphi=matrix(NA,nrow=n,ncol=d)
  partialmudf=matrix(NA,nrow=n,ncol=d)
  partialphiphi=matrix(NA,nrow=n,ncol=d)
  partialphidf=matrix(NA,nrow=n,ncol=d)
  partialdfdf=matrix(NA,nrow=n,ncol=d)
  
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    # second order derivatives of the log-margin with respect to alpha, mu and phi
    
    # alpha-alpha
    partialalphaalpha[,j]=(-2*alpha[j]^2+2*alpha[j]-1)/(alpha[j]^2*(1-alpha[j])^2)+
      1*(Z[,j]<=mu[j])*arg^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*arg^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                               (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    # alpha-mu
    partialalphamu[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                           arg*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                           (1-alpha[j])/phi[j]*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # alpha-phi
    partialalphaphi[,j]=1*(Z[,j]<=mu[j])*(-1/phi[j]*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                            arg^2*(1-alpha[j])/phi[j]*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                            (1-alpha[j])/phi[j]*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(-1/phi[j]*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         arg^2*alpha[j]/phi[j]*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         alpha[j]/phi[j]*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    # mu-mu
    partialmumu[,j]=1*(Z[,j]<=mu[j])*((1-alpha[j])/phi[j])^2*(ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j])-
                                                                (ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2)+
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j])^2*(ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j])-
                                             (ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2  )
    
    if(margfuncs[j]=='laplace'){
      # for the laplace this always gives zero so the value is replaced with its true expected value
      partialmumu[,j]=-alpha[j]*(1-alpha[j])/phi[j]^2
    }
    
    
    # mu-phi
    partialmuphi[,j]=1*(Z[,j]<=mu[j])*(-(1-alpha[j])/phi[j]^2*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                                         arg*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                                         (1-alpha[j])^2/phi[j]^2*arg*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(alpha[j]/phi[j]^2*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # phi-phi
    partialphiphi[,j]=1/phi[j]^2 + 
      1*(Z[,j]<=mu[j])*(-2*(1-alpha[j])/phi[j]^2*arg*ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) +
                          arg^2*(1-alpha[j])^2/phi[j]^2*ratiosecond(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]) -
                          (1-alpha[j])^2/phi[j]^2*arg^2*(ratiofirst(x = -(1-alpha[j])*arg,basefunc = margfuncs[j],df = df[j]))^2) +
      1*(Z[,j]>mu[j])*(2*alpha[j]/phi[j]^2*arg*ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) +
                         arg^2*alpha[j]^2/phi[j]^2*ratiosecond(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]) -
                         alpha[j]^2/phi[j]^2*arg^2*(ratiofirst(x = alpha[j]*arg,basefunc = margfuncs[j],df = df[j]))^2)
    
    
    # related to degrees of freedom (df) when basefunc is student-t 
    if(margfuncs[j]=='t'){
      
      L=1+(1-alpha[j])^2/df[j]*arg^2
      R=1+alpha[j]^2/df[j]*arg^2
      
      # alpha-df
      partialalphadf[,j]=1*(Z[,j]<=mu[j])*( (1-alpha[j])/df[j]^2*arg^2/L + (1-alpha[j])^3*(df[j]+1)/df[j]^3*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( alpha[j]/df[j]^2*arg^2/R - alpha[j]^3*(df[j]+1)/df[j]^3*arg^4/R^2 )
      
      # mu-df
      partialmudf[,j]=1*(Z[,j]<=mu[j])*( (-1/phi[j]+(df[j]+1)/(phi[j]*df[j]))*(1-alpha[j])^2/df[j]*arg/L + (1-alpha[j])^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/L^2 ) +
        1*(Z[,j]>mu[j])*( (1/phi[j]-(df[j]+1)/(phi[j]*df[j]))*alpha[j]^2/df[j]*arg/R + alpha[j]^4/df[j]^3*(df[j]+1)/phi[j]*arg^3/R^2 )
      
      # phi-df
      partialphidf[,j]=1*(Z[,j]<=mu[j])*( -(1-alpha[j])^2/(df[j]^2*phi[j])*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( -alpha[j]^2/(df[j]^2*phi[j])*arg^2/R + alpha[j]^4*(df[j]+1)/(df[j]^3*phi[j])*arg^4/R^2 )
      
      # df-df
      partialdfdf[,j]=1/(2*df[j]^2) + 1/4*trigamma((df[j]+1)/2) - 1/4*trigamma(df[j]/2) +
        1*(Z[,j]<=mu[j])*( ( (1-alpha[j])^2/df[j]^2 - (1-alpha[j])^2*(df[j]+1)/df[j]^3 )*arg^2/L + (1-alpha[j])^4*(df[j]+1)/(2*df[j]^4)*arg^4/L^2 ) +
        1*(Z[,j]>mu[j])*( ( alpha[j]^2/df[j]^2 - alpha[j]^2*(df[j]+1)/df[j]^3 )*arg^2/R + alpha[j]^4*(df[j]+1)/(2*df[j]^4)*arg^4/R^2 )
      
    }
    
  }
  
  
  
  ### derivatives of pseudo-obvervations U ( F(Z_j;eta_j) ) w.r.t the margin parameters
  
  # additional function for the student-t case as it requires the calculation of an integral
  derivIBeta=function(Z,nu){
    
    # function to numerically approximate the derivate of the incomplete beta function
    # approximation is in the integral term, other terms are calculated exactly
    
    # number of observations
    n=length(Z)
    
    # function passed down to integrate
    intfun=function(s,nu){
      return(log(s)*s^(nu/2-1)*(1-s)^(-1/2))
    }
    
    # exact part of the derivative
    gz=nu/(nu+Z)
    out=(1-gz)^(-1/2)*gz^(nu/2-1)*Z/(nu+Z)^2
    
    # calculation of the integral part per observation
    for(j in 1:n){
      x=1/2*integrate(f = intfun,lower = 0,upper = gz[j],nu=nu,subdivisions = 1000)$value
      out[j]=out[j]+x
    }
    
    return(out)
  }
  
  Ualpha=matrix(NA,nrow=n,ncol=d)
  Umu=matrix(NA,nrow=n,ncol=d)
  Uphi=matrix(NA,nrow=n,ncol=d)
  Udf=matrix(NA,nrow=n,ncol=d)
  for(j in 1:d){
    # argument without the alpha part
    arg=(Z[,j]-mu[j])/phi[j]
    
    temp=switch(margfuncs[j],
                'normal'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pnorm(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dnorm(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*pnorm(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dnorm(x = alpha[j]*arg))
                  , -QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dAND(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'logistic'=cbind(
                  1*(Z[,j]<=mu[j])*(2*plogis(q = (1-alpha[j])*arg)-2*alpha[j]*arg*dlogis(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*plogis(q = alpha[j]*arg)+2*(1-alpha[j])*arg*dlogis(x = alpha[j]*arg))
                  , -QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALoD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
                ,'t'=cbind(
                  1*(Z[,j]<=mu[j])*(2*pt(q = (1-alpha[j])*arg,df = df[j])-2*alpha[j]*arg*dt(x = (1-alpha[j])*arg,df = df[j])) +
                    1*(Z[,j]>mu[j])*(2-2*pt(q = alpha[j]*arg,df = df[j])+2*(1-alpha[j])*arg*dt(x = alpha[j]*arg,df = df[j]))
                  , -QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , -arg*QBAsyDist::dATD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j],nu = df[j])
                  , 1*(Z[,j]<=mu[j])*alpha[j]*( beta(df[j]/2,1/2)*derivIBeta(Z=((1-alpha[j])*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+((1-alpha[j])*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2 -
                    1*(Z[,j]>mu[j])*(1-alpha[j])*( beta(df[j]/2,1/2)*derivIBeta(Z=(alpha[j]*arg)^2,nu=df[j]) - 1/2*zipfR::Ibeta(x = df[j]/(df[j]+(alpha[j]*arg)^2),a = df[j]/2,b = 1/2)*beta(df[j]/2,1/2)*(digamma(df[j]/2)-digamma((df[j]+1)/2)) )/beta(df[j]/2,1/2)^2
                )
                ,'laplace'=cbind(
                  1*(Z[,j]<=mu[j])*(2*LaplacesDemon::plaplace(q = (1-alpha[j])*arg)-2*alpha[j]*arg*LaplacesDemon::dlaplace(x = (1-alpha[j])*arg)) +
                    1*(Z[,j]>mu[j])*(2-2*LaplacesDemon::plaplace(q = alpha[j]*arg)+2*(1-alpha[j])*arg*LaplacesDemon::dlaplace(x = alpha[j]*arg))
                  , -QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                  , -arg*QBAsyDist::dALaD(y = Z[,j],mu = mu[j],phi = phi[j],alpha = alpha[j])
                )
    )
    
    Ualpha[,j]=temp[,1]
    Umu[,j]=temp[,2]
    Uphi[,j]=temp[,3]
    if(margfuncs[j]=='t'){
      Udf[,j]=temp[,4]
    }
    
  }
  
  
  ### derivatives of the log-copula density w.r.t. all parameters ###
  ### three different cases: archimedean, normal or t             ###
  
  if(copula=="normal" | copula=="t"){
    
    if(copula=="normal"){
      nc=d*(d-1)/2
      nm=3*d+lt
      derivatives=derivsgaussian(theta = theta,U = U)
    } else {
      nc=d*(d-1)/2+1
      nm=3*d+lt
      derivatives=derivst(theta = theta,U = U,nu = nu)
    }
    
    
    # first derivative of the log-copula density w.r.t. the copula parameter
    # rows for observations
    # columns for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}  
    partialcop=derivatives$partialcop
    
    # second derivative of the log-copula density w.r.t. the copula parameter
    # first dimension is for observations
    # second and third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
    partialcopcop=derivatives$partialcopcop
    
    # mixed derivative of the log-copula density w.r.t. copula parameter and u_j
    # first dimension is for the observations
    # second dimension is for u_1,...,u_d
    # third dimension is for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
    partialcopu=derivatives$partialcopu
    
    # first derivative of the log-copula density w.r.t. the pseudo observations
    partialu=derivatives$partialu
    
    # observations which do not give NA in the derivatives
    ii=which(!is.na(apply(partialcop,1,sum)))
    
    
    ### first derivative of the log-copula density with respect to the margin parameters ###
    partialmarginpars=matrix(NA,nrow=length(ii),ncol=nm)
    counter=1
    for(j in 1:d){
      
      cc=2
      
      # Derivative of the log-copula wrt the argument u_j
      repeating=partialu[ii,j]
      
      
      # copula parameter-alpha
      partialcopalpha=repeating*Ualpha[ii,j]
      
      # copula parameter-mu
      partialcopmu=repeating*Umu[ii,j]
      
      # copula parameter-phi
      partialcopphi=repeating*Uphi[ii,j]
      
      # combining the previous three
      temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
      
      # when basefunc is student-t
      if(margfuncs[j]=='t'){
        cc=3
        
        # copula parameter-df
        partialcopdf=repeating*Udf[ii,j]
        temp=cbind(temp,partialcopdf)
      }
      
      
      
      partialmarginpars[,counter:(counter+cc)]=temp
      counter=counter+cc+1
    }
    
    
    # second derivative of the copula with respect to the copula parameters and margin parameters
    # first dimension for observations
    # second dimension for alpha_1,..,kappa_1,alpha_2,...,kappa_d
    # third dimension for Sigma_{2,1},...,Sigma_{d,1},Sigma_{3,2},....,Sigma_{d,d-1}
    partialcopeta=array(NA,dim =c(length(ii),nm,nc))
    for(k in 1:nc){
      counter=1
      for(j in 1:d){
        cc=2
        
        # this term is a common denominator in all partial derivatives, 
        # the derivate of the log-copula w.r.t u_j
        repeating=partialcopu[ii,j,k]
        
        
        # copula parameter-alpha
        partialcopalpha=repeating*Ualpha[ii,j]
        
        # copula parameter-mu
        partialcopmu=repeating*Umu[ii,j]
        
        # copula parameter-phi
        partialcopphi=repeating*Uphi[ii,j]
        
        # combining the previous three
        temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
        
        # when basefunc is student-t
        if(margfuncs[j]=='t'){
          cc=3
          
          # copula parameter-df
          partialcopdf=repeating*Udf[ii,j]
          temp=cbind(temp,partialcopdf)
        }
        
        
        
        partialcopeta[,counter:(counter+cc),k]=temp
        counter=counter+cc+1
      }
      
    }
    
    
    ### penalty term coming from the margin contribution ###
    PenaltyMargin=rep(NA,d)
    
    
    ### construction of I_eta ###
    
    # counter for creation of the full matrix
    counter=1
    
    I_eta=matrix(0,ncol=nm,nrow=nm)
    
    # minus average of second order derivatives
    paa=-1/length(ii)*apply(partialalphaalpha[ii,],2,sum)
    pam=-1/length(ii)*apply(partialalphamu[ii,],2,sum)
    pap=-1/length(ii)*apply(partialalphaphi[ii,],2,sum)
    pad=-1/length(ii)*apply(partialalphadf[ii,],2,sum)
    pmm=-1/length(ii)*apply(partialmumu[ii,],2,sum)
    pmp=-1/length(ii)*apply(partialmuphi[ii,],2,sum)
    pmd=-1/length(ii)*apply(partialmudf[ii,],2,sum)
    ppp=-1/length(ii)*apply(partialphiphi[ii,],2,sum)
    ppd=-1/length(ii)*apply(partialphidf[ii,],2,sum)
    pdd=-1/length(ii)*apply(partialdfdf[ii,],2,sum)
    
    # for putting all derivatives of the margins wrt their parameters next to eachother
    d_margins=c()
    
    counter=1
    for(j in 1:d){
      
      if(margfuncs[j]=='t'){
        
        cc=3
        
        # matrix I_eta_j
        I=matrix(NA,nrow=4,ncol=4)
        I[1,4]=I[4,1]=pad[j]
        I[2,4]=I[4,2]=pmd[j]
        I[3,4]=I[4,3]=ppd[j]
        I[4,4]=pdd[j]
        
        # for K^0_eta
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
        
      } else {
        
        cc=2
        
        # matrix I_eta_j
        I=matrix(NA,nrow=3,ncol=3)
        
        # for K^0_eta
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
        
      }
      
      I[1,1]=paa[j]
      I[1,2]=I[2,1]=pam[j]
      I[1,3]=I[3,1]=pap[j]
      I[2,2]=pmm[j]
      I[2,3]=I[3,2]=pmp[j]
      I[3,3]=ppp[j]
      
      # for the contribution of the copula to the penalty
      I_eta[counter:(counter+cc),counter:(counter+cc)]=I
      counter=counter+cc+1
      
      # for K^0_eta
      d_margins=cbind(d_margins,temp)
      
      # for the penalty
      K=1/length(ii)*t(temp)%*%temp
      PenaltyMargin[j]=sum(diag(solve(I)%*%K))
      
    }
    
    ### construction of K^0_eta ###
    K0eta=1/length(ii)*t(d_margins)%*%partialmarginpars
    
    
    ### construction of K_eta_thetaC ###
    
    # full matrix
    K_eta_thetaC=matrix(0,nrow=nm,ncol=nc)
    
    for(k in 1:nc){
      
      # counter for creation of the full matrix
      counter=1
      
      # block by block corresponding to the j-th margin
      for(j in 1:d){
        if(margfuncs[j]=='t'){
          cc=3
          temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
        } else {
          cc=2
          temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
        }
        K_eta_thetaC[counter:(counter+cc),k]=1/(length(ii))*t(temp)%*%partialcop[ii,k]
        counter=counter+cc+1
      }
      
    }
    
    ### construction of K_thetaC ###
    K_thetaC=1/(length(ii))*t(partialcop[ii,])%*%partialcop[ii,]
    
    
    ### construction of I_eta_thetaC ###
    if(d==2 & copula=="normal"){
      I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2),sum)
    } else {
      I_eta_thetaC=-1/(length(ii))*apply(partialcopeta[ii,,],c(2,3),sum)
    }
    
    ### construction of I_thetaC ###
    if(d==2 & copula=="normal"){
      I_thetaC=-1/(length(ii))*sum(partialcopcop[ii,,])
    } else {
      I_thetaC=-1/(length(ii))*apply(partialcopcop[ii,,],c(2,3),sum)
    }
    
  } else {
    
    # number of margin parameters
    nm=3*d+lt
    
    # derivatives of the log-copula density
    
    derivatives=switch(copula,"clayton"=derivsClayton(theta = theta,U = U)
                       ,"joe"=derivsJoe(theta = theta,U = U)
                       ,"frank"=derivsFrank(theta = theta,U = U)
                       ,"gumbel"=derivsGumbel(theta = theta,U = U)
                       ,"AMH"=derivsAMH(theta = theta,U = U)
    )
    
    # first derivative of the log-copula density w.r.t. the copula parameter
    partialcop=derivatives$partialcop
    
    # second derivative of the log-copula density w.r.t. the copula parameter
    partialcopcop=derivatives$partialcopcop
    
    # first derivative of the log-copula density w.r.t u_j (in the j-th column)
    partialu=derivatives$partialu
    
    # mixed derivative of the log-copula density w.r.t. copula parameter and u_j (in the j-th column)
    partialcopu=derivatives$partialcopu
    
    # observations which do not give NA in the derivatives
    ii=which(!is.na(partialcop))
    
    
    # first derivative of the log-copula density with restpect to the margin parameters
    partialmarginpars=matrix(NA,nrow=length(ii),ncol=nm)
    counter=1
    for(j in 1:d){
      
      cc=2
      
      # Derivative of the log-copula wrt the argument u_j
      repeating=partialu[ii,j]
      
      
      # copula parameter-alpha
      partialcopalpha=repeating*Ualpha[ii,j]
      
      # copula parameter-mu
      partialcopmu=repeating*Umu[ii,j]
      
      # copula parameter-phi
      partialcopphi=repeating*Uphi[ii,j]
      
      # combining the previous three
      temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
      
      # when basefunc is student-t
      if(margfuncs[j]=='t'){
        cc=3
        
        # copula parameter-df
        partialcopdf=repeating*Udf[ii,j]
        temp=cbind(temp,partialcopdf)
      }
      
      
      
      partialmarginpars[,counter:(counter+cc)]=temp
      counter=counter+cc+1
    }
    
    
    # mixed derivative of the log-copula with respect to the copula parameters and margin parameters
    partialcopeta=matrix(NA,nrow=length(ii),ncol=nm)
    counter=1
    for(j in 1:d){
      
      cc=2
      
      # mixed partial derivative wrt to the copula parameter and u_j of the log-copula density
      repeating=partialcopu[ii,j]
      
      
      # copula parameter-alpha
      partialcopalpha=repeating*Ualpha[ii,j]
      
      # copula parameter-mu
      partialcopmu=repeating*Umu[ii,j]
      
      # copula parameter-phi
      partialcopphi=repeating*Uphi[ii,j]
      
      # combining the previous three
      temp=cbind(partialcopalpha,partialcopmu,partialcopphi)
      
      # when basefunc is student-t
      if(margfuncs[j]=='t'){
        cc=3
        
        # copula parameter-df
        partialcopdf=repeating*Udf[ii,j]
        temp=cbind(temp,partialcopdf)
      }
      
      
      
      partialcopeta[,counter:(counter+cc)]=temp
      counter=counter+cc+1
    }
    
    
    ### penalty term coming from the margin contribution ###
    PenaltyMargin=rep(NA,d)
    
    ### construction of I_eta and K_eta ###
    
    I_eta=matrix(0,ncol=nm,nrow=nm)
    
    # minus average of second order derivatives
    paa=-1/length(ii)*apply(partialalphaalpha[ii,],2,sum)
    pam=-1/length(ii)*apply(partialalphamu[ii,],2,sum)
    pap=-1/length(ii)*apply(partialalphaphi[ii,],2,sum)
    pad=-1/length(ii)*apply(partialalphadf[ii,],2,sum)
    pmm=-1/length(ii)*apply(partialmumu[ii,],2,sum)
    pmp=-1/length(ii)*apply(partialmuphi[ii,],2,sum)
    pmd=-1/length(ii)*apply(partialmudf[ii,],2,sum)
    ppp=-1/length(ii)*apply(partialphiphi[ii,],2,sum)
    ppd=-1/length(ii)*apply(partialphidf[ii,],2,sum)
    pdd=-1/length(ii)*apply(partialdfdf[ii,],2,sum)
    
    # for putting all derivatives of the margins wrt their parameters next to eachother
    d_margins=c()
    
    counter=1
    for(j in 1:d){
      
      if(margfuncs[j]=='t'){
        
        cc=3
        
        # matrix I_eta_j
        I=matrix(NA,nrow=4,ncol=4)
        I[1,4]=I[4,1]=pad[j]
        I[2,4]=I[4,2]=pmd[j]
        I[3,4]=I[4,3]=ppd[j]
        I[4,4]=pdd[j]
        
        # for K^0_eta
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j],partialdf[ii,j])
        
      } else {
        
        cc=2
        
        # matrix I_eta_j
        I=matrix(NA,nrow=3,ncol=3)
        
        # for K^0_eta
        temp=cbind(partialalpha[ii,j],partialmu[ii,j],partialphi[ii,j])
        
      }
      
      I[1,1]=paa[j]
      I[1,2]=I[2,1]=pam[j]
      I[1,3]=I[3,1]=pap[j]
      I[2,2]=pmm[j]
      I[2,3]=I[3,2]=pmp[j]
      I[3,3]=ppp[j]
      
      # for the contribution of the copula to the penalty
      I_eta[counter:(counter+cc),counter:(counter+cc)]=I
      counter=counter+cc+1
      
      # for K^0_eta
      d_margins=cbind(d_margins,temp)
      
      # for the margin penalty
      K=1/length(ii)*t(temp)%*%temp
      PenaltyMargin[j]=sum(diag(solve(I)%*%K))
      
    }
    
    ### construction of K^0_eta ###
    K0eta=1/length(ii)*t(d_margins)%*%partialmarginpars
    
    ### construction of I_thetaC ###
    I_thetaC=-1/length(ii)*sum(partialcopcop[ii])
    
    ### construction of I_etathetaC ###
    I_eta_thetaC=-1/length(ii)*apply(partialcopeta[ii,],2,sum)
    
    ### construction of K_etathetaC ###
    K_eta_thetaC=1/length(ii)*t(d_margins)%*%partialcop[ii]
    
    ### construction of K_thetaC ###
    K_thetaC=1/length(ii)*t(partialcop[ii])%*%partialcop[ii]
  }  
  
  check=F
  try({
    solve(I_eta)
    check=T
  },silent=T)
  if(check==F){
    I_eta=I_eta+diag(10^-6,nrow = nrow(I_eta),ncol = ncol(I_eta))
  }
  
  ### penalty term coming from the copula contribution ###
  PenaltyCop=sum(diag(solve(I_eta)%*%K0eta))+sum(diag(-solve(I_thetaC)%*%t(I_eta_thetaC)%*%solve(I_eta)%*%K_eta_thetaC+solve(I_thetaC)%*%K_thetaC))
  
  return(list("copCIC"=PenaltyCop,"marginCIC"=PenaltyMargin))
  
}



