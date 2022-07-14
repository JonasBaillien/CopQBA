########################################################################
########################################################################
### Code accompagnying the paper:                                    ###
###     "Inference for copulas with two-piece margins"               ###
### Authors: Baillien Jonas, Gijbels irène and Verhasselt Anneleen   ###
########################################################################
########################################################################

##########################
### Part 2: Simulation ###
##########################

### source the functions file ###
#################################

source("~/PhD/andere/paper copulas/ES/code revisie/Code part 1 Functions.R")

### Important Notice ###
########################

# In certain sections, the model needs to be specified. All used models are in comments
# and the desired one can be uncommented.

# Output is immediately written away, for this the correct directory where to write needs
# to be specified. The same holds when loading in the output in later stages.





###############################################################
### script for plotting a copula and the corresponding data ###
###############################################################

### general information on the model

# seed
seed=14589
# number of monte carlo simulations
reps=500
# dimensions of the problem
d=2
# number of samples to be generated
sampsize=50
# number of starting points in the fitting of the marginals
nstart=10


### information on the margins

# type of univariate QB distribution (of length d)
margfuncs=c("t","normal")#,"logistic","laplace","t","normal","logistic","laplace","t","normal")
# alpha parameters (of length d)
alpha=c(0.2,0.7)#,0.3,0.8,0.6,0.3,0.3,0.5,0.7,0.2)
# mu parameters (of length d)
mu=c(2,3)#,-1,-1,4,-2,1,-2,-3,4)
# phi parameters (of length d)
phi=c(0.2,0.4)#,0.1,2,1,1,10,0.5,5,3)
# degrees of freedom (of length the number of "t" margins)
df=c(4)#,15,25)



### defining the copula object to simulate from
copulamodel=gumbelCopula(param = 2,dim=d) 


# generate a sample from the copula for the scale:
sample=samplecopula(n = sampsize,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed=seed)
X=sample$X
x11()
pairs(X)
x11()
pairs(sample$U)


### plot in case 2D

# data density plot
points=100

xax=seq(min(X[,1]),max(X[,1]),length.out=points)
yax=seq(min(X[,2]),max(X[,2]),length.out=points)
ax=cbind(xax,yax)


densityax=matrix(NA,nrow=points,ncol=d)
U=matrix(NA,nrow=points,ncol=d)

for(i in 1:d){
  if(margfuncs[i]=="normal"){
    U[,i]=pAND(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    densityax[,i]=dAND(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
  } else if(margfuncs[i]=="logistic"){
    U[,i]=pALoD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    densityax[,i]=dALoD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
  } else if(margfuncs[i]=="laplace"){
    U[,i]=pALaD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
    densityax[,i]=dALaD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i])
  } else {
    U[,i]=pATD(q = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i],nu=df[i])
    densityax[,i]=dATD(y = ax[,i],mu=mu[i],phi=phi[i],alpha=alpha[i],nu=df[i])
  }
}

Ugrid=as.matrix(expand.grid(U[,1],U[,2]))
distributioncopula=matrix(pCopula(Ugrid,copulamodel),nrow=points) # copula distribution on the [0,1]^2 grid
densitycopula=matrix(dCopula(Ugrid,copulamodel),nrow=points) # copula density on the [0,1]^2 grid

dens=matrix(NA,nrow=points,ncol=points) # density on the data scale
for(i in 1:points){
  for(j in 1:points){
    dens[i,j]=densityax[i,1]*densityax[j,2]*densitycopula[i,j]
  }
}
x11()
contour(xax,yax,dens)

# copula distribution plot
ax=seq(0.001,0.999,length.out=points) # a single axis with uniformly distributed points in [0,1]
gridax=as.matrix(expand.grid(ax,ax)) # create a grid on [0,1]^2
distributioncop=matrix(pCopula(gridax,copulamodel),nrow=points) # distribution value of the grid according to the copula
x11()
contour(ax,ax,distributioncop,col=1) # contourplot of copula distribution





















#############################################
### Simulation of a specific              ###
### copula with given margins in parallel ###
#############################################

### Important notice ###
# Uncomment the correct model from the options below 

# Manually change the directory into which the output is written.
# This also has to be done within the foreach loop itself and in the bundling of the output




### general information concerning the models (parameters etc.)

# # seed
# seed=14589
# # number of monte carlo simulations
# reps=1000
# # dimensions of the problem
# d=2 # 2 6 10
# # number of samples to be generated
# sampsize=75 # 75 150 300 600
# # number of starting points in the fitting of the marginals
# nstart=20
# 
# 
# 
# ### information on the margins
# 
# # type of univariate QB distribution (of length d)
# margfuncs=c("t","normal","logistic","laplace","t","normal","logistic","laplace","t","normal")[1:d]
# # alpha parameters (of length d)
# alpha=c(0.2,0.7,0.3,0.8,0.6,0.3,0.3,0.5,0.7,0.2)[1:d]
# # mu parameters (of length d)
# mu=c(2,3,-1,-1,4,-2,1,-2,-3,4)[1:d]
# # phi parameters (of length d)
# phi=c(0.2,0.4,0.1,2,1,1,10,0.5,5,3)[1:d]
# # degrees of freedom (of length the number of "t" margins)
# df=c(3,6,12)[1:round(d/3,0)]
# 
# ### copula models used (uncomment the desired one)
# 
# # ### Frank copula ###
# # ####################
# # 
# # copula parameters
# coppars=2
# 
# ### defining the copula object to simulate from
# copulamodel=frankCopula(param = coppars, dim = d)
# 
# 
# ### defining the copula for the model fit
# copulafit=frankCopula(dim = d)
# method="L-BFGS-B"


### Gaussian copula ###
#######################

# # copula parameters
# coppars=0.2
# 
# ### defining the copula object to simulate from
# copulamodel=normalCopula(param = coppars, dim = d, dispstr = "un")
# 
# 
# ### defining the copula for the model fit
# copulafit=normalCopula(dim = 2,dispstr="un")
# method="Nelder-Mead"
# 

# ### t copula ###
# ################
# 
# # copula parameters
# coppars=c(0.2, -0.05,  0.08, -0.2,  0.17,
#           0.18, -0.22,  0.03, -0.13,
#           -0.08,  0.3, -0.34,
#           0.01,  0.17,
#           -0.21)
# 
# ### defining the copula object to simulate from
# copulamodel=tCopula(param = coppars, dim = 6, dispstr = "un",df = 3)
# method="BFGS"
# 
# ### defining the copula for the model fit
# copulafit=tCopula(dim = d,dispstr = "un")


#################
### setting 2 ###
#################

# seed
seed=14589
# number of monte carlo simulations
reps=1000
# dimensions of the problem
d=3
# number of samples to be generated
sampsize=100 # 25 50 100 250 1000
# number of starting points in the fitting of the marginals
nstart=20

# type of univariate QB distribution (of length d)
margfuncs=c("laplace","normal","logistic")
# alpha parameters (of length d)
alpha=c(0.4,0.7,0.3)
# mu parameters (of length d)
mu=c(1,2,3)
# phi parameters (of length d)
phi=c(0.5,2.9,0.1)
# degrees of freedom (of length the number of "t" margins)
df=NULL
# 
### gumbel copula ###
#####################

# copula parameters
coppars=1.3

### defining the copula object to simulate from
copulamodel=gumbelCopula(param = coppars, dim = 3)


### defining the copula for the model fit
copulafit=gumbelCopula(dim = 3)
method="L-BFGS-B"



### file with the simulation setting for replicability
settings=list("seed"=seed,"d"=d,"sampsize"=sampsize,"nstart"=nstart,"margfuncs"=margfuncs,"alpha"=alpha,"mu"=mu,"phi"=phi,"df"=df,"coppars"=coppars,"copula"=copulamodel)


### simulating data and fitting the model to it
cores=detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)


### create folder to store output in
dir.create(paste0("~/PhD/code/copulas/output/normalcop_dim",d,"_S",sampsize))

# voor probleemgevallen normale copula
# library(readr)
# niet <- unlist(read_csv("~/PhD/code/copulas/output/normalcop_dim2_S500/Welke niet erbij.txt"))
# method="BFGS"


### main loop for drawing sample from the model and refitting it to the sample
result=foreach(i=niet,.packages=c('QBAsyDist','copula','nloptr'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 ### required files:
                 source("Code part 1 Functions.R")
                 
                 # generate a sample from the copula:
                 sample=samplecopula(n = sampsize,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed=seed+i)
                 X=sample$X
                 
                 # timekeeping
                 start.time <- Sys.time()
                 
                 # fit the copula model to the generated sample
                 fit=fitcopula(data = X,margfuncs = margfuncs,nstart = nstart, cop = copulafit,seed = seed+i,method = method)
                 
                 # more timekeeping
                 end.time <- Sys.time()
                 time.taken = difftime(end.time,start.time,units = "secs")
                 
                 # write output
                 fit[["time"]]=time.taken
                 filename=paste0("~/PhD/code/copulas/output/normalcop_dim",d,"_S",sampsize,"/run",i,".Rdata")
                 save(fit,file=filename)
                 
               }
stopCluster(cl)

# Bundling the output
dataset=list()
ind=c()
for (i in 1:reps){
  try({
    load(paste0("~/PhD/code/copulas/output/normalcop_dim",d,"_S",sampsize,"/run",i,".Rdata"))
    dataset[[i]] <- fit
    ind=c(ind,i)
  },silent=T)
}
save(dataset,file=paste0("~/PhD/code/copulas/output/normalcop_dim",d,"_S",sampsize,".Rdata"))
save(settings,file=paste0("~/PhD/code/copulas/output/normalcop_dim",d,"_S",sampsize,"settings.Rdata"))






### plots for simulation setting 1 ###
######################################


### Gaussian copula ###
#######################

# dimension
d=2
# number of monte carlo simulations
reps=1000

# type of univariate QBA distribution (of length d)
margfuncs=c("t","normal")
# alpha parameters (of length d)
alpha=c(0.2,0.7)
# mu parameters (of length d)
mu=c(2,3)
# phi parameters (of length d)
phi=c(0.2,0.4)
# degrees of freedom for t distribution
df=c(4)
# copula parameter (rho)
coppars=c(0.2)


# different sample sizes
ss=c(75,150,300,600)
# length of copula parameter vector
cp=1

# holding matrices for the MSE for the different parameters
MSEa=matrix(NA,nrow=d,ncol=4)
MSEm=matrix(NA,nrow=d,ncol=4)
MSEp=matrix(NA,nrow=d,ncol=4)
MSEd=matrix(NA,nrow=round(d/3,0),ncol=4)
MSEc=matrix(NA,nrow=cp,ncol=4)

# run over the sample size
for(j in 1:4){
  
  sampsize=ss[j]
  
  
  ### load the fits
  load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/normal/normalcop_dim",d,"_S",sampsize,".Rdata"))
  fits=dataset
  
  
  ### create matrices with parameters grouped together
  alphaf=matrix(NA,nrow=reps,ncol=d)
  muf=matrix(NA,nrow=reps,ncol=d)
  phif=matrix(NA,nrow=reps,ncol=d)
  dff=matrix(NA,nrow=reps,ncol=d)
  loglf=rep(NA,reps)
  copparsf=matrix(NA,nrow=reps,ncol=cp)
  time.taken=rep(NA,reps)
  for(i in 1:reps){
    temp_dat=fits[[i]]
    alphaf[i,]=temp_dat$alpha
    muf[i,]=temp_dat$mu
    phif[i,]=temp_dat$phi
    dff[i,]=temp_dat$df
    loglf[i]=temp_dat$logl
    copparsf[i,]=temp_dat$copulapar
    time.taken[i]=temp_dat$time
  }
  
  
  ### bias results of estimates parameters
  biasalpha=apply(alphaf,2,mean)-alpha
  biasmu=apply(muf,2,mean)-mu
  biasphi=apply(phif,2,mean)-phi
  biasdf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,mean)-df
  biascoppars=apply(copparsf,2,mean)-coppars
  
  ### empirical variance of estimated parameters
  varalpha=apply(alphaf,2,var)
  varmu=apply(muf,2,var)
  varphi=apply(phif,2,var)
  vardf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,var)
  varcoppars=apply(copparsf,2,var)
  
  ### MSE
  MSEalpha=biasalpha^2+varalpha
  MSEmu=biasmu^2+varmu
  MSEphi=biasphi^2+varphi
  MSEdf=biasdf^2+vardf
  MSEcoppars=biascoppars^2+varcoppars
  
  # fill in holding matrices when calculating corresponding values per sample size
  MSEa[,j]=MSEalpha
  MSEm[,j]=MSEmu
  MSEp[,j]=MSEphi
  MSEd[,j]=MSEdf
  MSEc[,j]=MSEcoppars
  
}

# plotting the MSE
# for alpha
colnames(MSEa)=ss
Ma=melt(data = MSEa,varnames = c("d","n"))
x11()
ggplot(data = Ma,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("alpha",1),paste0("alpha",2)),values=c("black","red")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18)) # oorspronkelijk 12 14 14

# for mu
colnames(MSEm)=ss
Mm=melt(data = MSEm,varnames = c("d","n"))
x11()
ggplot(data = Mm,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("mu",1),paste0("mu",2)),values=c("black","red")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for phi
colnames(MSEp)=ss
Mp=melt(data = MSEp,varnames = c("d","n"))
x11()
ggplot(data = Mp,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("phi",1),paste0("phi",2)),values=c("black","red")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for degrees of freedom
colnames(MSEd)=ss
Md=melt(data = MSEd,varnames = c("d","n"))
x11()
ggplot(data = Md,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("nu",1)),values=c("black")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for copula parameters
colnames(MSEc)=ss
Mc=melt(data = MSEc,varnames = c("d","n"))
x11()
ggplot(data = Mc,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("Sigma",1,2)),values=c("black")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))


# drawing sample
S=samplecopula(n = 100000,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed = seed)
Z=S[[2]]
# calculating the approximate FI-matrix
mats=FI_IFM_Gaussian(Z=Z,alpha=alpha,mu=mu,phi=phi,df=df,coppars=coppars,margfuncs=margfuncs)
FIM=mats$FIM
vars=solve(mats$FIM)

vcm=matrix(0,7,7)
vars2=matrix(NA,2,4)

temp=asymptoticvariance(alpha = alpha[1],mu = mu[1],phi = phi[1],nu = NULL,margfunc = margfuncs[1])
vars2[1,]=temp$vars
vcm[1:4,1:4]=temp$cov
temp=asymptoticvariance(alpha = alpha[2],mu = mu[2],phi = phi[2],nu = NULL,margfunc = margfuncs[2])
vars2[2,]=temp$vars
vcm[5:7,5:7]=temp$cov

vars[1:7,1:7]=vcm

vars=vars/100
round(vars,4)


### Student's t copula ###
##########################

# dimensions
d=6
# number of monte carlo simulations
reps=1000


# type of univariate QB distribution (of length d)
margfuncs=c("t","normal","logistic","laplace","t","normal")
# alpha parameters (of length d)
alpha=c(0.2,0.7,0.3,0.8,0.6,0.3)
# mu parameters (of length d)
mu=c(2,3,-1,-1,4,-2)
# phi parameters (of length d)
phi=c(0.2,0.4,0.1,2,1,1)
# degrees of freedom for t distribution
df=c(4,15)
# copula parameter vector (Sigma and nu)
coppars=c(0.2, -0.05,  0.08, -0.2,  0.17,
          0.18, -0.22,  0.03, -0.13,
          -0.08,  0.3, -0.34,
          0.01,  0.17,
          -0.21)
coppars=c(coppars,3)

# sample sizes considered
ss=c(75,150,300,600)
# length of copula parameter vector
cp=16

# holding matrices for MSE
MSEa=matrix(NA,nrow=d,ncol=4)
MSEm=matrix(NA,nrow=d,ncol=4)
MSEp=matrix(NA,nrow=d,ncol=4)
MSEd=matrix(NA,nrow=round(d/3,0),ncol=4)
MSEc=matrix(NA,nrow=cp,ncol=4)

# run over sample sizes
for(j in 1:4){
  
  sampsize=ss[j]
  
  
  ### load the fits
  load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/t/tcop_dim",d,"_S",sampsize,".Rdata"))
  fits=dataset
  
  
  ### create matrices with parameters grouped together
  alphaf=matrix(NA,nrow=reps,ncol=d)
  muf=matrix(NA,nrow=reps,ncol=d)
  phif=matrix(NA,nrow=reps,ncol=d)
  dff=matrix(NA,nrow=reps,ncol=d)
  loglf=rep(NA,reps)
  copparsf=matrix(NA,nrow=reps,ncol=cp)
  time.taken=rep(NA,reps)
  for(i in 1:reps){
    temp_dat=fits[[i]]
    alphaf[i,]=temp_dat$alpha
    muf[i,]=temp_dat$mu
    phif[i,]=temp_dat$phi
    dff[i,]=temp_dat$df
    loglf[i]=temp_dat$logl
    copparsf[i,]=temp_dat$copulapar
    time.taken[i]=temp_dat$time
  }
  
  
  ### bias results
  biasalpha=apply(alphaf,2,mean)-alpha
  biasmu=apply(muf,2,mean)-mu
  biasphi=apply(phif,2,mean)-phi
  biasdf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,mean)-df
  biascoppars=apply(copparsf,2,mean)-coppars
  
  ### empirical variance of parameters
  varalpha=apply(alphaf,2,var)
  varmu=apply(muf,2,var)
  varphi=apply(phif,2,var)
  vardf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,var)
  varcoppars=apply(copparsf,2,var)
  
  ### MSE
  MSEalpha=biasalpha^2+varalpha
  MSEmu=biasmu^2+varmu
  MSEphi=biasphi^2+varphi
  MSEdf=biasdf^2+vardf
  MSEcoppars=biascoppars^2+varcoppars
  
  # fill in when calculating corresponding values per sample size
  MSEa[,j]=MSEalpha
  MSEm[,j]=MSEmu
  MSEp[,j]=MSEphi
  MSEd[,j]=MSEdf
  MSEc[,j]=MSEcoppars
  
}

# plot of MSE
# for alpha
colnames(MSEa)=ss
Ma=melt(data = MSEa,varnames = c("d","n"))
x11()
ggplot(data = Ma,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(pch=19,size=2) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("alpha",1),paste0("alpha",2),
                                      paste0("alpha",3),paste0("alpha",4),
                                      paste0("alpha",5),paste0("alpha",6))
                     ,values=c("black","red","green","blue","cyan","magenta")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for mu
colnames(MSEm)=ss
Mm=melt(data = MSEm,varnames = c("d","n"))
x11()
ggplot(data = Mm,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("mu",1),paste0("mu",2),
                                      paste0("mu",3),paste0("mu",4),
                                      paste0("mu",5),paste0("mu",6)),
                     values=c("black","red","green","blue","cyan","magenta")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for phi
colnames(MSEp)=ss
Mp=melt(data = MSEp,varnames = c("d","n"))
x11()
ggplot(data = Mp,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("phi",1),paste0("phi",2),
                                      paste0("phi",3),paste0("phi",4),
                                      paste0("phi",5),paste0("phi",6)),
                     values=c("black","red","green","blue","cyan","magenta")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for degrees of freedom (both margin and copula parameters)
colnames(MSEd)=ss
Md=melt(data = rbind(MSEd,MSEc[16,]),varnames = c("d","n"))
x11()
ggplot(data = Md,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("nu",1),paste0("nu",5),
                                      paste0("nu ","copula")),values=c("black","cyan","orange")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for copula parameters (Sigma)
colnames(MSEc)=ss
Mc=melt(data = MSEc[1:5,],varnames = c("d","n"))
x11()
ggplot(data = Mc,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("Sigma",1,2),paste0("Sigma",1,3),
                                      paste0("Sigma",1,4),paste0("Sigma",1,5),
                                      paste0("Sigma",1,6)),values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

colnames(MSEc)=ss
Mc=melt(data = MSEc[6:10,],varnames = c("d","n"))
x11()
ggplot(data = Mc,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("Sigma",2,3),paste0("Sigma",2,4),
                                      paste0("Sigma",2,5),paste0("Sigma",2,6),
                                      paste0("Sigma",3,4)),values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

colnames(MSEc)=ss
Mc=melt(data = MSEc[11:15,],varnames = c("d","n"))
x11()
ggplot(data = Mc,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("Sigma",3,5),paste0("Sigma",3,6),
                                      paste0("Sigma",4,5),paste0("Sigma",4,6),
                                      paste0("Sigma",5,6)),values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))



### frank copula ###
####################

# dimension
d=10
# number of monte carlo simulations
reps=1000


# type of univariate QB distribution (of length d)
margfuncs=c("t","normal","logistic","laplace","t","normal","logistic","laplace","t","normal")
# alpha parameters (of length d)
alpha=c(0.2,0.7,0.3,0.8,0.6,0.3,0.3,0.5,0.7,0.2)
# mu parameters (of length d)
mu=c(2,3,-1,-1,4,-2,1,-2,-3,4)
# phi parameters (of length d)
phi=c(0.2,0.4,0.1,2,1,1,10,0.5,5,3)
# degrees of freedom for t distribution
df=c(4,15,25)
# theta
coppars=2

# sample sizes considered
ss=c(75,150,300,600)
# length of copula parameter vector
cp=1

# holding matrices for MSE
MSEa=matrix(NA,nrow=d,ncol=4)
MSEm=matrix(NA,nrow=d,ncol=4)
MSEp=matrix(NA,nrow=d,ncol=4)
MSEd=matrix(NA,nrow=round(d/3,0),ncol=4)
MSEc=matrix(NA,nrow=cp,ncol=4)

# run over sample size
for(j in 1:4){
  
  sampsize=ss[j]
  
  
  ### load the fits
  load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/frank/frankcop_dim",d,"_S",sampsize,".Rdata"))
  fits=dataset
  
  
  ### create matrices with parameters grouped together
  alphaf=matrix(NA,nrow=reps,ncol=d)
  muf=matrix(NA,nrow=reps,ncol=d)
  phif=matrix(NA,nrow=reps,ncol=d)
  dff=matrix(NA,nrow=reps,ncol=d)
  loglf=rep(NA,reps)
  copparsf=matrix(NA,nrow=reps,ncol=cp)
  time.taken=rep(NA,reps)
  for(i in 1:reps){
    temp_dat=fits[[i]]
    alphaf[i,]=temp_dat$alpha
    muf[i,]=temp_dat$mu
    phif[i,]=temp_dat$phi
    dff[i,]=temp_dat$df
    loglf[i]=temp_dat$logl
    copparsf[i,]=temp_dat$copulapar
    time.taken[i]=temp_dat$time
  }
  
  
  ### bias results
  biasalpha=apply(alphaf,2,mean)-alpha
  biasmu=apply(muf,2,mean)-mu
  biasphi=apply(phif,2,mean)-phi
  biasdf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,mean)-df
  biascoppars=apply(copparsf,2,mean)-coppars
  
  ### empirical variance of parameters
  varalpha=apply(alphaf,2,var)
  varmu=apply(muf,2,var)
  varphi=apply(phif,2,var)
  vardf=apply(matrix(dff[,which(margfuncs=="t")],ncol=length(which(margfuncs=="t"))),2,var)
  varcoppars=apply(copparsf,2,var)
  
  ### MSE
  MSEalpha=biasalpha^2+varalpha
  MSEmu=biasmu^2+varmu
  MSEphi=biasphi^2+varphi
  MSEdf=biasdf^2+vardf
  MSEcoppars=biascoppars^2+varcoppars
  
  # fill in when calculating corresponding values per sample size
  MSEa[,j]=MSEalpha
  MSEm[,j]=MSEmu
  MSEp[,j]=MSEphi
  MSEd[,j]=MSEdf
  MSEc[,j]=MSEcoppars
  
}

# plot of MSE
# for alpha
colnames(MSEa)=ss
Ma=melt(data = MSEa[1:5,],varnames = c("d","n"))
x11()
ggplot(data = Ma,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(pch=19,size=2) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("alpha",1),paste0("alpha",2),
                                      paste0("alpha",3),paste0("alpha",4),
                                      paste0("alpha",5))
                     ,values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

Ma=melt(data = MSEa[6:10,],varnames = c("d","n"))
x11()
ggplot(data = Ma,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(pch=19,size=2) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("alpha",6),paste0("alpha",7),
                                      paste0("alpha",8),paste0("alpha",9),
                                      paste0("alpha",10))
                     ,values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for mu
colnames(MSEm)=ss
Mm=melt(data = MSEm[1:5,],varnames = c("d","n"))
x11()
ggplot(data = Mm,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("mu",1),paste0("mu",2),
                                      paste0("mu",3),paste0("mu",4),
                                      paste0("mu",5),paste0("mu",6)),
                     values=c("black","red","green","blue","cyan","magenta")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

Mm=melt(data = MSEm[6:10,],varnames = c("d","n"))
x11()
ggplot(data = Mm,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("mu",6),paste0("mu",7),
                                      paste0("mu",8),paste0("mu",9),
                                      paste0("mu",10)),
                     values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for phi
colnames(MSEp)=ss
Mp=melt(data = MSEp[1:5,],varnames = c("d","n"))
x11()
ggplot(data = Mp,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("phi",1),paste0("phi",2),
                                      paste0("phi",3),paste0("phi",4),
                                      paste0("phi",5)),
                     values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

Mp=melt(data = MSEp[6:10,],varnames = c("d","n"))
x11()
ggplot(data = Mp,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) +
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("phi",6),paste0("phi",7),
                                      paste0("phi",8),paste0("phi",9),
                                      paste0("phi",10)),
                     values=c("black","red","green","blue","cyan")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for degrees of freedom
colnames(MSEd)=ss
Md=melt(data = MSEd,varnames = c("d","n"))
x11()
ggplot(data = Md,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("nu",1),paste0("nu",5),
                                      paste0("nu",9)),values=c("black","cyan","orange")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

# for copula parameter
colnames(MSEc)=ss
Mc=melt(data = MSEc,varnames = c("d","n"))
x11()
ggplot(data = Mc,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c(paste0("theta")),values=c("black")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))

### coparison over d=2,6,10 for frank
reps=1000 # number of monte carlo simulations
ds=c(2,6,10) # dimensions
ss=c(75,150,300,600) # sample sizes
cp=1 # length of copula parameter vector
MSEf=matrix(NA,nrow=3,ncol=4) # holding matrix

# run over dimensions
for(k in 1:3){
  
  d=ds[k]
  MSEc=matrix(NA,nrow=cp,ncol=4)
  
  # type of univariate QB distribution (of length d)
  margfuncs=c("t","normal","logistic","laplace","t","normal","logistic","laplace","t","normal")[1:d]
  # alpha parameters (of length d)
  alpha=c(0.2,0.7,0.3,0.8,0.6,0.3,0.3,0.5,0.7,0.2)[1:d]
  # mu parameters (of length d)
  mu=c(2,3,-1,-1,4,-2,1,-2,-3,4)[1:d]
  # phi parameters (of length d)
  phi=c(0.2,0.4,0.1,2,1,1,10,0.5,5,3)[1:d]
  # degrees of freedom (of length the number of "t" margins)
  df=c(4,15,25)[1:round(d/3,0)]
  # copula parameters
  coppars=2
  
  # run over sample sizes
  for(j in 1:4){
    
    sampsize=ss[j]
    
    
    ### load the fits
    load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/frank/frankcop_dim",d,"_S",sampsize,".Rdata"))
    fits=dataset
    
    
    ### create matrices with parameters grouped together
    copparsf=matrix(NA,nrow=reps,ncol=cp)
    for(i in 1:reps){
      temp_dat=fits[[i]]
      copparsf[i,]=temp_dat$copulapar
    }
    
    ### bias results
    biascoppars=apply(copparsf,2,mean)-coppars
    
    ### empirical variance of parameters
    varcoppars=apply(copparsf,2,var)
    
    ### MSE
    MSEcoppars=biascoppars^2+varcoppars
    
    # fill in when calculating corresponding values per sample size
    MSEf[k,j]=MSEcoppars
    
  }
}

# plot
colnames(MSEf)=ss
Mf=melt(data = MSEf,varnames = c("d","n"))
x11()
ggplot(data = Mf,aes(x=n,y=value,group=factor(d),color=factor(d))) +
  geom_line(lwd=1.2) + 
  geom_point(size=2,pch=19) +
  xlab("n") +
  ylab("AMSE") +
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="",labels=c("d=2","d=6","d=10"),values=c("black","cyan","orange")) +
  scale_x_continuous(labels=c("75","150","300","600"),breaks = c(75,150,300,600)) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(size=18))




### simulation setting 2 boxplots ###
#####################################

### model

# number of monte carlo simulations
reps=1000
# dimensions of the problem
d=3
# number of samples to be generated
sampsize=25  #25 50 100 250 1000
# number of starting points in the fitting of the marginals
nstart=10



### information on the margins

# type of univariate QB distribution (of length d)
margfuncs=c("laplace","normal","logistic")
# alpha parameters (of length d)
alpha=c(0.4,0.7,0.3)
# mu parameters (of length d)
mu=c(1,2,3)
# phi parameters (of length d)
phi=c(0.5,2.9,0.1)
# degrees of freedom (of length the number of "t" margins)
df=NULL
# copula parameters
coppars=1.3

### defining the copula object to simulate from
copulamodel=gumbelCopula(param = coppars, dim = d) 


### data of fits

load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",25,".Rdata"))
fits1=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",50,".Rdata"))
fits2=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",100,".Rdata"))
fits3=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",250,".Rdata"))
fits4=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",1000,".Rdata"))
fits5=dataset

# grouping parameters
alphaf=array(data=NA,dim=c(reps,d,5))
muf=array(data=NA,dim=c(reps,d,5))
phif=array(data=NA,dim=c(reps,d,5))
dff=array(data=NA,dim=c(reps,d,5))
loglf=matrix(NA,nrow=reps,ncol=5)
copparsf=matrix(data=NA,nrow=reps,ncol=5)
time.taken=matrix(NA,nrow=reps,ncol=5)
for(i in 1:reps){
  temp_dat=fits1[[i]]
  alphaf[i,,1]=temp_dat$alpha
  muf[i,,1]=temp_dat$mu
  phif[i,,1]=temp_dat$phi
  dff[i,,1]=temp_dat$df
  loglf[i,1]=temp_dat$logl
  copparsf[i,1]=temp_dat$copulapar
  time.taken[i,1]=temp_dat$time
  
  temp_dat=fits2[[i]]
  alphaf[i,,2]=temp_dat$alpha
  muf[i,,2]=temp_dat$mu
  phif[i,,2]=temp_dat$phi
  dff[i,,2]=temp_dat$df
  loglf[i,2]=temp_dat$logl
  copparsf[i,2]=temp_dat$copulapar
  time.taken[i,2]=temp_dat$time
  
  temp_dat=fits3[[i]]
  alphaf[i,,3]=temp_dat$alpha
  muf[i,,3]=temp_dat$mu
  phif[i,,3]=temp_dat$phi
  dff[i,,3]=temp_dat$df
  loglf[i,3]=temp_dat$logl
  copparsf[i,3]=temp_dat$copulapar
  time.taken[i,3]=temp_dat$time
  
  temp_dat=fits4[[i]]
  alphaf[i,,4]=temp_dat$alpha
  muf[i,,4]=temp_dat$mu
  phif[i,,4]=temp_dat$phi
  dff[i,,4]=temp_dat$df
  loglf[i,4]=temp_dat$logl
  copparsf[i,4]=temp_dat$copulapar
  time.taken[i,4]=temp_dat$time
  
  temp_dat=fits5[[i]]
  alphaf[i,,5]=temp_dat$alpha
  muf[i,,5]=temp_dat$mu
  phif[i,,5]=temp_dat$phi
  dff[i,,5]=temp_dat$df
  loglf[i,5]=temp_dat$logl
  copparsf[i,5]=temp_dat$copulapar
  time.taken[i,5]=temp_dat$time
}

# bias
ba=matrix(NA,5,3)
bm=matrix(NA,5,3)
bp=matrix(NA,5,3)
bc=matrix(NA,5,1)
for(j in 1:5){
  ba[j,]=colMeans(alphaf[,,j])-alpha
  bm[j,]=colMeans(muf[,,j])-mu
  bp[j,]=colMeans(phif[,,j])-phi
  bc[j]=mean(copparsf[,j])-coppars
}
round(cbind(ba[,1],bm[,1],bp[,1],ba[,2],bm[,2],bp[,2],ba[,3],bm[,3],bp[,3],bc),4)



### theoretical variance of parameter estimates
# seed
seed=1234
# number of observations to approximate the theoretical variance
n=100000

# drawing sample
S=samplecopula(n = n,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed = seed)
Z=S[[2]]
# calculating the approximate FI-matrix
mats=FI_IFM(Z=Z,alpha=alpha,mu=mu,phi=phi,df=NULL,coppars=coppars,margfuncs=margfuncs,copula="gumbel")
FIM=mats$FIM
vars=solve(mats$FIM)

# replace the margin part with a more accurate approximation. This is theoretically justified
vcm=matrix(0,9,9)
vars2=matrix(NA,3,3)
for(i in 1:3){
  temp=asymptoticvariance(alpha = alpha[i],mu = mu[i],phi = phi[i],nu = NULL,margfunc = margfuncs[i])
  vars2[i,]=temp$vars
  vcm[(3*i-2):(3*i),(3*i-2):(3*i)]=temp$cov
}

vars[1:9,1:9]=vcm

library(xtable)
a=xtable(x = vars,digits = 4)
print(a,file="~/PhD/code/copulas/output/covgumbel.txt")

# holding matrices for approximate theoretical and empirical variances and their ratio
thvars=matrix(NA,nrow=5,ncol=10)
empvars=matrix(NA,nrow=5,ncol=10)
ratio=matrix(NA,nrow=5,ncol=10)
# sample sizes considered
ss=c(25,50,100,250,1000)

for(j in 1:5){
  # thvars[j,]=c(vars2[1,],vars2[2,],vars2[3,],vars[10,10])/ss[j]
  # thvars[j,]=diag(vars)/ss[j]
  empvars[j,c(1,4,7)]=apply(alphaf[,,j],2,var)
  empvars[j,c(2,5,8)]=apply(muf[,,j],2,var)
  empvars[j,c(3,6,9)]=apply(phif[,,j],2,var)
  empvars[j,10]=var(copparsf[,j])
  # ratio[j,]=empvars[j,]/thvars[j,]
}

round(empvars,4)
round(ratio,4)


### boxplots test sample size
### alpha
alphaf.d1.df=data.frame(alphaf[,1,])
colnames(alphaf.d1.df)=c("25","50","100","250","1000")
alphaf.d1.df=melt(alphaf.d1.df)
x11()
ggplot(alphaf.d1.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = alpha[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

alphaf.d2.df=data.frame(alphaf[,2,])
colnames(alphaf.d2.df)=c("25","50","100","250","1000")
alphaf.d2.df=melt(alphaf.d2.df)
x11()
ggplot(alphaf.d2.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = alpha[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

alphaf.d3.df=data.frame(alphaf[,3,])
colnames(alphaf.d3.df)=c("25","50","100","250","1000")
alphaf.d3.df=melt(alphaf.d3.df)
x11()
ggplot(alphaf.d3.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = alpha[3],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

### mu
muf.d1.df=data.frame(muf[,1,])
colnames(muf.d1.df)=c("25","50","100","250","1000")
muf.d1.df=melt(muf.d1.df)
x11()
ggplot(muf.d1.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = mu[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

muf.d2.df=data.frame(muf[,2,])
colnames(muf.d2.df)=c("25","50","100","250","1000")
muf.d2.df=melt(muf.d2.df)
x11()
ggplot(muf.d2.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = mu[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

muf.d3.df=data.frame(muf[,3,])
colnames(muf.d3.df)=c("25","50","100","250","1000")
muf.d3.df=melt(muf.d3.df)
x11()
ggplot(muf.d3.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = mu[3],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

### phi
phif.d1.df=data.frame(phif[,1,])
colnames(phif.d1.df)=c("25","50","100","250","1000")
phif.d1.df=melt(phif.d1.df)
x11()
ggplot(phif.d1.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = phi[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

phif.d2.df=data.frame(phif[,2,])
colnames(phif.d2.df)=c("25","50","100","250","1000")
phif.d2.df=melt(phif.d2.df)
x11()
ggplot(phif.d2.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = phi[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

phif.d3.df=data.frame(phif[,3,])
colnames(phif.d3.df)=c("25","50","100","250","1000")
phif.d3.df=melt(phif.d3.df)
x11()
ggplot(phif.d3.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = phi[3],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

### copula parameter
copparsf.df=data.frame(copparsf)
colnames(copparsf.df)=c("25","50","100","250","1000")
copparsf.df=melt(copparsf.df)
x11()
ggplot(copparsf.df, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("sample size") + 
  ylab("") +
  geom_abline(slope = 0,intercept = coppars[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))



### overlay with the theoretical distribution ss=50 and 250
# alpha 1
a11.m <- melt(alphaf[,1,2])
xa1=c(0.1,0.7)
xa11=seq(xa1[1],xa1[2],by=0.0005)
ya11=dnorm(x = xa11,mean = alpha[1],sd = sqrt(thvars[2,1]))
linea11=data.frame(xa11,ya11)
a11 <- ggplot(a11.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xa1[1],xa1[2]) +
  ylim(0,13) +
  geom_line(data = linea11, aes(x=xa11,y=ya11), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

a12.m <- melt(alphaf[,1,4])
xa12=seq(xa1[1],xa1[2],by=0.001)
ya12=dnorm(x = xa12,mean = alpha[1],sd = sqrt(thvars[4,1]))
linea12=data.frame(xa12,ya12)
a12 <- ggplot(a12.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xa1[1],xa1[2]) +
  ylim(0,13) +
  geom_line(data = linea12, aes(x=xa12,y=ya12), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(a11, a12, nrow = 1)

# mu 1
m11.m <- melt(muf[,1,2])
xm1=c(0,2)
xm11=seq(xm1[1],xm1[2],by=0.0005)
ym11=dnorm(x = xm11,mean = mu[1],sd = sqrt(thvars[2,2]))
linem11=data.frame(xm11,ym11)
m11 <- ggplot(m11.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xm1[1],xm1[2]) +
  ylim(0,4.5) +
  geom_line(data = linem11, aes(x=xm11,y=ym11), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

m12.m <- melt(muf[,1,4])
xm12=seq(xm1[1],xm1[2],by=0.0005)
ym12=dnorm(x = xm12,mean = mu[1],sd = sqrt(thvars[4,2]))
linem12=data.frame(xm12,ym12)
m12 <- ggplot(m12.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) +  
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xm1[1],xm1[2]) +
  ylim(0,4.5) +
  geom_line(data = linem12, aes(x=xm12,y=ym12), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(m11, m12, nrow = 1)

# phi 1
p11.m <- melt(phif[,1,2])
pm1=c(0,1)
xp11=seq(pm1[1],pm1[2],by=0.0005)
yp11=dnorm(x = xp11,mean = phi[1],sd = sqrt(thvars[2,3]))
linep11=data.frame(xp11,yp11)
p11 <- ggplot(p11.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(pm1[1],pm1[2]) +
  ylim(0,12) +
  geom_line(data = linep11, aes(x=xp11,y=yp11), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

p12.m <- melt(phif[,1,4])
xp12=seq(pm1[1],pm1[2],by=0.002)
yp12=dnorm(x = xp12,mean = phi[1],sd = sqrt(thvars[4,3]))
linep12=data.frame(xp12,yp12)
p12 <- ggplot(p12.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(pm1[1],pm1[2]) +
  ylim(0,12) +
  geom_line(data = linep12, aes(x=xp12,y=yp12), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(p11, p12, nrow = 1)

# alpha 2
a21.m <- melt(alphaf[,2,2])
xa2=c(0.4,1)
xa21=seq(xa2[1],xa2[2],by=0.0005)
ya21=dnorm(x = xa21,mean = alpha[2],sd = sqrt(thvars[2,4]))
linea21=data.frame(xa21,ya21)
a21 <- ggplot(a21.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xa2[1],xa2[2]) +
  ylim(0,9.5) +
  geom_line(data = linea21, aes(x=xa21,y=ya21), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

a22.m <- melt(alphaf[,2,4])
xa22=seq(xa2[1],xa2[2],by=0.001)
ya22=dnorm(x = xa22,mean = alpha[2],sd = sqrt(thvars[4,4]))
linea22=data.frame(xa22,ya22)
a22 <- ggplot(a22.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xa2[1],xa2[2]) +
  ylim(0,9.5) +
  geom_line(data = linea22, aes(x=xa22,y=ya22), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(a21, a22, nrow = 1)

# mu 2
m21.m <- melt(muf[,2,2])
xm2=c(-10,10)
xm21=seq(xm2[1],xm2[2],by=0.0005)
ym21=dnorm(x = xm21,mean = mu[2],sd = sqrt(thvars[2,5]))
linem21=data.frame(xm21,ym21)
m21 <- ggplot(m21.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xm2[1],xm2[2]) +
  ylim(0,0.4) +
  geom_line(data = linem21, aes(x=xm21,y=ym21), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

m22.m <- melt(muf[,2,4])
xm22=seq(xm2[1],xm2[2],by=0.0005)
ym22=dnorm(x = xm22,mean = mu[2],sd = sqrt(thvars[4,5]))
linem22=data.frame(xm22,ym22)
m22 <- ggplot(m22.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) +  
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xm2[1],xm2[2]) +
  ylim(0,0.4) +
  geom_line(data = linem22, aes(x=xm22,y=ym22), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(m21, m22, nrow = 1)

# phi 2
p21.m <- melt(phif[,2,2])
pm2=c(0,4.5)
xp21=seq(pm2[1],pm2[2],by=0.0005)
yp21=dnorm(x = xp21,mean = phi[2],sd = sqrt(thvars[2,6]))
linep21=data.frame(xp21,yp21)
p21 <- ggplot(p21.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(pm2[1],pm2[2]) +
  ylim(0,1.5) +
  geom_line(data = linep21, aes(x=xp21,y=yp21), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

p22.m <- melt(phif[,2,4])
xp22=seq(pm2[1],pm2[2],by=0.002)
yp22=dnorm(x = xp22,mean = phi[2],sd = sqrt(thvars[4,6]))
linep22=data.frame(xp22,yp22)
p22 <- ggplot(p22.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(pm2[1],pm2[2]) +
  ylim(0,1.5) +
  geom_line(data = linep22, aes(x=xp22,y=yp22), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(p21, p22, nrow = 1)


# alpha 3
a31.m <- melt(alphaf[,3,2])
xa3=c(0.0,0.75)
xa31=seq(xa3[1],xa3[2],by=0.0005)
ya31=dnorm(x = xa31,mean = alpha[3],sd = sqrt(thvars[2,7]))
linea31=data.frame(xa31,ya31)
a31 <- ggplot(a31.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xa3[1],xa3[2]) +
  ylim(0,10.5) +
  geom_line(data = linea31, aes(x=xa31,y=ya31), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

a32.m <- melt(alphaf[,3,4])
xa32=seq(xa3[1],xa3[2],by=0.001)
ya32=dnorm(x = xa32,mean = alpha[3],sd = sqrt(thvars[4,7]))
linea32=data.frame(xa32,ya32)
a32 <- ggplot(a32.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xa3[1],xa3[2]) +
  ylim(0,10.5) +
  geom_line(data = linea32, aes(x=xa32,y=ya32), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(a31, a32, nrow = 1)

# mu 3
m31.m <- melt(muf[,3,2])
xm3=c(2.4,3.6)
xm31=seq(xm3[1],xm3[2],by=0.0005)
ym31=dnorm(x = xm31,mean = mu[3],sd = sqrt(thvars[2,8]))
linem31=data.frame(xm31,ym31)
m31 <- ggplot(m31.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(xm3[1],xm3[2]) +
  ylim(0,8) +
  geom_line(data = linem31, aes(x=xm31,y=ym31), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

m32.m <- melt(muf[,3,4])
xm32=seq(xm3[1],xm3[2],by=0.0005)
ym32=dnorm(x = xm32,mean = mu[3],sd = sqrt(thvars[4,8]))
linem32=data.frame(xm32,ym32)
m32 <- ggplot(m32.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) +  
  xlab("sample size = 250") + 
  ylab("") +
  xlim(xm3[1],xm3[2]) +
  ylim(0,8) +
  geom_line(data = linem32, aes(x=xm32,y=ym32), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(m31, m32, nrow = 1)

# phi 3
p31.m <- melt(phif[,3,2])
pm3=c(0,0.15)
xp31=seq(pm3[1],pm3[2],by=0.0005)
yp31=dnorm(x = xp31,mean = phi[3],sd = sqrt(thvars[2,9]))
linep31=data.frame(xp31,yp31)
p31 <- ggplot(p31.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(pm3[1],pm3[2]) +
  ylim(0,45) +
  geom_line(data = linep31, aes(x=xp31,y=yp31), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

p32.m <- melt(phif[,3,4])
xp32=seq(pm3[1],pm3[2],by=0.002)
yp32=dnorm(x = xp32,mean = phi[3],sd = sqrt(thvars[4,9]))
linep32=data.frame(xp32,yp32)
p32 <- ggplot(p32.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(pm3[1],pm3[2]) +
  ylim(0,45) +
  geom_line(data = linep32, aes(x=xp32,y=yp32), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(p31, p32, nrow = 1)

# copula parameter
c1.m <- melt(copparsf[,2])
pc=c(1,1.7)
xc1=seq(pc[1],pc[2],by=0.0005)
yc1=dnorm(x = xc1,mean = coppars,sd = sqrt(thvars[2,10]))
linec1=data.frame(xc1,yc1)
c1 <- ggplot(c1.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 50") + 
  ylab("") +
  xlim(pc[1],pc[2]) +
  ylim(0,9) +
  geom_line(data = linec1, aes(x=xc1,y=yc1), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

c2.m <- melt(copparsf[,4])
xc2=seq(pc[1],pc[2],by=0.002)
yc2=dnorm(x = xc2,mean = coppars,sd = sqrt(thvars[4,10]))
linec2=data.frame(xc2,yc2)
c2 <- ggplot(c2.m, aes(x=value,y = ..density..)) +
  geom_histogram(bins=30) + 
  xlab("sample size = 250") + 
  ylab("") +
  xlim(pc[1],pc[2]) +
  ylim(0,9) +
  geom_line(data = linec2, aes(x=xc2,y=yc2), color = "red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11(width=7,height=3)
grid.arrange(c1, c2, nrow = 1)





#########################################################
### average/median variance of all fitted repetitions ###
#########################################################

### Gaussian copula ###
#######################

# dimensions
d=2
# number of monte carlo simulations
reps=1000
# seed
seed=1245

# type of univariate QB distribution (of length d)
margfuncs=c("t","normal")
# alpha parameters (of length d)
alpha=c(0.2,0.7)
# mu parameters (of length d)
mu=c(2,3)
# phi parameters (of length d)
phi=c(0.2,0.4)
# degrees of freedom (of length the number of "t" margins)
df=3

# copula parameters
coppars=0.2

### defining the copula object to simulate from
copulamodel=normalCopula(param = coppars, dim = d, dispstr = "un")

# sample sizes considered
ss=c(75,600)
# length of copula parameter vector
cp=1

load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/normal/normalcop_dim",d,"_S",75,".Rdata"))
fits1=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 1 (revisie)/normal/normalcop_dim",d,"_S",600,".Rdata"))
fits2=dataset

# grouping parameters
alphaf=array(data=NA,dim=c(reps,d,2))
muf=array(data=NA,dim=c(reps,d,2))
phif=array(data=NA,dim=c(reps,d,2))
dff=array(data=NA,dim=c(reps,d,2))
loglf=matrix(NA,nrow=reps,ncol=2)
copparsf=matrix(data=NA,nrow=reps,ncol=2)
time.taken=matrix(NA,nrow=reps,ncol=2)
for(i in 1:reps){
  temp_dat=fits1[[i]]
  alphaf[i,,1]=temp_dat$alpha
  muf[i,,1]=temp_dat$mu
  phif[i,,1]=temp_dat$phi
  dff[i,,1]=temp_dat$df
  loglf[i,1]=temp_dat$logl
  copparsf[i,1]=temp_dat$copulapar
  time.taken[i,1]=temp_dat$time
  
  temp_dat=fits2[[i]]
  alphaf[i,,2]=temp_dat$alpha
  muf[i,,2]=temp_dat$mu
  phif[i,,2]=temp_dat$phi
  dff[i,,2]=temp_dat$df
  loglf[i,2]=temp_dat$logl
  copparsf[i,2]=temp_dat$copulapar
  time.taken[i,2]=temp_dat$time
}

### theoretical variance of parameter estimates
# seed
seed=1234
# number of observations to approximate the theoretical variance
n=25000


# holding matrix
covpar=array(NA,dim=c(reps,8,8,2))
for(j in 1:reps){
  for(i in 1:2){
    
    alpha=alphaf[j,,i]
    mu=muf[j,,i]
    phi=phif[j,,i]
    df=dff[j,,i]
    coppars=copparsf[j,i]
    copulamodel=normalCopula(param = coppars,dim = d,dispstr = "un")
    
    # drawing sample
    S=samplecopula(n = n,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed = seed)
    Z=S[[2]]
    try({
      # calculating the approximate FI-matrix
      mats=FI_IFM_Gaussian(Z=Z,alpha=alpha,mu=mu,phi=phi,df=df,coppars=coppars,margfuncs=margfuncs)
      FIM=mats$FIM
      vars=solve(mats$FIM)
      
      # replace the margin part with a more accurate approximation. This is theoretically justified
      vcm=matrix(0,7,7)
      count=1
      for(k in 1:2){
        temp=asymptoticvariance(alpha = alpha[k],mu = mu[k],phi = phi[k],nu = df[k],margfunc = margfuncs[k])
        if(k==1){
          vcm[count:(count+3),count:(count+3)]=temp$cov
          count=count+4
        } else {
          vcm[count:(count+2),count:(count+2)]=temp$cov
          count=count+3
        }
      }  
      
      vars[1:7,1:7]=vcm
      covpar[j,,,i]=vars
      
    },silent=T)
  }
  print(j)
}

write.table(round(apply(covpar[,,,1],c(2,3),median,na.rm=T),4),
            file="~/PhD/code/copulas/output/covgaussian1.txt", sep = "$ & $",
            row.names=F,col.names=F)

write.table(round(apply(covpar[,,,2],c(2,3),median,na.rm=T),4),
            file="~/PhD/code/copulas/output/covgaussian2.txt", sep = "$ & $",
            row.names=F,col.names=F)



### Gumbel copula ###
#####################


### model

# number of monte carlo simulations
reps=1000
# dimensions of the problem
d=3
# number of samples to be generated
sampsize=25  #25 50 100 250 1000
# number of starting points in the fitting of the marginals
nstart=10



### information on the margins

# type of univariate QB distribution (of length d)
margfuncs=c("laplace","normal","logistic")
# alpha parameters (of length d)
alpha=c(0.4,0.7,0.3)
# mu parameters (of length d)
mu=c(1,2,3)
# phi parameters (of length d)
phi=c(0.5,2.9,0.1)
# degrees of freedom (of length the number of "t" margins)
df=NULL
# copula parameters
coppars=1.3

### defining the copula object to simulate from
copulamodel=gumbelCopula(param = coppars, dim = d) 

### data of fits

load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",25,".Rdata"))
fits1=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",50,".Rdata"))
fits2=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",100,".Rdata"))
fits3=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",250,".Rdata"))
fits4=dataset
load(paste0("~/PhD/code/copulas/output/paper copula setting 2 (revisie)/gumbelcop_dim",d,"_S",1000,".Rdata"))
fits5=dataset

# grouping parameters
alphaf=array(data=NA,dim=c(reps,d,5))
muf=array(data=NA,dim=c(reps,d,5))
phif=array(data=NA,dim=c(reps,d,5))
dff=array(data=NA,dim=c(reps,d,5))
loglf=matrix(NA,nrow=reps,ncol=5)
copparsf=matrix(data=NA,nrow=reps,ncol=5)
time.taken=matrix(NA,nrow=reps,ncol=5)
for(i in 1:reps){
  temp_dat=fits1[[i]]
  alphaf[i,,1]=temp_dat$alpha
  muf[i,,1]=temp_dat$mu
  phif[i,,1]=temp_dat$phi
  dff[i,,1]=temp_dat$df
  loglf[i,1]=temp_dat$logl
  copparsf[i,1]=temp_dat$copulapar
  time.taken[i,1]=temp_dat$time
  
  temp_dat=fits2[[i]]
  alphaf[i,,2]=temp_dat$alpha
  muf[i,,2]=temp_dat$mu
  phif[i,,2]=temp_dat$phi
  dff[i,,2]=temp_dat$df
  loglf[i,2]=temp_dat$logl
  copparsf[i,2]=temp_dat$copulapar
  time.taken[i,2]=temp_dat$time
  
  temp_dat=fits3[[i]]
  alphaf[i,,3]=temp_dat$alpha
  muf[i,,3]=temp_dat$mu
  phif[i,,3]=temp_dat$phi
  dff[i,,3]=temp_dat$df
  loglf[i,3]=temp_dat$logl
  copparsf[i,3]=temp_dat$copulapar
  time.taken[i,3]=temp_dat$time
  
  temp_dat=fits4[[i]]
  alphaf[i,,4]=temp_dat$alpha
  muf[i,,4]=temp_dat$mu
  phif[i,,4]=temp_dat$phi
  dff[i,,4]=temp_dat$df
  loglf[i,4]=temp_dat$logl
  copparsf[i,4]=temp_dat$copulapar
  time.taken[i,4]=temp_dat$time
  
  temp_dat=fits5[[i]]
  alphaf[i,,5]=temp_dat$alpha
  muf[i,,5]=temp_dat$mu
  phif[i,,5]=temp_dat$phi
  dff[i,,5]=temp_dat$df
  loglf[i,5]=temp_dat$logl
  copparsf[i,5]=temp_dat$copulapar
  time.taken[i,5]=temp_dat$time
}



### theoretical variance of parameter estimates
# seed
seed=1234
# number of observations to approximate the theoretical variance
n=25000

# holding matrix
varpar=array(NA,dim=c(reps,10,10,5))
for(j in 1:reps){
  for(i in 1:5){
    
    alpha=alphaf[j,,i]
    mu=muf[j,,i]
    phi=phif[j,,i]
    df=NULL
    coppars=copparsf[j,i]
    copulamodel=gumbelCopula(param = coppars,dim = d)
    
    # drawing sample
    S=samplecopula(n = n,alpha = alpha,mu = mu,phi = phi,df = df,margfuncs = margfuncs,cop = copulamodel,seed = seed)
    Z=S[[2]]
    # calculating the approximate FI-matrix
    mats=FI_IFM(Z=Z,alpha=alpha,mu=mu,phi=phi,df=NULL,coppars=coppars,margfuncs=margfuncs,copula="gumbel")
    FIM=mats$FIM
    vars=solve(mats$FIM)
    
    # replace the margin part with a more accurate approximation. This is theoretically justified
    vcm=matrix(0,9,9)
    vars2=matrix(NA,3,3)
    
    for(k in 1:3){
      temp=asymptoticvariance(alpha = alpha[k],mu = mu[k],phi = phi[k],nu = NULL,margfunc = margfuncs[k])
      vars2[k,]=temp$vars
      vcm[(3*k-2):(3*k),(3*k-2):(3*k)]=temp$cov
      vars[1:9,1:9]=vcm
    }
    varpar[j,,,i]=vars
  }
  print(j)
}





#############################################
### Code for testing the impact of the    ###
### weighting matrix on the GMM estimates ###
#############################################

# function to compute g(theta,X) in a sample of size n
# returns a nxp matrix with p the number of moment equations
gId=function(par,x){
  u1=pAND(q = x[,1],mu = par[2],phi = par[3],alpha = par[1])
  u2=pAND(q = x[,2],mu = par[5],phi = par[6],alpha = par[4])  
  x1=qnorm(u1)
  x2=qnorm(u2)
  g1=(1-2*par[1])/(par[1]*(1-par[1]))+1*(x[,1]<=par[2])*(1-par[1])*(par[2]-x[,1])^2/par[3]^2-1*(x[,1]>par[2])*par[1]*(par[2]-x[,1])^2/par[3]^2
  g2=-1*(x[,1]<=par[2])*(1-par[1])^2/par[3]^2*(par[2]-x[,1])+1*(x[,1]>par[2])*(par[1])^2/par[3]^2*(-par[2]+x[,1])
  g3=-1/par[3]+1*(x[,1]<=par[2])*(1-par[1])^2*(par[2]-x[,1])^2/par[3]^3+1*(x[,1]>par[2])*(par[1])^2*(par[2]-x[,1])^2/par[3]^3
  g4=(1-2*par[4])/(par[4]*(1-par[4]))+1*(x[,2]<=par[5])*(1-par[4])*(par[5]-x[,2])^2/par[6]^2-1*(x[,2]>par[5])*(par[4])*(par[5]-x[,2])^2/par[6]^2
  g5=-1*(x[,2]<=par[5])*(1-par[4])^2/par[6]^2*(par[5]-x[,2])+1*(x[,2]>par[5])*(par[4])^2/par[6]^2*(-par[5]+x[,2])
  g6=-1/par[6]+1*(x[,2]<=par[5])*(1-par[4])^2*(par[5]-x[,2])^2/par[6]^3+1*(x[,2]>par[5])*(par[4])^2*(par[5]-x[,2])^2/par[6]^3
  g7=par[7]/(1-par[7]^2)-(par[7]*(x1^2+x2^2)-x1*x2)/(1-par[7]^2)-(par[7]^3*(x1^2+x2^2)-par[7]^2*x1*x2)/(1-par[7]^2)^2
  
  crit=cbind(g1,g2,g3,g4,g5,g6,g7)
  return(sum(colMeans(crit)^2))
}

gOp=function(par,x){
  u1=pAND(q = x[,1],mu = par[2],phi = par[3],alpha = par[1])
  u2=pAND(q = x[,2],mu = par[5],phi = par[6],alpha = par[4])  
  x1=qnorm(u1)
  x2=qnorm(u2)
  
  g1=(1-2*par[1])/(par[1]*(1-par[1]))+1*(x[,1]<=par[2])*(1-par[1])*(par[2]-x[,1])^2/par[3]^2-1*(x[,1]>par[2])*par[1]*(par[2]-x[,1])^2/par[3]^2
  g2=-1*(x[,1]<=par[2])*(1-par[1])^2/par[3]^2*(par[2]-x[,1])+1*(x[,1]>par[2])*(par[1])^2/par[3]^2*(-par[2]+x[,1])
  g3=-1/par[3]+1*(x[,1]<=par[2])*(1-par[1])^2*(par[2]-x[,1])^2/par[3]^3+1*(x[,1]>par[2])*(par[1])^2*(par[2]-x[,1])^2/par[3]^3
  g4=(1-2*par[4])/(par[4]*(1-par[4]))+1*(x[,2]<=par[5])*(1-par[4])*(par[5]-x[,2])^2/par[6]^2-1*(x[,2]>par[5])*(par[4])*(par[5]-x[,2])^2/par[6]^2
  g5=-1*(x[,2]<=par[5])*(1-par[4])^2/par[6]^2*(par[5]-x[,2])+1*(x[,2]>par[5])*(par[4])^2/par[6]^2*(-par[5]+x[,2])
  g6=-1/par[6]+1*(x[,2]<=par[5])*(1-par[4])^2*(par[5]-x[,2])^2/par[6]^3+1*(x[,2]>par[5])*(par[4])^2*(par[5]-x[,2])^2/par[6]^3
  g7=par[7]/(1-par[7]^2)-(par[7]*(x1^2+x2^2)-x1*x2)/(1-par[7]^2)-(par[7]^3*(x1^2+x2^2)-par[7]^2*x1*x2)/(1-par[7]^2)^2
  
  g=cbind(g1,g2,g3,g4,g5,g6,g7)
  w=t(g)%*%g
  mg=matrix(colMeans(g),nrow=1,ncol=7)
  n=length(x[,1])
  
  return(mg%*%w%*%t(mg)/n)
}


alpha=c(0.3,0.8)
mu=c(0.5,2)
phi=c(0.5,2)
rho=0.4
cop=normalCopula(param = rho,dim = 2,dispstr = "un")
ss=50
seed=2901

m=1000
paropt=matrix(NA,nrow=m,ncol=7)
objopt=rep(NA,m)
parid=matrix(NA,nrow=m,ncol=7)
objid=rep(NA,m)
parifm=matrix(NA,nrow=m,ncol=7)
k=10
for(i in 1:m){
  seed=seed+i
  start=cbind(runif(k,min=0.15,max=0.45),runif(k,min=0.25,max=0.75),runif(k,min=0.25,max=0.75),runif(k,min=0.65,max=0.95),runif(k,min=1.5,max=2.5),runif(k,min=1.5,max=2.5),runif(k,min=0.15,max=0.65))
  parsid=matrix(NA,nrow=k,ncol=7)
  valueid=rep(NA,k)
  parsopt=matrix(NA,nrow=k,ncol=7)
  valueopt=rep(NA,k)
  X=samplecopula(n = ss,alpha = alpha,mu = mu,phi = phi,df = NULL,margfuncs = c("normal","normal"),cop = cop,seed = seed)
  for(j in 1:k){
    try({
      resopt=bobyqa(x0 = start[j,],fn = gOp,lower = c(0.1,-Inf,0.1,0.1,-Inf,0.1,-0.9),upper = c(0.9,Inf,Inf,0.9,Inf,Inf,0.9),x=X$X)
      parsopt[j,]=resopt$par
      valueopt[j]=resopt$value},silent=T)
    try({
      resid=bobyqa(x0 = start[j,],fn = gId,lower = c(0.1,-Inf,0.1,0.1,-Inf,0.1,-0.9),upper = c(0.9,Inf,Inf,0.9,Inf,Inf,0.9),x=X$X)
      parsid[j,]=resid$par
      valueid[j]=resid$value},silent=T)
  }
  try({
    paropt[i,]=parsopt[which.min(valueopt),]
    objopt[i]=valueopt[which.min(valueopt)]
  },silent = T)
  try({
    parid[i,]=parsid[which.min(valueid),]
    objid[i]=valueid[which.min(valueid)]
  },silent = T)
  try({
    resifm=fitcopula(data = X$X,margfuncs = c("normal","normal"),cop = normalCopula(dim=2,dispstr = "un"),nstart = 10,seed = seed)
    parifm[i,]=c(resifm$alpha[1],resifm$mu[1],resifm$phi[1],resifm$alpha[2],resifm$mu[2],resifm$phi[2],resifm$copulapar)
  },silent = T)
  print(i)
}
example=list("parWId"=parid,"parWOpt"=paropt,"parIFM"=parifm)
save(example,file = "~/PhD/code/copulas/output/impactW/example.Rdata")

boxplot(cbind(parid,paropt,parifm))
colMeans(parid)
colMeans(paropt)
colMeans(parifm)


alpha1=data.frame(cbind(parid[,1],paropt[,1],parifm[,1]))
colnames(alpha1)=c("W=Id","W=W_opt","IFM")
alpha1=melt(alpha1)
x11()
ggplot(alpha1, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(alpha[1])) +
  geom_abline(slope = 0,intercept = alpha[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

mu1=data.frame(cbind(parid[,2],paropt[,2],parifm[,2]))
colnames(mu1)=c("W=Id","W=W_opt","IFM")
mu1=melt(mu1)
x11()
ggplot(mu1, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(mu[1])) +
  geom_abline(slope = 0,intercept = mu[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


phi1=data.frame(cbind(parid[,3],paropt[,3],parifm[,3]))
colnames(phi1)=c("W=Id","W=W_opt","IFM")
phi1=melt(phi1)
x11()
ggplot(phi1, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(phi[1])) +
  geom_abline(slope = 0,intercept = phi[1],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

alpha2=data.frame(cbind(parid[,4],paropt[,4],parifm[,4]))
colnames(alpha2)=c("W=Id","W=W_opt","IFM")
alpha2=melt(alpha2)
x11()
ggplot(alpha2, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(alpha[2])) +
  geom_abline(slope = 0,intercept = alpha[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

mu2=data.frame(cbind(parid[,5],paropt[,5],parifm[,5]))
colnames(mu2)=c("W=Id","W=W_opt","IFM")
mu2=melt(mu2)
x11()
ggplot(mu2, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(mu[2])) +
  geom_abline(slope = 0,intercept = mu[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


phi2=data.frame(cbind(parid[,6],paropt[,6],parifm[,6]))
colnames(phi2)=c("W=Id","W=W_opt","IFM")
phi2=melt(phi2)
x11()
ggplot(phi2, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(phi[2])) +
  geom_abline(slope = 0,intercept = phi[2],color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


rho1=data.frame(cbind(parid[,7],paropt[,7],parifm[,7]))
colnames(rho1)=c("W=Id","W=W_opt","IFM")
rho1=melt(rho1)
x11()
ggplot(rho1, aes(x=variable, y=value)) + 
  geom_boxplot(fill="white") +
  xlab("Choice W") + 
  ylab(expression(rho)) +
  geom_abline(slope = 0,intercept = rho,color="red",size=1) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

