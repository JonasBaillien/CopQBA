########################################################################
########################################################################
### Code accompagnying the paper:                                    ###
###     "Inference for copulas with two-piece margins"               ###
### Authors: Baillien Jonas, Gijbels irène and Verhasselt Anneleen   ###
########################################################################
########################################################################

##################################
### Part 3: Real data examples ###
##################################

### Source the functions file ###
#################################

source("Code part 1 Functions.R")


### Important notice ###
########################
# Please read comments as certain elements can change depending on the setting
# Pay special attention to lines 59-100 as the dataset needs to be changed there

# Data are loaded from a certain directory, you need to change this to where you saved the data
# Certain values are calculated and saved for easier processing, also make sure you set the correct 
# directories there.
# This occurs on lines: 37, 46, 407, 527, 699, ,708, 719, 727, 738, 746, 773, 774, 
# 1505, 1515,1525, 1547, 1548, 1549, 1564, 1565, 1566, 1580, 1581, 1582 




### Wine data ### 
#################
# White wine data
# Can be downloaded from https://www.kaggle.com/maitree/wine-quality-selection
winequality.white <- read.csv("~/PhD/code/datasets/winequality-white.csv", sep=";") # change directory
whitewinedata=winequality.white[winequality.white$quality==7,c(2,9,10)]
pairs(whitewinedata)



### Pokemon data ###
####################
# Can be downloaded from https://www.kaggle.com/mlomuscio/pokemon
Pokemon <- read.csv("~/PhD/code/datasets/Pokemon.csv", header=FALSE, comment.char="#") # change directory
pdata=as.matrix(Pokemon[,6:11])

set.seed(147)
pdata2=pdata+runif(n = 6*800,min = -2,max = 2) # add noise for visual representation

x11()
pairs(pdata2)

#####################################
### fitting the model to the data ###
#####################################

Z=whitewinedata   # insert correct data here
# Z=pdata

# dimensions of the data
d=length(Z[1,])
n=length(Z[,1])


### margin fitting
U=matrix(NA,nrow=n,ncol=d)
basefunc=rep(NA,d)
pars=matrix(NA,nrow=4,ncol=d)
marginlogl=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Z[,i],crit = "AIC",nstart = 10,all = F,seed = 1489)
  U[,i]=tempout$U
  basefunc[i]=tempout$margin
  pars[,i]=tempout$parameters
  marginlogl[i]=tempout$logl
}


### copula fitting
tempout=bestcopulafit(U = U,all = T,crit = "AIC")
tempout$info

copula=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulapars=unlist(tempout$parameters[copula])
copulalogl=as.numeric(tempout$info[which.min(tempout$info[,3]),2])


### CIC penalty for margins and copula
# change the copula, theta (reach: 1:3 for wine, 1:15 for pokemon) and nu (omit when not "t", 4 for wine, 16 for pokemon)
penalty=CICpenalty(Z = Z,margfuncs = basefunc,alpha = pars[1,],mu = pars[2,],phi = pars[3,],df = pars[4,],copula = "t",theta = tempout$parameters$t[1:3],nu = tempout$parameters$t[4])
penalty

# copula contribution
-2*as.numeric(tempout$info[6,2])+2*penalty$copCIC    # 1=gumbel, 2=clayton, 3=frank, 4=joe, 5=normal, 6=t

# margin contribution
-2*marginlogl+2*penalty$marginCIC



#############
### Plots ###
#############

### histogram margin pokemon with overlay fitted distributions ###
##################################################################

Z=pdata
variables=c("Hitpoints","Attack","Defence","Special Attack","Special Defence","Speed")
d=length(Z[1,])
n=length(Z[,1])

# i to select the variable index for plotting
i=3
# fitting all 4 QBA-margin for variable i
bf=bestmargin(data = Z[,i],crit = "AIC",nstart = 10,all = T,seed = 1489)
bf
pars=matrix(as.numeric(bf[,3:6]),nrow=4)

# plot of the data and 4 fitted QBA-margins
def=data.frame(Defence=Z[,3])
xx=seq(-10,260,length.out = 1000)
yN = QBAsyDist::dAND(y = xx,mu = pars[1,2],phi = pars[1,3],alpha = pars[1,1])
yLo = QBAsyDist::dALoD(y = xx,mu = pars[3,2],phi = pars[3,3],alpha = pars[3,1])
yLa = QBAsyDist::dALaD(y = xx,mu = pars[2,2],phi = pars[2,3],alpha = pars[2,1])
yT = QBAsyDist::dATD(y = xx,mu = pars[4,2],phi = pars[4,3],alpha = pars[4,1],nu = pars[4,4])
data=data.frame(cbind(yN,yLa,yLo,yT))
x11()
ggplot(data = def,aes(x=Defence)) + 
  geom_histogram(aes(y = ..density..),breaks=seq(0,260,by=20),fill="gray",color="black") + 
  ylab("Density") + 
  xlim(-20,270) +
  geom_line(data = data,aes(x=xx,y=yN,color="1"),size=1) + 
  geom_line(data = data,aes(x=xx,y=yLa,color="2"),size=1) + 
  geom_line(data = data,aes(x=xx,y=yLo,color="3"),size=1) + 
  geom_line(data = data,aes(x=xx,y=yT,color="4"),size=1) +
  scale_color_manual(name = "", labels = c("1"="QBA-normal", "2"="QBA-Laplace","3"="QBA-logistic","4"="QBA-Student's t"), values=c("black","red","green","blue")) +
  theme(legend.position="bottom")  +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14))



### fitted margins and scatterplot white wine quality 7 ###
###########################################################

# white wine of quality 7
ww7=winequality.white[winequality.white$quality==7,c(2,9,10)] 
# histogram of each variable of white wine quality 7
h1 <- ggplot(data = ww7,aes(x=volatile.acidity)) + 
  geom_histogram(aes(y = ..density..),breaks=seq(0,0.8,by=0.04),fill="gray",color="black") +
  xlab("Volatile Acidity") +
  ylab("Density") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
h2 <- ggplot(data = ww7,aes(x=pH)) + 
  geom_histogram(aes(y = ..density..),breaks=seq(2.8,3.9,by=0.1),fill="gray",color="black") +
  xlab("pH") +
  ylab("Density") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
h3 <- ggplot(data = ww7,aes(x=sulphates)) + 
  geom_histogram(aes(y = ..density..),breaks=seq(0.2,1.1,by=0.1),fill="gray",color="black") +
  xlab("Sulphates") +
  ylab("Density") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
# pairwise scatterplots of white wine quality 7
s1 <- ggplot(data = ww7,aes(x=volatile.acidity,y=pH)) +
  geom_point(size=2) +
  xlab("Volatile Acidity") +
  ylab("pH") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
s2 <- ggplot(data = ww7,aes(x=volatile.acidity,y=sulphates)) +
  geom_point(size=2) +
  xlab("Volatile Acidity") +
  ylab("Sulphates") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
s3 <- ggplot(data = ww7,aes(x=pH,y=sulphates)) +
  geom_point(size=2) +
  xlab("pH") +
  ylab("Sulphates") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11()
grid.arrange(h1,h2,h3,s1,s2,s3,nrow=2,ncol=3)




###################################
### Boothstrap smoothed KS test ### 
###################################

### Pokemon data ###
####################

### fitting the model to the data ###
#####################################

Z=pdata

# dimensions of the data
d=length(Z[1,])
n=length(Z[,1])

# holding matrices for the type of margin and its parameters
basefunc=rep(NA,d)
pars=matrix(NA,nrow=4,ncol=d)


# determining the best fitting margin of the QBA-family
for(i in 1:d){
  tempout=bestmargin(data = Z[,i],crit = "AIC",nstart = 10,all = F,seed = 1489)
  basefunc[i]=tempout$margin
  pars[,i]=tempout$parameters
  
}

# smoothed Kolmogorov-Smirnov test on the fitted margins
KSS1=KS.test(x=Z[,1],y = "pATD",hx=2,alternative = "two.sided",mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1],nu = pars[4,1],kernel = "unif")
KSS2=KS.test(x=Z[,2],y = "pAND",alternative = "two.sided",mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2],kernel = "unif")
KSS3=KS.test(x=Z[,3],y = "pALoD",alternative = "two.sided",mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3],kernel = "unif")
KSS4=KS.test(x=Z[,4],y = "pAND",alternative = "two.sided",mu = pars[2,4],phi = pars[3,4],alpha = pars[1,4],kernel = "unif")
KSS5=KS.test(x=Z[,5],y = "pATD",alternative = "two.sided",mu = pars[2,5],phi = pars[3,5],alpha = pars[1,5],nu = pars[4,5],kernel = "unif")
KSS6=KS.test(x=Z[,6],y = "pAND",alternative = "two.sided",mu = pars[2,6],phi = pars[3,6],alpha = pars[1,6],kernel = "unif")

# number of replicates
N=1000

### First variable: HP
# test statistic from the original data
P1.ks.teststat=KSS1$statistic

# seed for sample generation
set.seed(147)

# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rATD(n = N*n,mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1],nu = pars[4,1]),nrow=N,ncol=n)

# holding vector for teststatistics of the bootstrap samples
P1.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitATD(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pATD",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,nu = fit$nu,alternative = "two.sided",kernel = "unif")
  P1.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P1.ks.stats)
abline(v=P1.ks.teststat)
# approximate p-value
P1.ks.pvalue=sum(P1.ks.stats>P1.ks.teststat)/N


### Second variable: Attack
# test statistic from the original data
P2.ks.teststat=KSS2$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
P2.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  P2.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P2.ks.stats)
abline(v=P2.ks.teststat)
# approximate p-value
P2.ks.pvalue=sum(P2.ks.stats>P2.ks.teststat)/N


### Third variable: Defense
# test statistic from the original data
P3.ks.teststat=KSS3$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rALoD(n = N*n,mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
P3.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitALoD(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pALoD",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  P3.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P3.ks.stats)
abline(v=P3.ks.teststat)
# approximate p-value
P3.ks.pvalue=sum(P3.ks.stats>P3.ks.teststat)/N


### Fourth variable: Special Attack
# test statistic from the original data
P4.ks.teststat=KSS4$statistic


# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = pars[2,4],phi = pars[3,4],alpha = pars[1,4]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
P4.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  P4.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P4.ks.stats)
abline(v=P4.ks.teststat)
# approximate p-value
P4.ks.pvalue=sum(P4.ks.stats>P4.ks.teststat)/N


### Fifth variable: Special Defense
# test statistic from the original data
P5.ks.teststat=KSS5$statistic


# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rATD(n = N*n,mu = pars[2,5],phi = pars[3,5],alpha = pars[1,5],nu = pars[4,5]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
P5.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitATD(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pATD",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,nu = fit$nu,alternative = "two.sided",kernel = "unif")
  P5.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P5.ks.stats)
abline(v=P5.ks.teststat)
# approximate p-value
P5.ks.pvalue=sum(P5.ks.stats>P5.ks.teststat)/N


### Sixt variable: Speed
# test statistic from the original data
P6.ks.teststat=KSS6$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = pars[2,6],phi = pars[3,6],alpha = pars[1,6]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
P6.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  P6.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(P6.ks.stats)
abline(v=P6.ks.teststat)
# approximate p-value
P6.ks.pvalue=sum(P6.ks.stats>P6.ks.teststat)/N


# saving for later use
BSvaluesP=cbind(P1.ks.stats,P2.ks.stats,P3.ks.stats,P4.ks.stats,P5.ks.stats,P6.ks.stats)
teststatsP=c(P1.ks.teststat,P2.ks.teststat,P3.ks.teststat,P4.ks.teststat,P5.ks.teststat,P6.ks.teststat)
output=list("bootstrapvalues"=BSvaluesP,"teststat"=teststatsP)
save(output,file = "~/PhD/code/copulas/output/paper copula pokemon gof smoothed.Rdata")




### wine data ###
#################

# dimensions of the data
d=length(ww7[1,])
n=length(ww7[,1])


### holding matrices for margins and their parameters
basefuncw7=rep(NA,d)
parsw7=matrix(NA,nrow=4,ncol=d)

### determining the best fitting QBA margin
for(i in 1:d){
  tempout=bestmargin(data = ww7[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  basefuncw7[i]=tempout$margin
  parsw7[,i]=tempout$parameters
}

### smoothed Kolmogorov-Smirnov test on the fitted margins
KSSW1=KS.test(x=ww7[,1],y = "pAND",alternative = "two.sided",mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1],kernel = "unif")
KSSW2=KS.test(x=ww7[,2],y = "pAND",alternative = "two.sided",mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2],kernel = "unif")
KSSW3=KS.test(x=ww7[,3],y = "pALoD",alternative = "two.sided",mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3],kernel = "unif")

### number of replicates
N=1000

### First variable: Volatile Acidity
# test statistic from the original data
W1.ks.teststat=KSSW1$statistic


# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W1.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  W1.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W1.ks.stats)
abline(v=W1.ks.teststat)
# approximate p-value
W1.ks.pvalue=sum(W1.ks.stats>W1.ks.teststat)/N


### Second variable: pH
# test statistic from the original data
W2.ks.teststat=KSSW2$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W2.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  W2.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W2.ks.stats)
abline(v=W2.ks.teststat)
# approximate p-value
W2.ks.pvalue=sum(W2.ks.stats>W2.ks.teststat)/N


### Third variable: Sulphates
# test statistic from the original data
W3.ks.teststat=KSSW3$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rALoD(n = N*n,mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W3.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitALoD(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=KS.test(x = bootstrapsamples[i,],y = "pALoD",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided",kernel = "unif")
  W3.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W3.ks.stats)
abline(v=W3.ks.teststat)
# approximate p-value
W3.ks.pvalue=sum(W3.ks.stats>W3.ks.teststat)/N

# saving results for later use
BSvaluesW=cbind(W1.ks.stats,W2.ks.stats,W3.ks.stats)
teststatsW=c(W1.ks.teststat,W2.ks.teststat,W3.ks.teststat)
output=list("bootstrapvalues"=BSvaluesW,"teststat"=teststatsW)
save(output,file = "~/PhD/code/copulas/output/paper copula white wine gof smooth.Rdata")










##########################################
### Further investigation of wine data ###
##########################################


### Creating contourplots for a 3D Student's t copula used in the white wine of quality 4 ###
#############################################################################################

### Integrating the first margin out of a 3D Student's t-copula
### used for obtaining the bivariate margins
tcop1=function(x,y,nu,Sigma,pars){
  # x is the representative of dimension d, the variable over which to integrate
  # the 3D Student-copula with degrees of freedom nu and covariance Sigma
  # pars are the parameters for the margins (4x3 matrix).
  z=rep(NA,3)
  z[1]=x
  z[c(2,3)]=y
  
  # Student's t-quantile function evaluated in the suitable pseudo observations
  s=rep(NA,3)
  s[1]=qt(p = pALaD(q = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]),df = nu)
  s[2]=qt(p = pAND(q = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]),df = nu)
  s[3]=qt(p = pAND(q = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]),df = nu)
  
  # product of densities from univariate margins
  f=rep(NA,3)
  f[1]=dALaD(y = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
  f[2]=dAND(y = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
  f[3]=dAND(y = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])
  
  # density
  dens=( gamma((nu+3)/2)*gamma(nu/2)^(2) )/( sqrt(det(Sigma))*gamma((nu+1)/2)^3 )*(1+t(s)%*%solve(Sigma)%*%s/nu)^(-(nu+3)/2)/( prod(1+s^2/nu)^(-(nu+1)/2) )*prod(f) 
  
  return(dens)
}

### Integrating the second margin out of a 3D Student's t-copula
### used for obtaining the bivariate margins
tcop2=function(x,y,nu,Sigma,pars){
  # x is the representative of dimension d, the variable over which to integrate
  # the 3D Student-copula with degrees of freedom nu and covariance Sigma
  # pars are the parameters for the margins (4x3 matrix).
  z=rep(NA,3)
  z[2]=x
  z[c(1,3)]=y
  
  # Student's t-quantile function evaluated in the suitable pseudo observations
  s=rep(NA,3)
  s[1]=qt(p = pALaD(q = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]),df = nu)
  s[2]=qt(p = pAND(q = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]),df = nu)
  s[3]=qt(p = pAND(q = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]),df = nu)
  
  # product of densities from univariate margins
  f=rep(NA,3)
  f[1]=dALaD(y = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
  f[2]=dAND(y = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
  f[3]=dAND(y = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])
  
  # density
  dens=( gamma((nu+3)/2)*gamma(nu/2)^(2) )/( sqrt(det(Sigma))*gamma((nu+1)/2)^3 )*(1+t(s)%*%solve(Sigma)%*%s/nu)^(-(nu+3)/2)/( prod(1+s^2/nu)^(-(nu+1)/2) )*prod(f) 
  
  return(dens)
}

### function for integrating the third margin out of a 3D Student's t-copula
### used for obtaining the bivariate margins
tcop3=function(x,y,nu,Sigma,pars){
  # x is the representative of dimension d, the variable over which to integrate
  # the 3D Student-copula with degrees of freedom nu and covariance Sigma
  # pars are the parameters for the margins (4x3 matrix).
  z=rep(NA,3)
  z[3]=x
  z[c(1,2)]=y
  
  # Student's t-quantile function evaluated in the suitable pseudo observations
  s=rep(NA,3)
  s[1]=qt(p = pALaD(q = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]),df = nu)
  s[2]=qt(p = pAND(q = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]),df = nu)
  s[3]=qt(p = pAND(q = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]),df = nu)
  
  # product of densities from univariate margins
  f=rep(NA,3)
  f[1]=dALaD(y = z[1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
  f[2]=dAND(y = z[2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
  f[3]=dAND(y = z[3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])
  
  # density
  dens=( gamma((nu+3)/2)*gamma(nu/2)^(2) )/( sqrt(det(Sigma))*gamma((nu+1)/2)^3 )*(1+t(s)%*%solve(Sigma)%*%s/nu)^(-(nu+3)/2)/( prod(1+s^2/nu)^(-(nu+1)/2) )*prod(f) 
  
  return(dens)
}

### function for integrating out the first margin
int_tcop1=function(y,nu,Sigma,pars,lb,ub){
  return(trapzfun(f = tcop1,a = lb,b = ub,y=y,pars=pars,Sigma=Sigma,nu=nu)$value)
}

### function for integrating out the second margin
int_tcop2=function(y,nu,Sigma,pars,lb,ub){
  return(trapzfun(f = tcop2,a = lb,b = ub,y=y,pars=pars,Sigma=Sigma,nu=nu)$value)
}

### function for integrating out the third margin
int_tcop3=function(y,nu,Sigma,pars,lb,ub){
  return(trapzfun(f = tcop3,a = lb,b = ub,y=y,pars=pars,Sigma=Sigma,nu=nu)$value)
}



### best fitting parameters for white wine quality 4 ###
########################################################
# Sigma
param=c(0.17182445, -0.09956279,  0.22662979)
Sigma=matrix(1,nrow=3,ncol=3)
Sigma[upper.tri(Sigma)]=param
Sigma[lower.tri(Sigma)]=param
# degrees of freedom
nu=9.53413073
# margin parameters sorted per margin (alpha, mu, phi, degrees of freedom)
parsw4=matrix(c(0.24160029,0.26999942,0.03989715,NA,0.30703972,3.08886916,0.06662737,NA,0.21856928,0.37357839,0.03816755,NA),nrow=4,ncol=3)

# defining the Student's t-copula object
tc=tCopula(param = param,dim = 3,dispstr = "un",df = nu)
# generating data from the Student's t-copula to assess the limits for plotting
U=rCopula(copula = tc,n = 2000)
X=matrix(NA,nrow=2000,ncol=3)
X[,1]=qALaD(beta = U[,1],mu = parsw4[2,1],phi = parsw4[3,1],alpha = parsw4[1,1])
X[,2]=qAND(beta = U[,2],mu = parsw4[2,2],phi = parsw4[3,2],alpha = parsw4[1,2])
X[,3]=qAND(beta = U[,3],mu = parsw4[2,3],phi = parsw4[3,3],alpha = parsw4[1,3])

# number of grid points in each dimension
n=100
x1=seq(0,1,length.out = n) # Volatile Acidity
x2=seq(2.7,4,length.out = n) # pH
x3=seq(0.1,1.5,length.out = n) # Sulphates

# defining the bivariate grids
gr12=expand.grid(x1,x2)
gr13=expand.grid(x1,x3)
gr23=expand.grid(x2,x3)

# bounds for each variable used in integration
lb1=-0.1
ub1=1.3
lb2=2.7
ub2=3.9
lb3=0.1
ub3=1.1

# parallel computing of the bivariate margins to speed up the process
cores=detectCores()
cl <- makeCluster(4)
registerDoParallel(cl)

# for contourplot pH-Sulphates
result=foreach(i=1:(n^2),.packages=c('QBAsyDist','pracma'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 value=int_tcop1(y=as.numeric(gr23[i,]),nu=nu,Sigma=Sigma,pars=parsw4,lb=lb1,ub=ub1)
                 
                 # write output
                 filename=paste0("~/PhD/code/copulas/output/contour/contourww4_yz_run",i,".Rdata")
                 save(value,file=filename)
                 
               }

# Bundling the output
f23=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_yz_run",i,".Rdata"))
    f23[i] <- value
  },silent=T)
}

# for contourplot Volatile Acidity-Sulphates
result=foreach(i=1:(n^2),.packages=c('QBAsyDist','pracma'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 value=int_tcop2(y=as.numeric(gr13[i,]),nu=nu,Sigma=Sigma,pars=parsw4,lb=lb2,ub=ub2)
                 # write output
                 filename=paste0("~/PhD/code/copulas/output/contour/contourww4_xz_run",i,".Rdata")
                 save(value,file=filename)
                 
               }
# Bundling the output
f13=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_xz_run",i,".Rdata"))
    f13[i] <- value
  },silent=T)
}

# for contourplot Volatile Acidity-pH
result=foreach(i=1:(n^2),.packages=c('QBAsyDist','pracma'),
               .combine = "cbind",.verbose = T,.errorhandling="remove") %dopar% {
                 
                 value=int_tcop3(y=as.numeric(gr12[i,]),nu=nu,Sigma=Sigma,pars=parsw4,lb=lb3,ub=ub3)
                 # write output
                 filename=paste0("~/PhD/code/copulas/output/contour/contourww4_xy_run",i,".Rdata")
                 save(value,file=filename)
                 
               }
# Bundling the output
f12=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_xy_run",i,".Rdata"))
    f12[i] <- value
  },silent=T)
}

### generating contourplots
f23=matrix(f23,nrow=n,ncol=n)
x11()
contour(x2,x3,f23,xlab="x2",ylab="x3")

f13=matrix(f13,nrow=n,ncol=n)
x11()
contour(x1,x3,f13,xlab="x1",ylab="x3")

f12=matrix(f12,nrow=n,ncol=n)
x11()
contour(x1,x2,f12,xlab="x1",ylab="x2")






### Fitting and plotting of the wine data investigation ###
###########################################################

# load in the wine data
winequality.red <- read.csv("~/PhD/code/datasets/winequality-red.csv", sep=";") # change directory
winequality.white <- read.csv("~/PhD/code/datasets/winequality-white.csv", sep=";") # change directory

# white wine of quality 4
ww4=winequality.white[winequality.white$quality==4,c(2,9,10)] 
# red wine of quality 4
rw4=winequality.red[winequality.red$quality==4,c(2,9,10)] 
# white wine of quality 7
ww7=winequality.white[winequality.white$quality==7,c(2,9,10)] 
# red wine of quality 7
rw7=winequality.red[winequality.red$quality==7,c(2,9,10)] 

### white whine of quality 4 
#############################

# data
Zw4=ww4 

# dimensions of data
d=length(Zw4[1,])
n=length(Zw4[,1])


### margin fitting
Uw4=matrix(NA,nrow=n,ncol=d)
basefuncw4=rep(NA,d)
parsw4=matrix(NA,nrow=4,ncol=d)
marginloglw4=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Zw4[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  Uw4[,i]=tempout$U
  basefuncw4[i]=tempout$margin
  parsw4[,i]=tempout$parameters
  marginloglw4[i]=tempout$logl
}

m1=bestmargin(data = Zw4[,1],crit = "AIC",nstart = 20,all = T,seed = 1489)
m2=bestmargin(data = Zw4[,2],crit = "AIC",nstart = 20,all = T,seed = 1489)
m3=bestmargin(data = Zw4[,3],crit = "AIC",nstart = 20,all = T,seed = 1489)

### plotting of fitted margins
# Volatile Acidity
hist(Zw4[,1],freq=F)
xx=seq(min(Zw4[,1]),max(Zw4[,1]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m1[1,4]),phi = as.numeric(m1[1,5]),alpha = as.numeric(m1[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m1[2,4]),phi = as.numeric(m1[2,5]),alpha = as.numeric(m1[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m1[3,4]),phi = as.numeric(m1[3,5]),alpha = as.numeric(m1[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m1[4,4]),phi = as.numeric(m1[4,5]),alpha = as.numeric(m1[4,3]),nu = as.numeric(m1[4,6])),col=4)

# pH
hist(Zw4[,2],freq=F)
xx=seq(min(Zw4[,2]),max(Zw4[,2]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m2[1,4]),phi = as.numeric(m2[1,5]),alpha = as.numeric(m2[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m2[2,4]),phi = as.numeric(m2[2,5]),alpha = as.numeric(m2[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m2[3,4]),phi = as.numeric(m2[3,5]),alpha = as.numeric(m2[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m2[4,4]),phi = as.numeric(m2[4,5]),alpha = as.numeric(m2[4,3]),nu = as.numeric(m2[4,6])),col=4)

# Sulphates
hist(Zw4[,3],freq=F)
xx=seq(min(Zw4[,3]),max(Zw4[,3]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m3[1,4]),phi = as.numeric(m3[1,5]),alpha = as.numeric(m3[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m3[2,4]),phi = as.numeric(m3[2,5]),alpha = as.numeric(m3[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m3[3,4]),phi = as.numeric(m3[3,5]),alpha = as.numeric(m3[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m3[4,4]),phi = as.numeric(m3[4,5]),alpha = as.numeric(m3[4,3]),nu = as.numeric(m3[4,6])),col=4)



### copula fitting
tempout=bestcopulafit(U = Uw4,all = T,crit = "AIC")
tempout$info

copulaw4=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulaparsw4=unlist(tempout$parameters[copulaw4])
copulaloglw4=as.numeric(tempout$info[which.min(tempout$info[,3]),2])


### white whine of quality 7 
#############################

# data
Zw7=ww7   

# dimensions of data
d=length(Zw7[1,])
n=length(Zw7[,1])


### margin fitting
Uw7=matrix(NA,nrow=n,ncol=d)
basefuncw7=rep(NA,d)
parsw7=matrix(NA,nrow=4,ncol=d)
marginloglw7=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Zw7[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  Uw7[,i]=tempout$U
  basefuncw7[i]=tempout$margin
  parsw7[,i]=tempout$parameters
  marginloglw7[i]=tempout$logl
}

m1=bestmargin(data = Zw7[,1],crit = "AIC",nstart = 20,all = T,seed = 1489)
m2=bestmargin(data = Zw7[,2],crit = "AIC",nstart = 20,all = T,seed = 1489)
m3=bestmargin(data = Zw7[,3],crit = "AIC",nstart = 20,all = T,seed = 1489)

### plotting of fitted margins
# Volatile Acidity
hist(Zw7[,1],freq=F)
xx=seq(min(Zw7[,1]),max(Zw7[,1]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m1[1,4]),phi = as.numeric(m1[1,5]),alpha = as.numeric(m1[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m1[2,4]),phi = as.numeric(m1[2,5]),alpha = as.numeric(m1[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m1[3,4]),phi = as.numeric(m1[3,5]),alpha = as.numeric(m1[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m1[4,4]),phi = as.numeric(m1[4,5]),alpha = as.numeric(m1[4,3]),nu = as.numeric(m1[4,6])),col=4)

# pH
hist(Zw7[,2],freq=F)
xx=seq(min(Zw7[,2]),max(Zw7[,2]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m2[1,4]),phi = as.numeric(m2[1,5]),alpha = as.numeric(m2[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m2[2,4]),phi = as.numeric(m2[2,5]),alpha = as.numeric(m2[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m2[3,4]),phi = as.numeric(m2[3,5]),alpha = as.numeric(m2[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m2[4,4]),phi = as.numeric(m2[4,5]),alpha = as.numeric(m2[4,3]),nu = as.numeric(m2[4,6])),col=4)

# Sulphates
hist(Zw7[,3],freq=F)
xx=seq(min(Zw7[,3]),max(Zw7[,3]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m3[1,4]),phi = as.numeric(m3[1,5]),alpha = as.numeric(m3[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m3[2,4]),phi = as.numeric(m3[2,5]),alpha = as.numeric(m3[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m3[3,4]),phi = as.numeric(m3[3,5]),alpha = as.numeric(m3[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m3[4,4]),phi = as.numeric(m3[4,5]),alpha = as.numeric(m3[4,3]),nu = as.numeric(m3[4,6])),col=4)



### copula fitting
tempout=bestcopulafit(U = Uw7,all = T,crit = "AIC")
tempout$info

copulaw7=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulaparsw7=unlist(tempout$parameters[copulaw7])
copulaloglw7=as.numeric(tempout$info[which.min(tempout$info[,3]),2])



### red whine of quality 4 
###########################

# data
Zr4=rw4   

# dimensions of data
d=length(Zr4[1,])
n=length(Zr4[,1])


### margin fitting
Ur4=matrix(NA,nrow=n,ncol=d)
basefuncr4=rep(NA,d)
parsr4=matrix(NA,nrow=4,ncol=d)
marginloglr4=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Zr4[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  Ur4[,i]=tempout$U
  basefuncr4[i]=tempout$margin
  parsr4[,i]=tempout$parameters
  marginloglr4[i]=tempout$logl
}

m1=bestmargin(data = Zr4[,1],crit = "AIC",nstart = 20,all = T,seed = 1489)
m2=bestmargin(data = Zr4[,2],crit = "AIC",nstart = 20,all = T,seed = 1489)
m3=bestmargin(data = Zr4[,3],crit = "AIC",nstart = 20,all = T,seed = 1489)

### plotting of fitted margins
# Volatile Acidity
hist(Zr4[,1],freq=F)
xx=seq(min(Zr4[,1]),max(Zr4[,1]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m1[1,4]),phi = as.numeric(m1[1,5]),alpha = as.numeric(m1[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m1[2,4]),phi = as.numeric(m1[2,5]),alpha = as.numeric(m1[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m1[3,4]),phi = as.numeric(m1[3,5]),alpha = as.numeric(m1[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m1[4,4]),phi = as.numeric(m1[4,5]),alpha = as.numeric(m1[4,3]),nu = as.numeric(m1[4,6])),col=4)

# pH
hist(Zr4[,2],freq=F)
xx=seq(min(Zr4[,2]),max(Zr4[,2]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m2[1,4]),phi = as.numeric(m2[1,5]),alpha = as.numeric(m2[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m2[2,4]),phi = as.numeric(m2[2,5]),alpha = as.numeric(m2[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m2[3,4]),phi = as.numeric(m2[3,5]),alpha = as.numeric(m2[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m2[4,4]),phi = as.numeric(m2[4,5]),alpha = as.numeric(m2[4,3]),nu = as.numeric(m2[4,6])),col=4)

# Sulphates
hist(Zr4[,3],freq=F)
xx=seq(min(Zr4[,3]),max(Zr4[,3]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m3[1,4]),phi = as.numeric(m3[1,5]),alpha = as.numeric(m3[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m3[2,4]),phi = as.numeric(m3[2,5]),alpha = as.numeric(m3[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m3[3,4]),phi = as.numeric(m3[3,5]),alpha = as.numeric(m3[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m3[4,4]),phi = as.numeric(m3[4,5]),alpha = as.numeric(m3[4,3]),nu = as.numeric(m3[4,6])),col=4)



### copula fitting
tempout=bestcopulafit(U = Ur4,all = T,crit = "AIC")
tempout$info

copular4=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulaparsr4=unlist(tempout$parameters[copular4])
copulaloglr4=as.numeric(tempout$info[which.min(tempout$info[,3]),2])




### red whine of quality 7 
###########################

# data
Zr7=rw7  

# dimensions of data
d=length(Zr7[1,])
n=length(Zr7[,1])


### margin fitting
Ur7=matrix(NA,nrow=n,ncol=d)
basefuncr7=rep(NA,d)
parsr7=matrix(NA,nrow=4,ncol=d)
marginloglr7=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Zr7[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  Ur7[,i]=tempout$U
  basefuncr7[i]=tempout$margin
  parsr7[,i]=tempout$parameters
  marginloglr7[i]=tempout$logl
}

m1=bestmargin(data = Zr7[,1],crit = "AIC",nstart = 20,all = T,seed = 1489)
m2=bestmargin(data = Zr7[,2],crit = "AIC",nstart = 20,all = T,seed = 1489)
m3=bestmargin(data = Zr7[,3],crit = "AIC",nstart = 20,all = T,seed = 1489)

### plotting of fitted margins
# Volatile Acidity
hist(Zr7[,1],freq=F)
xx=seq(min(Zr7[,1]),max(Zr7[,1]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m1[1,4]),phi = as.numeric(m1[1,5]),alpha = as.numeric(m1[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m1[2,4]),phi = as.numeric(m1[2,5]),alpha = as.numeric(m1[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m1[3,4]),phi = as.numeric(m1[3,5]),alpha = as.numeric(m1[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m1[4,4]),phi = as.numeric(m1[4,5]),alpha = as.numeric(m1[4,3]),nu = as.numeric(m1[4,6])),col=4)

# pH
hist(Zr7[,2],freq=F)
xx=seq(min(Zr7[,2]),max(Zr7[,2]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m2[1,4]),phi = as.numeric(m2[1,5]),alpha = as.numeric(m2[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m2[2,4]),phi = as.numeric(m2[2,5]),alpha = as.numeric(m2[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m2[3,4]),phi = as.numeric(m2[3,5]),alpha = as.numeric(m2[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m2[4,4]),phi = as.numeric(m2[4,5]),alpha = as.numeric(m2[4,3]),nu = as.numeric(m2[4,6])),col=4)

# Sulphates
hist(Zr7[,3],freq=F)
xx=seq(min(Zr7[,3]),max(Zr7[,3]),length.out = 1000)
lines(xx,dAND(y = xx,mu = as.numeric(m3[1,4]),phi = as.numeric(m3[1,5]),alpha = as.numeric(m3[1,3])),col=1)
lines(xx,dALaD(y = xx,mu = as.numeric(m3[2,4]),phi = as.numeric(m3[2,5]),alpha = as.numeric(m3[2,3])),col=2)
lines(xx,dALoD(y = xx,mu = as.numeric(m3[3,4]),phi = as.numeric(m3[3,5]),alpha = as.numeric(m3[3,3])),col=3)
lines(xx,dATD(y = xx,mu = as.numeric(m3[4,4]),phi = as.numeric(m3[4,5]),alpha = as.numeric(m3[4,3]),nu = as.numeric(m3[4,6])),col=4)



### copula fitting
tempout=bestcopulafit(U = Ur7,all = T,crit = "AIC")
tempout$info

copular7=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulaparsr7=unlist(tempout$parameters[copular7])
copulaloglr7=as.numeric(tempout$info[which.min(tempout$info[,3]),2])





### contourplots of the fitted densities for white wine of quality 7 
#####################################################################

# number of grid points in each dimension
m=100
x1=seq(0,1,length.out = m) # Volatile Acidity
x2=seq(2.7,4,length.out = m) # pH
x3=seq(0.1,1.5,length.out = m) # Sulphates

# pseudo observations of the margins based on best fitting one obtained above
u1w7=pAND(q = x1,mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1])
u2w7=pAND(q = x2,mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2])
u3w7=pALoD(q = x3,mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3])

# density of margins based on best fitting one obtained above
d1w7=dAND(y = x1,mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1])
d2w7=dAND(y = x2,mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2])
d3w7=dALoD(y = x3,mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3])

# combining different margin densities/pseudo observations  
p12w7=expand.grid(u1w7,u2w7)
p13w7=expand.grid(u1w7,u3w7)
p23w7=expand.grid(u2w7,u3w7)

d12w7=expand.grid(d1w7,d2w7)
d13w7=expand.grid(d1w7,d3w7)
d23w7=expand.grid(d2w7,d3w7)


# calculating the bivariate marginal densities for the gaussian copula and plotting them
r12w7=copulaparsw7[1]
f12w7=1/sqrt(1-r12w7^2)*exp(-1/(2*(1-r12w7^2))*(r12w7^2*qnorm(p = p12w7[,1])+r12w7^2*qnorm(p = p12w7[,2])-2*r12w7*qnorm(p = p12w7[,1])*qnorm(p = p12w7[,2])))*d12w7[,1]*d12w7[,2]
f12w7=matrix(f12w7,nrow=m,ncol=m)
contour(x = x1,y = x2,z = f12w7)
points(ww7[,c(1,2)],col=2)

r13w7=copulaparsw7[2]
f13w7=1/sqrt(1-r13w7^2)*exp(-1/(2*(1-r13w7^2))*(r13w7^2*qnorm(p = p13w7[,1])+r13w7^2*qnorm(p = p13w7[,2])-2*r13w7*qnorm(p = p13w7[,1])*qnorm(p = p13w7[,2])))*d13w7[,1]*d13w7[,2]
f13w7=matrix(f13w7,nrow=m,ncol=m)
contour(x = x1,y = x3,z = f13w7)
points(ww7[,c(1,3)],col=2)

r23w7=copulaparsw7[3]
f23w7=1/sqrt(1-r23w7^2)*exp(-1/(2*(1-r23w7^2))*(r23w7^2*qnorm(p = p23w7[,1])+r23w7^2*qnorm(p = p23w7[,2])-2*r23w7*qnorm(p = p23w7[,1])*qnorm(p = p23w7[,2])))*d23w7[,1]*d23w7[,2]
f23w7=matrix(f23w7,nrow=m,ncol=m)
contour(x = x2,y = x3,z = f23w7)
points(ww7[,c(2,3)],col=2)


### contourplots of the fitted densities for red wine of quality 7 
###################################################################

x1=seq(0,1,length.out = m) # Volatile Acidity
x2=seq(2.7,4,length.out = m) # pH
x3=seq(0.1,1.5,length.out = m) # Sulphates

# pseudo observations of the margins based on best fitting one obtained above
u1r7=pALaD(q = x1,mu = parsr7[2,1],phi = parsr7[3,1],alpha = parsr7[1,1])
u2r7=pALoD(q = x2,mu = parsr7[2,2],phi = parsr7[3,2],alpha = parsr7[1,2])
u3r7=pALoD(q = x3,mu = parsr7[2,3],phi = parsr7[3,3],alpha = parsr7[1,3])

# density of margins based on best fitting one obtained above
d1r7=dALaD(y = x1,mu = parsr7[2,1],phi = parsr7[3,1],alpha = parsr7[1,1])
d2r7=dALoD(y = x2,mu = parsr7[2,2],phi = parsr7[3,2],alpha = parsr7[1,2])
d3r7=dALoD(y = x3,mu = parsr7[2,3],phi = parsr7[3,3],alpha = parsr7[1,3])

# combining different margin densities/pseudo observations  
p12r7=expand.grid(u1r7,u2r7)
p13r7=expand.grid(u1r7,u3r7)
p23r7=expand.grid(u2r7,u3r7)

d12r7=expand.grid(d1r7,d2r7)
d13r7=expand.grid(d1r7,d3r7)
d23r7=expand.grid(d2r7,d3r7)


# calculating the bivariate marginal densities for the gaussian copula and plotting them
r12r7=copulaparsr7[1]
f12r7=1/sqrt(1-r12r7^2)*exp(-1/(2*(1-r12r7^2))*(r12r7^2*qnorm(p = p12r7[,1])+r12r7^2*qnorm(p = p12r7[,2])-2*r12r7*qnorm(p = p12r7[,1])*qnorm(p = p12r7[,2])))*d12r7[,1]*d12r7[,2]
f12r7=matrix(f12r7,nrow=m,ncol=m)
contour(x = x1,y = x2,z = f12r7)
points(rw7[,c(1,2)],col=2)

r13r7=copulaparsr7[2]
f13r7=1/sqrt(1-r13r7^2)*exp(-1/(2*(1-r13r7^2))*(r13r7^2*qnorm(p = p13r7[,1])+r13r7^2*qnorm(p = p13r7[,2])-2*r13r7*qnorm(p = p13r7[,1])*qnorm(p = p13r7[,2])))*d13r7[,1]*d13r7[,2]
f13r7=matrix(f13r7,nrow=m,ncol=m)
contour(x = x1,y = x3,z = f13r7)
points(rw7[,c(1,3)],col=2)

r23r7=copulaparsr7[3]
f23r7=1/sqrt(1-r23r7^2)*exp(-1/(2*(1-r23r7^2))*(r23r7^2*qnorm(p = p23r7[,1])+r23r7^2*qnorm(p = p23r7[,2])-2*r23r7*qnorm(p = p23r7[,1])*qnorm(p = p23r7[,2])))*d23r7[,1]*d23r7[,2]
f23r7=matrix(f23r7,nrow=m,ncol=m)
contour(x = x2,y = x3,z = f23r7)
points(rw7[,c(2,3)],col=2)


### contourplots of the fitted densities for red wine of quality 4 
###################################################################

x1=seq(0,1.2,length.out = m) # Volatile Acidity
x2=seq(2.7,4,length.out = m) # pH
x3=seq(0.1,1.5,length.out = m) # Sulphates

# pseudo observations of the margins based on best fitting one obtained above
u1r4=pAND(q = x1,mu = parsr4[2,1],phi = parsr4[3,1],alpha = parsr4[1,1])
u2r4=pALaD(q = x2,mu = parsr4[2,2],phi = parsr4[3,2],alpha = parsr4[1,2])
u3r4=pATD(q = x3,mu = parsr4[2,3],phi = parsr4[3,3],alpha = parsr4[1,3],nu = parsr4[4,3])

# density of margins based on best fitting one obtained above
d1r4=dAND(y = x1,mu = parsr4[2,1],phi = parsr4[3,1],alpha = parsr4[1,1])
d2r4=dALaD(y = x2,mu = parsr4[2,2],phi = parsr4[3,2],alpha = parsr4[1,2])
d3r4=dATD(y = x3,mu = parsr4[2,3],phi = parsr4[3,3],alpha = parsr4[1,3],nu = parsr4[4,3])

# combining different margin densities/pseudo observations  
p12r4=expand.grid(u1r4,u2r4)
p13r4=expand.grid(u1r4,u3r4)
p23r4=expand.grid(u2r4,u3r4)

d12r4=expand.grid(d1r4,d2r4)
d13r4=expand.grid(d1r4,d3r4)
d23r4=expand.grid(d2r4,d3r4)

# calculating the bivariate marginal densities for the gaussian copula and plotting them
r12r4=copulaparsr4[1]
f12r4=1/sqrt(1-r12r4^2)*exp(-1/(2*(1-r12r4^2))*(r12r4^2*qnorm(p = p12r4[,1])+r12r4^2*qnorm(p = p12r4[,2])-2*r12r4*qnorm(p = p12r4[,1])*qnorm(p = p12r4[,2])))*d12r4[,1]*d12r4[,2]
f12r4=matrix(f12r4,nrow=m,ncol=m)
contour(x = x1,y = x2,z = f12r4)
points(rw4[,c(1,2)],col=2)

r13r4=copulaparsr4[2]
f13r4=1/sqrt(1-r13r4^2)*exp(-1/(2*(1-r13r4^2))*(r13r4^2*qnorm(p = p13r4[,1])+r13r4^2*qnorm(p = p13r4[,2])-2*r13r4*qnorm(p = p13r4[,1])*qnorm(p = p13r4[,2])))*d13r4[,1]*d13r4[,2]
f13r4=matrix(f13r4,nrow=m,ncol=m)
contour(x = x1,y = x3,z = f13r4)
points(rw4[,c(1,3)],col=2)

r23r4=copulaparsr7[3]
f23r4=1/sqrt(1-r23r4^2)*exp(-1/(2*(1-r23r4^2))*(r23r4^2*qnorm(p = p23r4[,1])+r23r4^2*qnorm(p = p23r4[,2])-2*r23r4*qnorm(p = p23r4[,1])*qnorm(p = p23r4[,2])))*d23r4[,1]*d23r4[,2]
f23r4=matrix(f23r4,nrow=m,ncol=m)
contour(x = x2,y = x3,z = f23r4)
points(rw4[,c(2,3)],col=2)


### for white wine quality 4, see separate file. 
### Saving these datasets for the contours is used in the 
### file plots paper copulas


### comparison red wine with white wine of quality 7 
#####################################################

x1=seq(0,1,length.out = m) # Volatile Acidity
x2=seq(2.7,4,length.out = m) # pH
x3=seq(0.1,1.5,length.out = m) # Sulphates


x11()
plot(x1,d1w7,type="l",col=1,xlab="Volatile Acidity",ylab="Density",cex.axis=1.5,cex.lab=1.5,lwd=2,ylim=c(-0.2,5.5))
lines(x1,d1r7,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)

x11()
plot(x2,d2w7,type="l",col=1,xlab="pH",ylab="Density",cex.axis=1.5,cex.lab=1.5,lwd=2,ylim=c(-0.2,3.2))
lines(x2,d2r7,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)

x11()
plot(x3,d3w7,type="l",col=1,xlab="Sulphates",ylab="Density",cex.axis=1.5,cex.lab=1.5,lwd=2,ylim=c(-0.2,4))
lines(x3,d3r7,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)


x11()
contour(x2,x1,t(f12w7),levels=c(0.25,0.5,0.75,1,2,4,6,9),lwd=2,ylab="Volatile Acidity",xlab="pH",cex.lab=1.5,cex.axis=1.5)
contour(x2,x1,t(f12r7),levels=c(0.25,0.5,0.75,1,2,4,6,9),add=T,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)

x11()
contour(x3,x1,t(f13w7),levels=c(0.25,0.5,0.75,1,2,4,6,9),lwd=2,ylab="Volatile Acidity",xlab="Sulphates",cex.lab=1.5,cex.axis=1.5)
contour(x3,x1,t(f13r7),levels=c(0.25,0.5,0.75,1,2,4,6,9),add=T,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)

x11()
contour(x3,x2,t(f23w7),levels=c(0.25,0.5,0.75,1,2,4,6,9),lwd=2,ylab="pH",xlab="Sulphates",cex.lab=1.5,cex.axis=1.5)
contour(x3,x2,t(f23r7),levels=c(0.25,0.5,0.75,1,2,4,6,9),add=T,col=2,lwd=2)
#legend("topright",legend = c("White wine of quality 7","Red wine of quality 7"),col=c(1,2),lwd=c(2,2),lty=c(1,1),cex=1.5)







### 3D plots for the wine data ###
##################################

# range of the data for VA (xr), pH (yr) and S (zr) 
xr=range(ww4[,1],ww7[,1],rw4[,1],rw7[,1])
yr=range(ww4[,2],ww7[,2],rw4[,2],rw7[,2])
zr=range(ww4[,3],ww7[,3],rw4[,3],rw7[,3])

# grid on which to determine the density 
n=50 # number of grid-points in each dimension
x1=seq(xr[1],xr[2],length.out=n)
x2=seq(yr[1],yr[2],length.out=n)
x3=seq(zr[1],1.2,length.out=n)
points=expand.grid(x1,x2,x3)

# number of levels for the 3D contour plots and values of the density at these levels
nlevel=7
level=seq(8,36,length.out = nlevel)

# limits of the plotting area
xlim=c(0.1,1.1)
ylim=c(2.9,3.9)
zlim=c(0.2,0.8)



### white wine quality 4
#########################

# copula parameters from the best fitting model (Student's t-copula)
param=c(0.17182445, -0.09956279,  0.22662979)
Sigma=matrix(1,nrow=3,ncol=3)
Sigma[upper.tri(Sigma)]=param
Sigma[lower.tri(Sigma)]=param
nu=9.53413073

# margin parameters ordened by margin consisting of alpha, mu, phi, nu
pars=matrix(c(0.24160029,0.26999942,0.03989715,NA,0.30703972,3.08886916,0.06662737,NA,0.21856928,0.37357839,0.03816755,NA),nrow=4,ncol=3)


# Student's t-quantile function evaluated in margins CDF for use in the copula density
s=matrix(NA,nrow=n^3,ncol=3)
s[,1]=qt(p = pALaD(q = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]),df = nu)
s[,2]=qt(p = pAND(q = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]),df = nu)
s[,3]=qt(p = pAND(q = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]),df = nu)

# densities of univariate margins
f=matrix(NA,nrow=n^3,ncol=3)
f[,1]=dALaD(y = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
f[,2]=dAND(y = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
f[,3]=dAND(y = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])

# full density calculation
dfunc=function(X,nu,Sigma){
  f=X[1:3]
  s=X[4:6]
  return(( gamma((nu+3)/2)*gamma(nu/2)^(2) )/( sqrt(det(Sigma))*gamma((nu+1)/2)^3 )*(1+t(s)%*%solve(Sigma)%*%s/nu)^(-(nu+3)/2)/( prod(1+s^2/nu)^(-(nu+1)/2) )*prod(f) )
}
X=cbind(f,s)
dens=apply(X = X,MARGIN = 1,FUN = dfunc,nu = nu,Sigma = Sigma)
value=array(dens,dim=c(n,n,n))


# 3D plot of density
m <- array(TRUE,dim = c(n,n,n))
m[x1>pars[2,1],x2<pars[2,2],]=FALSE
v <- contour3d(value,level=level,x1,x2,x3,mask=m,color=jet.col(n = nlevel),alpha=seq(0.4,0.9,length.out = nlevel),draw = F)
x11()
par(mai = c(1.02,0.82,0.82,0.82))
M <- persp(xlim, ylim, matrix(zlim[1], 2, 2), theta = 15, phi = 30,
           col = "lightgray", zlim = zlim, ticktype = "detailed",
           scale = FALSE, d = 4, xlab = "Volatile Acidity",
           ylab = "pH", zlab = "Sulphates")
colkey(col=jet.col(),clim = c(8,36) ,add=TRUE,clab="Density")
drawScene(v, screen = NULL, R.mat = t(M), add = TRUE, scale = FALSE,
          light = c(.5, 0, 1))




### white wine quality 7
#########################

# copula parameters from the best fitting model (Gaussian copula)
param=c(0.05346992, -0.04266732,  0.18098453)
Sigma=matrix(1,nrow=3,ncol=3)
Sigma[upper.tri(Sigma)]=param
Sigma[lower.tri(Sigma)]=param

# margin parameters ordened by margin consisting of alpha, mu, phi, nu
pars=matrix(c(0.23329731, 0.19050495, 0.03054464, NA, 0.40117890, 3.16364162, 0.07550313, NA, 0.24414791, 0.40681863, 0.02500304, NA),nrow=4,ncol=3)


# Gaussian quantile function evaluated in margins CDF for use in the copula density
s=matrix(NA,nrow=n^3,ncol=3)
s[,1]=qnorm(p = pAND(q = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]))
s[,2]=qnorm(p = pAND(q = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]))
s[,3]=qnorm(p = pALoD(q = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]))


# densities of univariate margins
f=matrix(NA,nrow=n^3,ncol=3)
f[,1]=dAND(y = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
f[,2]=dAND(y = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
f[,3]=dALoD(y = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])

# full density calculation
dfunc=function(X,Sigma){
  f=X[1:3]
  s=X[4:6]
  return(1/sqrt(abs(det(Sigma)))*exp(-1/2*t(s)%*%(solve(Sigma)-diag(3))%*%s)*prod(f) )
}
X=cbind(f,s)
dens=apply(X = X,MARGIN = 1,FUN = dfunc,Sigma = Sigma)
value=array(dens,dim=c(n,n,n))


# 3D plot of density
m <- array(TRUE,dim = c(n,n,n))
m[x1>pars[2,1],x2<pars[2,2],]=FALSE
v <- contour3d(value,level=level,x1,x2,x3,mask=m,color=jet.col(n = nlevel),alpha=seq(0.4,0.9,length.out = nlevel),draw = F)
x11()
par(mai = c(1.02,0.82,0.82,0.82))
M <- persp(xlim, ylim, matrix(zlim[1], 2, 2), theta = 15, phi = 30,
           col = "lightgray", zlim = zlim, ticktype = "detailed",
           scale = FALSE, d = 4, xlab = "Volatile Acidity",
           ylab = "pH", zlab = "Sulphates")
colkey(col=jet.col(),clim = c(8,36) ,add=TRUE,clab="Density")
drawScene(v, screen = NULL, R.mat = t(M), add = TRUE, scale = FALSE,
          light = c(.5, 0, 1))


### red wine quality 4 
#######################

# copula parameters from the best fitting model (Gaussian copula)
param=c(0.3414465, -0.2071128, -0.1849579)
Sigma=matrix(1,nrow=3,ncol=3)
Sigma[upper.tri(Sigma)]=param
Sigma[lower.tri(Sigma)]=param

# margin parameters ordened by margin consisting of alpha, mu, phi, nu
pars=matrix(c(0.39866438, 0.62305269, 0.10371903, NA, 0.41472991, 3.33999858, 0.06063362, NA, 0.33003812, 0.51175935, 0.02937946, 2.00000000),nrow=4,ncol=3)


# Gaussian quantile function evaluated in margins CDF for use in the copula density
s=matrix(NA,nrow=n^3,ncol=3)
s[,1]=qnorm(p = pAND(q = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]))
s[,2]=qnorm(p = pALaD(q = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]))
s[,3]=qnorm(p = pATD(q = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3],nu = pars[4,3]))

# densities of univariate margins
f=matrix(NA,nrow=n^3,ncol=3)
f[,1]=dAND(y = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
f[,2]=dALaD(y = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
f[,3]=dATD(y = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3],nu = pars[4,3])


# full density calculation
dfunc=function(X,Sigma){
  f=X[1:3]
  s=X[4:6]
  return(1/sqrt(abs(det(Sigma)))*exp(-1/2*t(s)%*%(solve(Sigma)-diag(3))%*%s)*prod(f) )
}
X=cbind(f,s)
dens=apply(X = X,MARGIN = 1,FUN = dfunc,Sigma = Sigma)
value=array(dens,dim=c(n,n,n))


# 3D plot of density
m <- array(TRUE,dim = c(n,n,n))
m[x1>pars[2,1],x2<pars[2,2],]=FALSE
v <- contour3d(value,level=level,x1,x2,x3,mask=m,color=jet.col(n = nlevel),alpha=seq(0.4,0.9,length.out = nlevel),draw = F)
x11()
par(mai = c(1.02,0.82,0.82,0.82))
M <- persp(xlim, ylim, matrix(zlim[1], 2, 2), theta = 15, phi = 30,
           col = "lightgray", zlim = zlim, ticktype = "detailed",
           scale = FALSE, d = 4, xlab = "Volatile Acidity",
           ylab = "pH", zlab = "Sulphates")
colkey(col=jet.col(),clim = c(8,36) ,add=TRUE,clab="Density")
drawScene(v, screen = NULL, R.mat = t(M), add = TRUE, scale = FALSE,
          light = c(.5, 0, 1))




### red wine quality 7 
#######################

# copula parameters from the best fitting model (Gaussian copula)
param=c(0.26637875 -0.23561507 -0.02660246)
Sigma=matrix(1,nrow=3,ncol=3)
Sigma[upper.tri(Sigma)]=param
Sigma[lower.tri(Sigma)]=param

# margin parameters ordened by margin consisting of alpha, mu, phi, nu
pars=matrix(c(0.23603460, 0.29999975, 0.03565239, NA, 0.44533251, 3.26529219, 0.04079342, NA, 0.45843424, 0.72300644, 0.03670445, NA),nrow=4,ncol=3)


# Gaussian quantile function evaluated in margins CDF for use in the copula density
s=matrix(NA,nrow=n^3,ncol=3)
s[,1]=qnorm(p = pALaD(q = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1]))
s[,2]=qnorm(p = pALoD(q = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2]))
s[,3]=qnorm(p = pALoD(q = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3]))

# densities of univariate margins
f=matrix(NA,nrow=n^3,ncol=3)
f[,1]=dALaD(y = points[,1],mu = pars[2,1],phi = pars[3,1],alpha = pars[1,1])
f[,2]=dALoD(y = points[,2],mu = pars[2,2],phi = pars[3,2],alpha = pars[1,2])
f[,3]=dALoD(y = points[,3],mu = pars[2,3],phi = pars[3,3],alpha = pars[1,3])


# full density calculation
dfunc=function(X,Sigma){
  f=X[1:3]
  s=X[4:6]
  return(1/sqrt(abs(det(Sigma)))*exp(-1/2*t(s)%*%(solve(Sigma)-diag(3))%*%s)*prod(f) )
}
X=cbind(f,s)
dens=apply(X = X,MARGIN = 1,FUN = dfunc,Sigma = Sigma)
value=array(dens,dim=c(n,n,n))


# 3D plot of density
m <- array(TRUE,dim = c(n,n,n))
m[x1>pars[2,1],x2<pars[2,2],]=FALSE
v <- contour3d(value,level=level,x1,x2,x3,mask=m,color=jet.col(n = nlevel),alpha=seq(0.4,0.9,length.out = nlevel),draw = F)
x11()
par(mai = c(1.02,0.82,0.82,0.82))
M <- persp(xlim, ylim, matrix(zlim[1], 2, 2), theta = 15, phi = 30,
           col = "lightgray", zlim = zlim, ticktype = "detailed",
           scale = FALSE, d = 4, xlab = "Volatile Acidity",
           ylab = "pH", zlab = "Sulphates")
colkey(col=jet.col(),clim = c(8,36) ,add=TRUE,clab="Density")
drawScene(v, screen = NULL, R.mat = t(M), add = TRUE, scale = FALSE,
          light = c(.5, 0, 1))





### comparing contour plots wine data ###
#########################################

# setting up a grid. x1 (VA), x2 (pH) and x3 (S) and bivariate combinations of both
n=100 # grid points in each dimension
x1=seq(0,1,length.out = n)
x2=seq(2.7,4,length.out = n)
x3=seq(0.1,1.5,length.out = n)
gr12=expand.grid(x1,x2)
gr13=expand.grid(x1,x3)
gr23=expand.grid(x2,x3)



### white wine quality 4

# bivariate margins loaded from earlier computation
# VA-pH
f12w4=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_xy_run",i,".Rdata"))
    f12w4[i] <- value
  },silent=T)
}
f12w4=matrix(f12w4,nrow=n,ncol=n)

# VA-S
f13w4=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_xz_run",i,".Rdata"))
    f13w4[i] <- value
  },silent=T)
}
f13w4=matrix(f13w4,nrow=n,ncol=n)

# pH-S
f23w4=c()
for (i in 1:(n^2)){
  try({
    load(paste0("~/PhD/code/copulas/output/contour/contourww4_yz_run",i,".Rdata"))
    f23w4[i] <- value
  },silent=T)
}
f23w4=matrix(f23w4,nrow=n,ncol=n)



# univariate margins
# margin parameters ordened per margin
parsw4=matrix(c(0.24160029,0.26999942,0.03989715,NA,0.30703972,3.08886916,0.06662737,NA,0.21856928,0.37357839,0.03816755,NA),nrow=4,ncol=3)

fw4=matrix(NA,nrow=n,ncol=3)
fw4[,1]=dALaD(y = x1,mu = parsw4[2,1],phi = parsw4[3,1],alpha = parsw4[1,1])
fw4[,2]=dAND(y = x2,mu = parsw4[2,2],phi = parsw4[3,2],alpha = parsw4[1,2])
fw4[,3]=dAND(y = x3,mu = parsw4[2,3],phi = parsw4[3,3],alpha = parsw4[1,3])



### white wine quality 7

# bivariate margins from earlier computation
load("~/PhD/code/copulas/output/contour/contourww7_xy.Rdata")
load("~/PhD/code/copulas/output/contour/contourww7_xz.Rdata")
load("~/PhD/code/copulas/output/contour/contourww7_yz.Rdata")

# univariate margins
# margin parameters ordened per margin
parsw7=matrix(c(0.23329731, 0.19050495, 0.03054464, NA, 0.40117890, 3.16364162, 0.07550313, NA, 0.24414791, 0.40681863, 0.02500304, NA),nrow=4,ncol=3)
fw7=matrix(NA,nrow=n,ncol=3)
fw7[,1]=dAND(y = x1,mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1])
fw7[,2]=dAND(y = x2,mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2])
fw7[,3]=dALoD(y = x3,mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3])



### red wine quality 4

# bivariate margins from earlier computation
load("~/PhD/code/copulas/output/contour/contourrw4_xy.Rdata")
load("~/PhD/code/copulas/output/contour/contourrw4_xz.Rdata")
load("~/PhD/code/copulas/output/contour/contourrw4_yz.Rdata")

# univariate margins
# margin parameters ordened per margin
parsr4=matrix(c(0.39866438, 0.62305269, 0.10371903, NA, 0.41472991, 3.33999858, 0.06063362, NA, 0.33003812, 0.51175935, 0.02937946, 2.00000000),nrow=4,ncol=3)
fr4=matrix(NA,nrow=n,ncol=3)
fr4[,1]=dAND(y = x1,mu = parsr4[2,1],phi = parsr4[3,1],alpha = parsr4[1,1])
fr4[,2]=dALaD(y = x2,mu = parsr4[2,2],phi = parsr4[3,2],alpha = parsr4[1,2])
fr4[,3]=dATD(y = x3,mu = parsr4[2,3],phi = parsr4[3,3],alpha = parsr4[1,3],nu = parsr4[4,3])


### red wine quality 7

# bivariate margins from earlier computation
load("~/PhD/code/copulas/output/contour/contourrw7_xy.Rdata")
load("~/PhD/code/copulas/output/contour/contourrw7_xz.Rdata")
load("~/PhD/code/copulas/output/contour/contourrw7_yz.Rdata") 

# univariate margins
# margin parameters ordened per margin
parsr7=matrix(c(0.23603460, 0.29999975, 0.03565239, NA, 0.44533251, 3.26529219, 0.04079342, NA, 0.45843424, 0.72300644, 0.03670445, NA),nrow=4,ncol=3)
fr7=matrix(NA,nrow=n,ncol=3)
fr7[,1]=dALaD(y = x1,mu = parsr7[2,1],phi = parsr7[3,1],alpha = parsr7[1,1])
fr7[,2]=dALoD(y = x2,mu = parsr7[2,2],phi = parsr7[3,2],alpha = parsr7[1,2])
fr7[,3]=dALoD(y = x3,mu = parsr7[2,3],phi = parsr7[3,3],alpha = parsr7[1,3])


# density plot volatile acidity
VA=rbind(cbind(x1,1,fw4[,1]),cbind(x1,2,fw7[,1]),cbind(x1,3,fr4[,1]),cbind(x1,4,fr7[,1]))
VA.d=data.frame(VA)
colnames(VA.d)=c("x","type","density")
m1 <- ggplot(data = VA.d,aes(x = x,y = density,group=factor(type),color=factor(type))) +
  geom_line(lwd=1) +
  xlab("Volatile Acidity") +
  ylab("Density") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quality 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="bottom",legend.text = element_text(size=14)) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) 
legend=get_legend(m1)
m1 <- m1 + theme(legend.position="none")

# density plot pH
pH=rbind(cbind(x2,1,fw4[,2]),cbind(x2,2,fw7[,2]),cbind(x2,3,fr4[,2]),cbind(x2,4,fr7[,2]))
pH.d=data.frame(pH)
colnames(pH.d)=c("x","type","density")
m2 <- ggplot(data = pH.d,aes(x = x,y = density,group=factor(type),color=factor(type))) +
  geom_line(lwd=1) +
  xlab("pH") +
  ylab("Density") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quatity 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

# density plot Sulphates
S=rbind(cbind(x3,1,fw4[,3]),cbind(x3,2,fw7[,3]),cbind(x3,3,fr4[,3]),cbind(x3,4,fr7[,3]))
S.d=data.frame(S)
colnames(S.d)=c("x","type","density")
m3 <- ggplot(data = S.d,aes(x = x,y = density,group=factor(type),color=factor(type))) +
  geom_line(lwd=1) +
  xlab("Sulphates") +
  ylab("Density") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quatity 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

# contourplot Volatile Acidity - pH
X12=rbind(cbind(as.matrix(gr12),1,c(f12w4)),cbind(as.matrix(gr12),2,c(f12w7)),
          cbind(as.matrix(gr12),3,c(f12r4)),cbind(as.matrix(gr12),4,c(f12r7)))
colnames(X12)=c("x","y","type","density")
X12.d=data.frame(X12)
c12 <- ggplot(data = X12.d,aes(x=x,y=y,z=density,group=factor(type),color=factor(type))) +
  geom_contour(lwd=1,breaks=c(2,4,6,8)) +
  xlab("Volatile Acidity") + 
  ylab("pH") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quatity 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


# contourplot Volatile Acidity - Sulphates
X13=rbind(cbind(as.matrix(gr13),1,c(f13w4)),cbind(as.matrix(gr13),2,c(f13w7)),
          cbind(as.matrix(gr13),3,c(f13r4)),cbind(as.matrix(gr13),4,c(f13r7)))
colnames(X13)=c("x","y","type","density")
X13.d=data.frame(X13)
c13 <- ggplot(data = X13.d,aes(x=x,y=y,z=density,group=factor(type),color=factor(type))) +
  geom_contour(lwd=1,breaks=c(2,4,6,8)) +
  xlab("Volatile Acidity") + 
  ylab("Sulphates") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quatity 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))


# contourplot Volatile pH - Sulphates 
X23=rbind(cbind(as.matrix(gr23),1,c(f23w4)),cbind(as.matrix(gr23),2,c(f23w7)),
          cbind(as.matrix(gr23),3,c(f23r4)),cbind(as.matrix(gr23),4,c(f23r7)))
colnames(X23)=c("x","y","type","density")
X23.d=data.frame(X23)
c23 <- ggplot(data = X23.d,aes(x=x,y=y,z=density,group=factor(type),color=factor(type))) +
  geom_contour(lwd=1,breaks=c(2,4,6,8)) +
  xlab("pH") + 
  ylab("Sulphates") +
  scale_color_manual(labels = c("1"="White quality 4","2"="White quatity 7","3"="Red quality 4","4"="Red quality 7"),
                     values = c("black","red","green","blue"), name="") + 
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))

x11()
grid.arrange(m1,m2,m3,c12,c13,c23,legend,nrow=3,ncol=3,layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,7,7)),heights=c(1.2,1.2,0.2))



#########################################
### White wine data with added jitter ###
#########################################

### wine data ### 
#################
winequality.white <- read.csv("~/PhD/code/datasets/winequality-white.csv", sep=";")
whitewinedata=winequality.white[winequality.white$quality==7,c(2,9,10)]

Z=whitewinedata
table(table(Z))

# dimensions of the data
d=length(Z[1,])
n=length(Z[,1])

# adding uniform jitter to the data
set.seed(124578)
Z$volatile.acidity=Z$volatile.acidity+runif(n,min = -0.0025,max = 0.0025)
Z$pH=Z$pH+runif(n,min = -0.005,max = 0.005)
Z$sulphates=Z$sulphates+runif(n = n,min = -0.005,max = 0.005)

pairs(Z)

### margin fitting
U=matrix(NA,nrow=n,ncol=d)
basefunc=rep(NA,d)
pars=matrix(NA,nrow=4,ncol=d)
marginlogl=rep(NA,d)

for(i in 1:d){
  tempout=bestmargin(data = Z[,i],crit = "AIC",nstart = 10,all = F,seed = 1489)
  U[,i]=tempout$U
  basefunc[i]=tempout$margin
  pars[,i]=tempout$parameters
  marginlogl[i]=tempout$logl
}


### copula fitting
tempout=bestcopulafit(U = U,all = T,crit = "AIC")
tempout$info

copula=tempout$info[which.min(tempout$info[,3]),1] # 3 for AIC, 4 for BIC
copulapars=unlist(tempout$parameters[copula])
copulalogl=as.numeric(tempout$info[which.min(tempout$info[,3]),2])


### CIC penalty for margins and copula
penalty=CICpenalty(Z = Z,margfuncs = basefunc,alpha = pars[1,],mu = pars[2,],phi = pars[3,],df = pars[4,],
                   copula = "gumbel",theta = tempout$parameters$Gumbel[1])#,nu = tempout$parameters$t[4])
penalty

# copula contribution
-2*as.numeric(tempout$info[1,2])+2*penalty$copCIC    # 1=gumbel, 2=clayton, 3=frank, 4=joe, 5=normal, 6=t


# bootstrapped KS test


library(QBAsyDist)
library(goftest)
library(snpar)
source("~/PhD/code/copulas/bestmarginfit.R")

### loading in the wine data. Can be downloaded from 
ww7= Z

# dimensions of the data
d=length(ww7[1,])
n=length(ww7[,1])


### holding matrices for margins and their parameters
basefuncw7=rep(NA,d)
parsw7=matrix(NA,nrow=4,ncol=d)

### determining the best fitting QBA margin
for(i in 1:d){
  tempout=bestmargin(data = ww7[,i],crit = "AIC",nstart = 20,all = F,seed = 1489)
  basefuncw7[i]=tempout$margin
  parsw7[,i]=tempout$parameters
}

### Kolmogorov-Smirnov test on the fitted margins
KSSW1=ks.test(x=ww7[,1],y = "pAND",alternative = "two.sided",mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1])
KSSW2=KS.test(x=ww7[,2],y = "pAND",alternative = "two.sided",mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2])
KSSW3=KS.test(x=ww7[,3],y = "pALoD",alternative = "two.sided",mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3])

### number of replicates
N=1000
start.time <- Sys.time()

### First variable: Volatile Acidity
# test statistic from the original data
W1.ks.teststat=KSSW1$statistic


# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = parsw7[2,1],phi = parsw7[3,1],alpha = parsw7[1,1]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W1.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=ks.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided")
  W1.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W1.ks.stats)
abline(v=W1.ks.teststat)
# approximate p-value
W1.ks.pvalue=sum(W1.ks.stats>W1.ks.teststat)/N


### Second variable: pH
# test statistic from the original data
W2.ks.teststat=KSSW2$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rAND(n = N*n,mu = parsw7[2,2],phi = parsw7[3,2],alpha = parsw7[1,2]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W2.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitAND(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=ks.test(x = bootstrapsamples[i,],y = "pAND",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided")
  W2.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W2.ks.stats)
abline(v=W2.ks.teststat)
# approximate p-value
W2.ks.pvalue=sum(W2.ks.stats>W2.ks.teststat)/N


### Third variable: Sulphates
# test statistic from the original data
W3.ks.teststat=KSSW3$statistic

# seed for sample generation
set.seed(147)
# generate samples of the same size as the data from the fitted distribution
bootstrapsamples=matrix(rALoD(n = N*n,mu = parsw7[2,3],phi = parsw7[3,3],alpha = parsw7[1,3]),nrow=N,ncol=n)
# holding vector for teststatistics of the bootstrap samples
W3.ks.stats=rep(NA,N)

# fitting
seed=148
for(i in 1:N){
  fit=fitALoD(data = bootstrapsamples[i,],start = NULL,nstart = 20,seed = seed+i)
  t1=ks.test(x = bootstrapsamples[i,],y = "pALoD",mu = fit$mu,phi = fit$phi,alpha = fit$alpha,alternative = "two.sided")
  W3.ks.stats[i]=t1$statistic
  print(i)
}

# plot of sample teststatistics and data version
x11()
hist(W3.ks.stats)
abline(v=W3.ks.teststat)
# approximate p-value
W3.ks.pvalue=sum(W3.ks.stats>W3.ks.teststat)/N

# time-keeping
stop.time <- Sys.time()
time.taken = difftime(stop.time,start.time,units = "secs")
