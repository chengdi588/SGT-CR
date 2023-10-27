### ======================================================================= ###
###                             APPLICATION                                 ###  
### ----------------------------------------------------------------------- ###
###                Linear censored regression models with                   ###
###              skew scale mixtures of normal distributions                ###
###  Authors: Cl?cio S. Ferreira, Camila B. Zeller and Daniel C. F. Guzman  ###
### ======================================================================= ###

rm(list=ls(all=TRUE))
source('functions.r') 

#######################################################################################################################
# 1. Artificial data ---- simulation   

beta0=matrix(c(5,1))
p=length(beta0)
sigma2=1 
lambda=3

n=50
x=runif(n,0,10)
X=matrix(1,n,2)
X[,2]=x
mu=X%*%beta0
nu=3


erro=rstn(n, location=0, scale=sigma2, shape=lambda,nu=nu) # StN model 
y=mu+erro 
# incluindo censura
perc     <- 0.10 #level censored
aa       <- sort(y,decreasing=TRUE)
cutof    <- aa[ceiling(perc*n)]
cc       <- matrix(1,n,1)*(y>=cutof)  # right censoring
#cc      <- matrix(1,n,1)*(y<cutof) #Left censoring
y[cc==1] <- cutof
# encontrando as posicoes com censura
aux=1:n
ind=aux[cc==1]
cens="right"
cutof=cutof*rep(1,length(ind))# cutof shoul be a vector; if ci=c, so cutof=c*rep(1,length(ind))

theta=stn.ecmeC(y,X,ind,cens,cutof)
theta

theta1=ssl.ecmeC(y,X,ind,cens,cutof)
theta1 # 此处 2 min   
time2 <- Sys.time()  



########################################################################################################################

# 2. Real Data: Stellar abundances dataset  ------ real data  

rm(list=ls(all=TRUE))    
source('functions.r')  

dado <-read.table("D:/R_code_2023.4/7.1_code/SSMN-CR-master/SSMN-CR-master/stellar.rData",header=T) 
attach(dado)

# A table containing 68 rows and 8 columns with header row. Column 2 is a sample indicator (1 = star with planet; 2 = star without planet). Columns 4 and 7 are censoring indicators (1 = detected; 0 = undetected, left-censored).
attach(censor_Be) 
y             <- logN_Be 
n=length(y)
x             <- Teff/1000
X             <- cbind(1,x)
p=2 

ccaux=Ind_Be   # conforme o conjunto de dados censura =0
# encontrando as posicoes com censura
aux=1:n
ind=aux[ccaux==0]  # == small1;    

cc=rep(0,n)
cc[ind]=1

#cens="right"
cens="left"
cutof=y[ind]     # cutof shoul be a vector; if ci=c, so cutof=c*rep(1,length(ind))

##################### Fit of SSMN-CR
# SN
thetaSN=sn.emC(y,X,ind,cens,cutof)   
thetaSN

# StN

theta=stn.ecmeC(y,X,ind,cens,cutof)
theta
# log0 = -1.7803;  −1.8713, 0.5243, 0.0333, -1.9058, 2.0101    

# Skew-Slash

theta=ssl.ecmeC(y,X,ind,cens,cutof)
theta  
   

## Contaminated Normal
theta=scn.ecmeC(y,X,ind,cens,cutof)
theta
time2 <- Sys.time()   

# 1min 

plot( theta$like , type="l")  
# like 先迅速增大，后面稍微下降一点，再缓慢增大； 但整体呈现递增；
# 整体迭代 几百次； 
 
