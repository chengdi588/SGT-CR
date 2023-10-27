require(mnormt)
require(moments)
require(sn)

# Para Slash
qtgamma <-function (p, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(p >= 0 & p <= 1) & all(scale > 0) & all(shape > 
        0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    pNew <- p * (Fb - Fa) + Fa
    x <- qgamma(pNew, shape, scale = scale)
    return(x)
}

dtgamma <-function (x, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(shape > 0) & all(scale > 0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    y <- dgamma(x, shape, scale = scale)
    inda <- which(x < a)
    indb <- which(x > b)
    if (length(inda) > 0) 
        y[inda] <- 0
    if (length(indb) > 0) 
        y[indb] <- 0
    return(y/(Fb - Fa))
}


elogu<-function(u,nu,d){
 log(u)*dtgamma(u, shape=nu+1/2, scale=2/d, a=0, b=1)
}


#Referencia: Inference and diagnostics in skew scale mixtures
#of normal regression models (2015) - Clecio S. Ferreira, Victor H. Lachos & Heleno Bolfarine

##########################################################################################################################################
# SSMN distributions: Random samples, pdfs e cdfs

rsnn <- function(n, location=0, scale=1, shape=0, dp=NULL) 
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
  if(scale<=0) {
     stop("Parameter scale must be positive")
   } 
  delta<-shape/sqrt(1+shape^2)
  u1<-rnorm(n)
  u2<-rnorm(n)
  y<-location+sqrt(scale)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}


rstn <- function(n, location=0, scale=1, shape=0, nu=30, dp=NULL)  
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<=0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
  if (nu<1) warning('Nu < 1 can generate values tending to infinite',call. = FALSE) 
  u<-rgamma(n,nu/2,nu/2)
  ku<-1/u
  shape1<-shape*sqrt(ku)
  scale1<-scale*ku
  delta<-shape1/sqrt(1+shape1^2)
  u1 <- rnorm(n)
  u2 <- rnorm(n)
  y <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}

rssl <- function(n, location=0, scale=1, shape=0, nu=30, dp=NULL) 
{
 if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]                                                                 
     nu <- dp[4]
    }
   if (nu<=0) {
    stop("Parameter nu must be positive") 
   }
 if(scale<0) {
    stop("Parameter scale must be positive")
   }
 v<-runif(n,0,1)
 u<-v^(1/nu)
 ku<-1/u
 shape1<-shape*sqrt(ku)
 scale1<-scale*ku
 delta<-shape1/sqrt(1+shape1^2)
 u1 <- rnorm(n)
 u2 <- rnorm(n)
 y <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
 return(y)
}


rscn <- function(n, location=0, scale=1, shape=0, nu=1, gama=1, dp=NULL)
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
     gama <- dp[5]
    }
   if (nu<=0|nu>=1) {
    stop("Parameter nu must be between 0 and 1.0") 
   }
   if(gama<=0|gama>=1) {
    stop("Parameter gama must be between 0 and 1.0")
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
  uu<-rbinom(n,1,nu)
  u<-gama*uu+1-uu
  ku<-1/u
  shape1<-shape*sqrt(ku)
  scale1<-scale*ku
  delta<-shape1/sqrt(1+shape1^2)
  u1<-rnorm(n)
  u2<-rnorm(n)
  y<-location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
  return(y)
}

# Density functions of SSMN Distributions


dsnn <- function (x, location=0, scale=1, shape=0, dp=NULL)
{
  if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
    }
  if(scale<=0) {
     stop("Parameter scale must be positive")
   }
   if (any(is.na(x))) x<- x[-which(is.na(x))] 
  z<-(x-location)/sqrt(scale)
  y<-2*dnorm(z)*pnorm(shape*z)/sqrt(scale)
  return(y)
}


# densidade da t-Student normal assimetrica
# ver a formula 3.1 do Clecio； shape = lambda in 2020；  

dstn <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
   if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<=0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  cte<-gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi)
  pdft<-cte*(1+z^2/nu)^(-(nu+1)/2)
  y<-2*pdft*pnorm(shape*z)/sqrt(scale)
  return(y)
}


dssl <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
 if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
    }
   if (nu<=0) {
    stop("Parameter nu must be positive") 
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
 if (any(is.na(x))) x<- x[-which(is.na(x))] 
 z<-(x-location)/sqrt(scale)
 n<-length(x)
 cte=nu/((2*pi*scale)^(1/2))
 d=z^2
 fdps=cte*gamma(nu+1/2)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+1/2))
 fdps[which(z==0)]=cte/(nu+1/2)
 fsl=2*fdps*pnorm(shape*z)
 return(fsl)
}


dscn <- function(x, location=0, scale=1, shape=0, nu=1,gama=1,dp=NULL)
{ 
    if(!is.null(dp)) {
     if(!missing(shape)) 
        stop("You cannot set both component parameters and dp")
     location <- dp[1]
     scale <- dp[2]
     shape <- dp[3]
     nu <- dp[4]
     gama <- dp[5]
    }
   if (nu<=0|nu>=1) {
    stop("Parameter nu must be between 0 and 1.0") 
   }
   if(gama<=0|gama>=1) {
    stop("Parameter gama must be between 0 and 1.0")
   }
   if(scale<=0) {
    stop("Parameter scale must be positive")
   }
   if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  z2<-z*sqrt(gama)
  y<-2*(nu*dnorm(x,mean=location,sd=sqrt(scale/gama))+(1-nu)*dnorm(x,mean=location,sd=sqrt(scale)))*pnorm(shape*z)
  return(y)
}

# Cumulative Distributions functions of SSMN Distributions

#referencia trabalho Thalita - lema 2 pagina 16. (ideia semelhante..)

psnn<- function (x, location=0, scale=1, shape=0, dp=NULL)
{
if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
     }
    else
      { 
        if(scale<=0) {
          stop("Parameter scale must be positive")
        }
        delta=shape/sqrt(1+shape^2)
        z <- (x-location)/sqrt(scale)
        p <- numeric(length(x))
        for (i in 1:length(x)){          
            if(abs(x[i]) == Inf) p[i] <- (1+sign(x[i]))/2
            else{
                covi=matrix(c(1,-delta,-delta,1),2,2)
                Phi2=pmnorm(c(z[i],0),c(0,0),covi)
                p[i]=2*Phi2
        }
      }
   pmax(0,pmin(1,p))
  } 
}


pstn<- function (x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
     }
    else
      { 
        if (nu<=0) {
           stop("Parameter nu must be positive")
        } 
        if (scale<=0) { 
           stop("Parameter scale must be positive")
        }
        p <- numeric(length(x))
        z=(x-location)/sqrt(scale)
        for (i in 1:length(x)){
            if(abs(x[i]) == Inf) p[i] <- (1+sign(x[i]))/2
            else{
                 p[i]=Ephi3stn(z[i],shape,nu) # scale=sigma2
        }
      }
   pmax(0,pmin(1,p))
  }
}




Ephi3stn <- function(z,lambda,nu){
f<-function(u,z,lambda,nu){ 
   detu=lambda/sqrt(u+lambda^2)    # 1-lambda^2/(u+lambda^2)
   phin=rep(0,length(u))
   varcovu=matrix(1,2,2)
   for (i in 1:length(u)){
       varcovu[1,2]=-detu[i]
       varcovu[2,1]=-detu[i]
       phin[i]=pmnorm(z*sqrt(u[i])*c(1,0), mean = rep(0,2), varcov=varcovu)
      }
   fu=2*phin*dgamma(u,nu/2,rate=nu/2)
}
aa=integrate(f,qgamma(0.001,nu/2,rate=nu/2),qgamma(0.999,nu/2,rate=nu/2),z,lambda,nu,subdivisions =100)
return(aa$value)
} 


pssl<- function (x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
     }
    else
      { 
        if (nu<=0) {
           stop("Parameter nu must be positive")
        } 
        if (scale<=0) { 
           stop("Parameter scale must be positive")
        }
        p <- numeric(length(x))
        z=(x-location)/sqrt(scale)
        for (i in 1:length(x)){
            if(abs(x[i]) == Inf) p[i] <- (1+sign(x[i]))/2
            else{
                 p[i]=Ephissl(z[i],shape,nu) # scale=sigma2
        }
      }
   pmax(0,pmin(1,p))
  }
}




Ephissl <- function(z,lambda,nu){
f<-function(u,z,lambda,nu){ 
   detu=lambda/sqrt(u+lambda^2) 
   phin=rep(0,length(u))
   varcovu=matrix(1,2,2)
   for (i in 1:length(u)){
       varcovu[1,2]=-detu[i]
       varcovu[2,1]=-detu[i]
       phin[i]=pmnorm(z*sqrt(u[i])*c(1,0), mean = rep(0,2), varcov=varcovu)
      }
   fu=2*phin*dbeta(u,nu,1)
}
aa=integrate(f,qbeta(0.001,nu,1),qbeta(0.999,nu,1),z,lambda,nu,subdivisions =100)
return(aa$value)
} 


pscn <- function(x, location=0, scale=1, shape=0, nu=1, gama=1, dp=NULL)
{             
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
      gama <- dp[5]
     }
    else
      { 
        if (nu<=0) {
           stop("Parameter nu must be positive")
        } 
        if (scale<=0) { 
           stop("Parameter scale must be positive")
        }
      a=nu*psnn(x, location=location, scale=scale/gama, shape=shape/sqrt(gama))+(1-nu)*psnn(x, location=location, scale=scale, shape=shape)
  }
 return(a)
}

#######################################################################################################################################################################################################################################################################################################################               Truncated SSMN distributions              ##################################################

# Densities and random samples  
## cc 删失指标（0/1取值）        

dsnC <- function (cc, x, location=0, scale=1, shape=0,cens)
{

	densn        <- vector(mode = "numeric", length = length(x))
	densn[cc==0] <- dsnn(x[cc==0], location[cc==0], scale, shape)

if(cens=="left") 
  {
	densn[cc==1] <- psnn(x[cc==1], location[cc==1], scale, shape)
}

if(cens=="right") 
  {
	densn[cc==1] <- 1-psnn(x[cc==1], location[cc==1], scale, shape)
}

  	if(length(which(densn == 0)) > 0) densn[which(densn == 0)] <- .Machine$double.xmin
  	return(densn)
}


dstnC <- function (cc, x, location=0, scale=1, shape=0,nu=30,cens)
{

	densn        <- vector(mode = "numeric", length = length(x))
	densn[cc==0] <- dstn(x[cc==0], location[cc==0], scale, shape, nu)

if(cens=="left") 
  {
	densn[cc==1] <- pstn(x[cc==1], location[cc==1], scale, shape,nu)
}

if(cens=="right") 
  {
	densn[cc==1] <- 1-pstn(x[cc==1], location[cc==1], scale, shape,nu)
}

  	if(length(which(densn == 0)) > 0) densn[which(densn == 0)] <- .Machine$double.xmin
  	return(densn)
}



dsslC <- function (cc, x, location=0, scale=1, shape=0,nu=30,cens)
{

	densn        <- vector(mode = "numeric", length = length(x))
	densn[cc==0] <- dssl(x[cc==0], location[cc==0], scale, shape, nu)

if(cens=="left") 
  {
	densn[cc==1] <- pssl(x[cc==1], location[cc==1], scale, shape,nu)
}

if(cens=="right") 
  {
	densn[cc==1] <- 1-pssl(x[cc==1], location[cc==1], scale, shape,nu)
}


  	if(length(which(densn == 0)) > 0) densn[which(densn == 0)] <- .Machine$double.xmin
  	return(densn)
}


dscnC <- function (cc, x, location=0, scale=1, shape=0,nu=1,gama=1,cens)
{              #print(c(scale,gama))
  densn        <- vector(mode = "numeric", length = length(x))
	densn[cc==0] <- dscn(x[cc==0], location[cc==0], scale, shape,nu,gama)

if(cens=="left") 
  {
	densn[cc==1] <- pscn(x[cc==1], location[cc==1], scale, shape,nu,gama)
}
if(cens=="right") 
  {
	densn[cc==1] <- 1-pscn(x[cc==1], location[cc==1], scale, shape,nu,gama)
}


  	if(length(which(densn == 0)) > 0) densn[which(densn == 0)] <- .Machine$double.xmin
  	return(densn)

}

###############       TSSMN - Random samples         ################################

rTN <- function(n,mu,sigma2,a,b)
{
  #a: lower
  #b: upper
  u <- runif(n)
  sigma <- sqrt(sigma2)
  aux <- u*(pnorm(b,mean=mu,sd=sigma)-pnorm(a,mean=mu,sd=sigma))+pnorm(a,mean=mu,sd=sigma) 
  amostra.x <- qnorm(aux,mean=mu,sd=sigma)
  return(amostra.x)
}


rTSN  <- function(n,mu,sigma2,lambda,a,b)
{    # print(c(n,mu,a,b))
w <-rTN(n,0,1+lambda^2,0,Inf)	
a1 <- ((a-mu)/sqrt(sigma2)-lambda*w/(1+lambda^2))*sqrt(1+lambda^2)
b1 <- ((b-mu)/sqrt(sigma2)-lambda*w/(1+lambda^2))*sqrt(1+lambda^2)
t1 <- rTN(n,0,1,a1,b1)
y=mu+sqrt(sigma2)*( lambda*w/(1+lambda^2) + 1/sqrt(1+lambda^2)*t1)      
return(y)
}


rTST <- function(n,mu,sigma2,lambda,nu,a,b) 
{
  u<-rgamma(n,shape=nu/2,rate=nu/2)
  w <-rTN(n,0,(u+lambda^2)/u,0,Inf)
  #v <- as.numeric(rtmvnorm(1, 0, ((u[i]+lambda^2)/u[i]), 0, Inf))
  a1 <- ((a-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  b1 <- ((b-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  t1 <- rTN(n,0,1,a1,b1) 
  #z <- as.numeric(rtmvnorm(n=1,0,1,a1,b1))
  y=mu+sqrt(sigma2)*( lambda*w/(u+lambda^2) + 1/sqrt(u+lambda^2)*t1) 
return(y)
}


rTSSl <- function(n,mu,sigma2,lambda,nu,a,b)
{
  u <- rbeta(n,shape1=nu,shape2=1)
  w <-rTN(n,0,(u+lambda^2)/u,0,Inf)
  #v <- as.numeric(rtmvnorm(1, 0, ((u[i]+lambda^2)/u[i]), 0, Inf))
  a1 <- ((a-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  b1 <- ((b-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  t1 <- rTN(n,0,1,a1,b1) 
  #z <- as.numeric(rtmvnorm(n=1,0,1,a1,b1))
  y=mu+sqrt(sigma2)*( lambda*w/(u+lambda^2) + 1/sqrt(u+lambda^2)*t1) 
return(y)
}


rTSCN <- function(n,mu,sigma2,lambda,nu,gama,a,b)
{
  uu<-rbinom(n,1,nu)
  u<-gama*uu+1-uu
  w <-rTN(n,0,(u+lambda^2)/u,0,Inf)
  #v <- as.numeric(rtmvnorm(1, 0, ((u[i]+lambda^2)/u[i]), 0, Inf))
  a1 <- ((a-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  b1 <- ((b-mu)/sqrt(sigma2)-(lambda*w)/(u+lambda^2))*sqrt(u+lambda^2)
  t1 <- rTN(n,0,1,a1,b1) 
  #z <- as.numeric(rtmvnorm(n=1,0,1,a1,b1))
  y=mu+sqrt(sigma2)*( lambda*w/(u+lambda^2) + 1/sqrt(u+lambda^2)*t1) 
return(y)
}

##########################################################################################################################################
########## Log-likelihoods

logsnC<-function(theta,y,X,cens,cc){
p=ncol(X)
beta0=matrix(theta[1:p])
mu=X%*%beta0
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
dstk= dsnC(cc, y, location=mu, scale=sigma2, shape=lambda,cens)
lnu=sum(log(dstk))
return(lnu)
}

logstnC<-function(theta,y,X,cens,cc){
p=ncol(X)
beta0=matrix(theta[1:p])
mu=X%*%beta0
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])
dstk= dstnC(cc, y, location=mu, scale=sigma2, shape=lambda,nu=nu,cens)
lnu=sum(log(dstk))
return(lnu)
}


logsslC<-function(theta,y,X,cens,cc){
p=ncol(X)
beta0=matrix(theta[1:p])
mu=X%*%beta0
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])
dstk= dsslC(cc, y, location=mu, scale=sigma2, shape=lambda,nu=nu,cens)
lnu=sum(log(dstk))
return(lnu)
}

logscnC<-function(theta,y,X,cens,cc){
p=ncol(X)
beta0=matrix(theta[1:p])
mu=X%*%beta0
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])
gama=as.numeric(theta[p+4])
dstk= dscnC(cc, y, location=mu, scale=sigma2, shape=lambda,nu=nu,gama=gama,cens)
lnu=sum(log(dstk))
return(lnu)
}

###########################################################################################################################################
######### Information matrix for beta, sigma^2 and lambda


misnC<-function(theta,y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
betas=matrix(theta[1:p])
mu=X%*%betas
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
yc=y
yvetor=as.vector(y)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=rep(1,length(y))   # especifico SN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSN(1000,mu[pos],sigma2,lambda,a,b)     # especifico SN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=rep(1,length(y_ast))    # especifico SN    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
   
miSN=matrix(0,p+2,p+2)
for (i in 1:n){
mui=as.numeric(mu[i])
xi=matrix(X[i,])
dQidbeta=-(-uy[i]+u[i]*mui+lambda*t[i]-yc[i]+lambda^2*mui)/(sigma2)*xi
dQidsigma2=-1/sigma2+(uy2[i]-2*mui*uy[i]+mui^2*u[i]+t2[i]-2*lambda*ty[i]+2*lambda*mui*t[i]+lambda^2*(y2[i]-2*mui*yc[i]+mui^2))/(2*(sigma2^2))
dQidlambda=-(-ty[i]+t[i]*mui+lambda*y2[i]-2*lambda*yc[i]*mui+lambda*(mui^2))/sigma2
mii=matrix(c(dQidbeta,dQidsigma2,dQidlambda))
miSN=miSN+mii%*%t(mii)
}
return(miSN)
}



mistnC<-function(theta,y,X,ind,cens,cutof){ 
# Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
betas=matrix(theta[1:p])
mu=X%*%betas
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])

yc=y
yvetor=as.vector(y)
res=as.vector(y)-mu
d=res^2/sigma2 # vetor
sigma=sqrt(sigma2)
eta=lambda*res
aux=eta/sigma
aux1=pmax(aux,-37)
Wphi=dnorm(aux1)/pnorm(aux1)
t=eta+sigma*Wphi
t2=eta^2+sigma2+sigma*eta*Wphi
u=(nu+1)/(nu+d)   # especifico StN
y2=yvetor^2
uy=u*yvetor
uy2=u*y2
ty=t*yvetor
yc=yvetor

for (i in 1:lc){
        pos=ind[i]
        y_ast= rTST(1000,mu[pos],sigma2,lambda,nu,a[i],b[i])    # especifico StN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=(nu+1)/(nu+d_ast)    # especifico StN    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
}   
miSTN=matrix(0,p+2,p+2)
for (i in 1:n){
    mui=as.numeric(mu[i])
    xi=matrix(X[i,])
    dQidbeta=-(-uy[i]+u[i]*mui+lambda*t[i]-yc[i]+lambda^2*mui)/(sigma2)*xi
    dQidsigma2=-1/sigma2+(uy2[i]-2*mui*uy[i]+mui^2*u[i]+t2[i]-2*lambda*ty[i]+2*lambda*mui*t[i]+lambda^2*(y2[i]-2*mui*yc[i]+mui^2))/(2*(sigma2^2))
    dQidlambda=-(-ty[i]+t[i]*mui+lambda*y2[i]-2*lambda*yc[i]*mui+lambda*(mui^2))/sigma2
    mii=matrix(c(dQidbeta,dQidsigma2,dQidlambda))
    miSTN=miSTN+mii%*%t(mii)
#print(miSTN)
}

 return(miSTN)

 }




misslC<-function(theta,y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
betas=matrix(theta[1:p])
mu=X%*%betas
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])

 yc=y
 yvetor=as.vector(y)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=as.vector((2*nu+1)/d*(pgamma(1,nu+1.5,scale=2/d)/pgamma(1,nu+0.5,scale=2/d)))  # especifico SSL
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSSl(1000,mu[pos],sigma2,lambda,nu,a[i],b[i])    # especifico SSL
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((2*nu+1)/d_ast*(pgamma(1,nu+1.5,scale=2/d_ast)/pgamma(1,nu+0.5,scale=2/d_ast)))    # especifico SSL    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
   
miSSL=matrix(0,p+2,p+2)
for (i in 1:n){
mui=as.numeric(mu[i])
xi=matrix(X[i,])
dQidbeta=-(-uy[i]+u[i]*mui+lambda*t[i]-yc[i]+lambda^2*mui)/(sigma2)*xi
dQidsigma2=-1/sigma2+(uy2[i]-2*mui*uy[i]+mui^2*u[i]+t2[i]-2*lambda*ty[i]+2*lambda*mui*t[i]+lambda^2*(y2[i]-2*mui*yc[i]+mui^2))/(2*(sigma2^2))
dQidlambda=-(-ty[i]+t[i]*mui+lambda*y2[i]-2*lambda*yc[i]*mui+lambda*(mui^2))/sigma2
mii=matrix(c(dQidbeta,dQidsigma2,dQidlambda))
miSSL=miSSL+mii%*%t(mii)
}
return(miSSL)
}




miscnC<-function(theta,y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu,gama]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
betas=matrix(theta[1:p])
mu=X%*%betas
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
nu=as.numeric(theta[p+3])
gama=as.numeric(theta[p+4])

 yc=y
 yvetor=as.vector(y)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d/2)))    # especifico SCN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSCN(1000,mu[pos],sigma2,lambda,nu,gama,a[i],b[i])    # especifico SCN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d_ast/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d_ast/2)))   # especifico SCN    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
   
miSCN=matrix(0,p+2,p+2)
for (i in 1:n){
mui=as.numeric(mu[i])
xi=matrix(X[i,])
dQidbeta=-(-uy[i]+u[i]*mui+lambda*t[i]-yc[i]+lambda^2*mui)/(sigma2)*xi
dQidsigma2=-1/sigma2+(uy2[i]-2*mui*uy[i]+mui^2*u[i]+t2[i]-2*lambda*ty[i]+2*lambda*mui*t[i]+lambda^2*(y2[i]-2*mui*yc[i]+mui^2))/(2*(sigma2^2))
dQidlambda=-(-ty[i]+t[i]*mui+lambda*y2[i]-2*lambda*yc[i]*mui+lambda*(mui^2))/sigma2
mii=matrix(c(dQidbeta,dQidsigma2,dQidlambda))
miSCN=miSCN+mii%*%t(mii)
}
return(miSCN)
}

#########################################################################################################################################
############### EM algorithm, tau fixed

sn.emC<-function(y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
} 
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 theta0<-rbind(beta0,sigma2,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=as.vector(y)
 u=rep(1,n)   # especifico SN
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1   
    mu=as.vector(X%*%beta0)
    res=yc-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi    
    y2=yc^2
    ty=t*yc

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSN(1000,as.numeric(mu[pos]),sigma2,lambda,a[i],b)     # especifico SN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-mu[pos]
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)   
    }
    beta0=solve(t(X)%*%((1+lambda^2)*diag(1,n))%*%X)%*%t(X)%*%matrix((1+lambda^2)*yc-lambda*t)
    sigma2=sum(y2-2*mu*yc+mu^2+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,sigma2,lambda)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta     
 } 
mi= misnC(theta,y,X,ind,cens,cutof)
ep=sqrt(diag(solve(mi))) ## sd of thetahat  
logi=logsnC(theta,y,X,cens,cc)
return(list(theta=theta,logmax=logi,iter=cont,ep=ep))  # list: theta,log-like,cont,sd 
 }


sn.sigma2.emC<-function(sigma2,y,X,ind,cens,cutof){ 
 # Theta=[beta,lambda]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 theta0<-rbind(beta0,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=yc-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi    
    y2=yc^2
    ty=t*yc

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSN(1000,mu[pos],sigma2,lambda,a,b)     # especifico SN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-mu[pos]
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2) 
    }
    beta0=solve(t(X)%*%((1+lambda^2)*diag(1,n))%*%X)%*%t(X)%*%matrix((1+lambda^2)*yc-lambda*t)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,lambda)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta    
 } 

 return(theta)

 }


stn.emC<-function(y,X,ind,nu,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 theta0<-rbind(beta0,sigma2,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-6)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=(nu+1)/(nu+d)    # especifico StN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor
    for (i in 1:lc){
        pos=ind[i]
	      y_ast= rTST(1000,mu[pos],sigma2,lambda,nu,a[i],b[i])     # especifico StN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        #print(length(y_ast))
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast= (nu+1)/(nu+d_ast)    # especifico StN    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,sigma2,lambda)
    #print(theta)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta    
 } 

 return(theta)

 }


stn.sigma2.emC<-function(sigma2,y,X,ind,nu,cens,cutof){
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 #sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 theta0<-rbind(beta0,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=(nu+1)/(nu+d)    # especifico StN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor      
    #
    for (i in 1:lc){
        pos=ind[i]
	      y_ast= rTST(1000,mu[pos],sigma2,lambda,nu,a,b)     # especifico StN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        #print(length(y_ast))
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast= (nu+1)/(nu+d_ast)    # especifico StN    
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    #sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,lambda)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta     
 } 

 return(theta)

 }


ssl.emC<-function(y,X,ind,nu,cens,cutof){
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 theta0<-rbind(beta0,sigma2,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1 
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=as.vector((2*nu+1)/d*(pgamma(1,nu+1.5,scale=2/d)/pgamma(1,nu+0.5,scale=2/d)))    # especifico SSL
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSSl(1000,mu[pos],sigma2,lambda,nu,a,b)     # especifico SSL
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((2*nu+1)/d_ast*(pgamma(1,nu+1.5,scale=2/d_ast)/pgamma(1,nu+0.5,scale=2/d_ast)))   # especifico SSL   
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,sigma2,lambda)
#print(theta)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta     
 } 

 return(theta)

 }






scn.emC<-function(y,X,ind,nu1,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu]

nu=nu1[1]
gama=nu1[2]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 theta0<-rbind(beta0,sigma2,lambda)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1   
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u= as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d/2)))   # especifico SCN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSCN(1000,mu[pos],sigma2,lambda,nu,gama,a,b)    # especifico SCN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d_ast/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d_ast/2)))    # especifico SCN   
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    theta=rbind(beta0,sigma2,lambda)
#print(theta)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta    
 } 

 return(theta)

 }


###################################      ECME         ##################################################################################    




stn.ecmeC<-function(y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 nu=5
 theta0<-rbind(beta0,sigma2,lambda,nu)
 criterio=1
 cont=0
 W=1000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=(nu+1)/(nu+d)    # especifico StN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor      
    #
    for (i in 1:lc){
        pos=ind[i]
	      y_ast= rTST(1000,mu[pos],sigma2,lambda,nu,a[i],b[i])     # especifico StN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast= (nu+1)/(nu+d_ast)    # especifico StN  
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    thetanu=rbind(beta0,sigma2,lambda)
    nu<-optimize(ftnu,c(2.001,30),thetanu,y,X,cens,cc)$minimum
    theta=rbind(beta0,sigma2,lambda,nu)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta    
 } 
logi=logstnC(theta,y,X,cens,cc)
mi= mistnC(theta,y,X,ind,cens,cutof)
ep=sqrt(diag(solve(mi)))
return(list(theta=theta,logmax=logi,iter=cont,ep=ep))
 }
 
ftnu<-function(nu,theta,y,X,cens,cc){
p=ncol(X)
beta0=theta[1:p]
mu=as.vector(X%*%beta0)
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
ft<-sum(log(dstnC(cc, y, location=mu, scale=sigma2, shape=lambda,nu,cens)))
return(-ft)
}




ssl.ecmeC<-function(y,X,ind,cens,cutof){ 
 # Theta=[beta,sigma2,lambda,nu]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 nu=10
 theta0<-rbind(beta0,sigma2,lambda,10)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u=as.vector((2*nu+1)/d*(pgamma(1,nu+1.5,scale=2/d)/pgamma(1,nu+0.5,scale=2/d)))    # especifico SSL
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor    
    #
    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSSl(1000,mu[pos],sigma2,lambda,nu,a[i],b[i])     # especifico SSL
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((2*nu+1)/d_ast*(pgamma(1,nu+1.5,scale=2/d_ast)/pgamma(1,nu+0.5,scale=2/d_ast)))   # especifico SSL   
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    thetanu=rbind(beta0,sigma2,lambda)
    nu<-optimize(fSLnu,c(1.001,30),thetanu,y,X,cens,cc)$minimum
    theta=rbind(beta0,sigma2,lambda,nu)
    #print(theta)
    dif=theta-theta0
    criterio=(sum(dif^2))^0.5
    theta0=theta     
 } 
logi=logsslC(theta,y,X,cens,cc)
mi= misslC(theta,y,X,ind,cens,cutof)
ep=sqrt(diag(solve(mi)))
return(list(theta=theta,logmax=logi,iter=cont,ep=ep))
 }


fSLnu<-function(nu,theta,y,X,cens,cc){
p=ncol(X)
beta0=theta[1:p]
mu=as.vector(X%*%beta0)
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
ft<-sum(log(dsslC(cc, y, location=mu, scale=sigma2, shape=lambda,nu,cens)))
return(-ft)
}




scn.ecmeC<-function(y,X,ind,cens,cutof){  
 # Theta=[beta,sigma2,lambda,nu,gama]
lc=length(ind)
if(cens=="left"){
a=-Inf*rep(1,lc)
b=cutof
}
if(cens=="right"){
a=cutof
b=Inf*rep(1,lc)
}
n=length(y)
p=ncol(X)
 beta0<-solve(t(X)%*%X)%*%t(X)%*%y
 res0=y-X%*%beta0
 lambda<-as.numeric(skewness(res0))
 sigma2<-as.numeric(t((y-X%*%beta0))%*%(y-X%*%beta0)/(n-p))
 nu=0.5
 gama=0.5
 theta0<-rbind(beta0,sigma2,lambda,nu,gama)
 criterio=1
 cont=0
 W=3000 # numero maximo de iteracoes
 yc=y
 yvetor=as.vector(y)
 like <- c()  
 while ((criterio > 1e-4)&&(cont<W)) {
    cont=cont+1  
    mu=as.vector(X%*%beta0)
    res=as.vector(y)-mu
    d=res^2/sigma2 # vetor
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    t=eta+sigma*Wphi
    t2=eta^2+sigma2+sigma*eta*Wphi
    u= as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d/2)))   # especifico SCN
    y2=yvetor^2
    uy=u*yvetor
    uy2=u*y2
    ty=t*yvetor
    yc=yvetor

    for (i in 1:lc){
        pos=ind[i]
        y_ast= rTSCN(1000,mu[pos],sigma2,lambda,nu,gama,a[i],b[i])    # especifico SCN
	      y_ast=y_ast[is.finite(y_ast)==TRUE]
        resT=y_ast-as.numeric(mu[pos])
        etaT=lambda*resT
        aux=etaT/sigma
        aux1=pmax(aux,-37)
        Wphi_ast=dnorm(aux1)/pnorm(aux1)
        d_ast=resT^2/sigma2
        u_ast=as.vector((1-nu+nu*gama^(1.5)*exp((1-gama)*d_ast/2))/(1-nu+nu*gama^(0.5)*exp((1-gama)*d_ast/2)))    # especifico SCN   
        t_ast=lambda*resT+sigma*Wphi_ast
        t2_ast=(lambda*resT)^2+sigma2+sigma*lambda*resT*Wphi_ast
        t[pos]=mean(t_ast) 
        ty[pos]=mean(t_ast*y_ast)
        t2[pos]=mean(t2_ast)   
        u[pos]=mean(u_ast)
        uy[pos]=mean(u_ast*y_ast)
        uy2[pos]=mean(u_ast*(y_ast^2))
        yc[pos]=mean(y_ast)
        y2[pos]=mean(y_ast^2)    
    }
    beta0=solve(t(X)%*%(diag(u)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%matrix(uy-lambda*t+lambda^2*yc)
    sigma2=sum(uy2-2*mu*uy+mu^2*u+t2-2*lambda*ty+2*lambda*mu*t+lambda^2*(y2-2*mu*yc+mu^2))/(2*n)
    lambda=sum(ty-mu*t)/sum(y2-2*mu*yc+mu^2)
    thetanu=rbind(beta0,sigma2,lambda)
    nugama<-optim(c(nu,gama),fCNnugama,gr=NULL,thetanu,y,X,cens,cc,method="L-BFGS-B",lower=c(0.1,0.1),upper=c(0.99,0.99),control=list(maxit=50))$par
    nu=as.numeric(nugama[1])
    gama=as.numeric(nugama[2])
    theta=rbind(beta0,sigma2,lambda,nu,gama)
    #print(theta)
    dif=theta-theta0
    log0 <- logscnC(theta,y,X,cens,cc)  ## 
    criterio=(sum(dif^2))^0.5
    theta0=theta     
    like <- c(like,log0)  
 } 
mi=miscnC(theta,y,X,ind,cens,cutof)
ep=sqrt(diag(solve(mi)))
logi=logscnC(theta,y,X,cens,cc)
return(list(theta=theta,logmax=logi,iter=cont,ep=ep,like=like))
 }

fCNnugama<-function(nugama,theta,y,X,cens,cc){
nu=nugama[1]
gama=nugama[2]
p=ncol(X)
beta0=theta[1:p]
mu=as.vector(X%*%beta0)
sigma2=as.numeric(theta[p+1])
lambda=as.numeric(theta[p+2])
ft<-sum(log(dscnC(cc, y, location=mu, scale=sigma2, shape=lambda,nu,gama,cens)))
return(-ft)
}
    