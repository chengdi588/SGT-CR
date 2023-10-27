
 
# source(SGT.CR.basic.functions.R)   
Newton.optim <- function(Q.012, beta_init, max_iter = 100, tolerance = 1e-6) { 
  beta <- beta_init  
  for (iter in 1:max_iter) {  
    gradient <- Q.012(beta)$f1 
    hessian <- Q.012(beta)$f2  
    beta_new <- beta - 0.5 * solve(hessian) %*% gradient    
    update_norm <- sqrt(sum((beta_new - beta)^2)) 
    beta <- beta_new  
    
    if (update_norm < tolerance) {   
      break
    }
  }  
  return(beta)    
}  

        
## data ----------------------------------------------    

library(astrodatR)     
data(censor_Be)         
data_stel <-  censor_Be  
class(data_stel)       
colnames(data_stel)       
summary(data_stel)      

y <- data_stel$logN_Be   
x <- data_stel$Teff/1000   
X <- cbind(1,x)   
ratio <- 12/68    
n <-length(y)       
small1 <- which(data_stel$Ind_Be==0 )       
vi <- y   

 
initial_theta <- function(vi,X,small1){
  n <- length(vi)  ; p <- ncol(X) ;   big1 <- setdiff(seq(1,n,by=1), small1) ; Ta <- 0.46 ; Tb <- 1.76     
  
  betaLS <- solve(t(X[big1,])%*%X[big1,])%*%t(X[big1,])%*%vi[big1]    
  res0 <- (vi-X%*%betaLS)[big1]       
  sigma0<-as.numeric( sqrt(t(res0)%*%res0/(length(big1)-p )) )   
  mu_0 <- X%*%betaLS       
  vi_stand<- ( vi- mu_0  )[big1]/sigma0            
  b <- 2      
  r_ini <-ifelse( abs( 1-2*sum(vi_stand <= 0)/n ) < 0.99 , 1-2*sum(vi_stand <= 0)/n, sign( sum(vi_stand <= 0)/n -1 )*runif(0,1) )     
  a_ini <- try( optim( Ta,  
                      function(par_a){ -sgt.logL(c( betaLS, sigma0, r_ini , par_a, b ),y,X,vi,ratio,small1) } 
                      , method = "L-BFGS-B", lower = 0.001 , upper =Ta+2 ) $par )    
  if('try-error' %in% class(a_ini)){    
    a_ini <- optim( Ta+0.1,  
                      function(par_a){ -sgt.logL(c( betaLS, sigma0, r_ini , par_a, b ),y,X,vi,ratio,small1) } 
                      , method = "L-BFGS-B", lower = 0.001 , upper =Ta+2 ) $par       
  }   
  theta0 <- c(betaLS, sigma0, r_ini, a_ini, b)    
  return(theta0)
} 
 
sgt.em1.data <- function(y,X,vi, ratio, small1){     
    n <- length(y) ;  p <- ncol(X) ;  Ta <- 0.461 ; Tb <- 1.76  
    big1 <- setdiff(seq(1,n,by=1), small1)      
    
    rho <- rep(0,n)    
    rho[small1] <- 1      

    theta0 <- initial_theta(vi,X,small1)  
    log0 = sgt.logL(theta0,y,X,vi,ratio,small1 )             
    
    beta1<-signif(theta0[1:p],digits=7)   
    mu=X%*%beta1     
    sigma<-signif(theta0[p+1],digits=7) 
    r<-signif(theta0[p+2],digits=7)
    a<-signif(theta0[p+3],digits=7)
    b<-signif(theta0[p+4],digits=7)      
    
    criterio <- 1      
    cont<-0    
    M<- 500 
    M_print <- c() ; delta.like <- 1      
    while ( criterio >1e-5 && cont < M ){  
      M_print <- c(M_print, log0 )  
      cont<-cont+1           
      mu1 <- X%*%beta1        
      c1 <- (vi-mu1)/sigma                    
      
      Estepi=Estep1(y,X,vi,theta0,ratio,small1)           
      e1i_o = as.vector(Estepi$e1i_o)    
      e2i_sstar0 = as.vector(Estepi$e2i_sstar0)        
      e2i_star0=as.vector(Estepi$e2i_star0)     
      e3i_o=as.vector(Estepi$e3i_o)        
      fq012_beta  <- function(par_beta ){  
        qb_obs1<-t(X[small1,])%*%( e2i_sstar0[small1] )
        qb_obs2<-t(X[big1,])%*%( my.sign(y[big1]-X[big1,]%*%par_beta)/(1+r*my.sign(y[big1]-X[big1,]%*%par_beta))^b *
                                   my.abs(y[big1]-X[big1,]%*%par_beta)^(b-1)*e1i_o[big1] )
        f <-qb_obs1 + qb_obs2  
        f <-as.numeric(-b*f) 
        
        my.sign.dao <- 2*lambda*( exp(-lambda*(y[big1]-X[big1,]%*%par_beta) ) )/(1+ exp(-lambda*(y[big1]-X[big1,]%*%par_beta) ) )^2 
        com_j <- my.sign(y[big1]-X[big1,]%*%par_beta)^2* my.abs(y[big1]-X[big1,]%*%par_beta)^(b-2)/(1+r*my.sign(y[big1]-X[big1,]%*%par_beta))^(b+1)*(
          (b-1)*(1+r*my.sign(y[big1]-X[big1,]%*%par_beta))+
            my.sign(y[big1]-X[big1,]%*%par_beta)^2* my.sign.dao*(
              my.abs(y[big1]-X[big1,]%*%par_beta) - (b-1)*r*(y[big1]-X[big1,]%*%par_beta)
            ) )
        j <- -t(X[big1,])%*%diag(as.numeric(com_j*e1i_o[big1]))%*%(X[big1,])  
        j <- -b*j  
        return(list( f1=f, f2=j ))      
      }      
      beta1 <- Newton.optim(fq012_beta, beta1, max_iter = 100, tolerance = 1e-6)       
      
      qs_obs <- sum(abs(y[big1]-X[big1,]%*%matrix(beta1,nrow=p))^b*e1i_o[big1] / (1+r*my.sign(y[big1]-X[big1,]%*%matrix(beta1,nrow=p)) )^b ) 
      qs <- sum( e2i_star0[small1] ) + qs_obs      
      sigma <- (b/(2*n)*qs)^(1/b)         
      
      fq012_r<-function(par_r){   
        uc1<-(1+abs(c1[small1])^b/(2*a*(1+r*sign(c1[small1]))^b ) )^(-1)  
        Iuc1<- pbeta(uc1,a, 1+1/b)   
        zzi<- (y[big1]-X[big1,]%*%beta1) 
        f<- sigma^b *sum(Iuc1/Fx( c1[small1], r,a,b ) )-b* matrix(e1i_o[big1]* my.abs(zzi)^b* my.sign(zzi), nrow=1 )%*%( (1+rep(par_r, length(big1) )*my.sign(zzi) )^(-b-1) )
        j<- b*(b+1) *e1i_o[big1]%*%( my.abs( zzi)^(b) /(1+ rep(par_r, length(big1))*my.sign(zzi))^(b+2) )
        return( list(f1=f, f2=j)  )      
      }   
      r <- try( Newton.optim(fq012_r, r, max_iter = 100, tolerance = 1e-6)  )       
      
      if('try-error' %in% class(r) || abs(r)>1 ){
        r <- optim( r_ini, 
                    function(par_r){ -sgt.logL(c( beta1, sigma, par_r , a, b ), y, X,vi,ratio,small1 ) }
                    , method = "L-BFGS-B", lower = -0.98  , upper = 0.98)$par
      }   
      r <- as.numeric(r)   
      
      bifun_a<-function(par_a){ 
        f<- log(par_a) - digamma(par_a ) + 1 -mean( e1i_o ) + mean( e3i_o )    
        return( f^2 )      
      }       
      a <- try( optim(Ta, bifun_a, method = "L-BFGS-B", lower = max(a-1,0.02) , upper = a+4  )$par )   
      if('try-error' %in% class(a)){ 
        a <- optim( a-0.3, 
                    function(par_a){ -sgt.logL(c( beta1, sigma, r , par_a, b ),y,X,vi,ratio,small1 ) }
                    , method = "L-BFGS-B", lower = a-1, upper = a+1 ) $par  
      }   
      b <- try(optim( b,  
                      function(par_b){ -sgt.logL(c( beta1, sigma, r , a, par_b ),y,X,vi,ratio,small1 ) } 
                      , method = "L-BFGS-B", lower = 0.1, upper = 15 ) $par ) 
      if('try-error' %in% class(a)){ 
        b <- Newton.optim_b(b, c1, rho, a, r)  
      }    
      
      theta_new <- matrix(c(beta1,sigma, r, a, b ),p+4,1)         
      log1 = sgt.logL(theta_new,y,X,vi,ratio, small1 )     
      criterio<-(  abs( (log0-log1)/log0 )  ) 
      delta.like <- log1 - log0    
      if (delta.like<0) {
        break   
      } 
      theta0 <- theta_new      
      log0 <- log1 
      print(c(log0,cont))   
    }   
    #  
    data <- data.frame(x = 1:length(M_print), y = M_print)
    smoothed <- loess(y ~ x, data = data)
    smoothed_data <- data.frame(x = data$x, y = predict(smoothed))
    plot(smoothed_data$x, smoothed_data$y, type = "p", pch = 20, col = "red", 
         xlab = "Iterations", ylab = "log likelihood", cex.lab = 2.2, cex.axis = 2.2, cex.main =2.2) 
    lines(smoothed_data$x, smoothed_data$y, col = "red", lwd = 2) 
   
    MIO <- I_O(y, X, vi, theta0, small1)       
    se <- diag(solve(MIO))   
    return( list(thetaMLE=theta0, cont= cont, se = se, loglike=log0 ) )   
}     

## AIC BIC --------------------------------------------    

par(mar = c(5, 6, 5, 2) + 0.1)   
EM_SGT_result <- sgt.em1.data(y,X,vi, ratio, small1)   
thetaMLE <- theta0 <- EM_SGT_result$thetaMLE    
log0 <- EM_SGT_result$loglike   
se_hat <- EM_SGT_result$se    
round( cbind(thetaMLE , se_hat), 4)  

 

ln_l <- c(round(log0,4), -2.1278, -1.7802, -4.5144 )     
K_VEC <- c( 6, 5, 5, 3 )    

# AIC     
-2* ln_l + 2* K_VEC 

# BIC  
-2* ln_l + K_VEC*log(n)    

# EDC  
-2* ln_l + K_VEC*0.2*sqrt(n)  

# ABIC   
N <- (n +2)/24       
-2* ln_l + K_VEC*log(N)   

cbind( theta0, se_hat ) 
rbind( ln_l,
       -2* ln_l + 2* K_VEC ,
       -2* ln_l + K_VEC*log(n) ,
       -2* ln_l + K_VEC*0.2*sqrt(n) ,
       -2* ln_l + K_VEC*log(N)  ) 

# SGT - CR   
time1 <- Sys.time()   
EM_SGT_result <- sgt.em1.data(y,X,vi, ratio, small1)   
time2 <- Sys.time() 
EM_SGT_result
time2 - time1   # 2.54 s  

# T - CR   (2018)   
# code  not available   


# source(functions.R)    
# SSMN - CR (2020)   
data(censor_Be)   
ccaux=data_stel$Ind_Be     
aux=1:n
ind=aux[ccaux==0]    
cc=rep(0,n)
cc[ind]=1 

cens="left"
cutof=y[ind]      

time1 <- Sys.time()  
thetaSSMN =stn.ecmeC(y,X,ind,cens,cutof)
time2 <- Sys.time() 
thetaSSMN  

time2 - time1  #   2.0027 min    

## N-CR model   

library(CensRegMod)           

time1 <- Sys.time() 
ncr <- em.cens(cc, X, y, dist="Normal",  
               diagnostic="TRUE", typediag="1")   
time2 <- Sys.time() 
time2 - time1  # 51.6s     

ncr$influents   
ncr_GD <- ncr$measure        

## QD + GD ----------  
par(mar = c(5, 6, 5, 2) + 0.1)  

p <- ncol(X)     
my.QD <- QD_apro(y,X,vi, theta0, small1)   
plot( x=(1:n), y= my.QD,   
      xlab = 'index', ylab = 'QD', cex.lab = 2.2, cex.axis = 2.2, cex.main =2.2  ) 
text( x = which.max(my.QD)+ 1.5, y = max(my.QD),    
      labels= which.max(my.QD) ,cex =2 )     

GD_apro <- function(i, y, X, vi, theta, small1)    
{  
  p <- ncol(X)    
  n <- nrow(X)
  M <- new_HessQ(y,X, vi, theta, small1)   
  Qderi <- Q_grad(i, y, X, vi, ratio, theta, small1)   
  MatrizQ <- -M[[1]]      
  
  GDi <- t(Qderi) %*% solve(MatrizQ) %*% Qderi    
  return(abs(GDi))       
}       

my.GD <- numeric(length(y)) 
for (i in 1:length(y)) {
  my.GD[i] <- GD_apro(i, y, X, vi, theta0, small1) 
}   
plot( x=(1:n), y= my.GD,   
      xlab = 'index', ylab = 'GD' , cex.lab = 2.2, cex.axis = 2.2, cex.main =2.2  ) 
text( x = which.max(my.GD)+ 1.5, y = max(my.GD),    
      labels= which.max(my.GD) ,cex=2 )      



## ----------------------  local ---------------------------------------------------  
  
library(MASS)    
 
types <- c("case", "resp", "exp", "scale") 
   
for (type in types){  
  All_B_fq <- calculate_B_fq(y = y, X = X, vi = vi, ratio = ratio, theta = theta0, small1 = small1, type = type)
 
  plot( x=(1:n), y= All_B_fq , main=as.character(type), 
        xlab = 'index',ylab = 'M(0)'  , cex.lab = 2.2, cex.axis = 2.2, cex.main =2.0 )   
  c_star <- 3.5
  bench <- sd(All_B_fq) * c_star + 1/n 
  abline(h = bench, lty = 2) 
  text( x = which.max(All_B_fq)+ 2.5, y = max(All_B_fq),    
        labels= which.max(All_B_fq) ,cex= 2 )     
  
  cat("Type:", type, "\n" )
  print(which(All_B_fq > bench))   
  
  is_in_range <- All_B_fq >= 0 & All_B_fq <= 1 
  out_of_range_values <- All_B_fq[!is_in_range]  
  cat("Values not between 0 and 1：", out_of_range_values, "\n")   
}    

# TRC RC  
# remove 5  
y <- data_stel$logN_Be   
x <- data_stel$Teff/1000   
X <- cbind(1,x)   
ratio <- 12/68    
n <-length(y)       
small1 <- which(data_stel$Ind_Be==0 )       
vi <- y   

y <- y[-5] ; X <- X[-5,] ; vi <- vi[-5]     
EM_SGT_remove5 <- sgt.em1.data(y, X, vi, ratio, small1)
thetaMLE_remove5 <- EM_SGT_remove5$thetaMLE   
log0_remove5 <- EM_SGT_remove5$loglike   

# remove 26   
y <- data_stel$logN_Be    
x <- data_stel$Teff/1000   
X <- cbind(1,x)   
ratio <- 12/68    
n <-length(y)       
small1 <- which(data_stel$Ind_Be==0 )       
vi <- y  

y <- y[-26] ; X <- X[-26,] ; vi <- vi[-26]        

EM_SGT_remove26 <- sgt.em1.data(y, X, vi, ratio, small1)   
thetaMLE_remove26 <- EM_SGT_remove26$thetaMLE    
log0_remove26 <- EM_SGT_remove26$loglike   

# remove 5 + 26    
y <- data_stel$logN_Be    
x <- data_stel$Teff/1000   
X <- cbind(1,x)   
ratio <- 12/68    
n <-length(y)       
small1 <- which(data_stel$Ind_Be==0 )       
vi <- y  

y <- y[-c(5,26)] ; X <- X[-c(5,26),] ; vi <- vi[-c(5,26)]        

EM_SGT_remove5_26 <- sgt.em1.data(y, X, vi, ratio, small1)   
thetaMLE_remove5_26 <- EM_SGT_remove5_26$thetaMLE    
log0_remove5_26 <- EM_SGT_remove5_26$loglike   

# TRC MRC  

RC <- function(theta, theta_remove){
  abs( (theta-theta_remove) / theta )  
} 
RC_values <- t( cbind( RC(thetaMLE, thetaMLE_remove5),   
          RC(thetaMLE, thetaMLE_remove26) ,
          RC(thetaMLE, thetaMLE_remove5_26)  ) ) 
cbind( RC_values, 
       apply(RC_values, 1, sum) ,  
       apply(RC_values, 1, max) ) 

# original  data   
y <- data_stel$logN_Be           
x <- data_stel$Teff/1000   
X <- cbind(1,x)   
ratio <- 12/68    
n <-length(y)       
small1 <- which(data_stel$Ind_Be==0 )       
vi <- y  

###  N-CR diagnostic  ------------------------------------------------------------------------  

par(mar = c(5, 6, 5, 2) + 0.1)   
   
plot( x=(1:n), y= ncr_GD[,1] ,  
      xlab = 'index',ylab = 'GD', cex.lab = 2.2, cex.axis = 2.2, cex.main =2.0, main ="N-CR"  )  
text( x =  ncr$influents$GD  + 1.5, y = ncr_GD[ncr$influents$GD,1],
      labels= ncr$influents$GD ,cex= 2 )

plot( x=(1:n), y= ncr_GD[,2] , 
      xlab = 'index',ylab = expression(GD(beta)), cex.lab = 2.2, cex.axis = 2.2, cex.main =2.0,main ="N-CR"  )    
text( x =  ncr$influents$GDbeta + 1.2, y = ncr_GD[ncr$influents$GDbeta,2], 
      labels= ncr$influents$GDbeta , cex= 2 )   

plot( x=(1:n), y= ncr_GD[,3] ,   
      xlab = 'index',ylab = expression(GD(sigma^2)), cex.lab = 2.2, cex.axis = 2.2, cex.main =2.0,main ="N-CR"  )    
text( x =  ncr$influents$GDsigma2 + 1.99, y = ncr_GD[ncr$influents$GDsigma2,3],  
      labels= ncr$influents$GDsigma2 ,cex= 2 )   


## ------------ 4 residual -----------------------------------------------    

library(ggpubr)      

# SGT  
rMI <- function(y, vi, X, small1, theta){    
  n <- length(y) ;  p <- ncol(X)  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  beta1<-signif(theta[1:p],digits=7)   
  mu=X%*%beta1     
  sigma<-signif(theta[p+1],digits=7) 
  r<-signif(theta[p+2],digits=7)
  a<-signif(theta[p+3],digits=7)
  b<-signif(theta[p+4],digits=7)      
  
  zi <- (vi-mu) / sigma 
  rmi <- rep(NA, length(zi)) ;  big1 <- setdiff(seq(1,n,by=1), small1)  
  Syi <- 1- Fx( zi, r, a, b  )    
  rmiC <- log( Syi ) 
  rmiO <- 1 + log( Syi )  
  rmi[small1] <- rmiC[small1] 
  rmi[big1] <- rmiO[big1] 
  return(rmi)     
}

rMTI <- function(y, vi, X, small1, theta){ 
  n <- length(y) ;  p <- ncol(X)   
  big1 <- setdiff(seq(1,n,by=1), small1)      

  rdi <- rep(NA, length(y)) ;  big1 <- setdiff(seq(1,n,by=1), small1)   
  rdiO <- sign( rMI(y, vi, X, small1, theta)[big1] )*sqrt( -2*(rMI(y, vi, X, small1, theta)[big1] + log(1-rMI(y, vi, X, small1, theta)[big1] )) )  
  rdiC <- sign( rMI(y, vi, X, small1, theta)[small1] )*sqrt( -2*( rMI(y, vi, X, small1, theta)[small1] ) )  
  rdi[small1] <- rdiC 
  rdi[big1] <- rdiO   
  return(rdi)   
} 
     
sgt_residual <- rMTI(y, vi, X, small1, theta0)    
stand_residual <- (sgt_residual-mean(sgt_residual)) / sd(sgt_residual)
print( ggqqplot( stand_residual,font.x = 27,font.y = 27,font.main = 27,
                 font.xtickslab =27, font.ytickslab =27 ,main="SGT-CR") )   


# 2018 Model   
library(mvtnorm) 
theta2018 <- c( -2.2426, 0.5453, 0.0674, -6.4219, 3 ) # beta; sigma2,lambda,nu 

# def         
rMI <- function(y, vi, X, small1, theta){    
  n <- length(y) ;  p <- ncol(X)  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  beta1<-signif(theta[1:p],digits=7)   
  mu <- X%*%beta1      
  sigma <- sqrt(signif(theta[p+1],digits=7) ) 
  lambda<- signif(theta[p+2],digits=7)
  nu<-signif(theta[p+3],digits=7)    
  
  rmi <- rep(NA, length(vi)) ;  big1 <- setdiff(seq(1,n,by=1), small1)  
  results <- rep(NA, length(vi))
  for (i in 1:length(vi)) {  
    vii <- vi[i]  
    result <- try( pst(vii, xi = mu[i], omega = sigma, alpha = lambda, nu = nu ) , silent = TRUE )  
    if (class(result) == "numeric") {    
      results[i] <- result   
    } else {
      results[i] <- NA   
    }
  } 
  results <- na.omit(results)  
  
  Syi <- 1-  results
  rmiC <- log( Syi )  
  rmiO <- 1 + log( Syi )   
  rmi[small1] <- rmiC[small1]  
  rmi[big1] <- rmiO[big1] 
  return(rmi)     
} 

st2018_residual <- na.omit(rMTI(y, vi, X, small1, theta2018) )      
stand_residual <- (st2018_residual-mean(st2018_residual)) / sd(st2018_residual) 
print( ggqqplot( stand_residual,font.x = 27,font.y = 27,font.main = 27,
                 font.xtickslab =27, font.ytickslab =27 , main="ST-CR") )           


## SSMN:     
 
theta2020 <- thetaSSMN$theta          
        
rMI <- function(y, vi, X, small1, theta){    
  n <- length(y) ;  p <- ncol(X)  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  beta1<-signif(theta[1:p],digits=7)   
  mu=X%*%beta1      
  sigma<-signif(theta[p+1],digits=7) 
  lambda<-signif(theta[p+2],digits=7)
  nu<-signif(theta[p+3],digits=7)
    
  zi <- (vi-mu) / sigma 
  rmi <- rep(NA, length(zi)) ;  big1 <- setdiff(seq(1,n,by=1), small1)  
  results <- rep(NA, length(zi))
  for (i in 1:length(zi)) {  
    zii <- zi[i]  
    result <- try( pstn(zii, shape = lambda, nu = nu) , silent = TRUE )  
    if (class(result) == "numeric") {    
      results[i] <- result   
    } else {
      results[i] <- NA   
    }
  } 
  results <- na.omit(results) 
 
  Syi <- 1-  results 
  rmiC <- log( Syi )   
  rmiO <- 1 + log( Syi )   
  rmi[small1] <- rmiC[small1]  
  rmi[big1] <- rmiO[big1]  
  return(rmi)     
} 
  
ssmn_residual <- na.omit(rMTI(y, vi, X, small1, theta2020) )     
stand_residual <- (ssmn_residual-mean(ssmn_residual)) / sd(ssmn_residual)
print( ggqqplot( stand_residual,font.x = 27,font.y = 27,font.main = 27,
                 font.xtickslab =27, font.ytickslab =27 , main="SSMN-CR" ) )     

# N - CR model   
thetaNCR <- c(ncr$betas , ncr$sigma2)   

rMI <- function(y, vi, X, small1, theta){    
  n <- length(y) ;  p <- ncol(X)  
  beta1<-signif(theta[1:p],digits=7)   
  mu=X%*%beta1       
  sigma<-sqrt(signif(theta[p+1],digits=7) )     
   
  zi <- (vi-mu) / sigma      
  rmi <- rep(NA, length(zi)) ;  big1 <- setdiff(seq(1,n,by=1), small1)  
  Syi <- 1-  pnorm( zi, mean =0 ,sd = 1 )    
  rmiC <- log( Syi )    
  rmiO <- 1 + log( Syi )   
  rmi[small1] <- rmiC[small1]  
  rmi[big1] <- rmiO[big1] 
  return(rmi)     
} 

NCR_residual <- na.omit(rMTI(y, vi, X, small1, thetaNCR) )       
stand_residual <- (NCR_residual-mean(NCR_residual)) / sd(NCR_residual) 
print( ggqqplot( stand_residual ,font.x = 27,font.y = 27,font.main = 27,
                 font.xtickslab =27, font.ytickslab =27 , main ="N-CR" ) )        
 
 
 


