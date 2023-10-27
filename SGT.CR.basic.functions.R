
 
library(numDeriv)   
library(progress)    

rsgt <- function(n, location=0, scale=1, r=0, a=1,b=2, dp=NULL)  
{ 
  if (!is.null(dp)) {
    if (!missing(r)) {
      stop("You cannot set both component parameters and dp")
    }
    location <- dp[1]
    scale <- dp[2]
    r <- dp[3]
    a <- dp[4]
    b <- dp[5]
  }
  
  if (scale <= 0 || a <= 0 || b <= 0) {
    stop("Parameters 'scale', 'a', and 'b' must be positive")
  }
  
  if (abs(r) > 1) {
    stop("Parameter 'r' must be bounded between -1 and 1")
  }  
  p <- (1 + r) / 2
  w1 <- rbinom(n, 1, p) * 2 - 1 + r
  y1 <- rgamma(n, shape = 1 / b, rate = 1) 
  z1 <- rgamma(n, shape = a, rate = a)
  error <- 2^(1/b) * w1 * y1^(1/b) * z1^(-1/b)
  y <- location + scale * error 
  return(y)          
}       

my.sign <- function(x) {
  lam <- 50
  y <- 1 / (1 + exp(-lam * x))
  return(2 * y - 1)   
}

my.abs <- function(x) {
  lam1 <- 0.00001
  ret <- Vectorize(function(x) ifelse(abs(x) < lam1, x^2 / (2 * lam1), abs(x) - lam1 / 2))(x)   
  return(ret) 
}     

fx <- function(x, r, a, b) {
  constant <- b / (2^(1 + 1/b) * a^(1/b) * beta(a, 1/b))
  d <- (1 +  my.abs(x)^b / (2 * a * (1 + r *  my.sign(x))^b))^(-a - 1/b)
  return(constant * d) 
}

uc<-function(x,r,a,b){ 
  ucre <- (1 + abs(x)^b / (2 * a * (1 + r * sign(x))^b))^(-1)
  return(ucre)
}    

Fx <- function(x, r, a, b) {
  ux <- (1 + my.abs(x)^b / (2 * a * (1 + r * my.sign(x))^b))^(-1)
  Iu <- pbeta(ux, a, 1/b)  
  d <- ifelse(x < 0, (1 - r) / 2 * Iu, 1 - (1 + r) / 2 * Iu)
  return(d) 
} 

 
sgt.logL <-function(theta,y,X,vi, ratio,small1)  
{      
  n=length(y)            
  p=dim(X)[2]   
  beta<-signif(theta[1:p],digits=7)
  mu=X%*%beta  
  sigma<-signif(theta[p+1],digits=7) 
  r<-signif(theta[p+2],digits=7) 
  a<-signif(theta[p+3],digits=7)
  b<-signif(theta[p+4],digits=7)       
  z1<- (vi-mu)/sigma    
  
  if(ratio!=0){    
    big1 <- setdiff(seq(1,n,by=1), small1)      
    vero<-signif( c( Fx(z1[small1],r,a,b), fx(z1[big1],r,a,b)/sigma ),digits=7 )  
    return(sum(log(vero)))  
  }  
  if(ratio==0){ 
    vero<-signif( fx(z1,r,a,b)/sigma,digits=7 )  
    return(sum(log(vero)))   
  } 
}  


e2i_sstar<-function(c1,X,theta){     
  p<-ncol(X)   
  sigma<-signif(theta[p+1],digits=7)   
  r<-signif(theta[p+2],digits=7)
  a<-signif(theta[p+3],digits=7)
  b<-signif(theta[p+4],digits=7)
  
  e2i_sstar0<-c() 
  uc1<-uc(c1,r,a,b)   
  ret<- sigma^(b-1)/Fx(c1, r, a, b)/(2*a)^(1/b)/beta(a,1/b)*pbeta(uc1,a+1/b,1) 
  e2i_sstar0<-c(e2i_sstar0, ret) 
  return(-e2i_sstar0)      
}  

   
e2i_star <- function(c1, X, theta) {
  p <- ncol(X)
  beta <- theta[1:p]
  sigma <- theta[p + 1]
  r <- theta[p + 2]
  a <- theta[p + 3]
  b <- theta[p + 4]
  if (length(c1) == 0) {
    return(NULL)   
  }
  uc1 <- uc(c1, r, a, b) 
  e2i_sstar0 <- numeric(length(c1))
  
  e2i_star_negc1 <- sigma^b / b / Fx(c1, r, a, b) * (1 - r) * pbeta(uc1, a, 1/b + 1)
  e2i_star_posc1 <- sigma^b / b / Fx(c1, r, a, b) * (2 - (1 + r) * pbeta(uc1, a, 1/b + 1))

  e2i_sstar0[c1 <= 0] <- e2i_star_negc1[c1 <= 0]
  e2i_sstar0[c1 > 0] <- e2i_star_posc1[c1 > 0]   
  
  return(round(e2i_sstar0, 7))    
}


int.der.pbeta <- function(x, a, b) {
  der.pbeta <- function(u, a, b) {
    ret <- u^(a - 1) * (1 - u)^(b - 1) * log(u)
    return(ret)
  }
  
  if (is.vector(x)) {
    result <- sapply(x, function(xi) {
      integrate(function(u) der.pbeta(u, a, b), lower = 0, upper = xi)$value
    })
    return(result) 
  } else {
    return(integrate(function(u) der.pbeta(u, a, b), lower = 0, upper = x)$value)
  }
}


Estep1 <- function(y,X,vi,theta,ratio,small1){   
  n <- length(y)
  p <- dim(X)[2]
  beta <- signif(theta[1:p], digits = 7)
  mu <- X %*% beta
  sigma <- signif(theta[p + 1], digits = 7)
  r <- signif(theta[p + 2], digits = 7)
  a <- signif(theta[p + 3], digits = 7)
  b <- signif(theta[p + 4], digits = 7)
  
  e1i_o <- matrix(0, n, 1) 
  e3i_o <- matrix(0, n, 1)
  aux <- (vi - mu) / sigma 
  auxi <- pmax(aux, -50) 
  
  kappa <- a + abs(auxi)^b / (2 * (1 + r * sign(auxi))^b)
  e1i_o <- (a + 1 / b) / kappa   
  e3i_o <- digamma(a + 1 / b) - log(kappa)  
  
  if(ratio==0){     
    return( list(e1i_o=e1i_o,  e3i_o =e3i_o ) )      
  }  
  if(ratio!=0){      
    uc1 <- uc(aux[small1], r, a, b)  
    cond <- ( aux[small1] < 0 )   
    if(sum(cond)==length(small1) ){  
      e1i_o[small1[cond]] <- 0.5 * (1 - r) / Fx(aux[small1[cond]], r, a, b) * pbeta(uc1[cond], a + 1, 1 / b)
      pbetadao.a <- (int.der.pbeta(uc1[cond], a, 1 / b) - pbeta(uc1[cond], a, 1 / b) * beta(a, 1 / b) * (digamma(a) - digamma(a + 1 / b))) / beta(a, 1 / b)
      e3i_o[small1[cond]] <- digamma(a) - log(a) + 0.5 * (1 - r) * pbetadao.a / Fx(aux[small1[cond]], r, a, b)
    } 
    if(sum(cond)==0){
      e1i_o[small1[!cond]] <- (1 - 0.5 * (1 + r) * pbeta(uc1[!cond], a + 1, 1 / b)) / Fx(aux[small1[!cond]], r, a, b)
      pbetadao.a <- (int.der.pbeta(uc1[!cond],  a, 1 / b) - pbeta(uc1[!cond], a, 1 / b) * beta(a, 1 / b) * (digamma(a) - digamma(a + 1 / b))) / beta(a, 1 / b)
      e3i_o[small1[!cond]] <- digamma(a) - log(a) - 0.5 * (1 + r) * pbetadao.a / Fx(aux[small1[!cond]], r, a, b)  
    } 
    if(sum(cond)>0 && sum(cond)<length(small1) ){
      e1i_o[small1[cond]] <- 0.5 * (1 - r) / Fx(aux[small1[cond]], r, a, b) * pbeta(uc1[cond], a + 1, 1 / b)
      pbetadao.a <- (int.der.pbeta(uc1[cond], a, 1 / b) - pbeta(uc1[cond], a, 1 / b) * beta(a, 1 / b) * (digamma(a) - digamma(a + 1 / b))) / beta(a, 1 / b)
      e3i_o[small1[cond]] <- digamma(a) - log(a) + 0.5 * (1 - r) * pbetadao.a / Fx(aux[small1[cond]], r, a, b)
      e1i_o[small1[!cond]] <- (1 - 0.5 * (1 + r) * pbeta(uc1[!cond], a + 1, 1 / b)) / Fx(aux[small1[!cond]], r, a, b)
      pbetadao.a <- (int.der.pbeta(uc1[!cond],  a, 1 / b) - pbeta(uc1[!cond], a, 1 / b) * beta(a, 1 / b) * (digamma(a) - digamma(a + 1 / b))) / beta(a, 1 / b)
      e3i_o[small1[!cond]] <- digamma(a) - log(a) - 0.5 * (1 + r) * pbetadao.a / Fx(aux[small1[!cond]], r, a, b) 
    }
    
    e2i_sstar0<- e2i_sstar(aux,X,theta) 
    e2i_star0<-e2i_star(aux,X,theta)  
    return( list(e1i_o=e1i_o,  
                 e2i_sstar0=e2i_sstar0 ,e2i_star0=e2i_star0 , e3i_o =e3i_o  ))   
  }  
} 


library(moments)  
lambda<- 50       

Newton.optim <- function(Q.012, beta_init, max_iter = 100, tolerance = 1e-6) {
  beta <- beta_init
  
  for (iter in 1:max_iter) { 
    gradient <- Q.012(beta)$f1 
    hessian <- Q.012(beta)$f2  
    beta_new <- beta - solve(hessian) %*% gradient
    update_norm <- sqrt(sum((beta_new - beta)^2)) 
    beta <- beta_new 
    
    if (update_norm < tolerance) { 
      break
    }
  }  
  return(beta) 
}      

 
Newton.optim_b <- function(b0, c1, rho, a, r){  
  
  b0_value <- b0     
  
  epsilon <- 1e-5   
  max_iterations <- 100    
  
  n <- length(c1) 
  M <- abs(c1) / (1 + r * sign(c1)) 
  
  for (ii in 1:max_iterations) {
    
    derivative1 <- numeric(n)
    hessian_M1 <- numeric(n)
    fd1 <- numeric(n)
    fd2 <- numeric(n)
    
    
    for (jj in 1:n) {
      x <- c1[jj] 
      derivative1[jj] <- ifelse(x <= 0,   
                                grad(function(b) log(1 - (1 + r) / 2 * pbeta((1 + (-x)^b / (2 * a * (1 - r)^b))^(-1), a, 1 / b)), b0_value), 
                                grad(function(b) log(pbeta((1 + x^b / (2 * a * (1 + r)^b))^(-1), a, 1 / b)), b0_value))

      hessian_M1[jj] <- ifelse(x <= 0, 
                               hessian(function(b) log(1 - (1 + r) / 2 * pbeta((1 + (-x)^b / (2 * a * (1 - r)^b))^(-1), a, 1 / b)), b0_value), 
                               hessian(function(b) log(pbeta((1 + x^b / (2 * a * (1 + r)^b))^(-1), a, 1 / b)), b0_value))  
      
      M_val <- M[jj]
      fd1[jj] <- 1/b0_value + 1/b0_value^2 * (digamma(1/b0_value) - digamma(a + 1/b0_value) + log(2 * a) + log(1 + 1/2/a * M_val^b0_value)) -
        (a + 1/b0_value) * ((2 * a)^(-1) * M_val^b0_value * log(M_val)) / (1 + (2 * a)^(-1) * M_val^b0_value)
      fd2[jj] <- 1/b0_value^2 * (1/a * M_val^b0_value * log(M_val) / (1 + M_val^b0_value/2/a) - 1) -
        2/b0_value^3 * (digamma(1/b0_value) - digamma(a + 1/b0_value) + log(2 * a) + log(1 + M_val^b0_value/2/a)) -
        1/b0_value^4 * (trigamma(1/b0_value) - trigamma(a + 1/b0_value)) -
        (a + 1/b0_value) / (2 * a) * (log(M_val))^2 * M_val^b0_value / (1 + M_val^b0_value/2/a)^2
    }
    
    derivative1[abs(derivative1) < 1e-8] <- 1e-8
    hessian_M1[abs(hessian_M1) < 1e-8] <- 1e-8 
    
    f <- t(rho) %*% derivative1 + t(1 - rho) %*% fd1    
    j <- t(rho) %*% hessian_M1 + t(1 - rho) %*% fd2     
    f <- as.numeric(f)
    j <- as.numeric(j) 
 
    b1 <- b0_value - f / j  
    
    if (b1 < 0.1 || b1 > 20){  
      b1 <- Tb -0.2 
      next     
    }       

    delta_x <- abs(b1 - b0_value)   
    b0_value <- b1  
    if (all(delta_x < epsilon)) { 
      break
    }   
  } 
  return(b0_value)    
}  


initial_theta <- function(vi,X,small1){
  n <- length(vi)  ; p <- ncol(X) ;   big1 <- setdiff(seq(1,n,by=1), small1)     
  
  betaLS <- solve(t(X[big1,])%*%X[big1,])%*%t(X[big1,])%*%vi[big1]    
  res0 <- (vi-X%*%betaLS)[big1]       
  sigma0<-as.numeric( sqrt(t(res0)%*%res0/(length(big1)-p )) )    
  mu_0 <- X%*%betaLS       
  vi_stand<- ( vi- mu_0  )[big1]/sigma0           
  b <- 2   
  r_ini <-ifelse( abs( 1-2*sum(vi_stand <= 0)/n ) < 0.99 , 1-2*sum(vi_stand <= 0)/n, sign( sum(vi_stand <= 0)/n -1 )*runif(0,1) )     
  a_ini <- try( optim(Ta-0.2,  
                      function(par_a){ -sgt.logL(c( betaLS, sigma0, r_ini , par_a, b ),y,X,vi,ratio,small1) } 
                      , method = "L-BFGS-B", lower = max(Ta-2, 0.2), upper =Ta+2 ) $par,silent = TRUE  )     
  if('try-error' %in% class(a_ini)){     
    a_ini <- Ta-0.1       
  }   
  theta0 <- c(betaLS, sigma0, r_ini, a_ini, b)  
  return(theta0) 
}

Rbias <- function(est,theta){ 
  ret <- mean( abs(est-theta)/abs(theta) )     
  return( ret )   
}   
RMSE<-function(est,theta){
  d<-sqrt( mean((est-theta)^2)  )  
  return( d )  
}       

# simulation   
sgt.em1<-function(y,X,ratio){ 
  n <- length(y) ;  p <- ncol(X)            
  if(ratio!=0){     
    ncens <- floor(n*ratio)     
    small1 <- sort(sample(n,ncens) )           
    big1 <- setdiff(seq(1,n,by=1), small1)      
    
    ui1 <- runif(n) ; ui2 <- runif(n) 
    ci <- pmin( y+ui2,y-ui1+1 )     
    rho <- rep(0,n)   
    rho[small1] <- 1      
    vi <- y       
    vi[small1] <- ci[small1]        
    
    theta0 <- matrix( initial_theta(vi,X,small1), p+4, 1)  
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
      cont<-cont+1           
      mu1 <- X%*%beta1        
      c1 <- (vi-mu1)/sigma                
      
      Estepi=Estep1(y,X,vi,theta0,ratio,small1)           
      e1i_o = as.vector(Estepi$e1i_o)    
      e2i_sstar0 = as.vector(Estepi$e2i_sstar0)        
      e2i_star0=as.vector(Estepi$e2i_star0)    
      e3i_o=as.vector(Estepi$e3i_o)    
      ##  
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
      Nbeta1 <- try(Newton.optim(fq012_beta, beta1, max_iter = 100, tolerance = 1e-6), silent = TRUE)    
      beta1 <- Nbeta1  
      qs_obs <- sum(abs(y[big1]-X[big1,]%*%matrix(beta1,nrow=p))^b*e1i_o[big1] / (1+r*my.sign(y[big1]-X[big1,]%*%matrix(beta1,nrow=p)) )^b ) 
      qs <- sum( e2i_star0[small1] ) + qs_obs      
      sigma <- (b/(2*n)*qs)^(1/b)         
      
      #   
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
      a <- try( optim(Ta-0.3, bifun_a, method = "L-BFGS-B", lower = 1.1 , upper = Ta+4  )$par )      
      if('try-error' %in% class(a)){
        a <- optim( Ta-0.3, 
                    function(par_a){ -sgt.logL(c( beta1, sigma, r , par_a, b ),y,X,vi,ratio,small1 ) }
                    , method = "L-BFGS-B", lower = max(Ta-1, 0.5), upper =Ta+1 )$par      
      }   
      b <- try(optim( b,  
                      function(par_b){ -sgt.logL(c( beta1, sigma, r , a, par_b ),y,X,vi,ratio,small1 ) } 
                      , method = "L-BFGS-B", lower = 0.1, upper = 15 ) $par ) 
      if('try-error' %in% class(a)){ 
        b <- Newton.optim_b(Tb-1, c1, rho, a, r)  
      }  
      theta_new <- matrix(c(beta1,sigma, r, a, b ),p+4,1)        
      log1 = sgt.logL(theta_new,y,X,vi,ratio, small1 )     
      criterio<-(  abs( (log0-log1)/log0 )  ) 
      theta0 <- theta_new      
      delta.like <- log1 - log0         
      log0 <- log1    
      M_print <- c(M_print, log0 ) 
      
    }   
  
  }    
  return( list(thetaMLE=theta0, cont= cont, loglike=log0, vi=vi, small1 = small1 ) )   
} 


## ------------------------- se -----------------------------------   

s_i <-function(i, y, X, vi, theta, small1){   
  n <- length(y)    
  p <- ncol(X)   
  beta <- signif(theta[1:p], digits = 7)     
  sigma <- signif(theta[p + 1], digits = 7)  
  r <- signif(theta[p + 2], digits = 7)
  a <- signif(theta[p + 3], digits = 7)
  b <- signif(theta[p + 4], digits = 7)
  mu <- X %*% beta  
  c1 <- (vi - mu) / sigma    
  rho <- rep(0, n)   
  rho[small1] <- 1  
  
  Estepi <- Estep1(y, X, vi, theta, ratio, small1)         
  e1i_o <- as.vector(Estepi$e1i_o)   
  e2i_sstar0 <- as.vector(Estepi$e2i_sstar0) 
  e2i_star0 <- as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  uc1 <- uc(c1[i], r, a, b)   
  cons <- e1i_o[i] * my.sign(y[i] - X[i,] %*% beta) * 
    my.abs(y[i] - X[i,] %*% beta)^(b - 1) /
    (1 + r * my.sign(y[i] - X[i,] %*% beta))^b 
  x_i <- X[i,]  
  
  d1 <- b / (2 * sigma^b) * ( rho[i] * x_i * e2i_sstar0[i] + (1 - rho[i]) * x_i * as.vector(cons))
  d2 <- -1 / sigma + b / (2 * sigma^(b + 1)) * (rho[i] * e2i_star0[i]  + (1 - rho[i]) * e1i_o[i] *
                                                  my.abs(y[i] - X[i,] %*% beta)^(b) /
                                                  (1 + r * my.sign(y[i] - X[i,] %*% beta))^b ) 
  d3 <- -rho[i] * pbeta(uc1, a, 1 + 1 / b) / (2 * Fx(c1[i], r, a, b)) + (1 - rho[i]) * b / (2 * sigma^(b)) * e1i_o[i] *
    my.sign(y[i] - X[i,] %*% beta) *
    my.abs(y[i] - X[i,] %*% beta)^(b) /
    as.vector( (1 + r * my.sign(y[i] - X[i,] %*% beta))^(b + 1) ) 
  d4 <- log(a) + 1 - digamma(a) - e1i_o[i] + e3i_o[i] 
  d5 <- 1 / b + 1 / (b^2) * (log(2) + digamma(1 / b) - e3i_o[i]) - (1 - rho[i]) * e1i_o[i] * 0.5 * (
    my.abs(y[i] - X[i,] %*% beta)^(b) /
      as.vector(sigma^b * (1 + r * my.sign(y[i] - X[i,] %*% beta))^(b) )  
  ) * ( (
    my.abs(y[i] - X[i,] %*% beta) /
      as.numeric(sigma * (1 + r * my.sign(y[i] - X[i,] %*% beta))  ) 
  ) ) - rho[i] * 0.5 * (a + 1 / b) / Fx(c1[i], r, a, b) * exp4i(c1[i], uc1, r, a, b)  
  d5 <- as.numeric(d5)  
  return(matrix(c(d1, d2, d3, d4, d5), nrow = p + 4 ))       
}   
 
exp4i <- function(c1, uc1, r, a, b){ 
  ret1 <- adaptIntegrate(function(u) { 
    u^(a - 1) * (1 - u)^(1 / b) * (log(2 * a) + log(1 - u) - log(u))
  }, lower = 0, upper = uc1)$integral
  
  ret2 <- 2 / (a * b^2 + b) * (log(2 * a) + digamma(1 / b + 1) - digamma(a)) - (1 + r) / (b * beta(a, 1 / b)) * ret1   
  if (c1 > 0) {    
    return(ret2) 
  } else { 
    return((1 - r) / (b * beta(a, 1 / b)) * ret1)    
  }
}   

I_O <- function(y, X, vi, theta, small1) {
  n <- length(y)
  p <- ncol(X)
  I <- matrix(0, nrow = p + 4, ncol = p + 4)
  for (i in 1:n) {
    ret1 <- s_i(i, y, X, vi, theta, small1)
    I <- I + ret1 %*% t(ret1)
  } 
  return(round(I, 4))
}   

library(cubature) 

## ------------------------- LOCAL -----------------------------------   

calculate_w_var <- function(y, X, beta, sigma, r, a, b, small1, e1i_o, c1) {
  n <- length(y) ; uc1<-uc(c1,r,a,b)  
  rho <- rep(0,n) ; rho[small1] <- 1    
  
  w1i <- function(i){  
    d <-  e1i_o[i]*my.sign(y[i]-X[i,]%*%beta)* 
      my.abs(y[i]-X[i,]%*%beta)^(b-1)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return(d) 
  }
  w2i <- function(i){
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return(d) 
  } 
  w3i <- function(i){  
    uc1<-uc(c1[i],r,a,b)   
    d <- pbeta(uc1, a, 1+1/b)/Fx(c1[i], r,a,b ) 
    return(d) 
  } 
  w4i <- function(i){
    d <-  e1i_o[i]*my.sign(y[i]-X[i,]%*%beta)* 
      my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^(b+1)   
    return(d) 
  }  
  w5i <- function(i){
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b-2)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^(b)   
    return(d) 
  }   
  w6i <- function(i){
    d <-  w2i(i)/my.abs(y[i]-X[i,]%*%beta)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))    
    return(d) 
  }   
  w7i <- function(i){  
    cont <- 1/Fx(c1[i], r,a,b )  
    if(c1[i]<0) { 
      d <-  (uc1[i])^(a+1/b)/(1-r) 
      return(cont*d) 
    }else{
      d <- 1/(1-r) + ( 1-(uc1[i])^(a+1/b) )/(1+r)  
      return(cont*d)    
    }
  }  
  w8i <- function(i){ 
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    d1 <- d/(1+r*my.sign(y[i]-X[i,]%*%beta))^2   
    return(d1)  
  }  
  w9i <- function(i){  
    cont <- 1/Fx(c1[i], r,a,b )  
    if(c1[i]<0) { 
      d <- pbeta(uc1[i],a,1+1/b)/(1-r) 
      return(cont*d) 
    }else{
      d <- 2/(1-r^2) - pbeta(uc1[i],a,1+1/b)/(1+r)  
      return(cont*d)   
    }
  } 
  w10i <- function(i){  
    d <- my.abs(y[i]-X[i,]%*%beta)/ (1+r*my.sign(y[i]-X[i,]%*%beta)) 
    return( log(d) ) 
  } 
  w11i <- function(i){    
    d <- my.abs(y[i]-X[i,]%*%beta)/ (1+r*my.sign(y[i]-X[i,]%*%beta)) / sigma 
    return( log(d)^2 * w2i(i) )    
  } 
  s1i <- function(i){     
    d <- my.abs(y[i]-X[i,]%*%beta)/ (1+r*my.sign(y[i]-X[i,]%*%beta)) / sigma 
    w2ii <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return( log(d) * w2ii )     
  } 
  s2i <- function(i){     
    d <- exp4i(c1[i], uc1[i], r, a, b)     
    return( d / Fx(c1[i], r,a,b )  )      
  }  
  s1 <- as.matrix( sapply(1:n, s1i),ncol=1 )
  s2 <- as.matrix( sapply(1:n, s2i),ncol=1 )  
  
  w1 <- sapply(1:n, w1i) 
  w2 <- as.matrix( sapply(1:n, w2i),ncol=1 )  
  w3 <- as.matrix( sapply(1:n, w3i),ncol=1 )  
  w4 <- as.matrix( sapply(1:n, w4i),ncol=1 )  
  W1 <- diag(w1)   
  w5 <- as.matrix( sapply(1:n, w5i),ncol=1 )  
  w6 <- as.matrix( sapply(1:n, w6i),ncol=1 )  
  w7 <- as.matrix( sapply(1:n, w7i),ncol=1 )  
  w8 <- as.matrix( sapply(1:n, w8i),ncol=1 )  
  w9 <- as.matrix( sapply(1:n, w9i),ncol=1 )  
  w10 <- as.matrix( sapply(1:n, w10i),ncol=1 ) 
  w11 <- as.matrix( sapply(1:n, w11i),ncol=1 )  
  return(list(w1=w1, w2=w2, w3=w3, w4=w4, w5=w5, w6=w6, w7=w7, w8=w8, w9=w9, w10=w10, w11=w11,
              s1=s1, s2=s2))  
}


Q_grad <- function(j, y, X, vi, ratio, theta, small1) {  
  n <- length(y)   
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)      
  sigma <- signif(theta[p+1],digits=7)   
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)
  mu <- X%*%beta   
  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  c1 <- (vi-mu)/sigma       
  uc1<-uc(c1,r,a,b)   
  rho <- rep(0,n) ; rho[small1] <- 1     
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)   
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0=as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  
  e2_sst <- diag( e2i_sstar0 )     
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )   
  
  w_s_vars <- calculate_w_var(y, X, beta, sigma, r, a, b, small1, e1i_o, c1) 
  w1 <- w_s_vars$w1
  w2 <- w_s_vars$w2   
  w3 <- w_s_vars$w3   
  w4 <- w_s_vars$w4   
  W1 <- diag(w1)   
  w5 <- w_s_vars$w5   
  w6 <- w_s_vars$w6   
  w7 <- w_s_vars$w7    
  w8 <- w_s_vars$w8  
  w9 <- w_s_vars$w9   
  w10 <- w_s_vars$w10  
  w11 <- w_s_vars$w11    
  s1 <- w_s_vars$s1  
  s2 <- w_s_vars$s2   
  
  if(j!=0){    
    y <- y[-j] ; X <- X[-j,] ; vi <- vi[-j]   
    if(j %in% small1){ small1 <-  small1[small1 != j]  }else{small1 <- small1} 
    
    e1i_o = as.vector(e1i_o)[-j]    
    E2_star = e2i_sstar0 = as.vector(e2i_sstar0)[-j]  
    e2i_star0 = as.vector(e2i_star0)[-j] 
    e3i_o <- as.vector(e3i_o)[-j]     
    big1 <- setdiff(seq(1,n,by=1), small1)     
    c1<- c1[-j]      
    uc1 <- uc1[-j]       
    rho <- rho[-j]    
    
    w1 <- w1[-j] ; w2 <- w2[-j,] ; w3 <- w3[-j,] ; w4 <- w4[-j,] ; w5 <- w5[-j,] 
    w6 <- w6[-j,] ; w7 <- w7[-j,] ; w8 <- w8[-j,]  ; w9 <- w9[-j,] ; w10 <- w10[-j,]  ; w11 <- w11[-j,]     
    W1 <- diag(w1) ; s1 <- s1[-j,] ; s2 <- s2[-j,]   
    
    e2_sst <- diag( e2i_sstar0 )     
    e2 <- as.matrix( e2i_star0,ncol=1 )   
    e1 <- as.matrix( e1i_o,ncol=1 )  
    e3 <- as.matrix( e3i_o,ncol=1 )    
  
    dQ_dbeta <- (b / (2 * sigma^b)) * t(X) %*% (E2_star * rho + w1 * (rep(1, n - 1) - rho))
    dQ_dsigma <- -((n - 1) / sigma) + (b / (2 * sigma^(b + 1))) * (t(e2) %*% rho + t(w2) %*% (rep(1, n - 1) - rho))
    dQ_dr <- (-1/2) * t(w3) %*% rho + (b / (2 * sigma^b)) * t(w4) %*% (rep(1, n - 1) - rho)
    dQ_da <- (n - 1) * (log(a) - digamma(a) + 1) - t(rep(1, n - 1)) %*% e1 + t(rep(1, n - 1)) %*% e3
    dQ_db <- (n - 1) * ((1/b) + (log(2) / (b^2)) + (digamma(1/b) / (b^2))) - (1 / (b^2)) * t(rep(1, n - 1)) %*% e3
    dQ_db <- dQ_db - (1 / (2 * sigma^b)) * t(rep(1, n - 1) - rho) %*% s1 - ((a + (1/b)) / 2) * t(rho) %*% s2 
    
     
    partial_derivatives <- rbind( dQ_dbeta, dQ_dsigma, dQ_dr,  dQ_da, dQ_db)  
    return(partial_derivatives)
  }   
  
  if(j==0){ 
    n <- length(y) ; p <- ncol(X) ;  Q0 <- matrix(0, nrow = p+4, ncol = 1) 
    for (i in 1:n) {
      ret1 <- s_i(i, y, X, vi, theta, small1) 
      Q0 <- Q0+ ret1 
    } 
    return(Q0) }      
}   


    
new_HessQ <- function(y,X,vi,theta,small1){ 
  n <- length(y)    
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)  
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7) 
  mu <- X%*%beta  
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)       # n*5        
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0 = as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)   
  
  E2_sst <- diag( e2i_sstar0 )   
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )  
  
  rho <- rep(0,n) ; rho[small1] <- 1   
  ncens   <- floor(n*ratio)     
  
  c1<-(vi-mu)/sigma      
  big1 <- setdiff(seq(1,n,by=1), small1)  
  uc1 <- uc(c1,r,a,b)
  
  w1i <- function(i){  
    d <-  e1i_o[i]*my.sign(y[i]-X[i,]%*%beta)* 
      my.abs(y[i]-X[i,]%*%beta)^(b-1)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return(d) 
  }
  w2i <- function(i){
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return(d) 
  } 
  w3i <- function(i){ 
    uc1<-uc(c1[i],r,a,b)   
    d <- pbeta(uc1, a, 1+1/b)/Fx(c1[i], r,a,b ) 
    return(d) 
  } 
  w4i <- function(i){
    d <-  e1i_o[i]*my.sign(y[i]-X[i,]%*%beta)* 
      my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^(b+1)   
    return(d) 
  }  
  
  e2i_ssstar <- function(i){
    cont <- sigma^(b-2)*(a+1/b)/(2*a)^(2/b)/ Fx(c1[i], r,a,b )*beta(a+2/b,-1/b+1)/beta(a,1/b)  
    if(c1[i]<0) { 
      d <- cont/(1-r)*pbeta(uc1[i],a+2/b,-1/b+1)  
    }else{
      d <- cont*(2/(1-r^2) - 1/(1+r)*pbeta(uc1[i],a+2/b,-1/b+1) )  
    }
    return(d)   
  }  
  w5i <- function(i){
    d <-  e1i_o[i]* 
      my.abs(y[i]-X[i,]%*%beta)^(b-2)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    return(d)  
  }
  w6i <- function(i){
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    d1 <- d / my.abs(y[i]-X[i,]%*%beta)/(1+r*my.sign(y[i]-X[i,]%*%beta)) 
    return(d1) 
  } 
  w7i <- function(i){  
    cont <- 1/Fx(c1[i], r,a,b )  
    if(c1[i]<0) { 
      d <-  uc1[i]^(a+1/b)/(1-r) 
      return(cont*d) 
    }else{
      d <- 1/(1-r) + ( 1-uc1[i]^(a+1/b) )/(1+r)  
      return(cont*d)    
    }
  }  
  w8i <- function(i){ 
    d <-  e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^b 
    d1 <- d/(1+r*my.sign(y[i]-X[i,]%*%beta))^2   
    return(d1)  
  } 
  
  w9i <- function(i){  
    cont <- 1/Fx(c1[i], r,a,b )  
    if(c1[i]<0) { 
      d <- pbeta(uc1[i],a,1+1/b)/(1-r) 
      return(cont*d) 
    }else{
      d <- 2/(1-r^2) - pbeta(uc1[i],a,1+1/b)/(1+r)  
      return(cont*d)   
    }
  } 
  e2i_ssstar0 <- sapply(1:n,e2i_ssstar) 
  w1 <- sapply(1:n, w1i) 
  w2 <- as.matrix( sapply(1:n, w2i),ncol=1 ) 
  w3 <- as.matrix( sapply(1:n, w3i),ncol=1 ) 
  w4 <- as.matrix( sapply(1:n, w4i),ncol=1 )  
  w5 <- sapply(1:n, w5i) 
  w6 <- sapply(1:n, w6i) 
  w7 <- sapply(1:n, w7i)  
  w8 <- sapply(1:n, w8i)   
  w9 <- sapply(1:n, w9i)    
  W1 <- diag(w1)   
  
  d_bt_bt <-  -b*(b-1)/(2*sigma^b)*( t(X)%*%diag(rho)%*%diag(e2i_ssstar0 )%*%X +
                                       t(X)%*%(diag(n)-diag(rho)) %*% diag(w5)%*%X ) 
  d_bt_sg <- -b/sigma* Q_grad(0, y, X,vi, ratio, theta, small1)[1:p] 
  d_bt_r <- -b^2/2/sigma^b*( t(X)%*% diag(w6)%*% (1-rho) )  
  d_sg_sg <- -n*b/sigma^2 - (b+1)/sigma*Q_grad(0, y, X, vi, ratio, theta, small1)[1+p] 
  d_sg_r <- - b^2/2/sigma^(b+1)*t(w4) %*% (1-rho)  
  d_r_r <- -b*(b+1)/2/sigma^b *t(w8) %*% (1-rho)  - 0.5/(a*b+1)*t(w9) %*%rho  
  d_a_a <-  n*( 1/a - trigamma(a) )     
  
   
  w10i <- function(i){ 
    d <-  my.abs(y[i]-X[i,]%*%beta)/ (1+r*my.sign(y[i]-X[i,]%*%beta)) 
    return(log(d))   
  } 
  w11i <- function(i){ 
    d <-  my.abs(y[i]-X[i,]%*%beta)/ (1+r*my.sign(y[i]-X[i,]%*%beta)) / sigma 
    return(w2i(i) * (log(d))^2 )   
  } 
  
  w10 <- as.matrix( sapply(1:n, w10i),ncol=1 ) 
  w11 <- as.matrix( sapply(1:n, w11i),ncol=1 ) 
  
  d_b_b <- -n/b^2*(1 + (2*log(2)+2*digamma(1/b)) / b + trigamma(1/b)/(b^2) ) + 2*n / b^3 *sum(e3) - 1/2/sigma^b*(1-rho) %*% w11       
  d_beta_b <- (1- b*log(sigma))/2/sigma^b *( t(X)%*%( E2_sst%*%rho + W1%*%(1-rho) ) ) + b/2/sigma^b * ( t(X)%*% W1 %*% diag( sapply(1:n, w10i) ) %*% (1-rho) ) 
  d_sig_b <-  (1- b*log(sigma))/2/sigma^(b+1) *( t(e2)%*%rho + t(w2)%*%(1-rho) ) + b/2/sigma^(b+1) * t(1-rho) %*% diag( sapply(1:n, w2i), ) %*% w10     
  d_r_b <- (1- b*log(sigma))/2/sigma^(b+1) * ( t(1-rho) %*%w4 ) +  b/2/sigma^b * t(1-rho) %*% diag( sapply(1:n, w4i)  ) %*% w10    
  
  MatrizQ<-matrix(0,nrow=(p+4),ncol=(p+4) )  
  MatrizQ[1:p,1:p] <-   d_bt_bt 
  MatrizQ[1:p,(1+p)] <- d_bt_sg 
  MatrizQ[1:p,(2+p)] <- d_bt_r 
  MatrizQ[1:p,(3+p)] <- 0  
  
  MatrizQ[1+p,1+p] <- d_sg_sg  
  MatrizQ[1+p,2+p] <- d_sg_r  
  MatrizQ[1+p,3+p] <- 0  
  MatrizQ[2+p,2+p] <- d_r_r 
  MatrizQ[2+p,3+p] <- 0 
  MatrizQ[3+p,3+p] <- d_a_a     
  
  MatrizQ[4+p,4+p] <- d_b_b  
  MatrizQ[1:p,(4+p)] <- d_beta_b 
  MatrizQ[1+p,(4+p)] <- d_sig_b 
  MatrizQ[2+p,(4+p)] <-  d_r_b
  MatrizQ[3+p,(4+p)] <- 0
  
  MatrizQ[lower.tri(MatrizQ)] <- t(MatrizQ)[lower.tri(MatrizQ)]  
  obj.out <- list( MatrizQ <- MatrizQ, 
                   MatrizQbeta <- d_bt_bt, 
                   MatrizQsigma <- d_sg_sg,
                   MatrizQr<- d_r_r,
                   MatrizQa <- d_a_a, 
                   MatrizQb <- d_b_b )   
  return(obj.out)   
}  


Q.function <- function(y,X, vi, theta,thetak,small1)   
{   
  n <- length(y)   
  p <- ncol(X)    
  theta <- as.vector(theta) ;thetak <- as.vector(thetak)  
  
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)  
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)  
  mu <- X%*%beta    
  big1 <- setdiff(seq(1,n,by=1), small1)  
  
  c1 <- (vi-mu)/sigma    
  Estepi=Estep1(y,X, vi, thetak, ratio,small1 )    
  e1i_o = as.vector(Estepi$e1i_o)    
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0)      
  e2i_star0 = as.vector(Estepi$e2i_star0)    
  e3i_o = as.vector(Estepi$e3i_o)     
  
  cont <- log( a^a*b/2^(1+1/b)/gamma(a)/gamma(1/b)/sigma )     
  e0i_obs <- abs(y[big1]-X[big1,]%*%beta)^(b)*e1i_o[big1]/(1+r*sign(y[big1]-X[big1,]%*%beta))^(b)   
  suma <- n*cont - a*sum(e1i_o) + (a+1/b-1)*sum(e3i_o ) -1/2/sigma^b *( sum(e0i_obs) + sum( e2i_star0[small1] ))   
  suma <- as.numeric( suma ) 
  return(suma)    
}   

QD_apro <- function(y, X, vi, thetak, small1)           
{          
  p <- ncol(X)    
  n <- nrow(X)        
  
  M <- new_HessQ(y,X,vi, thetak, small1)       
  MatrizQ <- -M[[1]]        
  
  Q.theta <- Q.function(y, X, vi, thetak, thetak, small1)   
  Q.theta <- rep(Q.theta,n)    
  
  apply_Qtheta_i <- function(i){    
    Qderi <- Q_grad(i, y, X, vi, ratio, thetak, small1)     
    theta.i <- thetak + solve(MatrizQ) %*%Qderi      
    
    theta.i[ c((p+1),(p+3), (p+4)) ] <- abs( theta.i[ c((p+1),(p+3), (p+4)) ])  
    Q.theta.i<- Q.function(y, X, vi, theta.i,thetak, small1) 
    return(Q.theta.i)  
    
  }    
  Q.theta.i.ret <- sapply(1:n, apply_Qtheta_i)   
  QD <- abs(2*(Q.theta-Q.theta.i.ret))  
  return(QD)  
}  

  
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

calculate_delta_case <- function(y, X, vi, ratio, theta, small1){   
  n <- length(y)  
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)   
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)
  mu <- X%*%beta   
  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  c1 <- (vi-mu)/sigma       
  uc1<-uc(c1,r,a,b)   
  rho <- rep(0,n) ; rho[small1] <- 1     
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)   
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0=as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  
  e2_sst <- diag( e2i_sstar0 )     
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )   
  
  w_s_vars <- calculate_w_var(y, X, beta, sigma, r, a, b, small1, e1i_o, c1) 
  w1 <- w_s_vars$w1
  w2 <- w_s_vars$w2   
  w3 <- w_s_vars$w3   
  w4 <- w_s_vars$w4   
  W1 <- diag(w1)   
  w5 <- w_s_vars$w5   
  w6 <- w_s_vars$w6   
  w7 <- w_s_vars$w7    
  w8 <- w_s_vars$w8  
  w9 <- w_s_vars$w9   
  w10 <- w_s_vars$w10  
  w11 <- w_s_vars$w11   
  s1 <- w_s_vars$s1
  s2 <- w_s_vars$s2 
  
  Delta_beta <- (b / (2 * sigma^b)) * (t(X) %*% (diag(rho) %*% e2_sst) + t(X) %*% diag((1 - rho) * w1) )   
  Delta_sigma <- - (1 / sigma) * rep(1, n) + (b / (2 * sigma^(b + 1))) * (t(rho) %*% diag(e2[,1]) + t(1 - rho) %*% diag(w2[,1]) )
  Delta_r <- - (1/2) * t(rho) %*% diag(w3[,1]) + (b / (2 * sigma^b)) * t(1 - rho) %*% diag(w4[,1])   
  Delta_a <- ((log(a) - digamma(a) + 1) * rep(1, n)) - t(e1) + t(e3)  
  Delta_b <- ((1 / b) + (log(2) / (b^2)) + (digamma(1/b) / (b^2))) * rep(1, n) - (1 / (b^2)) * t(e3) - 
    (1 / (2 * sigma^b)) * t(1 - rho) %*% diag(s1[,1])  - ((a + (1 / b)) / 2) * t(rho) %*% (diag(s2[,1]))  
  
  return(rbind(Delta_beta, Delta_sigma, Delta_r, Delta_a, Delta_b))   
}  
## delta resp   
calculate_delta_resp <- function(y, X, vi, ratio, theta, small1) {
  n <- length(y)  
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)   
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)
  mu <- X%*%beta  
  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  c1 <- (vi-mu)/sigma       
  uc1<-uc(c1,r,a,b)   
  rho <- rep(0,n) ; rho[small1] <- 1  
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)   
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0=as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  
  e2_sst <- diag( e2i_sstar0 )     
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )   
  
  w_s_vars <- calculate_w_var(y, X, beta, sigma, r, a, b, small1, e1i_o, c1)  
  w1 <- w_s_vars$w1
  w2 <- w_s_vars$w2   
  w3 <- w_s_vars$w3   
  w4 <- w_s_vars$w4   
  W1 <- diag(w1)   
  w5 <- w_s_vars$w5   
  w6 <- w_s_vars$w6   
  w7 <- w_s_vars$w7    
  w8 <- w_s_vars$w8  
  w9 <- w_s_vars$w9   
  w10 <- w_s_vars$w10  
  w11 <- w_s_vars$w11   
  s1 <- w_s_vars$s1
  s2 <- w_s_vars$s2 
  
   
  part1 <- ((b - 1) * t(X) %*% (diag(n) - diag(rho)) %*% diag(w5[,1])) 
  part2 <- b * t(1 - rho) %*% (diag(w1))
  part3 <- b * t(1 - rho) %*% (diag(w6[,1]))
  part4 <- matrix(0, nrow = 1, ncol = n)
  part5 <- t(1 - rho) %*% (diag((w1 * w10)[,1]) + ((1 - b * log(sigma) )/ b) * (diag(w1))) 
  
 
  Delta_resp <- (b / (2 * sigma^b)) * rbind(part1, part2, part3, part4, part5) 
  
  return(Delta_resp)
}


##   delta exp   
calculate_delta_exp <- function(y, X, vi, ratio, theta, small1) { 
  n <- length(y)   
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)   
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)
  mu <- X%*%beta   
  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  c1 <- (vi-mu)/sigma       
  uc1<-uc(c1,r,a,b)   
  rho <- rep(0,n) ; rho[small1] <- 1   
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)   
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0=as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  
  e2_sst <- diag( e2i_sstar0 )     
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )     
  S <-  (apply(X, 2, sd))         
  
  w1 <-   w6 <-   w5 <-   w10 <- rep(0, n)  
  
  for (i in 1:n) {
    w1[i] <- e1i_o[i] * my.sign(y[i] - X[i,] %*% beta) *
      my.abs(y[i] - X[i,] %*% beta)^(b - 1) /
      (1 + r * my.sign(y[i] - X[i,] %*% beta))^b 
    w6[i] <- e1i_o[i]*my.abs(y[i]-X[i,]%*%beta)^(b-1)/
      (1+r*my.sign(y[i]-X[i,]%*%beta))^(b+1)    
    w5[i] <- e1i_o[i] * my.abs(y[i] - X[i,] %*% beta)^(b - 2) /
      (1 + r * my.sign(y[i] - X[i,] %*% beta))^b  
    w10[i] <- log(my.abs(y[i] - X[i,] %*% beta) / (1 + r * my.sign(y[i] - X[i,] %*% beta)) ) 
  } 
  
  Delta_exp <- matrix(0, nrow = (p+4), ncol = n)    
  
 
  Delta_exp[1:p,] <- (S %*% t(w1 * (1 - rho)) - (b - 1) * as.numeric(t(S) %*% beta)  * t(X) %*% diag(w5 * (1 - rho)) )
 
  Delta_exp[p+1,] <- -b * as.numeric(t(S) %*% beta)  * t(w1 * (1 - rho)) / sigma  
 
  Delta_exp[p+2,] <- -b * as.numeric(t(S) %*% beta) * t(w6 * (1 - rho))  
   
  Delta_exp[p+4,] <- as.numeric(t(S) %*% beta)  * t(1 - rho) %*% diag(w1 * w10) 
  
  return((b / (2 * sigma ^b)) * Delta_exp)  
}

## delta scale    
calculate_delta_scale <- function(y, X, vi, ratio, theta, small1) {  
  n <- length(y)  
  p <- ncol(X)    
  beta <- signif(theta[1:p],digits=7)     
  sigma <- signif(theta[p+1],digits=7)   
  r <- signif(theta[p+2],digits=7)
  a <- signif(theta[p+3],digits=7)
  b <- signif(theta[p+4],digits=7)
  mu <- X%*%beta  
  
  big1 <- setdiff(seq(1,n,by=1), small1)      
  c1 <- (vi-mu)/sigma       
  uc1<-uc(c1,r,a,b)   
  rho <- rep(0,n) ; rho[small1] <- 1  
  
  Estepi = Estep1(y,X,vi,theta,ratio,small1)   
  e1i_o = as.vector(Estepi$e1i_o)   
  e2i_sstar0 = as.vector(Estepi$e2i_sstar0) 
  e2i_star0=as.vector(Estepi$e2i_star0)
  e3i_o <- as.vector(Estepi$e3i_o)    
  
  e2_sst <- diag( e2i_sstar0 )     
  e2 <- as.matrix( e2i_star0,ncol=1 )   
  e1 <- as.matrix( e1i_o,ncol=1 )  
  e3 <- as.matrix( e3i_o,ncol=1 )   
  
  w_s_vars <- calculate_w_var(y, X, beta, sigma, r, a, b, small1, e1i_o, c1)  
  w1 <- w_s_vars$w1
  w2 <- w_s_vars$w2   
  w3 <- w_s_vars$w3   
  w4 <- w_s_vars$w4   
  W1 <- diag(w1)   
  w5 <- w_s_vars$w5   
  w6 <- w_s_vars$w6   
  w7 <- w_s_vars$w7    
  w8 <- w_s_vars$w8  
  w9 <- w_s_vars$w9   
  w10 <- w_s_vars$w10  
  w11 <- w_s_vars$w11   
  s1 <- w_s_vars$s1
  s2 <- w_s_vars$s2 
  
  D_rho <- diag(rho)
  I_n <- diag(rep(1, n)) 
  
  Delta_beta <- (b^2 / (2 * sigma^b)) * (t(X) %*% (D_rho %*% e2_sst) + t(X) %*% ( (I_n - D_rho) %*% diag(w1)))  
  
  Delta_sigma <- (b^2 / (2 * (sigma^(b + 1)))) * (t(rho) %*% (diag(e2[,1])) + t(1 - rho) %*% (diag(w2[,1])))
  
  Delta_r <- -((b / 2) * t(rho) %*% (diag(w3[,1]))) + ((b^2 / (2 * sigma^b)) * t(1 - rho) %*% (diag(w4[,1])))
  
  Delta_a <- matrix(0, nrow = 1, ncol = n) 
  
  Delta_b <- (-1 / (2 * sigma^b)) * ((1 - b * log(sigma)) * (t(rho) %*% diag(e2[,1]) + t(1-rho) %*% diag(w2[,1]) )  +  
                                       b * t(1 - rho) %*% (diag(w2[,1] * w10[,1])))
  
  return(rbind(Delta_beta , Delta_sigma , Delta_r , Delta_a , Delta_b ))
} 


library(MASS)    

neg_Q_ddot <- function(y, X, vi, ratio, theta, small1, type) { 
  Delta <- switch(type,
                  "case" = calculate_delta_case(y, X, vi, ratio, theta, small1),
                  "resp" = calculate_delta_resp(y, X, vi, ratio, theta, small1),
                  "exp" = calculate_delta_exp(y, X, vi, ratio, theta, small1), 
                  "scale" = calculate_delta_scale(y, X, vi, ratio, theta, small1)
  )  
  HessQ_inv <- ginv(-new_HessQ(y, X, vi, theta, small1)[[1]])
  ret <- 2 * t(Delta) %*% HessQ_inv %*% Delta 
  return(ret) 
}   

calculate_B_fq <- function(y, X, vi, ratio, theta, small1, type) {
  n <- length(y)
  p <- ncol(X)
  key.M <- neg_Q_ddot(y, X, vi, ratio, theta, small1, type = type)
  
  B_fq_list <- vector("numeric", length = n)  
  for (l in 1:n) {
    dl <- rep(0, n)
    dl[l] <- 1
    C_fq <- t(dl) %*% key.M %*% dl
    trM <- sum(diag(key.M))
    B_fq <- C_fq / trM  
    B_fq_list[l] <- abs(B_fq) 
  }
  
  return(B_fq_list)
}

types <- c("case", "resp", "exp", "scale")   

MSE <- function(est,Ttheta){
  d<- sum( ((est-Ttheta)^2)  )    
  return( sqrt(d) )    
}     
 
# mix SGT  
mix_SGT <- function(n, weights, scale1, scale2, Tr, Ta, Tb){  
  mixed_data <- numeric(n)
  for (i in 1:n) {    
    if (runif(1) < weights[1]) { 
      mixed_data[i] <- rsgt(1, location=0, scale= scale1, r=Tr, a=Ta, b=Tb, dp=NULL)    
    } else {  
      mixed_data[i] <- rsgt(1, location=0, scale= scale2, r=Tr, a=Ta, b=Tb, dp=NULL)    
    }
  }    
  return(mixed_data) 
}


# mix N  
mix_N <- function(n, weights, scale1, Tr, Ta, Tb, Nmean, Nsd){   
  mixed_data <- numeric(n)
  for (i in 1:n) {     
    if (runif(1) < weights[1]) {
      mixed_data[i] <- rsgt(1, location=0, scale=scale1, r=Tr, a=Ta, b=Tb, dp=NULL)     
    } else {  
      mixed_data[i] <- rnorm(1, mean = Nmean, sd = Nsd)        
    }    
  }     
  return(mixed_data) 
}
# mix chisq  
mix_chisq <- function(n, weights, scale1, Tr, Ta, Tb, df1){        
  mixed_data <- numeric(n)
  for (i in 1:n) {     
    if (runif(1) < weights[1]) {
      mixed_data[i] <- rsgt(1, location=0, scale=scale1, r=Tr, a=Ta, b=Tb, dp=NULL)      
    } else {
      mixed_data[i] <- rchisq(1, df1)           
    }    
  }        
  return(mixed_data) 
}


    