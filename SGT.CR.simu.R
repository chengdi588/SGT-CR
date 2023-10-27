

graphics.off()   

library(progress)         

Tsigma <- 0.5  ; Tr <-  -0.5                                                        
Ta=  6 ;    Tb=  5            
# Tbeta0=matrix( c(2, 2, 1, 1), ncol=1)         

Tbeta0=matrix( c(1, -1, 2, -2), ncol=1)           
Ttheta <- c( Tbeta0, Tsigma, Tr, Ta, Tb )               

n = 50 ; p <- nrow(Tbeta0)                      
monin = 500                  

x1 = rnorm(n)  ; x2 = rnorm(n)  ; x3 <- rnorm(n)    
X =  cbind(1,x1,x2,x3)   
mu = X%*%Tbeta0        

ratios <- c( 0.05, 0.15, 0.25, 0.35)            


## ---------------- 1 -----------------------------  

ratio <- ratios[1] ;   Estimate<-matrix(data=NA, nrow =monin, ncol = length(Ttheta) )      
time1<-Sys.time()                   

pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = monin)

for(j in 1:monin){    
  repeat{ 
    erro=rsgt(n, location=0, scale=1, r=Tr, a=Ta, b=Tb,dp=NULL)
    y=mu+Tsigma*erro
    tetaj = try( sgt.em1(y,X,ratio)$thetaMLE,silent = TRUE)      
    
    if (!inherits(tetaj, "try-error")) { 
      break
    }   
  }  
  Estimate[j,]=  tetaj                    
  pb$tick()   
}  
pb$terminate()    
Estimate1 <- na.omit(Estimate)    

## ---------------- 2 -----------------------------        
ratio <- ratios[2] ;   Estimate<-matrix(data=NA, nrow =monin, ncol = length(Ttheta) )      
time1<-Sys.time()                   

pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = monin)

for(j in 1:monin){    
  repeat{
    erro=rsgt(n, location=0, scale=1, r=Tr, a=Ta, b=Tb,dp=NULL) 
    y=mu+Tsigma*erro
    tetaj = try( sgt.em1(y,X,ratio)$thetaMLE,silent = TRUE  )    
    
    if (!inherits(tetaj, "try-error")) { 
      break
    }  
  }  
  Estimate[j,]=  tetaj       
  pb$tick()   
} 
pb$terminate()     
Estimate2 <- na.omit(Estimate)  

## ----------------3 ----------------------------  
ratio <- ratios[3] ;   Estimate<-matrix(data=NA, nrow =monin, ncol = length(Ttheta) )      
time1<-Sys.time()                   

pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = monin)

for(j in 1:monin){    
  repeat{
    erro=rsgt(n, location=0, scale=1, r=Tr, a=Ta, b=Tb,dp=NULL)
    y=mu+Tsigma*erro
    tetaj = try( sgt.em1(y,X,ratio)$thetaMLE,silent = TRUE  )    
    
    if (!inherits(tetaj, "try-error")) { 
      break
    }  
  }  
  Estimate[j,]=  tetaj         
  pb$tick()  
} 
pb$terminate()  
Estimate3 <- na.omit(Estimate)  


## ----------------- 4 ----------------------------      
ratio <- ratios[4] ;   Estimate<-matrix(data=NA, nrow =monin, ncol = length(Ttheta) )      
time1<-Sys.time()                   

pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = monin)

for(j in 1:monin){    
  repeat{
    erro=rsgt(n, location=0, scale=1, r=Tr, a=Ta, b=Tb,dp=NULL)
    y=mu+Tsigma*erro
    tetaj = try( sgt.em1(y,X,ratio)$thetaMLE,silent = TRUE  )    
    
    if (!inherits(tetaj, "try-error")) { 
      break
    }  
  }  
  Estimate[j,]=  tetaj           
  pb$tick()  
} 
pb$terminate()     
Estimate4 <- na.omit(Estimate)  
time2 <- Sys.time()                            
time2                   
time2-time1    

 