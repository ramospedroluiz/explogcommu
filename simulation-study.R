# Clear all variables from the environment
rm(list=ls(all=TRUE))

# Load necessary packages
require(stats4)      
library(maxLik)      
require(rootSolve)   
require(VGAM)        
require(copula)      

# Set output precision for numerical results
old = options(digits=4)
options(warn=-1)
############ Declaring variables and some parameters ##############

# The following lines contain six different combinations of `pp` and `pbeta`.
# Uncomment the desired combination to execute the code with specific parameter values.
# Also, ensure to uncomment the matching result-saving section at the end of the code.

 pp <- 0.3; pbeta <- 0.5    # To be saved in df1
# pp <- 0.3; pbeta <- 1      # To be saved in df2
# pp <- 0.3; pbeta <- 2      # To be saved in df3
# pp <- 0.5; pbeta <- 0.5    # To be saved in df4
# pp <- 0.5; pbeta <- 1      # To be saved in df5
# pp <- 0.5; pbeta <- 2      # To be saved in df6

# Define number of iterations and sample parameters
B <- 50     # Number of iterations
NN <- 2      # Incremental sample size
ninitial<-50 # Initial sample size
lim <- 50     # Upper limit for parameters in estimations
set.seed(2023) # Set seed for reproducibility

# Matrices to store estimation results across methods (rows: samples, cols: different methods)
mrep <- matrix(nrow=NN, ncol=12)
msep <- matrix(nrow=NN, ncol=12)
mrebeta <- matrix(nrow=NN, ncol=12)
msebeta <- matrix(nrow=NN, ncol=12)
Dabs <- matrix(nrow=NN, ncol=12)   # For D_abs values
Dmax <- matrix(nrow=NN, ncol=12)   # For D_max values

# Initialize arrays to store D_abs and D_max values across methods
D_abs_value_mle <- numeric(B)
D_max_value_mle <- numeric(B)
D_abs_value_me <- numeric(B)
D_max_value_me <- numeric(B)
D_abs_value_mme <- numeric(B)
D_max_value_mme <- numeric(B)
D_abs_value_lse <- numeric(B)
D_max_value_lse <- numeric(B)
D_abs_value_wlse <- numeric(B)
D_max_value_wlse <- numeric(B)
D_abs_value_pce <- numeric(B)
D_max_value_pce <- numeric(B)
D_abs_value_mps1 <- numeric(B)
D_max_value_mps1 <- numeric(B)
D_abs_value_mps2 <- numeric(B)
D_max_value_mps2 <- numeric(B)
D_abs_value_mps3 <- numeric(B)
D_max_value_mps3 <- numeric(B)
D_abs_value_cvm <- numeric(B)
D_max_value_cvm <- numeric(B)
D_abs_value_mad <- numeric(B)
D_max_value_mad <- numeric(B)
D_abs_value_rmad <- numeric(B)
D_max_value_rmad <- numeric(B)

# Function to generate random samples from the exponential-logarithmic distribution
rEL <- function(n, p, beta) {              
  t = rexplog(n, scale = 1/beta, shape = p) 
  return(t)
}

# CDF of the exponential-logarithmic distribution
pEL <- function(t, p, beta){                      
  den <- 1 - log(1 - (1 - p) * exp(-beta * t)) / log(p)
  return(den)
}

###################### Likelihood Function for MLE ##################

# Log-likelihood function for Maximum Likelihood Estimation (MLE)
loglike <- function(theta) {
  p <- theta[1]
  beta <- theta[2]
  aux <- -n * log(-log(p)) + n * log(beta) + n * log(1 - p) - beta * sum(t) - sum(log(1 - (1 - p) * exp(-beta * t)))
  return(aux)
}

###################################### Moment Estimators

# Function for calculating moment estimators
fmoment <- function(p){
  aux <- ((2 * polylog((1 - p), 3, n.sum = 101) * (mean(t)^2) * log(p)) / ((polylog((1 - p), 2, n.sum = 101))^2)) + mean(t^2)
  return(aux)
}

###################################### Modified Moment Estimators

# Function for calculating modified moment estimators
fmmom <- function(p){
  aux <- sqrt((-2 * polylog((1 - p), 3, n.sum = 201) * log(p)) / ((polylog((1 - p), 2, n.sum = 201))^2) - 1) - sd(t) / mean(t)
  return(aux)
}

######################## Maximum Product Spacing Estimators #########################

# Function for Maximum Product Spacings (MPS)
fmps1 <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  D <- numeric(n + 1)
  D[1] <- pEL(t[1], p, beta)
  D[n+1] <- 1 - pEL(t[n], p, beta)
  for (i in 2:n){
    if(t[i] == t[i-1]){
      D[i] <- pEL(t[i], p, beta)
    } else {
      D[i] <- pEL(t[i], p, beta) - pEL(t[i-1], p, beta)
    }
  }
  aux <- 1/(n + 1) * sum(log(D))
  return(aux)
}

# Function for Minimum Product Spacings with |X-Y| form
fmps2 <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  D <- numeric(n + 1)
  D[1] <- pEL(t[1], p, beta)
  D[n+1] <- 1 - pEL(t[n], p, beta)
  for (i in 2:n){
    if(t[i] == t[i-1]){
      D[i] <- pEL(t[i], p, beta)
    } else {
      D[i] <- pEL(t[i], p, beta) - pEL(t[i-1], p, beta)
    }
  }
  aux <- sum(abs(D - 1 / (n + 1)))
  return(-aux)
}

# Function for Minimum Product Spacings with |log X - log Y| form
fmps3 <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  D <- numeric(n + 1)
  D[1] <- pEL(t[1], p, beta)
  D[n+1] <- 1 - pEL(t[n], p, beta)
  for (i in 2:n){
    if(t[i] == t[i-1]){
      D[i] <- pEL(t[i], p, beta)
    } else {
      D[i] <- pEL(t[i], p, beta) - pEL(t[i-1], p, beta)
    }
  }
  aux <- sum(abs(log(D) - log(1 / (n + 1))))
  return(-aux)
}

############################## Anderson Darling Estimator ##################

# Function for Anderson Darling Estimators 
fmad <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  D <- sapply(1:n, function(i) log(pEL(t[i], p, beta)) + log(1 - pEL(t[n + 1 - i], p, beta)))
  aux <- n + (1/n) * sum((2 * seq(1, n, 1) - 1) * D)
  return(aux)
}

# Function for Right-tailed Anderson Darling Estimators 
frmad <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  auxh <- -n/2 + 2 * sum(pEL(t, p, beta)) + (1/n) * sum((2 * seq(1, n, 1) - 1) * log(1 - pEL(t[order(-seq(1, n))], p, beta)))
  return(auxh)
}

############################ Calculation of D_abs and D_max ##############################

# Function to calculate D_abs and D_max for different estimators
calc_D_values <- function(t, est_p, est_beta) {
  F_true <- pEL(t, pp, pbeta)      # True distribution CDF
  F_est <- pEL(t, est_p, est_beta) # Estimated distribution CDF
  
  abs_diff <- abs(F_true - F_est)  # Absolute difference between true and estimated CDF
  
  D_abs <- mean(abs_diff)          # Mean of absolute differences (D_abs)
  D_max <- max(abs_diff)           # Maximum of absolute differences (D_max)
  
  return(list(D_abs = D_abs, D_max = D_max))
}

############################ Main Simulation Loop ###############################

n <- ninitial-10  # Initial sample size
for(k in 1:NN) {
  # Initialize result matrices for different estimators
  mle <- matrix(nrow=B, ncol=2)   # MLE
  me <- matrix(nrow=B, ncol=2)    # Moment estimators
  mme <- matrix(nrow=B, ncol=2)   # Modifiedmoments estimators
  lse <- matrix(nrow=B, ncol=2)   # Least squares estimators
  wlse <- matrix(nrow=B, ncol=2)  # Weighted least squares estimators
  pce <- matrix(nrow=B, ncol=2)   # Percentile-based estimators
  mps1 <- matrix(nrow=B, ncol=2)  # MPS 1
  mps2 <- matrix(nrow=B, ncol=2)  # MPS 2
  mps3 <- matrix(nrow=B, ncol=2)  # MPS 3
  cvm <- matrix(nrow=B, ncol=2)   # Cramer-von Mises estimators
  mad <- matrix(nrow=B, ncol=2)   # Anderson Darling estimators
  rmad <- matrix(nrow=B, ncol=2)  # Right-tailed Anderson Darling estimators
  
  n<-n+10
  y<-seq(1,n,1)/(n+1)                                ## FOR OLSE, PCE
  w1<-((n+1)^2)*(n+2)/(seq(1,n,1)*(n-seq(1,n,1)+1))  ### For WLSE
  y2<-(2*seq(1,n,1)-1)/(2*n)                         ### For CVM Estimators
  o<-0           # iteration counter
  ite<-0         # iteration counter
  
  ################## calculating #######################
  ### Try() command means that even if the program encounters an error it does not stop executing ###
  while(o<B)
  {
    out<-NULL
    t<-sort(rEL(n, pp, pbeta))  ##Generate values from a dist. generalized range
    
    
    #Maximum Likelihood Estimates
    res  <- try(maxBFGS(loglike, start=c(pp, pbeta)), silent=TRUE)
    out<- try(res$estimate, silent=TRUE)
    if(is.double(out[1]) & is.double(out[2]) & out[2]>0 & out[1]>0 & out[2]<lim & out[1]<1){
      
    ### Moment Estimators###
      moma<-NULL
      moma[1]<-try(try(uniroot(f=fmoment, interval=c(0,0.99),extendInt = "yes")$root, silent=TRUE))
      if(is.double(moma[1]) & moma[1]>0 & moma[1]<1){
        moma[2]<- try(-polylog((1-moma[1]),2)/(mean(t)*log(moma[1])), silent=TRUE)
        if(is.double(moma[1]) & is.double(moma[2]) & moma[2]>0 & moma[1]>0 & moma[2]<lim & moma[1]<1){
          
     ## modified moments
      mmoma<-NULL
      mmoma[1]<-try(try(uniroot(f=fmmom, interval=c(0,0.99),extendInt = "yes")$root, silent=TRUE))
      if(is.double(mmoma[1]) & mmoma[1]>0 & mmoma[1]<1){
      mmoma[2]<- try(-polylog((1-mmoma[1]),2)/(mean(t)*log(mmoma[1])), silent=TRUE)
      if(is.double(mmoma[1]) & is.double(mmoma[2]) & mmoma[2]>0 & mmoma[1]>0 & mmoma[2]<lim & mmoma[1]<1 ){
              
       ## OLSE
      res1<-NULL
      res1<-try(coef(nls(y~pEL(t,p,beta),start=list(p=pp,beta=pbeta))), silent=TRUE)
      if(is.double(res1[1]) & is.double(res1[2]) & res1[1]>0 & res1[2]>0 & res1[1]<1 & res1[2]<lim ){
                
      ## WLSE
      res2<-NULL
      res2<-try(coef(nls(y~pEL(t,p,beta),start=list(p=pp,beta=pbeta),weights=w1)), silent=TRUE)
      if(is.double(res2[1]) & is.double(res2[2]) & res2[1]>0 & res2[2]>0 & res2[1]<1 & res2[2]<lim ){
                  
      ## Percentile Estimators
      aper<-NULL
      aper<-try(coef(nls(t~((1/beta)*log((1-p)/(1-p^(1-y)))) ,start=list(p=pp,beta=pbeta), algorithm = "port")), silent=TRUE)
      if(is.double(aper[1]) & is.double(aper[2]) & aper[1]>0 & aper[2]>0  & aper[1]<1 & aper[2]<lim){
                    
                    
      ## Maximum Product Spacings ##
      amps<-try(maxBFGS(fmps1, start=c(pp,pbeta))$estimate, silent=TRUE)
      if(is.double(amps[1]) & is.double(amps[2]) & amps[1]>0 & amps[2]>0 & amps[1]<1 & amps[2]<lim ){
                      
      ## Minimum Product Spacings Type 1  |x-y|##
      amps1<-NULL
      amps1  <- try(maxBFGS(fmps2, start=c(pp,pbeta))$estimate, silent=TRUE)
      if(is.double(amps1[1]) & is.double(amps1[2]) & amps1[1]>0 & amps1[2]>0 &  amps1[1]<1 & amps1[2]<lim ){
                        
      ## Minimum Product Spacings Type 2  |log x-log y|##
      amps2<-NULL
      amps2  <- try(maxBFGS(fmps3, start=c(pp,pbeta))$estimate, silent=TRUE)
      if(is.double(amps2[1]) & is.double(amps2[2]) & amps2[1]>0 & amps2[2]>0 &  amps2[1]<1 & amps2[2]<lim ){
                          
      ### Cramer-Von-Mises Estimators
      res11<-NULL 
      res11<-try(coef(nls(y2~pEL(t,p,beta),start=list(p=pp,beta=pbeta))), silent=TRUE)
      if(is.double(res11[1]) & is.double(res11[2]) & res11[1]>0 & res11[2]>0 & res11[1]<1 & res11[2]<lim ){
                            
      ##Anderson Darling Estimators###
      amad  <- try(maxBFGS(fmad,  start=c(pp,pbeta))$estimate, silent=TRUE)
      if(is.double(amad[1]) & is.double(amad[2]) & amad[1]>0 & amad[2]>0 & amad[1]<1 & amad[2]<lim ){
                              
      ##Right Tailed Anderson Darling Estimators###
      amadr<-NULL 
      amadr  <- try(maxBFGS(frmad, start=c(pp,pbeta))$estimate, silent=TRUE)
      if(is.double(amadr[1]) & is.double(amadr[2]) & amadr[1]>0 & amadr[2]>0 & amadr[1]<1 & amadr[2]<lim ){
  # Increment the iteration counter
  o<-o+1  
  
  ##Extract results 
  mle[o,] <- out
  me[o,] <-  moma
  mme[o,] <- mmoma
  lse[o,] <- res1
  wlse[o,] <- res2
  pce[o,] <- aper
  mps1[o,] <-  amps
  mps2[o,] <-  amps1
  mps3[o,] <-  amps2
  cvm[o,] <-  res11  
  mad[o,] <- amad  
  rmad[o,] <- amadr
                                
  cat("Number of iterations:    ",o,"Sample Size:  ",n,   "\n")
                                
 # Calculate D_abs and D_max for the current estimates
 D_values_mle <- calc_D_values(t, out[1], out[2])
 D_abs_value_mle[o] <- D_values_mle$D_abs
 D_max_value_mle[o] <- D_values_mle$D_max  
                                
 D_values_me <- calc_D_values(t, moma[1], moma[2])
 D_abs_value_me[o] <- D_values_me$D_abs
 D_max_value_me[o] <- D_values_me$D_max  
                                
 D_values_mme <- calc_D_values(t, mmoma[1], mmoma[2])
 D_abs_value_mme[o] <- D_values_mme$D_abs
 D_max_value_mme[o] <- D_values_mme$D_max  
                                
 D_values_lse <- calc_D_values(t, res1[1], res1[2])
 D_abs_value_lse[o] <- D_values_lse$D_abs
 D_max_value_lse[o] <- D_values_lse$D_max  
                                
  D_values_wlse <- calc_D_values(t, res2[1], res2[2])
  D_abs_value_wlse[o] <- D_values_wlse$D_abs
  D_max_value_wlse[o] <- D_values_wlse$D_max  
                                
  D_values_pce <- calc_D_values(t, aper[1], aper[2])
  D_abs_value_pce[o] <- D_values_pce$D_abs
  D_max_value_pce[o] <- D_values_pce$D_max  
                                
  D_values_mps1 <- calc_D_values(t, amps[1], amps[2])
  D_abs_value_mps1[o] <- D_values_mps1$D_abs
  D_max_value_mps1[o] <- D_values_mps1$D_max  
                                
  D_values_mps2 <- calc_D_values(t, amps1[1], amps1[2])
  D_abs_value_mps2[o] <- D_values_mps2$D_abs
  D_max_value_mps2[o] <- D_values_mps2$D_max  
                                
  D_values_mps3 <- calc_D_values(t, amps2[1], amps2[2])
  D_abs_value_mps3[o] <- D_values_mps3$D_abs
  D_max_value_mps3[o] <- D_values_mps3$D_max  
                                
  D_values_cvm <- calc_D_values(t, res11[1], res11[2])
  D_abs_value_cvm[o] <- D_values_cvm$D_abs
  D_max_value_cvm[o] <- D_values_cvm$D_max  
                                
  D_values_mad <- calc_D_values(t, amad[1], amad[2])
  D_abs_value_mad[o] <- D_values_mad$D_abs
  D_max_value_mad[o] <- D_values_mad$D_max  
                                
  D_values_rmad <- calc_D_values(t, amadr[1], amadr[2])
  D_abs_value_rmad[o] <- D_values_rmad$D_abs
  D_max_value_rmad[o] <- D_values_rmad$D_max  
                                
                                
                              } } } } } } } } } } } } }}
    ite<-ite+1
  }
  #############Save MRE and MSE for parameter `p` and `beta` #######################
  
  
  mediamle<-c((sum(mle[,1]/pp)/B),(sum((mle[,1]-pp)^2)/B),(sum(mle[,2]/pbeta)/B),(sum((mle[,2]-pbeta)^2)/B))
  
  mediame<-c((sum(me[,1]/pp)/B),(sum((me[,1]-pp)^2)/B),(sum(me[,2]/pbeta)/B),(sum((me[,2]-pbeta)^2)/B))
  
  mediamme<-c((sum(mme[,1]/pp)/B),(sum((mme[,1]-pp)^2)/B),(sum(mme[,2]/pbeta)/B),(sum((mme[,2]-pbeta)^2)/B))
  
  medialse<-c((sum(lse[,1]/pp)/B),(sum((lse[,1]-pp)^2)/B),(sum(lse[,2]/pbeta)/B),(sum((lse[,2]-pbeta)^2)/B))
  
  mediawlse<-c((sum(wlse[,1]/pp)/B),(sum((wlse[,1]-pp)^2)/B),(sum(wlse[,2]/pbeta)/B),(sum((wlse[,2]-pbeta)^2)/B))
  
  mediapce<-c((sum(pce[,1]/pp)/B),(sum((pce[,1]-pp)^2)/B),(sum(pce[,2]/pbeta)/B),(sum((pce[,2]-pbeta)^2)/B))
  
  mediamps1<-c((sum(mps1[,1]/pp)/B),(sum((mps1[,1]-pp)^2)/B),(sum(mps1[,2]/pbeta)/B),(sum((mps1[,2]-pbeta)^2)/B))
  
  mediamps2<-c((sum(mps2[,1]/pp)/B),(sum((mps2[,1]-pp)^2)/B),(sum(mps2[,2]/pbeta)/B),(sum((mps2[,2]-pbeta)^2)/B))
  
  mediamps3<-c((sum(mps3[,1]/pp)/B),(sum((mps3[,1]-pp)^2)/B),(sum(mps3[,2]/pbeta)/B),(sum((mps3[,2]-pbeta)^2)/B))
  
  mediacvm<-c((sum(cvm[,1]/pp)/B),(sum((cvm[,1]-pp)^2)/B),(sum(cvm[,2]/pbeta)/B),(sum((cvm[,2]-pbeta)^2)/B))
  
  mediamad<-c((sum(mad[,1]/pp)/B),(sum((mad[,1]-pp)^2)/B),(sum(mad[,2]/pbeta)/B),(sum((mad[,2]-pbeta)^2)/B))
  
  mediarmad<-c((sum(rmad[,1]/pp)/B),(sum((rmad[,1]-pp)^2)/B),(sum(rmad[,2]/pbeta)/B),(sum((rmad[,2]-pbeta)^2)/B))
  
  
  mrep[k,]<-c(mediamle[1],mediame[1],mediamme[1],medialse[1],mediawlse[1],mediapce[1],
              mediamps1[1],mediamps2[1],mediamps3[1],mediacvm[1],mediamad[1],mediarmad[1])
  
  msep[k,]<-c(mediamle[2],mediame[2],mediamme[2],medialse[2],mediawlse[2],mediapce[2],
              mediamps1[2],mediamps2[2],mediamps3[2],mediacvm[2],mediamad[2],mediarmad[2])
  
  
  mrebeta[k,]<-c(mediamle[3],mediame[3],mediamme[3],medialse[3],mediawlse[3],mediapce[3],
                 mediamps1[3],mediamps2[3],mediamps3[3],mediacvm[3],mediamad[3],mediarmad[3])
  
  msebeta[k,]<-c(mediamle[4],mediame[4],mediamme[4],medialse[4],mediawlse[4],mediapce[4],
                 mediamps1[4],mediamps2[4],mediamps3[4],mediacvm[4],mediamad[4],mediarmad[4])
  
  # Save D_abs and D_max 
  Dabs[k,]<-c(mean(D_abs_value_mle), mean(D_abs_value_me), mean(D_abs_value_mme),mean(D_abs_value_lse),mean(D_abs_value_wlse),mean( D_abs_value_pce),mean(D_abs_value_mps1),mean(D_abs_value_mps2),mean(D_abs_value_mps3),mean(D_abs_value_cvm),mean(D_abs_value_mad), mean(D_abs_value_rmad))
  Dmax[k,]<-c(mean(D_max_value_mle), mean(D_max_value_me), mean(D_max_value_mme),mean(D_max_value_lse),mean(D_max_value_wlse),mean( D_max_value_pce),mean(D_max_value_mps1),mean(D_max_value_mps2),mean(D_max_value_mps3),mean(D_max_value_cvm),mean(D_max_value_mad), mean(D_max_value_rmad))
}
################################################################# Extracted results in text tables
n_values = seq(ninitial, (ninitial+10*NN), by = 10)


# Combine data into a data frame
df <- data.frame()
for (i in 1:NN) {
  df <- rbind(df,
              data.frame(n = n_values[i], Estimation = "MRE($\\alpha$)", t(mrep[i,])),
              data.frame(n = "", Estimation = "MSE($\\alpha$)", t(msep[i,])),
              data.frame(n = "", Estimation = "MRE($\\beta$)", t(mrebeta[i,])),
              data.frame(n = "", Estimation = "MSE($\\beta$)", t(msebeta[i,])),
              data.frame(n = "", Estimation = "Dabs", t(Dabs[i,])),
              data.frame(n = "", Estimation = "Dmax", t(Dmax[i,]))
  )
}

colnames<- c("n", "Estimates", "MLE", "ME", "MME", "LSE", "WLSE", "PCE", "MPSE", "MSADE", "MSALDE", "CVM", "ADE", "RTADE")
colnames(df) <- colnames

df

#### Save Results in CSV format, remove # in the same order as did in the top to save results 


#write.csv(df,file="df1.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 0.5  

#write.csv(df,file="df2.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 1 

#write.csv(df,file="df3.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 2 

#write.csv(df,file="df4.csv", row.names = TRUE)     #save result for pp <- 0.5; pbeta <- 0.5 

#write.csv(df,file="df5.csv", row.names = TRUE)     #save result forr pp <- 0.5; pbeta <- 1 

#write.csv(df,file="df6.csv", row.names = TRUE)     #save result for pp <- 0.5; pbeta <- 2 