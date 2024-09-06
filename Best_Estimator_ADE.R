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
B <- 1     # Number of iterations
NN <- 14      # Incremental sample size
lim <- 50     # Upper limit for parameters in estimations
set.seed(2023) # Set seed for reproducibility

# Matrices to store estimation results across methods (rows: samples, cols: different methods)
mrep <- matrix(nrow=NN, ncol=1)
msep <- matrix(nrow=NN, ncol=1)
mrebeta <- matrix(nrow=NN, ncol=1)
msebeta <- matrix(nrow=NN, ncol=1)
Dabs <- matrix(nrow=NN, ncol=1)   # For D_abs values
Dmax <- matrix(nrow=NN, ncol=1)   # For D_max values

# Initialize arrays to store D_abs and D_max values across methods
D_abs_value_mad <- numeric(B)
D_max_value_mad <- numeric(B)

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

############################## Anderson Darling Estimator ##################

# Function for Anderson Darling Estimators 
fmad <- function(theta){
  p <- theta[1]
  beta <- theta[2]
  D <- sapply(1:n, function(i) log(pEL(t[i], p, beta)) + log(1 - pEL(t[n + 1 - i], p, beta)))
  aux <- n + (1/n) * sum((2 * seq(1, n, 1) - 1) * D)
  return(aux)
}

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

n <- 10  # Initial sample size
for(k in 1:NN) {
  # Initialize result matrices for different estimators
  mad <- matrix(nrow=B, ncol=2)   # Anderson Darling estimators
  
  n <- n + 10
  y <- seq(1, n, 1) / (n + 1)                                ## FOR OLSE, PCE
  w1 <- ((n + 1)^2) * (n + 2) / (seq(1, n, 1) * (n - seq(1, n, 1) + 1))  ### For WLSE
  y2 <- (2 * seq(1, n, 1) - 1) / (2 * n)                         ### For CVM Estimators
  o <- 0           # iteration counter
  ite <- 0         # iteration counter
  
  ################## calculating #######################
  ### Try() command means that even if the program encounters an error it does not stop executing ###
  while(o < B) {
    out <- NULL
    t <- sort(rEL(n, pp, pbeta))  ## Generate values from a dist. generalized range
    
    ## Anderson Darling Estimators ###
    amad <- try(maxBFGS(fmad, start=c(pp, pbeta))$estimate)
    if(is.double(amad[1]) & is.double(amad[2]) & amad[1] > 0 & amad[2] > 0 & amad[1] < 1 & amad[2] < lim ) {
      o <- o + 1  # Increment the iteration counter                           
      mad[o,] <- amad  
      
      cat(o, "    ", ite + 1, "    ", n,   "\n")
      
      # Calculate D_abs and D_max for the current estimates
      D_values_mad <- calc_D_values(t, amad[1], amad[2])
      D_abs_value_mad[o] <- D_values_mad$D_abs
      D_max_value_mad[o] <- D_values_mad$D_max
    }
    ite <- ite + 1
  } # Closing bracket for while loop
  
  # Save results for each k
  mrep[k, ] <- sum(mad[, 1] / pp) / B
  msep[k, ] <- sum((mad[, 1] - pp) ^ 2) / B
  mrebeta[k, ] <- sum(mad[, 2] / pbeta) / B
  msebeta[k, ] <- sum((mad[, 2] - pbeta) ^ 2) / B
  
  # Save D_abs and D_max 
  Dabs[k, ] <- mean(D_abs_value_mad)   # Notice the [k, ] to store by iteration
  Dmax[k, ] <- mean(D_max_value_mad)   # Same here to store correctly
} # Closing bracket for for loop

################################################################# Extracted results in text tables
n_values = seq(20, 150, by = 10)

# Combine data into a data frame
df <- data.frame()

for (i in 1:NN) {
  df <- rbind(df,
              data.frame(n = n_values[i], Estimation = "MRE($\\alpha$)", ADE = mrep[i]),
              data.frame(n = n_values[i], Estimation = "MSE($\\alpha$)", ADE = msep[i]),
              data.frame(n = n_values[i], Estimation = "MRE($\\beta$)", ADE = mrebeta[i]),
              data.frame(n = n_values[i], Estimation = "MSE($\\beta$)", ADE = msebeta[i]),
              data.frame(n = n_values[i], Estimation = "Dabs", ADE = Dabs[i]),
              data.frame(n = n_values[i], Estimation = "Dmax", ADE = Dmax[i])
  )
}

# Set column names correctly
colnames(df) <- c("n", "Estimation", "ADE")

#### Save Results in CSV format, remove # in the same order as did in the top to save results 

write.csv(df, file="df1.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 0.5  
#write.csv(df, file="df2.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 1 
#write.csv(df, file="df3.csv", row.names = TRUE)     #save result for pp <- 0.3; pbeta <- 2 
#write.csv(df, file="df4.csv", row.names = TRUE)     #save result for pp <- 0.5; pbeta <- 0.5 
#write.csv(df, file="df5.csv", row.names = TRUE)     #save result for pp <- 0.5; pbeta <- 1 
#write.csv(df, file="df6.csv", row.names = TRUE)     #save result for pp <- 0.5; pbeta <- 2
