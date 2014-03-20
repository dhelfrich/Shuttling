# Gavin Medley
# Time Dependent Shuttling Model



#------------------Function Definitions-------------------
  
#==============Metabolic rate dependent on Tb==============
fNet_M <- function(Mp,Eb_p,Tp,P,lam,Tb){

  Net_M <- (Mp-lam*Eb_p) + P*(Tb-Tp)
  
  return(Net_M)
}


#================Operative Environmental Temp============
fTe <- function(Ta,Tr_bar,Qn,H,R){
  
  Te <- (H*Ta + R*Tr_bar + Qn)/(H + R) # operative environmental temp
  
  return(Te)
}


#==========Tb_t evaluates Tb into a timeseries=============
#for the given input conditions new and old 
# which are vectors of values for (Te_prime, Td_prime), (Te, Td)
fTb_t <- function(newTa,oldTa,Tb_0){
  
  #------------------Variable Definitions-------------------
  Tr_bar <- 25 # mean radiative surface temp (deg C)
  eps <- .6 # emissivity of beetle shell (decimal fraction)
  sig <- 5.670e-8 # Stefan-Boltzmann constant (W m^-2 K^-4)
  lam <- 2260000 # J/kg latent heat of vaporization 
  R <- 4*eps*sig*(Tr_bar+273.15)^3 # linearized radiation term (W m^-2)
  H <- 1e-6 # convection coefficient (W m^-2 C^-1)
  Qa <- 37500 # absorbed incoming radiation (W m^-2)
  Qn <- Qa - 3*eps*sig*(Tr_bar + 273.15)^4 - 273.15*R # net radiance (W m^-2)
  K <- .5 # conductance of beetle body (W m^-2 C^-1)
  C <- .6 # Heat capacity of beetle body  (J m^-2 C^-1)
  I <- (K + H + R)/(K*(H + R)) # insulation term (m^2 C W^-1)
  maxt <- 100 # number of timesteps
  
  Mp <- 5.432003 # Metabolic rate at preferred temp (W/m^2)
  Eb_p <- .97975 # Respiration rate at preferred temp
  Er <- 0
  Tp <- 25 # temp at which there is a linear relationship between Tb and Net_M (David's #)
  P <- .2797931 # slope of relationship between Net_M and Tb (W/m^2/degC)
  
  # Calculate terms for Tb equation
  Te_prime <- fTe(newTa,Tr_bar,Qn,H,R)
  Te <- fTe(oldTa,Tr_bar,Qn,H,R)
  
  # Initialize timeseries to fill
  Tb_series <- rep(NA,maxt)
  Tb_series[1] <- Tb_0
  
  # Loop to fill the vector of body temps
  for (t in (1:maxt)){
    # Get Net_M values to calculate Td and Td_prime
    Net_M_prime <- fNet_M(Mp,Eb_p,Tp,P,lam,Tb_series[t])
    Net_M <- fNet_M(Mp,Eb_p,Tp,P,lam,Tb_series[t])
    
    # Calculate Td from the old and new temperatures
    Td_prime <- (Net_M_prime*I - lam*Er*(1/(H+R)))
    Td <- (Net_M*I - lam*Er*(1/(H+R)))
    
    Tb_series[t+1] <- (Te_prime + Td_prime) + ((Te+Td)-(Te_prime+Td_prime))*exp(-t/(I*C))
    #if ((Tb_series(t) < 0) || (Tb_series(t) > 30)){break}
  }
  
  return(Tb_series)
}


time_track <- fTb_t(25,28,30)
plot(time_track)

#--------------------Begin Analysis-----------------------

# # Finds the sum of squares when comparing to our real data
# Sum_of_Squares <- function(){
#   n<-Tb_t()
#   
#   I_d <- n[[2]] - I
#   I_ssd <- sum(I_d^2)
#   R_d <- n[[3]] - R
#   R_ssd <- sum(R_d^2)
#   ssd <- I_ssd + R_ssd
#   
#   return(ssd)
# }
# 
# 
# #real_timeseries <- read.csv(file.choose())
# 
# S_real <- real_timeseries[,1]
# I_real <- real_timeseries[,2]
# R_real <- real_timeseries[,3]
# 
# S0 <- S_real[1]
# I0 <- I_real[1]
# R0 <- R_real[1]
# a <- .000113
# g <- .2207
# maxt <- 39
# 
# # Plot the deterministic generated data
# vect <- SIR(S0,I0,R0,a,g,maxt)
# #windows()
# plot(vect[[1]],col="blue",main="WOOOO!",type="n",xlim=c(0,40),ylim=c(0,100),xlab="Time",ylab="# Individuals")
# points(vect[[1]],col="blue")
# points(vect[[2]],col="red")
# points(vect[[3]],col="green")
# 
# points(S_real,col="grey30")
# points(I_real,col="grey60")
# points(R_real,col="grey90")
# 
# maxt <- 40
# p_best <- NA
# ss_best <- Inf
# 
# #alpha_range <- seq(9e-05,1e-4, by= 1e-04)
# #gamma_range <- seq(.05,.15, by=1e-03)
# #alpha <- 9.534617e-05,gamma=2.158401e-01,ss=248358.9
# alpha_range <- exp(seq(log(1e-08),log(1),length.out=10))
# gamma_range <- exp(seq(log(1e-08),log(1),length.out=10))
# 
# #Initialize Vectors
# n <- length(alpha_range)*length(gamma_range)
# results <- matrix(NA,n,3)
# p_best <- NA
# ss_best <- Inf
# y <- seq(log(.000000001),log(.99),length.out=100)
# 
# #Direct search for SS minimum
# j <-1
# 
# for(alpha in alpha_range){
#   for(gamma in gamma_range){
#     
#     p <- c(alpha,gamma)
#     ss <- SIR_ss(S0,I0,R0,a,g,maxt,S,I,R)
#     results[j,] <- c(alpha,gamma,ss)
#     j <- j+1
#     
#     if(ss < ss_best){
#       ss_best <- ss
#       p_best <- cbind(alpha,gamma)
#     }
#   }
#   print(paste(round(100*j/n),"%",sep=""),quote=FALSE)
# }
# 
# #Print Minimum values found from direct search
# ss_best
# p_best
# 
# #-------Plot sum of squares profiles for direct search-----------------------
# #Figure 2
# windows("Process Error Parameter Profile")
# par(mfrow=c(2,2))
# scale <- 100
# plot(results[,1],results[,3],xlab="alpha =",ylab="Sum of squares",ylim=c(ss_best,ss_best+(scale*ss_best)),col="blue")
# 
# # Look at sum of squares based on 'best' a and g
# sum_squares <- SIR_ss(S0,I0,R0,a,g,maxt,S_real,I_real,R_real)
# sum_squares

