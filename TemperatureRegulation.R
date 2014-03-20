# Gavin Medley, Daniel Helfrich
# Time Dependent Solution for Thermoregulation
# Key pages of book: 469-471
#===TO DO=======

# Add in a function for M-lE_b as a function of T_b
# Investigate why the it overshoots the 
# theoretical asymptote so much (probably just the values we are using)

cat("\014") # Sends ctrl + L to Rstudio to clear the console
graphics.off() # Closes all plot windows
rm(list=ls(all=TRUE)) # Clears all variables

#=====Define constants and put into vector======
K_s <- 205 # conductance of skin
K_f <- 10000000000000000000 # conductance of pelage (fat)
K_sf <- K_s*K_f/(K_s+K_f) # combined conductance term (see p 471 for details)

G <- 0 # conductance with ground
C <- 9947 # heat capacity of body
H <- 8.38 # convection coefficient
T_rbar <- 45 # mean radiative surface temp (distribution in future)
T_g <- 30 # ground temperature

R <- 7.29 # Eq p. 466 4*eps*sig*(T_rbar + 273)^3 linearized radiation term
K_0 <- ((K_sf + G)*(H + R) + G*K_sf)/(K_sf + H + R) # Eq (13.24) clumped insulation term
tau <- C/K_0

M <- 10.5 # basal metabolic rate
lE_b <- 6.28 # respirative energy loss
lE_r <- 0 # evap from pelage surface
lE_s <- 0 # evap from skin surface (lE_r and lE_s are equivalent with no pelage p. 471)

eps <- .000006 # emissivity of the radiative surface
sig <- 5.670e-8 # Stefan-Boltzmann constant (W m^-2 K^-4)
Q_a <- 2500
Q_n <- Q_a - 3*eps*sig*(T_rbar + 273)^4 - 273*R # Eq p 466

h <- 1
max_iter <- 10000

# to keep things straight, here's an addressed vector of coefficients
#        1   2    3    4  5    6      7   8   9  10   11    12    13   14  15  16  17   18
cc <- c(K_s,K_f, K_sf, G, H, T_rbar, T_g, R, Q_n, M, lE_b, lE_r, lE_s, C, K_0, tau, h, max_iter)


#=======ODE===========
# Eq (13.21) or (13.48)
ddt <- function(cc,T_b,T_a){
  T_delta <- fM_star(cc,T_b)/cc[15] #M_star/K_0
  
  dTbdt <- (fT_e(cc,T_a)+T_delta-T_b)/cc[16] #(T_e+T_delta+T_b)/tau
  
  return(dTbdt)
}


#======Calculate T_e=======
# Eq (13.25)
fT_e <- function(cc,T_a){
  # K_sf*(H*T_a + R*T_rbar + Q_n + T_g) + G*T_g*(H+R)
  numerator <- cc[3]*(cc[5]*T_a + cc[8]*cc[6] + cc[9] + cc[7]) + cc[4]*cc[7]*(cc[5]+cc[8])
  # (K_sf + G)*(H + R) + G*K_sf
  denominator <- (cc[3] + cc[4])*(cc[5] + cc[8]) + cc[4]*cc[3]
  T_e <- numerator/denominator
  return(T_e)
}

#=======Calculate M_star=====
# Eq (13.28)
fM_star <- function(cc,T_b){
  # (1 + ((H + R)/K_f)) / (1 + ((H + R)/K_sf))
  bracket1 <- (1 + ((cc[5] + cc[8])/cc[2]))/(1 + ((cc[5] + cc[8])/cc[3]))
  # 1/(1 + ((H + R)/K_sf))
  bracket2 <- 1/(1 + ((cc[5] + cc[8])/cc[3]))
  # (M - lE_b) - lE_s*bracket1 - lE_r*bracket2
  M_star <- fnet_M(cc,T_b) - cc[13]*bracket1 - cc[12]*bracket2
  # NOTE: (M-lE_b) is a function of T_b!!! See fnet_M()
  return(M_star)
}

#=======Calculate (M-lE_b)=========
fnet_M <- function(cc,T_b){
  P <- .973931 # slope of relationship between net_M and T_b
  T_p <- 25 # temp for basal metabolic rate M = cc[10]
  # (M-lE_b) = (M-lE_b)_p + P*(T_b - T_p)
  net_M <- (cc[10]-cc[11]) + P*(T_b-T_p)
  return(net_M)
}

#=======Integration=======
# Runge Kutta 4 Method
# Take values for air temp and initial body temp
integrate <- function (cc, T_b0, T_a, max_iter) {
  h <- cc[17] # Resolution of integration
  T_btrack <- c()
  T_btrack <- append(T_btrack,T_b0)
  
  time_track <- seq(0,max_iter,by=h)
  
  for (i in seq(1,length(time_track)-1,by=1)){
    
    # First we define the terms in the RK4 method k1, k2, k3, k4
    k1 <- ddt(cc,T_btrack[i],T_a)
    k2 <- ddt(cc,T_btrack[i]+h*k1/2,T_a)
    k3 <- ddt(cc,T_btrack[i]+h*k2/2,T_a)
    k4 <- ddt(cc,T_btrack[i]+h*k3,T_a)
    
    # The iterate for RK4 (see wikipedia, or ANY numerical analysis book)
    T_btrack <- append(T_btrack,T_btrack[i]+(1/6)*h*(k1+2*k2+2*k3+k4))
    
    #     if ((T_btrack[i+1]>100) | (T_btrack[i+1]<(-10))){
    #       print("broke")
    #       break}
  }
  
  # Our output from the function is two vectors in a list. One of relevant time series and
  # one of relevant temperatures
  
  series <- rbind(time_track[1:i],T_btrack[1:i]) # rowbinds the two vectors
  
  #print(series,digits=3)
  return(series)
}

#======MAIN DRIVER FUNCTION======
# Takes a matrix of temperature and time values, initial body temp, and constant vector
# returns a plot of the timeseries as the organism shuttles through all air temperatures
driver <- function(Tt_states,T_b0,cc){
  Temp_states <- Tt_states[1,]
  time_states <- Tt_states[2,]
  
  temps <- c(T_b0) # Declare a vector starting with T_b0
  
  # Call integrate() in a for() loop through our vector of air temperatures
  for (j in seq(1,length(Temp_states),by=1)){
    # Call functions above
    series_data <- integrate(cc,T_b0,Temp_states[j],time_states[j])  
    temps <- append(temps,series_data[2,])
    T_b0 <- temps[length(temps)]
  }
  
  times <- h*seq(0,length(temps)-1,by=1) # Just make a vector to plot against. 
  
  windows()
  plot(times,temps,main='Plot of Body Temp over Time',xlab='Time Steps',
       ylab=expression('Body Temperature '*' '*degree*'C'))
}

# Create a vector of states
Tt_states <- matrix(nrow=2,ncol=10)
Tt_states[1,] <- 10*sin(seq(1,10,by=1))-60
Tt_states[2,] <- sample(3:10,10,replace=TRUE)*1000

driver(Tt_states,25,cc)

