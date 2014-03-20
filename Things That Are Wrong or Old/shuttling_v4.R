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
T_rbar <- 45 # mean radiative surface temp
T_g <- 20 # ground temperature

R <- 7.29 # Eq p. 466 4*eps*sig*(T_rbar + 273)^3 linearized radiation term
K_0 <- ((K_sf + G)*(H + R) + G*K_sf)/(K_sf + H + R) # Eq (13.24) clumped insulation term
tau <- C/K_0

M <- 10.5 # basal metabolic rate
lE_b <- 6.28 # respirative energy loss
lE_r <- 0 # evap from pelage surface
lE_s <- 0 # evap from skin surface (lE_r and lE_s are equivalent with no pelage p. 471)

eps <- .6 # emissivity of the radiative surface
sig <- 5.670e-8 # Stefan-Boltzmann constant (W m^-2 K^-4)
Q_a <- 0
Q_n <- Q_a - 3*eps*sig*(T_rbar + 273)^4 - 273*R # Eq p 466

# to keep things straight, here's an addressed vector of coefficients
#        1   2    3    4  5    6      7   8   9  10   11    12    13   14  15  16
cc <- c(K_s,K_f, K_sf, G, H, T_rbar, T_g, R, Q_n, M, lE_b, lE_r, lE_s, C, K_0, tau)


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
integrate <- function (cc, T_b0, T_a) {
  h <- 2 # Resolution of integration
  maxt <- 10000
  T_btrack <- rep(NA,maxt)
  T_btrack[1] <- T_b0
  
  time_track <- seq(0,maxt,by=h)
  
  for (i in seq(1,length(time_track)-1,by=1)){
    
    # First we define the terms in the RK4 method k1, k2, k3, k4
    k1 <- ddt(cc,T_btrack[i],T_a)
    k2 <- ddt(cc,T_btrack[i]+h*k1/2,T_a)
    k3 <- ddt(cc,T_btrack[i]+h*k2/2,T_a)
    k4 <- ddt(cc,T_btrack[i]+h*k3,T_a)
    
    # The iterate for RK4 (see wikipedia, or ANY numerical analysis book)
    T_btrack[i+1] <- T_btrack[i]+(1/6)*h*(k1+k2+k3+k4)
    
#     if ((T_btrack[i+1]>100) | (T_btrack[i+1]<(-10))){
#       print("broke")
#       break}
  }
  
  # Our output from the function is two vectors in a list. One of relevant time series and
  # one of relevant temperatures
  series <- list()
  series[[1]] <- time_track[1:i]
  series[[2]] <- T_btrack[1:i]
  #print(series,digits=3)
  return(series)
}


#======Call Iteratively=========
# Declare a vector of temperatures to go through
T_states <- c(30,50,160,200)

# THE Starting Body Temperature
T_b0 <- 28

temps <- list() # Just declare this as a list ahead of time

# Call integrate() in a for() loop through our vector of air temperatures
for (T_a_now in T_states){
  # Call functions above
  series_data <- integrate(cc,T_b0,T_a_now)  
  temps <- append(temps,series_data[[2]])
  T_b0 <- as.numeric(temps[length(temps)])
}

times <- seq(0,length(temps)-1,by=1)

windows()
plot(times,temps,main='Plot of Body Temp over Time',xlab='Time Steps',
           ylab=expression('Body Temperature '*' '*degree*'C'))

