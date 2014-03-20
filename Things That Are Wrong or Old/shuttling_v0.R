# Gavin Medley
# Time Dependent Solution
# Key pages of book: 469-471
#===TO DO=======
# Figure out why it's straight. Check function outputs
# Add in a function for M-lE_b as a function of T_b

#=====Define constants and put into vector======

K_s <- 205 # conductance of skin
K_f <- 100000000000 # conductance of pelage (fat)
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

#        1   2    3    4  5    6      7   8   9  10   11    12    13   14  15  16
cc <- c(K_s,K_f, K_sf, G, H, T_rbar, T_g, R, Q_n, M, lE_b, lE_r, lE_s, C, K_0, tau)


#=======ODE===========
# Eq (13.21) or (13.48)
ddt <- function(cc,T_b,T_a){
  T_delta <- fM_star(cc)/cc[15] #M_star/K_0
  
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
fM_star <- function(cc){
  # (1 + ((H + R)/K_f)) / (1 + ((H + R)/K_sf))
  bracket1 <- (1 + ((cc[5] + cc[8])/cc[2]))/(1 + ((cc[5] + cc[8])/cc[3]))
  # 1/(1 + ((H + R)/K_sf))
  bracket2 <- 1/(1 + ((cc[5] + cc[8])/cc[3]))
  # M - lE_b - lE_s*bracket1 - lE_r*bracket2
  M_star <- cc[10] - cc[11] - cc[13]*bracket1 - cc[12]*bracket2
  return(M_star)
}

#=======Integration=======
# Euler's method
h <- 1
maxt <- 1000
track <- rep(NA,maxt)
track[1] <- 25
for (t in (1:maxt)){
  track[t+1] <- ddt(cc,track[t],30)*h + track[t]
}

plot(track)