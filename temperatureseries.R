
#=======Plotting temperatures over time=======
# Euler's method
# The following values create a nice looking graph for testing until we 
# are able to find realistic values.
# T_bInitial <- 35
# h <- 1
# tmax <- 100
# tsteps <- tmax/h
# T_e <- 40
# T_delta <- .3
# tau <- 10

#Initializing temperature list
tempList<- rep(NA,tsteps)
tempList[1]<-T_bInitial

#Discrete version of (13.48)
for (i in (2:tsteps-1)){
#  tempList[i+1] <- ((T_e + T_delta - tempList[i])/tau)*h + tempList[i]
  tempList[i+1] <- ((T_e + T_delta - tempList[i])/tau)*h + tempList[i]
}

plot(tempList)

tempList
