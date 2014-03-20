# Calculations for P, Tp and Eb,p
x <- c(15,25,37)

y <- c(181.55,275.31,490.7)
y <- y*6.9733333e-8
# 6.97333e-8 watts = 1 microcalorie/min
# so y is in W now

# assuming a beetle is a cylinder .5mm in diameter and 2mm long
beetle_SA <- .00025^2*pi*2+.0005*pi*.002 # 2*pi*r^2 + 2*pi*r*l
# SA is in m^2
y <- y/beetle_SA
# now y is in W/m^2

linfit <- lm(y ~ x)

plot(x,y)
lines(x,fitted(linfit))

P <- linfit[[1]][2]
# P is the slope in W/m^2/degC
P
