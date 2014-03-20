Version Control

v0: Uses Euler's method to integrate the equation. Values we use are loosely the iguana values from the book.

v1: changed to RK4 integration for accuracy

v2: Added the functionality of a T_b dependent net metabolic rate

v3: Functionalized so we can call "integrate()" and get back temperature values. It breaks once it gets past the tolerance temperatures that we can set. 

v4: We can now integrate over multiple different environments and see the output.

v5: We added h and max_iter to the cc vector of constants. Also fixed RK4 coefficients.

v6: Math is now all done in vectors and matrices rather than lists (ew). driver(Tt_states,T_b0,cc) is now a driver function where Tt_states is a matrix with a row of temperature states and a row of times over which to integrate. driver() returns a plot of the entire timeseries.

v7: IDEAS: in addition to passing air temps and times, also pass radiation levels, ground temperatures, conductances, anything else (look for them all) that might change from one evironment to another.