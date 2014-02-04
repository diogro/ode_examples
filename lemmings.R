#Something with delay
## libraries
library(deSolve) # to solve
library(ggplot2) # to plot

## parameters: growth rate, carrying capacity and time lag 
r <- 3.5; k <- 19; lag=0.74

## initial conditions: population size slightly
## larger than k at t=0  to start oscilations
yinit <- c(y = K + 0.001)

## Time vector: start sequence with minus (lag value);
## lag should be a multiple of time steps
## given in argument "by" of the sequence
times <- seq(-0.74, 40, by = 0.01)

## function to be derived
## contains a call for thd special function 'lagvalue' from ode library
## that returns the values of the state variable (population size)
## in the time= time - lag
derivs <- function(t, y, parms) {
    ## This part returns the state variable at the time-lag if t>0
        ## For time<0 it returns the carrying capacity
        ## Meaning that the population was at k before t=0
    if (t < 0)
        y.lag <- k
    else
        y.lag <- lagvalue(t - lag)
    ## and here goes the expression to be integrated   
    dy <- r * y * (1 - y.lag/m)
    ## Output is the populatio sizes and the derivatives at each time
    list(dy, dy = dy)
}

## solve the model with 'dede', the solver for models with delays
## see help page for details
out <- dede(y = yinit, times = times, func = derivs, parms = NULL, atol = 1e-10)

# plots
## Coerces output to a dataframe, necessary to ggplot
out.df = as.data.frame(out)
ggplot(data = out.df, aes(time, y)) + geom_line()
ggplot(data = out.df, aes(y, dy.y, color = time)) + geom_point()
