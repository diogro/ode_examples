## Basic Suscetible - Infected - Recovered ODE model
## Libraries
library(deSolve)
library(ggplot2)

## The function to be integrated
sir <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dS = -r*S*I
    dI = r*S*I - a*I
    dR = a*I
    list(c(dS, dI, dR))
  })
}

#time vector
time <- seq(0, 50, by = 0.01)

# parameters: transmission and recovering rates
parameters <- c(r = 0.01, a = 0.2)

# initial population sizes: start with a single infected individual to test invasibility
state <- c(S = 99, I = 1, R = 0)
## In this model Total population is constant, so N = S + I + R
N <- sum(state)

## Integrating
out <- ode(y = state, times = time, func = sir, parms = parameters)

## Plots: matplot is a quick way to plot many temporal series in the same frame
##
matplot(x = out[,1], y= out[,-1], xlab="Time", ylab="N individuals",
        col=1:3, type="l", lty=1)
## Adds a legend
legend(y="center", x="right", c("Suscetible", "Infected", "Recovered"), col=1:3, lty=1)
## Phase space
plot(out[,2:3], xlab="Suscetible", ylab="Infected",
     col=heat.colors(nrow(out)), ylim=c(0,N), xlim=c(0,N))
abline(a=N, b=-1)

# parameters: now try S0 < a/r
parameters <- c(r = 0.01, a = 1)
state <- c(S = 90, I = 10, R = 0)
## integrate
out <- ode(y = state, times = time, func = sir, parms = parameters)

## Plot: disease not spreading
##
matplot(x = out[,1], y= out[,-1], xlab="Time", ylab="N individuals",
        col=1:3, type="l", lty=1)
## Adds a legend
legend("topright", c("Suscetible", "Infected", "Recovered"), col=1:3, lty=1)
