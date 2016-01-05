## Time sequence
t <- seq(0, 10, by = 0.01)
# parameters
parameters <- c(r = 1.5, K = 10)

# initial condition
state <- c(x = 0.1)

# let's define the right-hand side of the differential equation
# It must be a function of the dependent variable (x) and of the 
# time (t), even if time does not appear explicitly
# this is how you define a function:
logistic <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dx <- r*x*(1-x/K)
    list(dx)
  })
}

## Necessary R packages
library(deSolve)
library(ggplot2)

## Run the numerical integration
out <- ode(y = state, times = t, func = logistic, parms = parameters)

## Plots
## Standard plot library
plot(out, lwd=5, col="lightblue")
curve(0.1*10*exp(1.5*x)/(10+0.1*(exp(1.5*x)-1)), add=TRUE, lty=2, lwd=1.5)
legend("topleft", c("Numerical", "Analytical"), lty=c(1,2), col=c("lightblue", "black"), cex=1.5)
## ggplot2 library
ggplot(data = as.data.frame(out), aes(t, x)) + geom_point()
