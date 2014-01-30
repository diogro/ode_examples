t <- seq(0, 10, by = 0.01)
# parameters
parameters <- c(r = 2, K = 10)

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

library(deSolve)
out <- ode(y = state, times = t, func = logistic, parms = parameters)

plot(out)
