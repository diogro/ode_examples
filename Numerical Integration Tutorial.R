## Part 1 - How numerical integration works ##
## time intervals: a sequence from zero to ten at 0.5 steps
time <- seq(0, 10, by = 0.5)
## initial condition
x0 <- 0.1
#### The function to be integrated (right-hand expression of the derivative above)
f <- function(x){x * (1.-x)}

#### An empty R vector to store the results
x <- c()
#### Store the initial condition in the first position of the vector
x[1] <- x0

## loop over time: approximate the function at each time step
for (i in 1:(length(time)-1)){
    x[i+1] = x[i] + 0.5 * f(x[i])
}

## plotting 
plot(x~time)
curve(0.1 * exp(x)/(1+0.1*(exp(x)-1.)), add=T)
legend("topleft", c("approximation", "analytical"), 
       pch=c(1,NA), lty=c(NA,1))

## plotting with ggplot2 
library(ggplot2)#load each library once per R session
p <- ggplot(data = data.frame(x = x, t = time), aes(t, x)) + geom_point()
analytic <- stat_function(fun=function(t){0.1 * exp(t)/(1+0.1*(exp(t)-1.))})
print(p+analytic)

## Part 2 - Using Scientific libraries  to do numnerical integration ##
library(deSolve)# loads the library

## 2. 1 First case: a single equation ##
## time sequence
time <- seq(from=0, to=10, by = 0.01)
# parameters: a named vector
parameters <- c(r = 1.5, K = 10)

# initial conditions: also a named vector
state <- c(x = 0.1)

# let's define the right-hand side of the differential equation.
# To be recognized by the integration routines of deSolve
# it must be an R function that computes the values
# of the derivative on a time t 
## There are many ways to do this, but the recommended format is:
# 1. Make a function with three arguments:
# time sequence, state variables and parameters, in this order.
# 2. The function returns a list with results of the function to be integrated
# 3. Inside the R function use 'with(as.list(c(state, parameters){ ... }'
# and include between brackets the function(s) to be integrated
# and then close returning the list of the calculated values
logistic <- function(t, state, parameters){
    with(
        as.list(c(state, parameters)),{
            dx <- r*x*(1-x/K)
            return(list(dx))
        }
        )
}

# Now call the R function 'ode', to perform the integration
# the basic arguments of 'ode' are
# y: the vector of initial conditions
# times: the vector with the time sequence
# func: the R function as described above
# parms: vector of parameters
out <- ode(y = state, times = time, func = logistic, parms = parameters)
# the resulting object has the values of the integration
# at each time point in the vector of times
print(head(out)) # first 6 lines, print function not necessary when pasting in R

## Ploting with ggplot2
p <- ggplot(data = as.data.frame(out), aes(time, x)) + geom_point()
analytic <- stat_function(fun=function(t){0.1*10*exp(1.5*t)/(10+0.1*(exp(1.5*t)-1))})
print(p+analytic)

## 2.2. A system of two equations ##
# time sequence 
time <- seq(0, 50, by = 0.)

# parameters: a named vector
parameters <- c(r = 2, k = 0.5, e = 0.1, d = 1)

# initial condition: a named vector
state <- c(V = 1, P = 3)

# R function to calculate the value of the derivatives at each time value
# Use the names of the variables as defined in the vectors above
lotkaVolterra <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dV = r * V - k * V * P
    dP = e * k * V * P - d * P
    return(list(c(dV, dP)))
  })
}
## Integration with 'ode'
out <- ode(y = state, times = time, func = lotkaVolterra, parms = parameters)

## Ploting
out.df = as.data.frame(out) # required by ggplot: data object must be a data frame
library(reshape2)
out.m = melt(out.df, id.vars='time') # this makes plotting easier by puting all variables in a single column

p <- ggplot(out.m, aes(time, value, color = variable)) + geom_point()
print(p)

p2 <- ggplot(data = out.df[1:567,], aes(x = P, V, color = time)) + geom_point()
print(p2)
