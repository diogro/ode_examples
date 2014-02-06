library(deSolve)
library(ggplot2)
library(reshape2)

#time
t <- seq(0, 50, by = 0.01)

# parameters
parameters <- c(r = 2, k = 0.5, e = 0.1, d = 1)

# initial condition: this is an array now!
state <- c(V = 1, P = 3)

lotkaVolterra <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dV = r * V - k * V * P
    dP = e * k * V * P - d * P
    list(c(dV, dP))
  })
}

out <- ode(y = state, times = t, func = lotkaVolterra, parms = parameters)
out_df = as.data.frame(out)
out_m = melt(out_df, id.vars='time')

ggplot(out_m, aes(time, value, color = variable)) + geom_point()
ggplot(data = out_df[1:567,], aes(x = P, V, color = time)) + geom_point()
