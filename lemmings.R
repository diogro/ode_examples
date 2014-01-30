#Something with delay

library(deSolve)
library(ggplot2)

derivs <- function(t, y, parms) {
  if (t < 0)
    lag <- 19
  else
    lag <- lagvalue(t - 0.74)
  
  dy <- r * y * (1 - lag/m)
  list(dy, dy = dy)
}

# parameters
r <- 3.5; m <- 19

# initial values and times
yinit <- c(y = 19.001)
times <- seq(-0.74, 40, by = 0.01)

# solve the model
out <- dede(y = yinit, times = times, func = derivs, parms = NULL, atol = 1e-10)
out_df = as.data.frame(out)

# plots
ggplot(data = out_df, aes(time, y)) + geom_line()
ggplot(data = out_df, aes(y, dy.y, color = time)) + geom_point()
