##Installing package if necessary
##install.packages("deSolve")

##Loading the necessary libraries
library(deSolve)

### Rosenzweig-MacArthur model###

## ODE system to be integrated
RM <- function(t, y, parms){
    with(as.list(c(y, parms)),
         {
             dR = r*R*(1-R/K) - a*R*C/(1+a*h*R)
             dC = e*a*R*C/(1 + a*h*R) - d*C
             return(list(c(R=dR,C=dC)))
         })
}

##initial parameters
parms1 = c(r=1, K=10, a=1, h=0.1, e=0.1, d=0.1)

##initial population sizes
y0 = c(R=1,C=1)

##running the numerical solution from time 1 to 1000
out = ode(y=y0, times = seq(from=1, to=1000, by=0.5), func = RM, parms = parms1)

##Looking at the out object. Notice it is a matrix
head(out)

##plotting the solution
matplot(x=out[,1], y=out[,2:3], type = "l", lwd=2, lty=1,
        col=c("blue", "darkgreen"), xlab="Time", ylab="Population size")
legend("topright", legend = c("Resource","Consumer"), lty = 1, lwd=2, col=c("blue", "darkgreen"))

## Phase space flow and fixed points

##first we need to calculate the coordinates for the vectors
xs = seq(from=0.9, to=1.3, length.out=7)
ys = seq(from=0.94, to=1.04, length.out=7)
coords = expand.grid(R=xs, C=ys)
## And then loop over all coordinates to calculate the derivatives at that points
dvs = matrix(NA, ncol=2, nrow=nrow(coords))
for(i in 1:nrow(coords))
    dvs[i,] = unlist(RM(t=1, y = coords[i,], parms = parms1))

##now we plot the trajectory
plot(x=out[,2], y=out[,3], type="l", lwd=2, col="blue",
     xlab="Resource",ylab="Consumer", xlim=c(0.9,1.3), ylim=c(0.94,1.04))
##and add the vector field
arrows(x0=coords[,1], y0=coords[,2], x1=coords[,1]+dvs[,1]*0.5,
       y1=coords[,2]+dvs[,2]*0.5, length=0.1, lwd=2)

##Messing a little with the parameters
parms2 = c(r=1, K=15, a=1, h=0.1, e=0.1, d=0.1)

out2 = ode(y=y0, times = seq(from = 1, to = 1000, by=0.5), func = RM, parms = parms2)

##plotting the solution
matplot(x=out2[,1], y=out2[,2:3], type = "l", lwd=2, lty=1,col=c("blue", "darkgreen"),
        xlab="Time", ylab="Population size")
legend("topleft", legend = c("Resource","Consumer"), lty = 1, lwd=2,
       col=c("blue", "darkgreen"))

##calculating the vectors again
xs = seq(from=0,to=6,length.out=8)
ys = seq(from=0,to=2,length.out=8)
coords = expand.grid(R=xs,C=ys)
dvs = matrix(NA, ncol=2, nrow=nrow(coords))
for(i in 1:nrow(coords))
    dvs[i,] = unlist(RM(t=1, y = coords[i,], parms = parms2))

##The trajectory
plot(x=out2[,2],y=out2[,3], type="l", lwd=2, col="blue",
     xlab="Resource", ylab="Consumer", xlim=c(0,6),ylim=c(0,2))
##and vectors
arrows(x0=coords[,1], y0=coords[,2], x1=coords[,1]+dvs[,1]*0.5,
       y1=coords[,2]+dvs[,2]*0.5, length=0.1, lwd=2)

##plotting minimum and maximum population sizes with different K values

##object with the line numbers we will use for the following plot
lines = (nrow(out)-500):nrow(out)

##creating an empty plot
plot(0.1, type="n", xlim=c(0,20), ylim=range(out2[lines,2:3]), log="y",
     xlab="K", ylab="Min and Max population")

##points for the K = 10
points(x = c(10,10), y = range(out[lines,3]), pch=21,bg="darkgreen", cex=3)
points(x = c(10,10), y = range(out[lines,2]), pch=21,bg="blue", cex=3)

##points for the K = 15
points(c(15,15), range(out2[lines,3]), pch=21,bg="darkgreen", cex=3)
points(c(15,15), range(out2[lines,2]), pch=21,bg="blue", cex=3)


##The paradox of enrichment, varying K values

KK = seq(from = 0.5, to=25, by=0.5)
rminmax = matrix(NA, ncol=2, nrow=length(KK))#resource minimum and maximum
cminmax = matrix(NA, ncol=2, nrow=length(KK))#consumer minimux ans maximum

## Loop over all values of K andd get min and max population sizes
for(i in 1:length(KK)){
    parmsi = c(r=1, K=KK[i], a=1, h=0.1, e=0.1, d=0.1)  
    y0 = c(R=1,C=1)
    out3 = ode(y=y0, times = seq(from = 1, to = 1000, by=0.5), func = RM, parms = parmsi)
    rminmax[i,] = range(out3[(nrow(out3)-500):nrow(out3),2])
    cminmax[i,] = range(out3[(nrow(out3)-500):nrow(out3),3])
}
## Plot of min and max population sizes
plot(x=KK, y=rminmax[,1], type="l", lwd=2, col="blue",ylim=range(rminmax), log="y",
     xlab="K", ylab="Min and Max population")
points(x=KK, y=rminmax[,2], type="l", lwd=2, col="blue")
points(x=KK, y=cminmax[,1], type="l", lwd=2, col="darkgreen",ylim=range(rminmax))
points(x=KK, y=cminmax[,2], type="l", lwd=2, col="darkgreen",ylim=range(rminmax))


### Consumer-resource dynamics in a seasonal environment ###
## time sequence 
time <- seq(0, 2000, by = .5)

## parameters: a named vector
parameters <- c(r0=1, alpha=0.1, T=80, K=10, a=1, h=0.1, e=0.1, d=0.1)

## initial conditions: a named vector
state <- c(R = 1, C = 1)

## R function to calculate the value of the derivatives at each time value
## Use the names of the variables as defined in the vectors above
RM2 <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
        r = r0 * (1 + alpha*sin(2*pi*t/T))
        dR = R * ( r*(1 - R/K) - a*C / (1 + a*h*R) )
        dC = e*a*R*C / (1 + a*h*R) - d*C
        return(list(c(dR, dC)))
    })
}

## Integration with 'ode'
out <- ode(y = state, times = time, func = RM2, parms = parameters)

## Ploting with matplot
matplot(x = out[,1], y = out[,2:3], type="l", lwd=2, lty = 1,
        col=c("blue", "darkgreen"), xlab = "Time", ylab = "Population size")
legend("topright", c("Resource", "Consumer"), lty=1, lwd=2, col=c("blue", "darkgreen"))

## A resonance diagram ##
## New time sequence 
time <- seq(0, 6000, by = 1)
## Sequence of values of T
TT <- seq(1, 80, by = 2)
## A matrix to store the results
results <- matrix(ncol=4, nrow=length(TT),
                  dimnames=list(NULL, c("R.min","R.max","C.min","C.max")))
## Loop over all values in TT
for(i in 1:length(TT)){
    parameters <- c(r0=1, alpha=0.1, T=TT[i], K=10, a=1, h=0.1, e=0.1, d=0.1)
    tmp1 <- ode(y = state, times = time, func = RM2, parms = parameters)
    results[i,1:2] <- range(tmp1[1001:nrow(tmp1), 2])
    results[i,3:4] <- range(tmp1[1001:nrow(tmp1), 3])
}

## Plot of resonance diagram
plot(R.min ~ TT , data=results, type="l", lwd=2, lty = 1,
     col="blue", xlab = "T", ylab = "Min / Max population size",
     log="y", ylim = range(results))
lines(R.max ~ TT, data=results,  type="l", lwd=2, lty = 1, col=c("blue"))
lines(C.min ~ TT, data=results,  type="l", lwd=2, lty = 1, col=c("darkgreen"))
lines(C.max ~ TT, data=results,  type="l", lwd=2, lty = 1, col=c("darkgreen"))   
legend("topright", c("Resource", "Consumer"), lty=1, lwd=2, col=c("blue", "darkgreen"))
