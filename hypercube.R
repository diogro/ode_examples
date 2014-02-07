## Parameter space exploration with latin hypercube ##
## Examples for ode models ##
## For details please see the vignette of the pse package

## Load libraries
library(deSolve) # numerical integration library
library(rootSolve) # for runsteady integration function
library(pse) # for latin hypercube


### First example: two competitors in a Lotka-Volterra System ###

## 1. Create a function that receive model parameters and time sequence
## and outputs the final size of populations
## Argument runsteady=TRUE use the runsteady function of package rootSolve
## to run numerical integration till a stationary point (hopefully, please see help of the function)
oneRun <- function(X0, Y0, a, b,  times, runsteady=TRUE){
    ## parameters: initial sizes and competition coeficients
    parameters <- c(a = a, b = b)
    ## initial population sizes: start with a single infected individual to test invasibility
    state <- c(X = X0, Y = Y0)
    ## The function to be integrated
    LV <- function(t, state, parameters){
        with(as.list(c(state, parameters)), {
            dX = X*(1 - X - a*Y)
            dY = Y*(1 - Y - b*X)
            list(c(dX, dY))
        })
    }
    ## Integrating: runsteady to run untill convergence (be carefull using that!)
    ## Or the usual integration with the function ode
    if(runsteady){
        out <- runsteady(y = state, times = times, func = LV, parms = parameters)
        return(out$y) # runsteady returns only final states
    }
    else{
        out <- ode(y = state, times = times, func = LV, parms = parameters)
        return(out[nrow(out),-1]) # indexing to get final state values
    }
}

## 2. Create a function with mapply that gets a matrix of parameter combinations
## (combinations are lines and parameters are columns)
## and returns the output of the function above for each parameter combination.
## Use argument 'MoreArgs' to input the values for additional arguments (runsteady and times)
modelRun <- function (my.pars) {
    return(mapply(oneRun, my.pars[,1], my.pars[,2], my.pars[,3], my.pars[,4],
                  MoreArgs=list(times=c(0,Inf), runsteady=TRUE)))
}

## 3. Now prepare the hypercube. You have to specify:
## A character vector to name the parameters
## In this case the initial conditions and the two competition coefficients
factors <- c("X0", "Y0", "a", "b")
## The probability quantile function used to draw values of each parameter.
## If you do not information on the distribution of parameters use uniform distribution
## to have equiprobable values in a range from min to max:
q <- c("qunif", "qunif", "qunif", "qunif")
## A list of the parameters of the above probability distributions
## for the uniform the parameters are minimum and maximun values
## So in this case each element of the list is another list
## with the max and min of the corresponding model parameter, in the same order as above
q.arg <- list(list(min=0, max=1), list(min=0, max=1), list(min=0, max=2), list(min=0, max=2))

## 4. Now create the hypercube and run the model for each combination of parameters
## with LHS, the working horse function of pse package.
## N=200 combinations of parameters, please see help for further information
myLHS <- LHS(model=modelRun, factors=factors, N=200, q=q, q.arg=q.arg, nboot=100) # use nboot=0 if do not need partial correlations

## 5. And now use exploratory statistics to identify interesting regions of the parameter space
## pse package provides some options like
## 5a. Cumulative distribution of outputs
plotecdf(myLHS, stack=TRUE) # about 40% of the population sizes are ~zero, so coexistence is not granted
## 5b. Scatterplots of response x parameter values
plotscatter(myLHS, add.lm=FALSE)
## Again we see that many runs resulted in exclusion of one of the species
## We see also that initial conditions do not affect the final result, but
## competition coefficients do.
## Partial rank correlations (run LHS with nboot to get confidence intervals)
## Which are the non-parametric correlation of each parameter with each output,
## with the effects of the other variables partialed out
## Confidence bars crossing the zero horizontal line indicate no-significant LINEAR partial correlation
plotprcc(myLHS)
## Final population sizes have a strong correlation with competition coefficients


## 6. You get the table of parameter combinations with get.data
hypercube <- get.data(myLHS)

## 7. and the output for each combination with get.results
hypercube <- cbind(hypercube, get.results(myLHS)) # appending columns of outputs to the hypercube matrix
names(hypercube)[5:6] <- c("X", "Y") # cosmetic

##7. And now you are free to use any exploratory tool to identify interesting regions of the parameter space
## There is not a single recipe, but a first guess is to try univariate and then bivariate and multivariate analyses.
## In some case you can also use linear models such as multiple regressions to identify patterns.

## We already know that there is coexistence in this parameter space.
## bu in how many runs?
sum(hypercube$X>1e-6&hypercube$Y>1e-6)

## Which parameter combinations ensue coexistence?
## Univariate exploration with boxplot
## Create a logical variable to flag coexistence
hypercube$coexistence <- hypercube$X>1e-6&hypercube$Y>1e-6
## boxplots
par(mfrow=c(2,2))
boxplot(X0~coexistence, data=hypercube, xlab="Coexistence")
boxplot(Y0~coexistence, data=hypercube, xlab="Coexistence")
boxplot(a~coexistence, data=hypercube, xlab="Coexistence")
boxplot(b~coexistence, data=hypercube, xlab="Coexistence")
par(mfrow=c(1,1))
## There is a upper bound of coefficients to heva coexistence

## Which combination of coefficients ensue coexistence?
## A bivariate scatterplot
plot(a~b, data=hypercube, type="n") ## to scale the plot for the whole parameter space
points(a~b, data=hypercube, subset=hypercube$X>1e-6&hypercube$Y>1e-6, col="blue")
points(a~b, data=hypercube, subset=!(hypercube$X>1e-6&hypercube$Y>1e-6), col="red")
## Coexistence possible only if both  coeficients < 1, roughly
## Further exploration can reveal more!


### A more complicated example: Disease spread in a Suscetible-Infected-Recovered epidemic model ###

## 1. Create a function that receive model parameters and time sequence
## and outputs a logical value: disease has spread or not
oneRun.sir <- function(S0, I0, R0, r, a, time=seq(0, 50, by = 0.01)){
    ## parameters: transmission and recovering rates
    parameters <- c(r = r, a = a)
    ## initial population sizes: start with a single infected individual to test invasibility
    state <- c(S = S0, I = I0, R = R0)
    ## The function to be integrated
    sir <- function(t, state, parameters){
        with(as.list(c(state, parameters)), {
            dS = -r*S*I
            dI = r*S*I - a*I
            dR = a*I
            list(c(dS, dI, dR))
        })
    }
    ## Integrating
    out <- ode(y = state, times = time, func = sir, parms = parameters)
    return(any(out[,3]>I0)) # returns a logical: any infected number larger than initial number?
}

## 2. Create a function that gets matrix of parameter combinations
## (combinations are lines and parameters are columns)
## and returns the output of the integration for each parameter combination.
## To do this, use mapply
## As we are interested in disease spread we set I0 to a small value and R0 to zero
modelRun.sir <- function (my.pars) {
    return(mapply(oneRun.sir, my.pars[,1], 1, 0, my.pars[,2], my.pars[,3]))
}

## 3. Now prepare the hypercube. You have to specify:
## A string vector to name the parameters
factors <- c("S0", "r", "a")
## The probability quantile function used to draw values of each parameter
## for equiprobable values use uniform distribution:
q <- c("qunif", "qunif", "qunif")
## A list of the parameters of the above probability distributions
## for the uniform the parameters are minimum and maximun values
## So in this case Each element of the list are is a list
## with the max and min of the corresponding model parameter
q.arg <- list(list(min=10, max=200), list(min=0.001, max=0.01), list(min=0, max=1))

## 4. Now create the hypercube and run the model for each combination of parameters
## 200 combinations of parameters
myLHS <- LHS(model=modelRun.sir, factors=factors, N=200, q=q, q.arg=q.arg)

## 5. you get the table of parameter combinations with get.data
hypercube <- get.data(myLHS)

## 6. and the output for each combination with get.results
hypercube$Spread <- get.results(myLHS)

## Now explore the relationships between parameters and outputs
## Start looking at the effect of each parameter
## In this case a boxplot of parameter values by output
par(mfrow=c(1,3))
boxplot(S0~Spread, data=hypercube, xlab="Disease spread", ylab="Initial number of susceptibles")
boxplot(r~Spread, data=hypercube, xlab="Disease spread", ylab="Transmission rate")
boxplot(a~Spread, data=hypercube, xlab="Disease spread", ylab="Recovering rate")
par(mfrow=c(1,1))

## You can also explora relatioships between parameters conditioned
## conditioned to the output
## Relationship a X r when disease spread or not
## only populations with initial size < median of all initial size
plot(a~r, data=hypercube, type="n") # to show all ranges of the parameters
points(a~r, data=hypercube, subset=S0<median(S0)&Spread==1, col="red")
points(a~r, data=hypercube, subset=S0<median(S0)&Spread==0, col="blue")

## And you can also explore rate among two parameters and an additional  parameters
plot(I(a/r)~S0, data=hypercube, type="n") # to show all ranges of the parameters
points(I(a/r)~S0, subset=Spread==1, data=hypercube, col="red")
points(I(a/r)~S0, subset=Spread==0, data=hypercube, col="blue")
abline(0,1) ## equivalence line
## BINGO! Disease spreads if S0 > a/r
