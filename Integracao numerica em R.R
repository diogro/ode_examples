## ----Exemplo_Euler------------------------------------------------------------
# intervalos de tempo: uma sequência de zero a dez em passos de 0,5
time <- seq(0, 10, by = 0.5)
# condição inicial
x0 <- 0.1
## A função a ser integrada (expressão à direita da derivada acima)
f <- function(x){x * (1.-x)}
## Um vetor vazio para armazenar os resultados
x <- c()
## Armazena a condição inicial na primeira posição do vetor
x[1] <- x0


# loop ao longo do tempo: aproxima a função a cada passo de tempo
for (i in 1:(length(time)-1)){
    x[i+1] = x[i] + 0.5 * f(x[i])
}

## plotando com ggplot2
library(ggplot2)#carrega cada biblioteca uma vez por sessão R
p <- ggplot(data = data.frame(x = x, t = time), aes(t, x)) + geom_point()
analytic <- stat_function(fun=function(t){0.1 * exp(t)/(1+0.1*(exp(t)-1.))})
print(p+analytic)


## -----------------------------------------------------------------------------
# paramentros: devem estar em um vetor
parameters <- c(r = 1.5, K = 10)

# condições iniciais: também devem estar em um vetor, mesmo que seja um só valor
state <- c(x = 0.1)


## -----------------------------------------------------------------------------
## seqüência de tempo
time <- seq(from=0, to=10, by = 0.01)


## -----------------------------------------------------------------------------
## A ODE logística a ser integrada
logistic <- function(t, state, parameters){
    with(
        as.list(c(state, parameters)),{
            dx <- r*x*(1-x/K)
            return(list(dx))
        }
        )
}


## -----------------------------------------------------------------------------
library(deSolve)# Carrega a biblioteca para integração, basta chamar uma vez por seção do R
## Executa a integração
out <- ode(y = state, times = time, func = logistic, parms = parameters)


## -----------------------------------------------------------------------------
head(out) # primeiras 6 linhas


## -----------------------------------------------------------------------------
#### Para o jupyter notebook apenas
options(jupyter.plot_mimetypes = 'image/png', repr.plot.height=5)
####

plot(out, lwd=6, col="lightblue", main="", ylab="N(t)")
curve(0.1*10*exp(1.5*x)/(10+0.1*(exp(1.5*x)-1)), add=TRUE)
legend("topleft", c("Numerical", "Analytical"), lty=1, col=c("lightblue", "black"), lwd=c(6,1))


## -----------------------------------------------------------------------------
## Plotando com o ggplot2
p <- ggplot(data = as.data.frame(out), aes(time, x)) + geom_point()
analytic <- stat_function(fun=function(t){0.1*10*exp(1.5*t)/(10+0.1*(exp(1.5*t)-1))})
print(p+analytic)


## -----------------------------------------------------------------------------
# Paramêtros: vetor
parameters <- c(r = 2, k = 0.5, e = 0.1, d = 1)

# condições iniciais: vetor
state <- c(V = 1, P = 3)


## -----------------------------------------------------------------------------
# sequência do tempo
time <- seq(0, 50, by = 0.01)

# Paramêtros: vetor
parameters <- c(r = 2, k = 0.5, e = 0.1, d = 1)

# condições iniciais: vetor
state <- c(V = 1, P = 3)


## -----------------------------------------------------------------------------
# Função R para calcular o valor das derivadas a cada valor de tempo
# Use os nomes das variáveis conforme definido nos vetores acima
lotkaVolterra <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    dV = r * V - k * V * P
    dP = e * k * V * P - d * P
    return(list(c(dV, dP)))
  })
}


## -----------------------------------------------------------------------------
## Integração com'ode'
out <- ode(y = state, times = time, func = lotkaVolterra, parms = parameters)

## Plotando
out.df = as.data.frame(out) # exigido pelo ggplot
library(reshape2)
out.m = melt(out.df, id.vars='time') # isso facilita a plotagem colocando todas as variáveis em uma única coluna

p <- ggplot(out.m, aes(time, value, color = variable)) + geom_point()
print(p)


## -----------------------------------------------------------------------------
p2 <- ggplot(data = out.df[1:567,], aes(x = P, V, color = time)) + geom_point()
print(p2)

