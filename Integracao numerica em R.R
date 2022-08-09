
## intervalos de tempo: uma sequência de zero a dez em passos de 0,5
time <- seq(0, 10, by = 0.5)
## condição inicial
x0 <- 0.1
#### A função a ser integrada (expressão à direita da derivada acima)
f <- function(x){x * (1.-x)}
#### Um vetor vazio para armazenar os resultados
x <- c()
#### Armazena a condição inicial na primeira posição do vetor
x[1] <- x0


## loop ao longo do tempo: aproxima a função a cada passo de tempo
for (i in 1:(length(time)-1)){
    x[i+1] = x[i] + 0.5 * f(x[i])
}

#### plotando com ggplot2
library(ggplot2)#carrega cada biblioteca uma vez por sessão R
p <- ggplot(data = data.frame(x = x, t = time), aes(t, x)) + geom_point()
analytic <- stat_function(fun=function(t){0.1 * exp(t)/(1+0.1*(exp(t)-1.))})
print(p+analytic)
