library(tidyverse)

teste <- "a = 21"

tibble::tibble("a = 21" = 2)


a <- c(1.5,-1.5,0.5,-0.5,1,-1) # constantes
mu_0 <- c(8/25,-4,4/9,0,0,3/4)
mu_1 <- c(1/5,-1,1/3,1,1/2,1/4)
y_0 <- 2                # condição inicial
y <- tibble("a=1.5, b0=1, b1=1/2" = rep(0,20),"a=-1.5, b0=1, b1=1/2" = rep(0,20),
            "a=0.5, b0=1, b1=1/2" = rep(0,20),"a=-0.5, b0=1, b1=1/2" = rep(0,20),
            "a=1, b0=1, b1=1/2" = rep(0,20),"a=-1, b0=1, b1=1/2" = rep(0,20),
            t = 1:20)

for (i in 1:5){
  for(t in 1:20){ 
    y[t,i] <- (-a[i])^t*y_0 + mu_0[i] + mu_1[i]*t
  }
}  

for(t in 1:20){ 
  y[t,6] <- (-a[i])^t*y_0 + mu_0[i]*t + mu_1[i]*t^2
}
 
y %>% 
  pivot_longer(cols = 1:6,
               names_to="funcao",
               values_to = "values") %>% 
  ggplot(aes(x = t, y = values))+
  facet_wrap(vars(factor(funcao,levels = unique(funcao))),scales = "free",ncol=2)+
  geom_line(color = "red")+
  geom_point(alpha = .7,size=.5,color = "red")+
  labs(y = NULL)+
  theme_bw()



## Preços

a <- 10
b <- -2.5
a1 <- 2
b1 <- 1.5 

p_eq <- (a1-a)/(b-b1)
t <- 0:10
p_0 <- 1

oferta <- data.frame(p = c(p_0,rep(NA, 10)),
                     s = c(0,rep(NA,10)),
                     d = c(a + b*p_0,rep(NA, 10)),
                     t = 0:10)
for (t in 2:11){
  oferta[t,"p"] <- (p_0 - p_eq)*(b1/b)^(t-1) + p_eq
  oferta[t,"s"] <- a + b*oferta[t-1,"p"]
  oferta[t,"d"] <- a1 + b1*oferta[t,"p"]
}

oferta %>% 
  ggplot(aes(t,p))+
  geom_line(colour = "blue")+
  geom_hline(yintercept = p_eq, linetype= "dashed")

fases <- 
oferta %>%
  mutate("qtd" = s) %>% 
  pivot_longer(cols = s:d,
               names_to = "tipo",
               values_to = "qtde") 

fases <- as.data.frame(fases)
fases$qtd <- c(fases[4:14,"qtd"],rep(NA,11))


fases %>% 
  ggplot(aes(qtd,p))+
  geom_path(aes(color = qtd),arrow = arrow(type = "closed",length = unit(0.15, "cm")))+
  geom_abline(slope = 1/b, intercept = -a/b, color = "green")+
  geom_text(aes(x = 11, y = .3, label = "Demanda"), color = "darkgreen")+
  geom_abline(slope = 1/b1, intercept = -a1/b1, color = "red")+
  geom_text(aes(x = 11, y = 5, label = "Oferta"), color = "darkred")+
  scale_x_continuous(limits = c(0,20))+
  scale_y_continuous(limits = c(0,6))+
  labs(x = "quantidade",
       y = "preço")+
  theme_bw()+
  theme(legend.position = "none")


# 2a ordem

a_1 <- c(1.5,-1.5,0.5,-0.5,1,-1) # constantes
a_2
lambda_1 <- c(0.5,0.5,0.5)
lambda_2 <- c(0.2,1.3,-1.3)
y <- tibble("lambda_1 = 0.5, lambda_2 = 0.2" = rep(0,20),
            "lambda_1 = 0.5, lambda_2 = 1.3" = rep(0,20),
            "lambda_1 = -0.5, lambda_2 = 1.3" = rep(0,20),
            t = 1:20)

for (i in 1:3){
  for(t in 1:20){ 
    y[t,i] <- lambda_1[i]^t + lambda_2[i]^t
  }
}  


y %>% 
  pivot_longer(cols = 1:3,
               names_to="funcao",
               values_to = "values") %>% 
  ggplot(aes(x = t, y = values))+
  facet_wrap(vars(factor(funcao,levels = unique(funcao))),
             scales = "free",
             ncol=1)+
  geom_line(color = "red")+
  geom_point(alpha = .7,size=1.5,color = "red")+
  labs(y = NULL,
       x = "tempo")+
  theme_bw()


# resolvendo característico
a_0 <- 1
a_1 <- -3
a_2 <- 2
coeficientes <- c(a_2,a_1,a_0)

raizes <- polyroot(coeficientes)
Re(raizes)

matriz <-  matrix(c(1,raizes[1],1,raizes[2]),ncol = 2)
iniciais <- c(2,3)

cte <- solve(Re(matriz))%*%iniciais
cte

for(t in 3:20){ 
  y[t,"y"] <- cte[1]*Re(raizes[1])^t + cte[2]*raizes[2]^t
}

# theta = arccos(c/r)

y <- tibble(y = rep(0,100),
            t = 1:100)
y[1,"y"] <- 0   # y_0
y[2,"y"] <- 100   # y_1

a_0 <- 1
a_1 <- -sqrt(3)
a_2 <- 1
coeficientes <- c(a_2,a_1,a_0)

raizes <- polyroot(coeficientes)
r <- Mod(raizes[1])
omega <- Arg(raizes[1])

A_1 <- y[1,"y"]
A_2 <- (y[2,"y"]-r*y[1,"y"]*cos(omega))/(r*sin(omega))

for(t in 3:100){ 
  y[t,"y"] <- r^t*(A_1*cos(omega*t) + A_2*sin(omega*t))
}

#Gráfico
y %>% 
  ggplot(aes(x = t, y = y))+
  geom_line(color = "red")+
  geom_point(alpha = .7,size=1.5,color = "red")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(y = NULL,
       x = "tempo")+
  theme_bw()

