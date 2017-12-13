library(downloader) ##use install.packages to install
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())
filename <- "mice_pheno.csv" 
url <- paste0("https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/",filename)
download(url, destfile=filename)
dat <- read.csv(filename)

controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>%  
  dplyr::select(Bodyweight) %>% unlist

hfPopulation <- filter(dat,Sex == "F" & Diet == "hf") %>%  
  dplyr::select(Bodyweight) %>% unlist

mu_hf <- mean(hfPopulation)
mu_control <- mean(controlPopulation)
print(mu_hf - mu_control)

reject <- function(N, alpha=0.05){
  hf <- sample(hfPopulation,N) 
  control <- sample(controlPopulation,N)
  pval <- t.test(hf,control)$p.value
  pval < alpha
}
B <- 2000
Ns <- seq(5, 50, 5)
power <- sapply(Ns,function(N){
  rejections <- replicate(B, reject(N))
  mean(rejections)
})
N <- 30
alphas <- c(0.1,0.05,0.01,0.001,0.0001)
power <- sapply(alphas,function(alpha){
  rejections <- replicate(B,reject(N,alpha=alpha))
  mean(rejections)
})
plot(alphas, power, xlab="alpha", type="b", log="x")


N <- 30
power30 <- sapply(c(alphas),function(alpha){
  rejections <- replicate(B,reject(N,alpha=alpha))
  mean(rejections)
})

N <- 15
power15 <- sapply(c(alphas),function(alpha){
  rejections <- replicate(B,reject(N,alpha=alpha))
  mean(rejections)
})

for (N in Ns){
  print (N)
}













crossing(x = alphas, N = Ns) %>%
  mutate(power <- sapply(c(alphas,N),function(alpha, N){
    rejections <- replicate(B,reject(N,alpha=alpha))
    mean(rejections)
  })) %>%
  mutate(N = factor(N)) %>%
  ggplot(aes(x, power, color = N, group = N )) +
  geom_line()+
  xlab("alpha")+
  ylab("power")
