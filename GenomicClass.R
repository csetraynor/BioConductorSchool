.libPaths("C:/RFolder/R-3.4.2/library")
#Try now with GSE3141 58 cases of lung cancer doi:10.1038/nature18268
source("http://www.bioconductor.org/biocLite.R")
library(tidyverse)



####################### Download Data (*Not optimized)
library(downloader)
filename <- "mice_pheno.csv"
url <- paste0("https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/", filename)
if (!file.exists(filename)) download(url,destfile=filename)


###################### 

library(GEOquery)


gse <- getGEO('GSE3526', GSEMatrix = TRUE)

e <- gse[[1]]

#######################

population <- read.csv(filename)
population <- unlist(population) # turn it into a numeric

control <- sample(population,12)
mean(control)

control <- sample(population,12)
mean(control)

##12 control mice
control <- sample(population,12)
##another 12 control mice that we act as if they were not
treatment <- sample(population,12)
print(mean(treatment) - mean(control))

n <- 10000
null <- vector("numeric",n)
for (i in 1:n) {
  control <- sample(population,12)
  treatment <- sample(population,12)
  null[i] <- mean(treatment) - mean(control)
}



#########################
dat <- read.csv(filename)
controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>%  
  select(Bodyweight) %>% unlist

ttestgenerator <- function(n){
  cases <- sample(controlPopulation, n)
  control <- sample(controlPopulation,n)
  tstat <- (mean(cases) - mean(control)) /
    sqrt( var(cases) / n + var(control)/ n)
  return(tstat)
}
ttests <- replicate(1000, ttestgenerator(10))

hist(ttests)
qqnorm(ttests)
abline(0,1)

ttests <- replicate(1000, ttestgenerator(3))
qqnorm(ttests)
abline(0,1)

# when the sample size is not large enough and the population values follow a normal distribution, then the t-distribution is a better approximation

ps <- (seq(0, 999) +0.5) / 1000
qqplot(qt(ps, df=2*3-2), ttests, xlim = c(-6,6), ylim=c(-6,6))
abline(0,1)

qqnorm(controlPopulation)
qqline(controlPopulation)

# Monte Carlo simulations in practice, it is much more typical to assume a parametric distribution and generate a population from this, which is called a parametric simulation. This means that we take parameters estimated from the real data (here the mean and the standard deviation), and plug these into a model (here the normal distribution). 

controls <- rnorm(5000, mean=24, sd=3.5)

ttestgenerator <- function(n, mean =24, sd = 3.5){
  cases <- rnorm(n, mean, sd)
  controls <- rnorm(n, mean, sd)
  tstat <- (mean(cases) - mean(controls)) /
    sqrt(var(cases) / n + var(controls)/n)
  return(tstat)
}

set.seed(1)

ttests <- replicate(1000, ttestgenerator(5, 0, 1))
100 * length(ttests[ttests > 2]) /length(ttests)
