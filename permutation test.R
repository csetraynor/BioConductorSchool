.libPaths("C:/RFolder/R-3.4.2/library")
source("http://www.bioconductor.org/biocLite.R")
library("Biobase")
library(devtools)
library(githubinstall)

# install_github("genomicsclass/dagdata")
# install_github("genomicsclass/ph525x")

library(ph525x)
library(dagdata)

dat=read.csv("femaleMiceWeights.csv")

library(dplyr)
library(downloader)

control <- filter(dat, Diet=="chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat, Diet=="hf") %>% select(Bodyweight) %>% unlist
obsdiff <- mean(treatment)-mean(control)

#Permutation tests take advantage of the fact that if we randomly shuffle the cases and control labels, then the null is true. 

N <- 12
avgdiff <- replicate(1000, {
  all <- sample(c(control, treatment))
  newcontrols <- all[1:N]
  newtreatments <- all[(N+1):(2*N)]
  return(mean(newtreatments) - mean(newcontrols))
})
hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)

# Permutation P-values should never be zero
#the proportion of permutations with larger difference
(sum(abs(avgdiff) > abs(obsdiff)) + 1) / (length(avgdiff) + 1)

#Now let’s repeat this experiment for a smaller dataset. 
N <- 5
control <- sample(control, N)
treatment <- sample(treatment, N)
obsdiff <- mean(treatment) - mean(control)


url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile = filename)
babies <- read.table("babies.txt", header=T)
bwt.nonsmoke <- filter(babies, smoke == 0) %>% select(bwt) %>% unlist
bwt.smoke <- filter(babies, smoke == 1) %>% select(bwt) %>% unlist

N=10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke ,N)
smokers <- sample(bwt.smoke, N)
obs <- mean(smokers) - mean(nonsmokers)

dat <- c(smokers, nonsmokers)
shuffle <- sample(dat)
smokersstar <- shuffle[1:N]
nonsmokersstar <-  shuffle[(N+1):(2*N)]
mean(smokersstar)-mean(nonsmokersstar)

avgdiff <- replicate(1000, {
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <-  shuffle[(N+1):(2*N)]
  return(mean(smokersstar)-mean(nonsmokersstar))
})

(sum(abs(avgdiff) > abs(obs)) +1) / (length(avgdiff) + 1)


#Repeat the above exercise, but instead of the differences in mean, consider the differences in median 
N=10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke ,N)
smokers <- sample(bwt.smoke, N)
obs <- median(smokers) - median(nonsmokers)

dat <- c(smokers, nonsmokers)
shuffle <- sample(dat)
smokersstar <- shuffle[1:N]
nonsmokersstar <-  shuffle[(N+1):(2*N)]
median(smokersstar)-median(nonsmokersstar)

avgdiff <- replicate(1000, {
  shuffle <- sample(dat)
  smokersstar <- shuffle[1:N]
  nonsmokersstar <-  shuffle[(N+1):(2*N)]
  return(median(smokersstar)-median(nonsmokersstar))
})

(sum(abs(avgdiff) > abs(obs)) +1) / (length(avgdiff) + 1)
hist(avgdiff)
abline(v=obs, col="red", lwd=2)
