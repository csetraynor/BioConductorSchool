.libPaths("C:/RFolder/R-3.4.2/library")
#source("http://www.bioconductor.org/biocLite.R")
#library("Biobase")
library(devtools)
#library(githubinstall)
# install_github("genomicsclass/dagdata")
# install_github("genomicsclass/ph525x")
library(ph525x)
library(dagdata)

#Lady Tasting Tea

#Two By Two Tables
tab <- matrix(c(3,1,1,3),2,2)
rownames(tab) <- c('Poured Before','Poured After')
colnames(tab) <- c('Guessed Before', 'Guessed After')
tab
fisher.test(tab, alternative = "greater")

#Chi-sqaure Test
#Manhattan plot 
