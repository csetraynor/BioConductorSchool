#Cheat to read the GEO microarray SOFT files
#Check out this information for a comprehensive review http://genomicsclass.github.io/book/
#Install GEOquery

source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
update.packages(repos = biocinstallRepos(), ask = FALSE) #update packages



#Access GEO Series Data

#Try now with GSE3141 58 cases of lung cancer doi:10.1038/nature18268
library(GEOquery)


gse <- getGEO('GSE3141', GSEMatrix = TRUE)

#aCCES RAW DATA FROM geo

#If raw data such as .CEL files exist on GEO
filePaths = getGEOSuppFiles("GSE21653")
filePaths

#Acces GSE Data Tables from GEO, access the phenotypic info about the samples

dim(pData(gse[[1]]))
head(pData(gse[[1]])[, 1:3])

#Sometimes GSE include separate data tables with the sample information
df1 <- getGSEDataTables("GSE3494")
lapply(df1, head)


#######################################################################################################
#######################################################################################################


#Loading a GEO Dataset (GDS File)
library(Biobase)
library(GEOquery)

gds858 <- getGEO('GDS858', destdir = )

#We can extract two main things from a GDS object : meta data and table of expression data

#Meta(gds858)
Meta(gds858)$description

colnames(Table(gds858))
#Table(gds858)
Table(gds858)[1:10,1:16]

#we can turn a GDS object into an expression set object using base 2 log 

eset <- GDS2eSet(gds858, do.log2 = TRUE)

eset

sampleNames(eset)
featureNames(eset)
pData(eset)

#Loading a GPL Annotation file

Meta(gds858)$platform

library(Biobase)
library(GEOquery)

#Download a GEOPlatform (GPL) file, put it in the current directory, and load it:
gpl96 <- getGEO('GPL96', destdir=".")

#Or, open an existing GPL file:
#gpl96 <- getGEO(filename='GPL96.soft')

#extract meta and table data

Meta(gpl96)$title
Table(gpl96)[1:10,1:4]
