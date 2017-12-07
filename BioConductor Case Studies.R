#Examples and exercises from Bioconductor Case studies
#Florian Hahne, Wolfgang Huber, Robert Gentleman, and Seth Falcon




source("http://www.bioconductor.org/biocLite.R")
library("Biobase")
library("ALL") #load ALL package and attach the data
library("genefilter")
data("ALL")

#select sample originating from B-cell tunours by searching the BT variable which distinguishes B cell and T cell tumours

bcell = grep("^B", as.character(ALL$BT))

#Next we want to know which samples are molecular types BCR/ABL or NEG

types = c("NEG", "BCR/ABL")
moltyp = which(as.character(ALL$mol.bio) %in% types)

ALL_bcrneg = ALL[, intersect(bcell, moltyp)]

#reduce the subset of levels of sample annotation which is kep in var factor
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
ALL_bcrneg$BT = factor(ALL_bcrneg$BT)

#Filter genes that were not expressed or the data didn't show any relevant variation to extract any information, feature.exclude="^AFFX" removes the control probes, as a measure of dispersion can be used the IQR

varCut = 0.5
filt_bcrneg = nsFilter(ALL_bcrneg, require.entrez=TRUE,
                       require.GOBP=TRUE, remove.dupEntrez=TRUE,
                       var.func=IQR, var.cutoff=varCut,
                       feature.exclude="^AFFX")
filt_bcrneg$filter.log

################# R review
apropos("mean")
addstr <- function(x) {x = paste0("^",x)
 return(x)
}

#apply functin eaply function apply function to data stored in evironments
library('hgu95av2.db')
hgu95av2MAP$"1001_at"


#extract all of the map location for a particular chromosome

myPos = eapply(hgu95av2MAP, function(x) grep("^17p", x, 
                                             value = TRUE))#we map the chromosome 17 the caret ^ , means that we should match the start of the word
myPos= unlist(myPos)
length(myPos)

#map the probes to any chromosome

mapchromo <- function(chrom) {myPos=eapply(hgu95av2MAP, function(x) grep( addstr(chrom), x, 
                                                  value = TRUE))
}
mapchromo(14)
myPos= unlist(myPos)
length(myPos)

myFindMap = function(mapEnv, which){
  myg = addstr(which)
  a1 = eapply(mapEnv, function(x)
    grep(myg, x, value =TRUE))
  unlist(a1)
}
myFindMap(hgu95av2MAP, 14)

#review environment

e1 = new.env(hash = TRUE)
e1$a = rnorm(10)
e1$b = runif(20)
ls(e1)
xx = as.list(e1)
names(xx)
rm(a, envir = e1)

#
envchrom = new.env(has =TRUE)
envchrom$locations = myFindMap(hgu95av2MAP, 18)

envchrom$strip = function(x) gsub("18", "", x)

myExtract = function(env) env$strip(env$locations)
myExtract(envchrom)[1:5]
