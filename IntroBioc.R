.libPaths("C:/RFolder/R-3.4.2/library")
source("http://www.bioconductor.org/biocLite.R")
library("Biobase")
library(devtools)
library(githubinstall)

# install_github("genomicsclass/dagdata")
# install_github("genomicsclass/ph525x")

library(ph525x)
library(dagdata)


dir <- system.file(package="dagdata") #extracts the location of package
list.files(dir)
list.files(file.path(dir,"extdata")) 
list.files(file.path(dir, "script"))

filename <- file.path(dir,"extdata/femaleMiceWeights.csv")
dat <- read.csv(filename)


BiocStyle::markdown()
library(knitcitations)
library(bibtex)
allbib = read.bibtex("allbib.bib")

##########################################
library(GEOquery)

x = getGEO("GSE5859")

#how are the entities related
all.equal(sampleInfo$filename, colnames((geneExpression)))

#Likewise, the rownames of geneExpression coincide exactly with the PROBEID field of geneAnnotation.
options(digits = 2)
cbind(sampleInfo[1:3,], colnames(geneExpression)[1:3],
      t(geneExpression)[1:3,1:4])
#Binding the tables together in an ExpressionSet
rownames(sampleInfo) = sampleInfo$filename
rownames(geneAnnotation) = geneAnnotation$PROBEID

#Now we make the ExpressionSet.

library(Biobase)

es5859 = ExpressionSet(assayData = geneExpression)
pData(es5859) = sampleInfo
fData(es5859) = geneAnnotation
es5859

#One of the nice things about this arrangement is that we can easily select features using higher level concepts annotated in the fData and pData components. For example to obtain expression data for genes on the Y chromosome only:

es5859[which(fData(es5859)$CHR=="chrY"),]

#The full set of methods to which ExpressionSet instances respond can be seen using

methods(class = "ExpressionSet")

#The most important methods are
# exprs(): get the numerical expression values
# pData(): get the sample-level data
# fData(): get feature-level data
# annotation(): get a tag that identifies nomenclature for feature names
# experimentData(): get a MIAME-compliant metadata structure

annotation(es5859) = "hgfocus.db" # need to look at GSE record in GEO, and know .db

#Second, acquire a MIAME-compliant document of metadata about the experiment.

library(annotate)

mi = pmid2MIAME("17206142")
experimentData(es5859) = mi
es5859

nchar(abstract(es5859))
substr(abstract(es5859),1,95)

#GEO, GEOquery, ArrayExpress for expression array archives

library(GEOquery)
glioMA = getGEO("GSE78703")[[1]]

glioMA

#ArrayExpress: searching and harvesting from EMBL-EBI ArrayExpress

library(ArrayExpress)
sets = queryAE(keywords = "glioblastoma", species = "homo+sapiens")
dim(sets)

initdir = dir()
if (!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")

#SummarizedExperiment: accommodating more diverse feature concepts

methods(class="SummarizedExperiment")


######IRanges and GRanges

library(IRanges)
ir <- IRanges(5,10)
ir
#A single IRanges object can hold more than one range. 
IRanges(start=c(3,5,17), end=c(10,8,20))
#Intra-range operations
ir
shift(ir, -2)

narrow(ir, end=5)
flank(ir, width=3, start = TRUE, both = FALSE)
ir

#Inter-range operations
(ir <- IRanges(start = c(3,5,17), end=c(10,8,20)))
range(ir)
reduce(ir)
gaps(ir)
disjoin(ir)


#GRanges
# GRanges are objects which contain IRanges and two more important pieces of information:
# the chromosome we are referring to (called seqnames in Bioconductor)
# the strand of DNA we are referring to

library(GenomicRanges)
#Let’s create a set of two ranges on a made-up chromosome, chrZ. 
gr <- GRanges("chrZ", IRanges(start = c(5,10), end=c(35,45)),
              strand = "+", seqlengths = c(chrZ=100L))

genome(gr) <- "hg19"
gr
seqnames(gr)
seqlengths(gr)

mcols(gr)
mcols(gr)$value <- c(-1,4)
gr

mcols(gr)$value <- NULL

#GRangesList
gr2 <- GRanges("chrZ", IRanges(11:13,51:53))
grl <- GRangesList(gr, gr2 )
grl

length(grl)

elementNROWS(grl)
grl[[1]]

width(grl)
sum(width(grl))

# We can add metadata columns as before, now one row of metadata for each GRanges object, not for each range. It doesn’t show up when we print the GRangesList, but it is still stored and accessible with mcols.

mcols(grl)$value <- c(5,7)
grl
mcols(grl)

# findOverlaps and %over%

(gr1 <- GRanges("chrZ", IRanges(c(1,11,21,31,41), width = 5),strand = "*"))
(gr2 <- GRanges("chrZ", IRanges(c(19,33), c(38,35)), strand="*"))
fo <- findOverlaps(gr1,gr2)
queryHits(fo)   

#Rle and Views
(r <- Rle(c(1,1,1,0,0,-2,-2,-2,rep(-1,20))))
str(r)
as.numeric(r)

#A Views object can be thought of as “windows” looking into a sequence.
(v <- Views(r, start = c(4,2), end= c(7,6)))
str(v)


#Applications with genomic elements: strand-aware operations

ir <- IRanges(c(3,8,14,15,19,34,40),
              width = c(12,6,6,15,6,2,7))
par(mfrow=c(4,1),mar=c(4,2,2,2))
plotRanges(ir, xlim=c(0,60))
plotRanges(reduce(ir), xlim=c(0,60))
plotRanges(disjoin(ir), xlim=c(0,60))
plotRanges(gaps(ir), xlim=c(0,60))
# reduce(x) produces a set of nonoverlapping ranges that cover all positions covered by x. This can be used to reduce complexity of a gene model with many transcripts, where we may just want the addresses of intervals known to be transcribed, regardless of transcript of residence.
# 
# disjoin(x) produces a set of ranges that cover all positions covered by x, such that none of the ranges in the disjoin output overlaps any end points of intervals in x. This gives us the largest possible collection of contiguous intervals that are separated wherever the original set of intervals had an endpoint.
# 
# gaps(x) produces a set of ranges covering the positions in [start(x), end(x)] that are not covered by any range in x. Given coding sequence addresses and exon intervals, this can be used to enumerate introns.

library(GenomicRanges)
library(ggbio)
library(trackViewer)
gir=GRanges(seqnames = "chr1", ir, strand = c(rep("+",4), rep("-", 3)))
par(mfrow=c(4,1), mar=c(4,2,2,2))
plotGRanges(gir, xlim=c(0,60))
plotGRanges(resize(gir,1), xlim=c(0,60), col="green")
plotGRanges(flank(gir,3), xlim=c(0,60), col="purple")
plotGRanges(flank(gir,2,start=F), xlim=c(0,60), col="brown")


#Applications to visualization of methylation array data

library(ArrayExpress)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
if(!file.exists("E-MTAB-5797.sdrf.txt")) nano = getAE("E-MTAB-5797")
library(minfi)
pref = unique(substr(dir(patt="idat"),1,17))
raw = read.metharray(pref)
glioMeth = preprocessQuantile(raw) # generate SummarizedExperiment

MbyGene = function(mset, symbol="TP53", rad=5000) {
  # phase 1: annotated GRanges for the gene
  require(erma)
  require(Gviz)
  gmod = suppressMessages(genemodel(symbol))     # erma utility
  gseq = as.character(seqnames(gmod)[1])
  gmod$transcript = symbol
  # phase 2: filter down to the region of interest
  mlim = mset[which(seqnames(mset)==gseq),] # restrict to chromosome
  # now focus the methylation data to vicinity of gene
  d1 = subsetByOverlaps(GRanges(rowRanges(mlim),,, getM(mlim)), 
                        range(gmod)+rad)
  # phase 3: use Gviz
  plotTracks(list(DataTrack(d1), 
                  GeneRegionTrack(gmod, 
                                  transcriptAnnotation="transcript", name=gseq), 
                  GenomeAxisTrack(name=gseq, showTitle=TRUE)))
}


pdf("MethArrayData.pdf", 7, 5)
MbyGene(glioMeth, symbol = "TERT")
dev.off()