## Nicolas Servant
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

rm(list=ls())
options(mc.cores=6)

## input = .matrix file of raw Hi-C data
## xgi = BED file of Hi-C data
## cnv = BED file of Hi-C data with CNV information
## odir = output folder
cnv <- NULL

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## Add copynumber data                                                                                                                                   
setGenomeCn <- function(x, cnv){
  stopifnot(inherits(x, "HTCexp"))
  xgi <- x_intervals(x)
  stopifnot(length(xgi)==length(cnv))
  elementMetadata(xgi)[,"cnv"]<- as.numeric(cnv$name)
  x_intervals(x) <- xgi
  y_intervals(x) <- xgi
  x
}

require(HiTC)
require(rtracklayer)
require(Matrix)
require(RColorBrewer)
require(ggplot2)
require(gridExtra)

## Load data
d <- importC(input, xgi, ygi=xgi)
d <- forceSymmetric(d)
## Expect to have half of trans + cis

##################################
##
## Annotation Hi-C data from know bias
##
##################################

## Generate MboI_hg19.RData file
##map_hg19<- import("wgEncodeCrgMapabilityAlign100mer_hg19.bigWig", format="BigWig", asRangedData=FALSE)
##cutSites <- getAnnotatedRestrictionSites(resSite="AAGCTT", overhangs5=1, chromosomes=seqlevels(hic_cis), genomePack="BSgenome.Hsapiens.UCSC.hg19", wingc=200, mappability=map_hg19, winmap=500)
##cutSites_MboI <- getAnnotatedRestrictionSites(resSite="GATC", overhangs5=0, chromosomes=paste0("chr", c(1:22, "X", "Y", "M")), genomePack="BSgenome.Hsapiens.UCSC.hg19", wingc=200, mappability=map_hg19, winmap=500)
##save(cutSites_MboI, file="MboI_hg19.RData")

load("./files/MboI_hg19.RData")
cutSites <- GRangesList(lapply(cutSites_MboI, sort))

## Set genomic features based on restriction sites
rd <- file.path(odir, gsub(".mat(rix)*", ".RData", basename(input)))
d.annot <- mclapply(d[grep("chrM", names(d), invert=TRUE)], setGenomicFeatures, cutSites=cutSites)
d.annot <- HTClist(d.annot)

if (!is.null(cnv)){
  d.cnv <- import(cnv)
  d.annot[isIntraChrom(d.annot)] <- HTClist(lapply(d.annot[isIntraChrom(d.annot)], function(x, d.cnv){
    setGenomeCn(x, d.cnv[which(as.vector(seqnames(d.cnv)==seqlevels(x)))])}, d.cnv=d.cnv))  
}

save(d.annot, file=rd)



