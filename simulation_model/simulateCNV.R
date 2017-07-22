## Nicolas Servant
## HiC-Norm
## Copyleft 2015 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.

#rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)

if (la > 0){
  for (i in 1:la){
    eval(parse(text=args[[i]]))
  }
}

message("## input=", input)
message("## xgi=", xgi)
message("## cnv=", cnv)
message("## output_prefix=", output_prefix)

require(HiTC)
require(rtracklayer)
require(Matrix)
require(ggplot2)
require(reshape2)

## Load cancer lib
source("lib_simulateCNV.R")

## Set seed for reproducible results
set.seed(seed=123456789)

## Load Hi-C data
d <- importC(input, xgi, ygi=xgi)
d <- d[grep("chrY|chrM",names(d), invert=TRUE)]
d <- forcePairwise(forceSymmetric(d))
d.sim <- d

## Load CNV data
d.cnv <- import(cnv)

## export BED file
exportCNVBED(d, d.cnv, file=paste0(output_prefix, ".bed"))

## Calculate C'ij
d.adj <- getCNVadjustedMaps(d, d.cnv)

## Downsampling

##mediane per block
#d.sim <- poisson_downsampling(d, d.adj, d.cnv, smooth.method="median")
## no smoothing - just rm outliers
d.sim <- poisson_downsampling(d, d.adj, d.cnv, smooth.method="filter.high", filter.high.th=0.95)
## take the median per distance
#d.sim <- poisson_downsampling(d, d.adj, d.cnv, smooth.method="median.per.dist")


## Export
exportC(d.sim$x.sim, output_prefix)

## Make p_ij plots for low resolution
plotpij <- function(x){
  chrs <- c("chr1", "chr2", "chr3", "chr4")
  
  dat <- getCombinedContacts(reduce(x, chr=chrs))
  datint <- getCombinedIntervals(reduce(x, chr=chrs))
  xl <- cumsum(runLength(seqnames(datint$ygi)))
  
  mdat <- melt(as.matrix(dat))
  mdat$Y=cut(mdat$value, breaks = c(0,0.5,0.99,1.01,2:6,Inf),right = FALSE)
  reds <- colorRampPalette(c("#fdae61", "#a50026"))(6)
  greens <- c("#006837", "#a6d96a")
  
  ggplot(mdat, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=Y)) +
    scale_fill_manual(breaks=c("[0,0.5)","[0.5,0.99)","[0.99,1.01)","[1.01,2)","[2,3)","[3,4)","[4,5)","[5,6)","[6,Inf)"), values=c(greens, "white", reds)) +
      scale_y_reverse() + geom_vline(xintercept=xl) + geom_hline(yintercept = xl)+
        xlab("") + ylab("") + theme(legend.title=element_blank())
}

png(paste0(output_prefix,"_pij_map.png"), width=1500, height=1400)
plotpij(d.sim$p_ij)
dev.off()

