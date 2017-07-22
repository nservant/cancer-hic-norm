################################################################################################################################
## Nicolas Servant
##
## Hi-C Analysis
## Segment Hi-C data
##
################################################################################################################################
rm(list=ls())

require(GLAD)
require(HiTC)
require(rtracklayer)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
require(Matrix)

options(mc.cores=6)
source("lib_cnv_hic.R")

##output_prefix
##input.annot
##glad.stringency = normal/low/high

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
    for (i in 1:la)
          eval(parse(text=args[[i]]))
  }

chrs <- c("chr1", "chr2", "chr3", "chr4")

## Get chromosome size
genomePack <- "BSgenome.Hsapiens.UCSC.hg19"
stopifnot(require(genomePack, character.only=TRUE))
genome <- eval(as.name(genomePack))
chrs.all <- paste0("chr", c(1:22,"X"))##seqlevels(d)
chrlen <- seqlengths(genome)[chrs.all]
chrstart <- c(1, cumsum(as.numeric(chrlen))[-length(chrlen)]+1)
names(chrstart) <- names(chrlen)

######
##
## Get copy number from Hi-C data
##
######

load(file=input.annot)
gr <-  HiTC:::getCombinedIntervals(d.annot, merge=TRUE)
xpos <- start(gr) + as.vector(Rle(values=chrstart, lengths=runLength(seqnames(gr))))

## Correct for GC/map/len
rs.corrected <- correctFromHiCBias(d.annot, model="poisson")

## Plots
dat <- data.frame(chr=as.vector(seqnames(rs.corrected)), pos = xpos, counts = rs.corrected$rcounts, counts.cor = rs.corrected$counts.cor)

p <- ggplot(dat) + geom_point(aes(x=pos, y=counts), col="gray25", size=3) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  ylab("") + geom_vline(xintercept = chrstart, col="gray80", size=.4,linetype="twodash") + ylim(c(0, 150000)) ## 0.8e5 MCF7 = 150e3 Rao = 250e3

## segmentation with cghseg
rs.seg <- runSegPerChromosome(rs.corrected, rm.outliers=TRUE, calling=FALSE)
rs.seg <- smoothGLAD(rs.seg, bkp.stringency = glad.stringency)

## plots
rs.seg.gr <- unlist(rs.seg)
xpos <- start(rs.seg.gr) + as.vector(Rle(values=chrstart, lengths=runLength(seqnames(rs.seg.gr))))
dat <- data.frame(chr=as.vector(seqnames(rs.seg.gr)), pos = xpos, counts.cor = rs.seg.gr$counts.cor, smt = rs.seg.gr$smt, cn=rs.seg.gr$cnv)

p2 <- ggplot(dat) + geom_point(aes(x=pos, y=counts.cor), col="gray25", size=3) +
  geom_line(aes(x=pos, y=smt), col="red", size=2) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  ylab("") + geom_vline(xintercept = chrstart, col="gray80", size=.4,linetype="twodash") + geom_hline(yintercept=1, size=.4,linetype="twodash") + ylim(c(0, 5))

pdf(paste0(output_prefix, "_correctFromHiCBias_segmented.pdf"), width=30, height=10)
grid.arrange(p, p2, ncol=1, nrow=2)
dev.off()

## export results
segperchr<-sapply(rs.seg, function(xx){xx$smt})
obs_smt <- lapply(segperchr, estimateMissing)
w.seg <- lapply(obs_smt, function(x){x[which(x==0)]=10e-4; round(x, 2)})
out <- data.frame(chr = as.vector(seqnames(rs.seg.gr)), start = start(rs.seg.gr), end = end(rs.seg.gr), score = unlist(w.seg))

write.table(out, file=paste0(output_prefix, "_seg.bed"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

