## Nicolas Servant
## Copyleft 2016 Institut Curie
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.


##
## Make CNV plots
##
plotCountsCNV <- function(x, x.cn){
  x.intra <- x[isIntraChrom(x)]
  
  ## normPerExpected
  res <- lapply(x.intra, function(xx){
    xx <- normPerExpected(xx)
    chr <- seqlevels(xx)
    counts <- intdata(xx)
    
    cn <- rep(2, dim(intdata(xx))[1])
    ov <- findOverlaps(x_intervals(xx), x.cn)
    cn[queryHits(ov)] <- x.cn$name[subjectHits(ov)]
    cn.mat <- matrix(cn, ncol=1) %*% matrix(cn, nrow=1)
    a <- split(counts, cn.mat)
    sapply(a, mean, na.rm=TRUE)
  })
  names(res) <- NULL
  res<-unlist(res)
  
  mp <- data.frame(Cn=log2(as.numeric(names(res))), Contacts=log2(as.vector(res)))
  p <- ggplot(na.omit(mp), aes(x=Cn, y=Contacts)) + geom_point(colour="black" , size=3) + ylab("log2(O/E)") + xlab("Copy Number\n(log2CiCj)")
  p
}

make1Dplot <- function(x, x.cnv=NULL, chrs=NULL, ylab="", xlab="", use.cis=TRUE, use.trans=TRUE, yl=NULL, ploidy=2){

  stopifnot(use.cis || use.trans)
  
  if (!use.trans){
    ## cis
    x <- x[isIntraChrom(x)]
  }
  if (!use.cis){
    ## trans
    x <- x[!isIntraChrom(x)]
  }
  
  counts <- log2(rowSums(getCombinedContacts(x), na.rm=TRUE))
  counts.mean <- mean(counts)
  
  ## genomic position
  pos <- getCombinedIntervals(x)$xgi
  pos <- pos[order(as.numeric(pos$name))]
  sl <-as.character(runValue(seqnames(pos)))

  chr.len <- sapply(sl, function(x){end(pos[max(which(as.vector(seqnames(pos))==x))])})   
  chr.cumsum <- c(0, cumsum(as.numeric(chr.len))[1:length(chr.len)-1])
  names(chr.cumsum) <- names(chr.len)
  xpos <- start(pos) + rep.int(chr.cumsum, runLength(seqnames(pos)))
  chrvec <- rep(names(chr.cumsum),  runLength(seqnames(pos)))
  
  dat <- data.frame(chr=chrvec, pos = xpos, counts=counts)

  if (!is.null(chrs)){
    dat <- dat[which(dat$chr %in% chrs),]
    chr.cumsum <- chr.cumsum[chrs]
  }
  p <- ggplot(dat) + geom_point(aes(x=pos, y=counts), col="gray25", size=3) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
      ylab(ylab) + geom_vline(xintercept = chr.cumsum, col="gray80", size=.4,linetype="twodash")

  if (!is.null(yl)){
    p <- p + ylim(yl)
  }
  
  ## Add CNV regions
  if (!is.null(x.cnv)){
    if (!is.null(chrs))
      x.cnv <- x.cnv[IRanges::which(as.vector(seqnames(x.cnv)%in%chrs))]
    if (ploidy==4){
        sel.col <- c("#238B45", "#A1D99B", "#A1D99B", "#FFFFFF","#FFFFFF","#FB6A4A","#EF3B2C", "#CB181D","#CB181D","#A50F15","#67000D")
        names(sel.col) <- as.character(0:10)
    }else{
        sel.col <- c("#238B45", "#A1D99B","#FFFFFF","#FB6A4A","#EF3B2C", "#CB181D","#CB181D","#A50F15","#67000D","#CB181D","#CB181D")
        names(sel.col) <- as.character(0:10)
    }

    
    if (is.null(yl)){
      yr <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
    }else{
      yr <- yl
    }
    #rs <- start(x.cnv) + rep.int(chr.cumsum[levels(seqnames(x.cnv))], runLength(seqnames(x.cnv)))
    #re <- end(x.cnv) + rep.int(chr.cumsum[levels(seqnames(x.cnv))], runLength(seqnames(x.cnv)))
    rs <- start(x.cnv) + rep.int(chr.cumsum[as.character(runValue(seqnames(x.cnv)))], runLength(seqnames(x.cnv)))
    re <- end(x.cnv) + rep.int(chr.cumsum[as.character(runValue(seqnames(x.cnv)))], runLength(seqnames(x.cnv)))
 
    rects <- data.frame(start=rs, end=re, ystart=yr[1], yend=yr[2], group=as.character(x.cnv$name))

    p <-  p + geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=ystart, ymax=yend, fill=as.character(group)),  colour="transparent", alpha=.3) +
      scale_fill_manual(values=sel.col, labels=names(sel.col)) +  theme(legend.position="none")
  }
  p
}


## exportCNVBED
## Export BED file required by ice for the normalization
## TODO : deal with sexual chromosomes !
exportCNVBED <- function(d, d.cnv, file){
    z <- unlist(GRangesList(getCombinedIntervals(d)$ygi))
    z$cnv <- 2
    ov <- findOverlaps(z, d.cnv)
    z[queryHits(ov)]$cnv = d.cnv[subjectHits(ov)]$name
    ## BED start coordinate
    start(z) <-  start(z) -1
    write.table(as.data.frame(z)[,c("seqnames", "start", "end", "cnv")], file = file, quote=FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}
    
## estimate_transH
## Calculate Trans Homologous contact frequency
##
## x : HTClist ; Hi-C data
## chrom : character ; if specified the trans contact is calculated for a given chromosome. Otherwise, it will be calculated genome-wide
## use.zero : logical ; if true, all values including zeros are used to compute the mean of trans contact

estimate_transH <- function(x, chrom=NULL, use.zero=TRUE, N=2){

  alltrans <- x[which(!isIntraChrom(x))]
  
  if (!is.null(chrom)){
    ntrans <- t(data.frame(HiTC:::pair.chrom(seqlevels(x), rm.cis=TRUE)))
    ntrans <- names(which(ntrans[,1]==chrom | ntrans[,2]==chrom))
    alltrans <- x[intersect(ntrans, names(x))]
  }

  ## Remove sexual chromosome and chrM
  alltrans <- alltrans[grep("chrM|chrX|chrY", names(alltrans), invert = TRUE)]
  if (length(alltrans)>0){
      if (use.zero){
        vtrans <- lapply(alltrans, function(xx){as.vector(intdata(xx))})
        mtrans <- median(unlist(vtrans), na.rm=TRUE)
        ##mtrans <- median(sapply(alltrans, function(x){xd <- as.matrix(intdata(x)); median(xd, na.rm=TRUE)}), na.rm=TRUE)    
      }else{
        vtrans <- lapply(alltrans, function(xx){xd=as.vector(intdata(xx)); xd[which(xd>0)]})
        mtrans <- median(unlist(vtrans), na.rm=TRUE)
        ##mtrans <- median(sapply(alltrans, function(x){xd <- as.matrix(intdata(x)); median(xd[which(xd>0)], na.rm=TRUE)}), na.rm=TRUE)
      }
  }else{
      mtrans <- NA
  }

  ## But what we observed is actually Ni * Nj * transh
  ## transh is therefore equal to C_ij_trans / (Ni*Nj)
  mtrans/(N*N)
}#estimate_transH


## The goal here is to estimate cis from an observed contact frequency under the assumption that there is no cnv
## So usually N=2
estimate_cis <- function(c_ij, transh, N){

    ## If i and j belong to the same segment, the observed counts C_ij is expected
    ## to be equal to C_ij = N*cis + N(N-1)transh
    cis = (c_ij - N * (N-1) * transh) / N
    #cis[which(cis<0)] <- 0
    cis
}##estimate_cis


check_cnv_to_simulate <- function(x, res=NA){
    if (!is.na(res)){
        torm <- x[which(width(x)<res)]
        warning("Remove ", length(torm), " alteration smaller than ", res, "bp in ", paste0(unique(seqnames(torm)), sep=","))
        x <- x[which(width(x)>=res)]
    }
    if (length(which(countOverlaps(x, x)>1))>0){
        stop("some alterations are overlapping. Please remove them")
    }
    x
}

smooth <- function(x, smooth.method=c("median", "median.per.dist", "filer.high"), smooth.filter.high.th=1, xdist=NULL){
    x <- as.matrix(x)

    if (smooth.method=="filter.high" && smooth.filter.high.th<1){
      th <- quantile(x, na.rm=TRUE, probs=smooth.filter.high.th)
      x[which(x>th)] <- th
      out <- x
    }
    else if(smooth.method=="median"){
      out <- median(x, na.rm=TRUE)
      out <- matrix(out, ncol=ncol(x), nrow=nrow(x))
    }else if (smooth.method=="median.per.dist" && !is.null(xdist)){
      th <- quantile(x, na.rm=TRUE, probs=.99)
      x[which(x>th)] <- th
            
      mi <- split(as.vector(x), as.character(xdist))
      milen <- lapply(mi, length)
      mimed <- lapply(mi, median, na.rm = TRUE)
      miexp <- lapply(1:length(milen), function(i) {
        rep(mimed[[i]], milen[[i]])
      })
      names(miexp) <- names(mi)
      
      out <- matrix(unsplit(miexp, as.character(xdist)), nrow = nrow(x), 
                    ncol = ncol(x))
    }else{
      stop("Unknown method")
    }
     
    out[which(is.na(x))] <- NA
    #z[which(x==0)] <- 0
    out
}

## cis_cnv
## Simulate CNV from a cis contact map
##
## x : HTCexp ; object of cis map
## seg : GRanges ; genomic coordinate of altered segment(s) from a given chromosome
## transh : numeric ; trans homologuous contact probability
## ploidy : numeric ; ploidy of original data
##

updateCIS <- function(x, seg, transh=NA, do.method="estimate.cnv", smooth.filter.high.th=1, verbose=TRUE, upa=TRUE){

    x <- forceSymmetric(x)
    ## convert to count matrix in dense format for speed ... at the expense of the RAM
    idata <- as.matrix(intdata(x))
    if (do.method == "median.per.dist"){
        idistance <- intervalsDist(x)
    }else{
        idistance <- NULL
    }

    ## get bins ranges
    xgi <- x_intervals(x)
    if (upa){
        start(xgi) <- end(xgi) <- start(xgi) + (width(xgi)/2)
    }
        
    ## get all bins falling in all altered segments
    ov <- findOverlaps(xgi, seg)
    ploidy <- 2
    
    ## 1- update interaction counts within each segment
    for (i in 1:length(seg)){
        cn <- as.numeric(seg$name)[i]
        if (verbose)
            message("-- contact within ", as.character(seqnames(seg[i])), ":", start(seg)[i], "-", end(seg)[i], "[", cn, "]")
        ## get bins of the matrix overlapping with that segment
        idx <- queryHits(ov)[which(subjectHits(ov)==i)]
        ## sub-matrix related to the fragment
        xx <-  idata[idx, idx]
        xd <- idistance[idx, idx]

        if (do.method=="estimate.cnv"){
          ## estimate cis from the observed data
          cis <-  estimate_cis(c_ij=xx, transh=transh, N=ploidy)          
          ## Calculate new expected interaction
          ## N_i * cis + N_i(N_j-1) * transH
          l <- cn * cis  + (cn * (cn -1)) * transh
          ## In theory the loci that have NO interactions, cannot interact more, even after simulation ?
          l[which(xx==0)] <- 0
        }else {
          l <- smooth(xx, xdist=xd, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
        }
        idata[idx, idx] <- l
      }
    
    ## 2 - update interaction between altered segments
    if (length(seg)>1){
        ## for all segments
        for (i in 1:length(seg)){
            idx.x <- queryHits(ov)[which(subjectHits(ov)==i)]
            cn.x <- as.numeric(seg$name)[i]
            
            ## and its other partners
            op <- setdiff(1:length(seg), i)
            for (j in op){
                idx.y <- queryHits(ov)[which(subjectHits(ov)==j)]
                cn.y <- as.numeric(seg$name)[j]
                ##print(idata[c(142, 23), c(142, 23)])
                xx <- idata[idx.x, idx.y]
                xd <- idistance[idx.x, idx.y]
                if (verbose)
                  message("-- contact between ", as.character(seqnames(seg[i])), ":", start(seg)[i], "-", end(seg)[i], "[", cn.x, "]",
                          " and ", as.character(seqnames(seg[j])), ":", start(seg)[j], "-", end(seg)[j], "[", cn.y, "]")
                
                if (do.method=="estimate.cnv"){
                  cis = estimate_cis(c_ij=xx, transh=transh, N=ploidy)
                  ## p * cis + (N_i * N_j - p) * transh
                  if (cn.x < ploidy || cn.y < ploidy){
                    cn <- min(c(cn.x, cn.y))
                  }else{
                    cn <- ploidy
                  }
                  l  <- cn * cis + (cn.x * cn.y - cn) * transh
                  l[which(xx==0)] <- 0 ## force zeros
                }else{
                  l <- smooth(xx, xdist=xd, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
                }
                
                idata[idx.x, idx.y] <- l
            }
        }
    }
    
    ## 3- update interaction counts outside all segments
    for (i in 1:length(seg)){
        idx.allseg <- queryHits(ov)
        idx.curseg <- queryHits(ov)[which(subjectHits(ov)==i)]

        ## get all bins that do not overlap any altered segment
        ## and which therefore have a cn = ploidy
        idx.rev <- setdiff(1:length(xgi), idx.allseg)
        cn.x <- as.numeric(seg$name)[i]
        xx <- idata[idx.curseg, idx.rev]
        xd <- idistance[idx.curseg, idx.rev]
        if (verbose)
            message("-- contact outside ", as.character(seqnames(seg[i])), ":", start(seg)[i], "-", end(seg)[i], "[", cn.x, "]")

        if (do.method=="estimate.cnv"){
            cis = estimate_cis(c_ij=xx, transh=transh, N=ploidy)
            ## p * cis + (N_i * N_j - p) * transh
            if (cn.x < ploidy){
                cn <- min(c(cn.x, ploidy))
            }else{
                cn <- ploidy
            }
            l <- cn * cis + (cn * cn.x - cn) * transh
            l[which(xx==0)] <- 0 ## force zeros
        }else{
            l <- smooth(xx, xdist=xd, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
        }
        idata[idx.curseg, idx.rev] <- l
        idata[idx.rev, idx.curseg] <- t(l)
    }
    ## Return final object
    out.mat <- as(Matrix(idata), "sparseMatrix")
    rownames(out.mat) <- rownames(idata)
    colnames(out.mat) <- colnames(idata)
    intdata(x) <- out.mat
    x
}#updateCIS


## updateTRANS
## Simulate CNV from trans contact maps
##
## x : HTClist ; object of trans maps
## seg : GRanges ; genomic coordinate of altered segments
## transh : numeric ; trans homologuous contact
## ploidy : numeric ; ploidy of original data

updateTRANS <- function(x, allseg, transh=NA, do.method="estimate.cnv", ploidy=2, verbose=TRUE, upa=TRUE){
    
    ## get only the upper part of trans maps
    if (isComplete(x))
        x <- forceSymmetric(x)

    ## restrict cnv to chromosome in the dataset
    allseg <- allseg[IRanges::which(seqnames(allseg) %in% seqlevels(x))]

    
    ## For all chromosomes
    for (chr in seqlevels(x)){
        if (verbose)
            message("update trans interactions from chromosome ", chr, "...")
        
        ## Altered segments from the current chromosome
        seg <- allseg[which(as.vector(seqnames(allseg))==chr)]
        ## Altered segments from the other chromosomes
        oseg <- allseg[which(as.vector(seqnames(allseg))!=chr)]
    
        ## for all trans map
        x <- HTClist(lapply(x, function(y){
          
            ## get matrix count
            idata <- as.matrix(intdata(y))

            xgi <- x_intervals(y)
            ygi <- y_intervals(y)
            
            if (upa) {
                start(xgi) <- end(xgi) <- start(xgi) + (width(xgi)/2)
                start(ygi) <- end(ygi) <- start(ygi) + (width(ygi)/2)
            }
            
            ## indices of current segments - working in rows
            if (seqlevels(ygi) != chr){
                return (y)
            }
            ## coordinates of segments from current chromosomes
            ov <- suppressWarnings(findOverlaps(ygi, seg))
            cn <- as.numeric(seg[unique(subjectHits(ov))]$name)
            idx.seg <- split(queryHits(ov), subjectHits(ov))
            idx.rev.seg <- setdiff(1:length(ygi), unlist(idx.seg, use.name=FALSE))
            
            ## Is there any other alterations in this trans map ?
            ov.oseg <- suppressWarnings(findOverlaps(xgi, oseg))

            cn.oseg <- as.numeric(oseg[unique(subjectHits(ov.oseg))]$name)
            idx.oseg <- split(queryHits(ov.oseg), subjectHits(ov.oseg))
            idx.rev.oseg <- setdiff(1:length(xgi), unlist(idx.oseg, use.name=FALSE))
            
            ## 1. First update interaction in trans between altered segments
             if (length(idx.seg)>0 & length(idx.oseg)>0){
                 ## for segment in chr
                 for (i in 1:length(idx.seg)){
                    idx.s <- idx.seg[[i]]
                    c.s <- cn[i]
                    ## for other alterations
                    for (j in 1:length(idx.oseg)){
                        idx.o <- idx.oseg[[j]]
                        c.o <- cn.oseg[j]
                        if (verbose)
                            message("-- contact on ", paste0(seqlevels(y), collapse="-")," between ", 
                                    as.character(seqnames(seg[i])), ":", start(seg)[i], "-", end(seg)[i], "[", c.s, "]",
                                    " and ", as.character(seqnames(oseg[j])), ":", start(oseg)[j], "-", end(oseg)[j], "[", c.s, "]")

                        xx <- idata[idx.s, idx.o]
                        ##print(head(xx))
                        if (do.method=="estimate.cnv"){
                            ## Ni * Nj * trans
                            ##message("c.s =", c.s, "[", paste0(idx.s, sep=","), "]", " c.o=", c.o, "[", paste0(idx.o, sep=","), "]")
                            l <- (c.s * c.o)/ploidy^2 * xx
                            l[which(xx==0)] <- 0
                        }else{
                            l <- smooth(xx, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
                        }
                        ##print(head(l))
                        idata[idx.s, idx.o] <- l
                    }
                }
            }
            
            ## 2. Then update intractions in trans outside of other alterations
            if (length(idx.oseg)>0){
                  for (i in 1:length(idx.oseg)){
                      idx.s <- idx.oseg[[i]]
                      c.s <- cn.oseg[i]
                      if (verbose)
                          message("-- contact on ", paste0(seqlevels(y), collapse="-")," outside ",
                                  as.character(seqnames(oseg[i])), ":", start(oseg)[i], "-", end(oseg)[i], "[", c.s, "]")
                      xx <- idata[idx.rev.seg, idx.s]
                      if (do.method=="estimate.cnv"){
                          ## Ni * Nj * trans
                          l <- (c.s * ploidy)/ploidy^2 * xx
                          l[which(xx==0)] <- 0
                      }else{
                          l <- smooth(xx, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
                      }
                      idata[idx.rev.seg, idx.s] <- l
                  }
              }
            
            if (length(idx.rev.oseg)>0 & length(idx.seg)>0){
                ## for segment in chr
                for (i in 1:length(idx.seg)){
                    idx.s <- idx.seg[[i]]
                    c.s <- cn[i]
                    if (verbose)
                        message("-- contact on ", paste0(seqlevels(y), collapse="-")," outside ",
                                as.character(seqnames(seg[i])), ":", start(seg)[i], "-", end(seg)[i], "[", c.s, "]")
                    xx <- idata[idx.s, idx.rev.oseg]
                    ## Ni * Nj * trans
                    if (do.method=="estimate.cnv"){
                        l <- (c.s * ploidy)/ploidy^2 * xx
                        l[which(xx==0)] <- 0
                    }else{
                        l <- smooth(xx, smooth.method=do.method, smooth.filter.high.th=smooth.filter.high.th)
                    }
                    idata[idx.s, idx.rev.oseg] <- l
                }
            }
            intdata(y) <- Matrix(idata)
            y
        }))
    }
    x
}##UpdateTRANS



getCNVadjustedMaps <- function(x, cnv, transh=NA){
    x.sim <- x

    cnv <- check_cnv_to_simulate(cnv, width(x_intervals(x[[1]]))[1])

    ## Estimate TransH from complete data
    ## if use.zero = TRUE -> transh < 1
    if (is.na(transh)){
      transh <- estimate_transH(x, use.zero=TRUE)
    }
    message("transh = ", transh)
        
    ## Intrachromosomal maps
    x.sim[isIntraChrom(x)] <- HTClist(lapply(x[isIntraChrom(x)], function(xintra){
        chr <- seqlevels(xintra)
        ## Simulate CNVs
        if (is.element(chr,  as.character(unique(seqnames(cnv))))){
            message("update cis interactions from chromosome ", chr, "...")
            xintra.sim <- updateCIS(xintra, do.method="estimate.cnv", cnv[which(as.vector(seqnames(cnv))==chr)], transh)
        }else{
            xintra.sim <- xintra
        }
        return(xintra.sim)
    }))

    ##Interchromosomal maps
    x.trans <-  x[!isIntraChrom(x)]
    if (length(x.trans)>0)
      x.sim[names(x.trans)] <- updateTRANS(x.trans, cnv, transh)
    x.sim
}


getSmoothPij <- function(x, cnv, smooth.method=c("median","median.per.dist","filter.high"), filter.high.th=1){
    x.smooth <- x

    cnv <- check_cnv_to_simulate(cnv, width(x_intervals(x[[1]]))[1])

    ## Intrachromosomal maps
    x.smooth[isIntraChrom(x)] <- HTClist(lapply(x[isIntraChrom(x)], function(xintra){
        chr <- seqlevels(xintra)
        ## Simulate CNV
        if (is.element(chr, as.character(unique(seqnames(cnv))))){
            xintra.sim <- updateCIS(xintra, cnv[which(as.vector(seqnames(cnv))==chr)], do.method=smooth.method, smooth.filter.high=filter.high.th, verbose=TRUE)
        }else{
            xintra.sim <- xintra
        }
        return(xintra.sim)
    }))

    ##Interchromosomal maps
    x.trans <-  x[!isIntraChrom(x)]
    if (length(x.trans)<0)
        x.smooth[names(x.trans)] <- updateTRANS(x.trans, cnv, do.method="median", smooth.filter.high=filter.high.th, verbose=TRUE)
    x.smooth
}

##smooth.method=c("none", "median", "median.dist")
poisson_downsampling <- function(x, x.adj, x.cnv=NULL, smooth.method=NA, filter.high.th=1){

    x.adj[isIntraChrom(x.adj)] <- lapply(x.adj[isIntraChrom(x.adj)], forceSymmetric)
    x[isIntraChrom(x)] <- lapply(x[isIntraChrom(x)], forceSymmetric)

    ## Generate cnv multiplicative factor p_ij
    p_ij <- HTClist(lapply(names(x), function(chrsname){
        a <- x.adj[[chrsname]]
        b <- x[[chrsname]]
        intdata(b) <- intdata(a) / intdata(b)
        b
    }))
    
    ## smooth pij
    if (!is.na(smooth.method)){
      p_ij <- getSmoothPij(p_ij, x.cnv, smooth.method=smooth.method, filter.high.th=filter.high.th)
    }
    ## proportion of q_ij = p_ij / max(p_ij)
    mv <- sapply(p_ij, function(x){max(intdata(x), na.rm=TRUE)})
    mx <- max(mv)
    message("max p_ij = ", mx, " [", names(mv)[which.max(mv)], "]")
    
    q_ij <- HTClist(lapply(p_ij, function(x){intdata(x) <- intdata(x) / mx; x}))

    ## back to list
    #q_ij <- lapply(q_ij, intdata)

    ## downsampling using rbinom
    cc <- lapply(names(x), function(chrsname){
        probs = intdata(q_ij[[chrsname]])
        counts = intdata(x[[chrsname]])
        l <- length(counts)
        ## counts = 0 => probs = NA will generate warnings
        suppressWarnings(matrix(rbinom(n=l, size=as.matrix(round(counts)), prob=as.matrix(probs)), ncol=ncol(counts)))
    })
    names(cc) <- names(x)
    
    ## create final object for export
    x.sim <- HTClist(lapply(names(x), function(chrsname){
        y <- x[[chrsname]] 
        z <- cc[[chrsname]]
        colnames(z) <- colnames(intdata(x[[chrsname]]))
        rownames(z) <- rownames(intdata(x[[chrsname]]))
        intdata(y) <- as(Matrix(z), "sparseMatrix")
        y
    }))

    ## forceSymmetric
    x.sim <- forceSymmetric(x.sim)
    #x.sim[isIntraChrom(x.sim)] <- lapply(x.sim[isIntraChrom(x.sim)], forceSymmetric)
    x.sim[isIntraChrom(x.sim)] <- lapply(x.sim[isIntraChrom(x.sim)], forceTriangular)

    list(x.sim=x.sim, p_ij=p_ij, q_ij=q_ij)
}

