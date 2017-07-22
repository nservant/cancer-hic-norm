##
## Keep Zero or not !?
##

plotHiCHeatmap <- function(x, log=FALSE, rm.diag=FALSE, forcePairwise=TRUE, title="", col=brewer.pal(9, "Blues")){
  require(ggplot2)
  require(reshape2)
    
  if (forcePairwise & !isPairwise(x))
    x <- forcePairwise(x)

  ##contact.sum <- mclapply(x, function(xx){z=intdata(xx);mean(z[which(!is.na(z))])})
  contact.sum <- mclapply(x, function(xx){
      im <- as.matrix(intdata(xx))
      ridx <- which(rowSums(im)>0)
      cidx <- which(colSums(im)>0)
      mean(im[ridx, cidx], na.rm=TRUE)
  })
  chrs <- sortSeqlevels(seqlevels(x))
  pchrs <- HiTC:::pair.chrom(chrs)
  
  contact.mat <- matrix(NA, nrow=length(chrs), ncol=length(chrs))
  colnames(contact.mat) <- rownames(contact.mat) <- chrs
  
  for (i in 1:length(contact.sum)){
    n <- names(contact.sum)[i]
    contact.mat[pchrs[[n]][1], pchrs[[n]][2]] <- contact.sum[[n]]
  }
  
  if (log)
    contact.mat <- log2(contact.mat)
 
  contact.plot <- t(contact.mat)
  if (rm.diag)
    diag(contact.plot) <- NA
  
  cm <- melt(contact.plot)
  cm$Var1<-factor(cm$Var1, levels=chrs)
  cm$Var2<-factor(cm$Var2, levels=rev(chrs))
  ggplot(cm, aes(Var1, Var2)) + geom_tile(aes(fill=value), colour="white")+
    #scale_fill_gradient(low=col[1], high=col[9], name="")+
    scale_fill_gradientn(colours=col)+
      labs(x="", y="")+ theme(axis.text=element_text(size=11), axis.title=element_text(face="bold"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
          ggtitle(title)
}



plotHiCgenome <- function(x, signal=NULL, ctrl=NULL, sf=1, ploidy=2, cn=NULL, ylim=c(5,20), genomePack="BSgenome.Hsapiens.UCSC.hg19", log=TRUE, ...){

  stopifnot(require(genomePack, character.only=TRUE))
  genome <- eval(as.name(genomePack))
  chrs <- sortSeqlevels(setdiff(seqlevels(x), "chrY"))
  pchr <- t(as.data.frame(HiTC:::pair.chrom(seqlevels(x))))
  
  chrlen <- seqlengths(genome)[chrs]
  chrstart <- c(1, cumsum(as.numeric(chrlen))[-length(chrlen)]+1)
  names(chrstart) <- names(chrlen)
  
  xlim <- sum(as.numeric(chrlen), na.rm=TRUE)##sum(sapply(x[isIntraChrom(x)], function(xx){dim(intdata(xx))[2]}))
  plot(x=NA, xlim=c(1,xlim), ylim=ylim, type="n", ylab="Contacts (log2)", xlab="Genome", frame=FALSE, font.lab=2, cex.lab=.7, cex.main=.8, xaxt="n", ...)

  if (!isPairwise(x) & isComplete(x))
      x <- forcePairwise(x)
 
  if (!is.null(cn)){
    require(RColorBrewer)
    col.g<-brewer.pal(6,"Reds")
    col.l<-brewer.pal(6,"Greens")
    col <- c(col.l,"#FFFF33", col.g)
    names(col) <-  c(as.character((ploidy-1):(ploidy-6)), as.character(ploidy), as.character((ploidy+1):(ploidy+6)))
  }
  
  xgr <- getCombinedIntervals(x, merge=TRUE)
  xgr <- split(xgr, seqnames(xgr))
  
  for (chr in chrs){
      xi <- xgr[[chr]]      
      xpos <- chrstart[chr] + start(xi)+round(width(xi)/2)
      if (is.null(signal)){
          allpairs <- intersect(names(x), rownames(pchr)[which(pchr[,1]==chr)])
          bigMat <- do.call(cBind,lapply(x[allpairs], function(xx){as.matrix(intdata(xx))}))
          ## bigMat <- do.call(cBind,lapply(x[allpairs], intdata))
          data <- rowSums(bigMat, na.rm=TRUE)
      }else{
         data <- elementMetadata(xi)[,signal]
      }
      if (!is.null(ctrl)){
          if (is.null(signal)){
              bigMat.ctrl <- do.call(cBind,lapply(ctrl[allpairs], intdata))
              data.ctrl <- rowSums(bigMat.ctrl, na.rm=TRUE)
              ## norm
              ##sf <- sum(data.ctrl, na.rm=TRUE)/sum(data, na.rm=TRUE)
              ##data <- data*sf
              data.ctrl <- data.ctrl*sf
          }else{
              data.ctrl <- elementMetadata(y_intervals(ctrl[[paste0(chr,chr)]]))[,signal]
          }

          ## log
          if (log){
              data <- log2(1+data) - log2(1+data.ctrl)
              data[which(data==0)] <- NA
          }
          else
              data <- data/data.ctrl
      }else{
          if(log)
              data <- log2(1+data)
      }
        
      if (!is.null(cn)){
          col.cn <- col[as.character(cn[[paste0(chr,chr)]])]
          points(x=xpos, y=data, cex=.5, pch=16, col=col.cn)
      }else{
          points(x=xpos, y=data, cex=.5, pch=16)
      }
      abline(v=max(xpos), lty=2)
      
      if(!is.null(xi$smt)){
          yseg <- xi$smt
          points(xpos, yseg, type="l", lwd=1.2, col="red", cex=.7)
      }
  }

  
}

##**********************************************************************************************************##
##
## Bias to normalize
##
##**********************************************************************************************************##


plotHiCBias <- function(x, l.out=20, use=c("cis", "trans", "all"), lim=NA, norm=FALSE, type="GC", colors= rev(brewer.pal(4, "RdYlGn")), ...){
  stopifnot()
  use <- match.arg(use)
  
  if (use=="cis"){
    x <- x[isIntraChrom(x)]
    if (norm)
      x <- normPerExpected(x)
    #x <- HTClist(lapply(x, function(xx){diag(intdata(xx)) <- NA_real_ ; xx}))
    main <- "Hi-C Intrachromosomal Contacts\n"
  }
  else if(use == "trans"){
    x <- x[!isIntraChrom(x)]
    main <- "Hi-C Interchromosomal Contacts\n"
  }else if(use == "all"){
    main <- "Hi-C Intra/Interchromosomal Contacts\n"
  }
  else {
    stop("cis and trans cannot be false")
  }

  if (length(lim)==1 || is.na(lim)){
    rgc <- sapply(x, function(xx){
      z <- elementMetadata(x_intervals(xx))[,type]
      range(z[which(z>0)], na.rm=TRUE)})
    index <- seq(from=round(min(rgc),1), to=round(max(rgc),1), length.out=l.out)
  }else{
    index <- seq(from=min(lim), to=max(lim), length.out=l.out)
  }
  
  allres <- mclapply(x, function(xx, index){
    xgi <- x_intervals(xx)
    elt.xgi <- elementMetadata(xgi)[,type]
    ygi <- y_intervals(xx)
    elt.ygi <- elementMetadata(ygi)[,type]
    idata <- intdata(xx)
    
    res <- matrix(0, ncol=length(index), nrow=length(index))
    for (i in 1:length(index)-1){
      i.st <- index[i]
      i.end <- index[i+1]
      iy <- which(elt.ygi > i.st & elt.ygi <= i.end)
      for (j in 1:length(index)-1){
        j.st <- index[j]
        j.end <- index[j+1]
        ix <- which(elt.xgi > j.st & elt.xgi <= j.end)
        if (length(ix)>0 && length(iy)>0)
          res[i,j] <- mean(as.matrix(idata)[id(ygi)[iy], id(xgi)[ix]], na.rm=TRUE)
      }
    }
    res
  }, index)
  
  res <- Reduce("+", allres)/length(allres)
  res <- log2(1+res)
  
  #par(cex.main=.8, mar=c(4,4,2.5,2.5), cex.axis=.7)
  #mypal <- colorRampPalette(colors)
  #filled.contour(res, nlevels=15, col = mypal(22), main=main, frame.plot=FALSE,
  #               plot.axes={
  #                 axis(1,cex.axis=.7, labels=round(index, 3), at=seq(0,1,length.out=l.out), las=2)
  #                 axis(2,cex.axis=.7, labels=round(index, 3), at=seq(0,1,length.out=l.out))
  #               },...)
  dat <- melt(as.matrix(res))
  p <- ggplot(dat, aes(x=Var1, y=Var2)) + geom_tile(aes(fill=as.numeric(value)), colour="white") +
  #  scale_fill_gradient(low="gray", high="red", limits=c(0,1)) +
    scale_fill_gradient2("log2(1+O/E)", low="blue", mid="white", high="red", midpoint=0.1)+
      xlab(type) + ylab(type) + scale_x_continuous(breaks = 1:length(index), labels = as.character(round(index, 3))) +
        scale_y_continuous(breaks = 1:length(index), labels = as.character(round(index, 3)))+
          theme(axis.text.x=element_text(angle=90, hjust=1))
  
  return(p)
}


##**********************************************************************************************************##
##
## Copy Number estimation
##
##**********************************************************************************************************##

correctFromHiCBias <- function(x, asGRanges=TRUE, model=c("poisson","loess","nb"), ...){
  model <- match.arg(model, c("poisson", "loess", "nb"))

  if (inherits(x, "HTClist")){
    idata <- HiTC:::getCombinedContacts(x)
    counts <- rowSums(idata, na.rm=TRUE)
    gr <-  HiTC:::getCombinedIntervals(x, merge=TRUE)
  }else if (inherits(x, "HTCexp")){
    idata <- intdata(x)
    counts <-  rowSums(idata, na.rm=TRUE)
    gr <- x_intervals(x)
  }
  gr$Smt <- gr$Ratio <- gr$CnGAP <- NULL
  gr$rcounts <- counts 
  map <- gr$map
  ## scale covariable
  gc <- scale(gr$GC)
  len <- scale(gr$len)

  
  d.counts <- data.frame(x=counts, gc=gc, map=map, len=len)
 
  ## all interactions
  if (model == "loess"){
    ## try with log to constraint the range of data (negative values). But will it change the dynamic of the signal ?
    l <- loess(x ~ gc + len + map, data=log2(1+d.counts)-log2(median(counts, na.rm=TRUE)))
    corrected <-exp(l$y - l$fitted)
    
    l <- loess(x ~ gc + len + map, data=d.counts, ...)
    corrected <-exp(l$y - l$fitted)
  }else if (model=="poisson"){
    l <- glm(x ~ gc + len + offset(map), data=d.counts, family="poisson")
    ## PG - sustraction -> center on zero
    #corrected <-l$y - l$fitted
    ## Try to add mu to rescale the data (mu=84 :( )
    #corrected <- l$y-l$fitted + predict(l, data.frame(gc=0, map=0, len=0), type="resp")
    ## why not doing a ratio
    ##corrected <- l$y/l$fitted
    l.cor <- (l$y - l$fitted + median(l$y, na.rm=TRUE))/median(l$y, na.rm=TRUE)
    ## exactly as fitted(l) 
    ## corrected <- round(l$y/exp(coeff[1] + coeff[2] * l$model$gc + 
    ##                           coeff[3] * l$model$len +  coeff[4] * l$model$map), 4)
  }else if (model=="nb"){
    require(MASS)
    l <-glm.nb(x ~ len + gc + offset(map), data=d.counts, ...)
  }

  ##corrected[l$fitted <= 0] <- 0
  x.cor <- rep(NA_real_, length(counts))
  x.cor[setdiff(1:length(counts), l$na.action)] <- l.cor
  
  if(asGRanges){
    gr$counts.cor <- round(x.cor,3)
    return(gr)
  }else{
    return(round(corrected,3))
  }
}


estimateCNV <- function(x, ctrl, genomePack="BSgenome.Hsapiens.UCSC.hg19"){

    stopifnot(require(genomePack, character.only=TRUE))
    genome <- eval(as.name(genomePack))
    
    if (!is.null(ctrl)){
        com <- intersect(names(x), names(ctrl))
        x <- x[com]
        ctrl <- ctrl[com]
    }
    
    chrs <- seqlevels(x)
    chrlen <- seqlengths(genome)[chrs]
    chrstart <- c(1, cumsum(as.numeric(chrlen))[-length(chrlen)]+1)
    names(chrstart) <- names(chrlen)
    
    ## Get rowSums corrected from Hi-C bias
    rs.cor <- correctFromHiCBias(x)
    ctl.cor <- correctFromHiCBias(ctrl)

    #sx <- summary(x)
    #sc <- summary(ctrl)

    gr <-  HiTC:::getCombinedIntervals(x, merge=TRUE)
    
    xpos <- start(gr) + as.vector(Rle(values=chrstart, lengths=runLength(seqnames(gr))))
    ratio <- rs.cor$cor/ctl.cor$cor
    par(mfrow=c(3,1))
    plot(x=xpos, y=rs.cor$cor*1/6.69, pch=16, frame=FALSE, xlab="Genome", ylim=c(0,2e4))
    plot(x=xpos, y=ctl.cor$cor, pch=16, frame=FALSE, xlab="Genome", ylim=c(0,2e4))
    plot(x=xpos, y=ratio, pch=16, frame=FALSE, xlab="Genome", ylab="log(x/ctrl)", ylim=c(0,20))

    seg <- runSeg.cghp(ratio, xpos)
}

##Segmentation genome
runSeg <- function(x, rm.outliers=TRUE, call.nlevels=3, genomePack="BSgenome.Hsapiens.UCSC.hg19"){
  require(cghseg)
  
  stopifnot(require(genomePack, character.only=TRUE))
  genome <- eval(as.name(genomePack))
  chrs <- seqlevels(rs.corrected)
  chrlen <- seqlengths(genome)[chrs]
  chrstart <- c(1, cumsum(as.numeric(chrlen))[-length(chrlen)]+1)
  names(chrstart) <- names(chrlen)
  chrstart.rle <- Rle(values=chrstart, lengths=runLength(seqnames(rs.corrected)))

  data <- x$fit
  gpos <- start(x) + width(x)/2 + as.vector(chrstart.rle)

  ## rm NA
  idx <- which(!is.na(data))
  data <- data[idx]
    
  CGHo = new('CGHoptions')
   nbprocs(CGHo)<-4
  calling(CGHo) <- TRUE
  nblevels(CGHo) <- call.nlevels
  Y <- as.data.frame(data)
  cgd <- new("CGHdata", Y)
  res <- uniseg(cgd, CGHo, uniKmax=NULL)
  y <- getsegprofiles(res)

  
  x$seg <- NA
  x$seg[idx]<- y
  split(x, seqnames(x))
}

##Segmentation per chromosome
runSegPerChromosome <- function(x, rm.outliers=TRUE, calling=TRUE, out.thres=.02, call.nlevels=3, genomePack="BSgenome.Hsapiens.UCSC.hg19"){
    require(cghseg)
    #require(copynumber)

    grl <- split(x, seqnames(x))

    lseg <- GRangesList(lapply(grl, function(xx, nlevels=3, rm.outl=TRUE, out.t=.02){
      message("\nStart ",as.character(unique(seqnames(xx)))," segmentation ...")
      data <- xx$counts.cor
      gpos <- start(xx)+width(xx)/2
      
      if (rm.outliers){
        require(DNAcopy)
        obj <- CNA(data, as.vector(seqnames(xx)), gpos, data.type="logratio")
        smoothed.obj <- smooth.CNA(obj, outlier.SD.scale=3, smooth.SD.scale=2)
        ## order results
        smoothed.gr <- GRanges(seqnames=as.character(smoothed.obj$chrom), IRanges(smoothed.obj$maploc, smoothed.obj$maploc), val=smoothed.obj$Sample.1)
        seqlevels(smoothed.gr)<-seqlevels(xx)
        smoothed.gr <- sort(smoothed.gr)
        data <- smoothed.gr$val
        ##qt <- quantile(data, na.rm=TRUE, probs=c(out.thres/2, 1-(out.thres/2)))
        ##idx <- which(data>qt[1] & data<qt[2])
      }
      idx <- which(!is.na(data))
      data <- data[idx]
      
      if(length(data)>10){
        
        #xx$used=FALSE
        #xx$used[idx] <- TRUE
        
        CGHo = new('CGHoptions')
        calling(CGHo) <- FALSE
        nblevels(CGHo)<- nlevels
        Y <- data.frame(data)
        cgd <- new("CGHdata", Y)
        res <- uniseg(cgd, CGHo)
        y <- getsegprofiles(res)
        nbseg <- length(unique(diff(y))) + 1
        message("# segments = ", nbseg)
        
        if (nbseg >= nlevels & calling){
          message("Calling ...")
          calling(CGHo) <- TRUE
          res.call <- uniseg(cgd, CGHo)
          y.call <- getsegprofiles(res.call)
          nbseg.call <- length(unique(diff(y.call))) + 1
          message("# segments after calling = ", nbseg.call)
          y <- y.call
        }
        
        xx$smt <- NA
        xx$smt[idx] <- y
      }else{
        xx$smt <- NA
      }
      xx
    }, nlevels=call.nlevels, rm.outl=rm.outliers, out.t=out.thres))
    lseg
}

smoothGLAD <- function(rs.seg, bkp.stringency="normal"){
    require(GLAD)
    ## glad
    temp <- unlist(rs.seg)
    prof <- as.profileCGH(data.frame(LogRatio=temp$counts.cor, PosOrder=1:length(temp), PosBase=start(temp), Chromosome=seqnames(temp)))
    prof$profileValues$Smoothing <- estimateMissing(temp$smt)[prof$profileValues$PosOrder]
    
    ##out.glad <- daglad(prof, OnlyOptimCall=TRUE, amplicon=10, deletion=-10, deltaN=4,  forceGL=c(-5, 5), lambdabreak=8, lambdaclusterGen=40, param=c(d=6))
    if (bkp.stringency == "normal"){
        out.glad <- daglad(prof, OnlyOptimCall=TRUE,
                           amplicon=5, deletion=0, deltaN=1, forceGL=c(-.5, 2), lambdabreak=65, lambdaclusterGen=40, param=c(d=8))
    }else if (bkp.stringency == "low"){
        out.glad <- daglad(prof, OnlyOptimCall=TRUE,
                           amplicon=5, deletion=0, deltaN=1, forceGL=c(-.5, 2), lambdabreak=8, lambdaclusterGen=40, param=c(d=6))
    }else if  (bkp.stringency == "high"){
        out.glad <- daglad(prof, OnlyOptimCall=TRUE,
                           amplicon=5, deletion=0, deltaN=1, forceGL=c(-.5, 2), lambdabreak=80, lambdaclusterGen=40, param=c(d=10))
    }
    
    out.glad.list <- split(out.glad$profileValues, out.glad$profileValues$ChromosomeChar)
    lm <- 0
    unlisted <- unlist(rs.seg, use.names=FALSE)
    values(unlisted)$smt <- NA
    rs.seg<-relist(unlisted, rs.seg)
    for (i in names(rs.seg)){
        print(i)
        rs.seg[[i]]$smt <- rep(NA, length(rs.seg[[i]]))
        rs.seg[[i]]$smt[out.glad.list[[i]]$PosOrder - lm] <- as.vector(out.glad.list[[i]]$Smoothing)
        lm <- lm + length(rs.seg[[i]])
    }
    rs.seg
}

estimateMissing <- function(x){

   nav <- which(is.na(x))
   s <- c(0, which(diff(nav)!=1))+1
   e <- c(which(diff(nav)!=1), length(nav))
   ir <- IRanges(nav[s], nav[e])

   stopifnot(sum(width(ir))==length(nav))

   xcor <- x
   for (i in 1:length(ir)){
       xx <- ir[i]
       vna <- start(xx):end(xx)
       dw <- start(xx)-1
       up <- end(xx)+1

       if (up > length(x)){
           narep <- rep(dw, length(vna))
       }else if (dw < 1){
           narep <- rep(up, length(vna))
       }else{
           narep <- ifelse(abs(vna-up)>=abs(vna-dw), dw, up)
       }
       xcor[start(xx): end(xx)] <- x[narep]
   }
   xcor
}
