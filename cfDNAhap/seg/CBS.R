###Seshan VE, Olshen A (2018). DNAcopy: DNA copy number data analysis. R package version 1.54.0.
###Modified from DNAcopy CBS segmentation

script.dir <- dirname(sys.frame(1)$ofile)
CBStoLoad <- file.path(script.dir, "src", "CBS")
CBSdatatoLoad <- file.path(script.dir, "data", "default.DNAcopy.bdry.R")
source(CBSdatatoLoad)

if (!is.loaded("CBS")) {
    dyn.load(CBStoLoad)
} else {
    warning("Fail to load CBS library")
}

cfFAR <- function(infh, data.type="FAR", sampleid=NULL, presorted=FALSE) {
    fh <- read.delim(infh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    dt <- fh[,c(1,2,11,12,13)]
    colnames(dt) <- c("chrom", "pos", "far", "siteinfer", "type")
    if (!is.character(dt$chrom)) {
      	dt$chrom <- as.character(dt$chrom)
    }
    if (!is.numeric(dt$pos)) {
	dt$pos <- as.numeric(dt$pos)
    }
    if (!is.numeric(dt$far)) {
	dt$far <- as.numeric(dt$far)
    }
    rownums <- dim(dt)[1]
    colnums <- dim(dt)[2]
    ina <- (!is.na(dt$chrom) & is.finite(dt$pos))
    if (sum(!ina)>0)
      warning("markers with missing chrom and/or pos removed\n")
    if (!presorted) {
      sortindex <- which(ina)[order(dt$chrom[ina], dt$pos[ina])]
    } else {
      sortindex <- which(ina)
    }
    dt <- dt[sortindex,]

# check for duplicate pos within a chromosome)
    if (length(ii <- which(diff(dt$pos)==0)) > 0) {
        if (any(dt$chrom[ii]==dt$chrom[ii+1])) warning("dt has repeated positions\n")
    }
    return(dt)
}

inflfact <- function(trim) {
    a <- qnorm(1-trim)
    x <- seq(-a,a,length.out=10001)
    x1 <- (x[-10001] + x[-1])/2
    1/(sum(x1^2*dnorm(x1)/(1-2*trim))*(2*a/10000))
}


trimmed.variance <- function(far, trim=0.025) {
    n <- length(far)
    n.keep <- round((1-2*trim)*(n-1))
    inflfact(trim)*sum((sort(abs(diff(far)))[1:n.keep])^2 / (2*n.keep))
}

smoothFAR <- function(dt, smooth.region=10, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025) {
    chrom <- dt$chrom
    chrlist <- unique(chrom)
    FAR <- dt[,"far"]
    ina <- which(is.finite(FAR))
    trimmed.SD <- sqrt(trimmed.variance(FAR[ina], trim))
    outlier.SD <- outlier.SD.scale*trimmed.SD
    smooth.SD <- smooth.SD.scale*trimmed.SD
    k <- smooth.region
    n <- length(FAR[ina])
    cfrq <- diff(c(which(!duplicated(chrom[ina])), n+1))
    nchr <- length(cfrq) # to allow for some chrom with all missing
    smoothed.data <- .Fortran("smoothLR",
                              as.integer(n),
                              as.double(FAR[ina]),
                              as.integer(nchr),
                              as.integer(cfrq),
                              sgdat=double(n),
                              as.integer(k),
                              as.double(outlier.SD),
                              as.double(smooth.SD))$sgdat
    dt[,"far"][ina] <- smoothed.data
    return(dt)							  
}

getbdry <- function(eta, nperm, max.ones, tol= 1e-2) {
    bdry <- rep(0, max.ones*(max.ones+1)/2)
    zz <- .Fortran("getbdry",
                   as.double(eta),
                   as.integer(max.ones),
                   as.integer(nperm),
                   as.integer(max.ones*(max.ones+1)/2),
                   bdry=as.integer(bdry),
                   etastr=double(max.ones),
                   as.double(tol))
#  list("eta.star"=zz$etastr, "boundary"=zz$bdry)
    zz$bdry
}

changepoints.prune <- function(FAR, lseg, change.cutoff=0.05) {
    n <- length(FAR)
    nseg <- length(lseg)
    ncpt <- nseg-1
    zzz <- .Fortran("prune",
                    as.integer(n),
                    as.double(FAR),
                    as.integer(nseg),
                    as.integer(lseg),
                    as.double(change.cutoff),
                    double(nseg),
                    as.integer(ncpt),
                    loc=integer(ncpt),
                    integer(2*ncpt),
                    pncpt=integer(1))
    pruned.ncpt <- zzz$pncpt
    pruned.cpts <- cumsum(lseg)[zzz$loc[1:pruned.ncpt]]
    pruned.lseg <- diff(c(0,pruned.cpts,n))
    pruned.lseg
}

changepoints.sdundo <- function(FAR, lseg, trimmed.SD, change.SD=3) {
    change.SD <- trimmed.SD*change.SD
    cpt.loc <- cumsum(lseg)
    sdundo <- TRUE
    while(sdundo) {
        k <- length(cpt.loc)
        if (k>1) {
            segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
            segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]])}, FAR)
            adsegmed <- abs(diff(segmed))
            if (min(adsegmed) < change.SD) {
                i <- which(adsegmed == min(adsegmed))
                cpt.loc <- cpt.loc[-i]
            } else {
                sdundo <- FALSE
            }
        } else {
            sdundo <- FALSE
        }
    }
    lseg.sdundo <- diff(c(0,cpt.loc))
    lseg.sdundo
}

changepoints <- function(FAR, data.type="FAR", alpha=0.01, weights=NULL, sbdry, sbn, nperm=10000, p.method="hybrid", min.width=2, kmax=25, nmin=200, trimmed.SD=NULL, undo.splits="none", undo.prune=0.05, undo.SD=3, verbose=1, ngrid=100, tol=1e-6) {
    n <- length(FAR)
    if (missing(trimmed.SD)) trimmed.SD <- mad(diff(FAR))/sqrt(2)
#   start with the whole 
    seg.end <- c(0,n)
    k <- length(seg.end)
    change.loc <- NULL
    weighted <- ifelse(is.null(weights), FALSE, TRUE) 
    while (k > 1) {
        current.n <- seg.end[k]-seg.end[k-1]
        if (verbose>=3) {
	    cat(".... current segment:",seg.end[k-1]+1,"-",seg.end[k],"\n")
	}
        if (current.n >= 2*min.width) {
            current.FAR <- FAR[(seg.end[k-1]+1):seg.end[k]]
#   check whether hybrid method needs to be used
            hybrid <- FALSE
            delta <- 0
            if ((p.method=="hybrid") & (nmin < current.n)) {
                hybrid <- TRUE
                delta <- (kmax+1)/current.n
            }
#   call the changepoint routine
            if (weighted) {
#   get the weights for the current set of probes
                current.wts <- weights[(seg.end[k-1]+1):seg.end[k]]
                current.rwts <- sqrt(current.wts)
                current.cwts <- cumsum(current.wts)/sqrt(sum(current.wts))
#   if all values of current.genomdat are the same don't segment
                if (isTRUE(all.equal(diff(range(current.FAR)), 0))) {
                    zzz <- list()
                    zzz$ncpt <- 0
                } else {
#   centering the current data will save a lot of computations later
                    current.avg <- sum(current.FAR*current.wts)/sum(current.wts)
                    current.FAR <- current.FAR - current.avg
#   need total sum of squares too
                    current.tss <- sum(current.wts*(current.FAR^2))
                    zzz <- .Fortran("wfindcpt",
                                    n=as.integer(current.n),
                                    x=as.double(current.FAR),
                                    tss=as.double(current.tss),
                                    wts=as.double(current.wts),
                                    rwts=as.double(current.rwts),
                                    cwts=as.double(current.cwts),
                                    px=double(current.n),
                                    sx=double(current.n),
                                    nperm=as.integer(nperm),
                                    cpval=as.double(alpha),
                                    ncpt=integer(1),
                                    icpt=integer(2),
                                    hybrid=as.logical(hybrid),
                                    al0=as.integer(min.width),
                                    hk=as.integer(kmax),
                                    mncwt=double(kmax),
                                    delta=as.double(delta),
                                    ngrid=as.integer(ngrid),
                                    sbn=as.integer(sbn),
                                    sbdry=as.integer(sbdry),
                                    tol= as.double(tol))
                }
            } else { 
#   if all values of current.genomdat are the same don't segment
                if (isTRUE(all.equal(diff(range(current.FAR)), 0))) {
                    zzz <- list()
                    zzz$ncpt <- 0
                } else {
#   centering the current data will save a lot of computations later
                    current.avg <- mean(current.FAR)
                    current.FAR <- current.FAR - current.avg
#   need total sum of squares too
                    current.tss <- sum(current.FAR^2)
                    zzz <- .Fortran("fndcpt",
                                    n=as.integer(current.n),
                                    x=as.double(current.FAR),
                                    tss=as.double(current.tss),
                                    px=double(current.n),
                                    sx=double(current.n),
                                    nperm=as.integer(nperm),
                                    cpval=as.double(alpha),
                                    ncpt=integer(1),
                                    icpt=integer(2),
                                    ibin=as.logical(data.type=="binary"),
                                    hybrid=as.logical(hybrid),
                                    al0=as.integer(min.width),
                                    hk=as.integer(kmax),
                                    delta=as.double(delta),
                                    ngrid=as.integer(ngrid),
                                    sbn=as.integer(sbn),
                                    sbdry=as.integer(sbdry),
                                    tol= as.double(tol))
                }
            }
        } else {
            zzz <- list()
            zzz$ncpt <- 0
        }
        if(zzz$ncpt==0) change.loc <- c(change.loc,seg.end[k])
        seg.end <- switch(1+zzz$ncpt,seg.end[-k],
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt[1],seg.end[k]),
                          c(seg.end[1:(k-1)],seg.end[k-1]+zzz$icpt,seg.end[k]))
        k <- length(seg.end)
        if(verbose>=3) cat(".... segments to go:",seg.end,"\n")
    }
    seg.ends <- rev(change.loc)
    nseg <- length(seg.ends)
    lseg <- diff(c(0,seg.ends))
    if (nseg > 1) {
        if (undo.splits == "prune") {
            lseg <- changepoints.prune(FAR, lseg, undo.prune)
        }
        if (undo.splits == "sdundo") {
            lseg <- changepoints.sdundo(FAR, lseg, trimmed.SD, undo.SD)
        }
    }
    segmeans <- 0*lseg
    ll <- uu <- 0
    for (i in 1:length(lseg)) {
        uu <- uu + lseg[i]
        if (weighted) {
            segmeans[i] <- sum(FAR[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
        } else {
            segmeans[i] <- mean(FAR[(ll+1):uu])
        }
        ll <- uu
    }
    list("lseg" = lseg, "segmeans" = segmeans)
}


segment <- function(dt, sampleid=NULL, weights=NULL, alpha=0.01, nperm=10000, p.method= c("hybrid","perm"), min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=NULL, trim = 0.025, undo.splits=c("none","prune", "sdundo"), undo.prune=0.05, undo.SD=3, verbose=1) {
    call <- match.call()
    if (min.width < 2 | min.width > 5) stop("minimum segment width should be between 2 and 5")
    if (nmin < 4*kmax) stop("nmin should be >= 4*kmax")
    if (missing(sbdry)) {
        if (nperm==10000 & alpha==0.01 & eta==0.05) {
            #source(CBSdatatoLoad)
	    sbdry <- get("default.DNAcopy.bdry", envir=environment())
        } else {
            max.ones <- floor(nperm*alpha) + 1
            sbdry <- getbdry(eta, nperm, max.ones)
        }
    }
    weighted <- ifelse(missing(weights), FALSE, TRUE)
#   rudimentary error checking for weights
    if (weighted) {
        if (length(weights) != nrow(x)) stop("length of weights should be the same as the number of pos")
        if (min(weights) <= 0) stop("all weights should be positive")
    }
    sbn <- length(sbdry)
    uchrom <- unique(dt$chrom)
    p.method <- match.arg(p.method)
    undo.splits <- match.arg(undo.splits)
    segres <- list()
    allsegs <- list()
    allsegs$ID <- NULL
    allsegs$chrom <- NULL
    allsegs$loc.start <- NULL
    allsegs$loc.end <- NULL
    allsegs$num.mark <- NULL
    allsegs$seg.mean <- NULL
    segRows <- list()
    segRows$startRow <- NULL
    segRows$endRow <- NULL
	
    if (verbose>=1) cat(paste("Analyzing: sample", sampleid,"\n"))
    FAR <- dt[,"far"]
    ina <- which(is.finite(FAR))
    FAR <- FAR[ina]
    trimmed.SD <- sqrt(trimmed.variance(FAR, trim))
    chromi <- dt$chrom[ina]
    if (weighted) {
        wghts <- weights[ina]
    } else {
        wghts <- NULL
    }
    sample.lsegs <- NULL
    sample.segmeans <- NULL
    for (ic in uchrom) {
        if (verbose>=2) cat(paste("  current chromosome:", ic, "\n"))
        segci <- changepoints(FAR[chromi==ic], data.type="FAR", alpha, wghts, sbdry, sbn, nperm, p.method, min.width, kmax, nmin, trimmed.SD, undo.splits, undo.prune, undo.SD, verbose)
        sample.lsegs <- c(sample.lsegs, segci$lseg)
        sample.segmeans <- c(sample.segmeans, segci$segmeans)
    }
    sample.nseg <- length(sample.lsegs)
    sample.segs.start <- ina[cumsum(c(1,sample.lsegs[-sample.nseg]))]
    sample.segs.end <- ina[cumsum(sample.lsegs)]
    allsegs$ID <- c(allsegs$ID, rep(sampleid,sample.nseg))
    allsegs$chrom <- c(allsegs$chrom, dt$chrom[sample.segs.end])
    allsegs$loc.start <- c(allsegs$loc.start, dt$pos[sample.segs.start])
    allsegs$loc.end <- c(allsegs$loc.end, dt$pos[sample.segs.end])
    allsegs$num.mark <- c(allsegs$num.mark, sample.lsegs)
    allsegs$seg.mean <- c(allsegs$seg.mean, sample.segmeans)
    segRows$startRow <- c(segRows$startRow, sample.segs.start)
    segRows$endRow <- c(segRows$endRow, sample.segs.end)
	
    allsegs$seg.mean <- round(allsegs$seg.mean, 5)
    allsegs <- as.data.frame(allsegs)
    allsegs$ID <- as.character(allsegs$ID)
    segres$output <- allsegs
    segres$segRows <- as.data.frame(segRows)
    segres$call <- call
    if (weighted) segres$weights <- weights
    return(segres)
}

FetalHapSeg <- function(infh, sampleid=NULL, ffest, ffsex) {
    ffest <- as.numeric(ffest)
    #threshold setting - normal threshold; trisomy threshold
    norm.T1.ffupper <- ffest*1.25
    norm.T1.fflower <- ffest*0.75
    norm.base.ffupper <- ffest*0.25
    norm.base.fflower <- ffest*(-0.25)
    norm.T2.ffupper <- -ffest*0.75
    norm.T2.fflower <- -ffest*1.25

    dt <- cfFAR(infh, data.type="FAR", sampleid, presorted=FALSE)
    #check whether data available in case when phased data is only for one parent
    dtoutfull <- data.frame(chr=character(), pos=numeric(), siteinfer=character(), FAR=numeric(), CBSmean=numeric(), type=character(), CBShomolog=character(), stringsAsFactors=FALSE)
    dtoutsum <- data.frame(ID=character(), chrom=character(), loc.start=numeric(), loc.end=numeric(), num.mark=numeric(), seg.mean=numeric(), stringsAsFactors=FALSE)
    if (nrow(dt)==0) {
        fo0 <- gsub(".tsv", ".CBSseg.tsv", infh)
        fo1 <- gsub(".tsv", ".CBSseg.sum.tsv", infh)
        write.table(dtoutfull, fo0, quote=FALSE, row.names=FALSE, sep="\t")
        write.table(dtoutsum, fo1, quote=FALSE, row.names=FALSE, sep="\t")

    } else {
	#if male fetus, then sex chromosome calculation needs modification to fit the ratio computation
	if (ffsex=="male") {
	    replen1 = nrow(dt[(dt$chrom=="Y" & dt$type=="P1"),])
	    if (replen1 > 0) {
		dt$type[dt$chrom=="Y" & dt$type=="P1"] = "P2"
	    }

	    #replen2 = nrow(dt[(dt$chrom=="X" & dt$type=="M2"),])
	    #if (replen2 > 0) {
	    #	dt$type[dt$chrom=="X" &dt$type=="M2"] = "M1"
	    dt[dt$chrom=="X",]$far = 2*dt[dt$chrom=="X",]$far
	    #}
	}

	dtsub1 <- subset(dt, type=="M1" | type=="P1")
	dtsub2 <- subset(dt, type=="M2" | type=="P2")
	if (unique(dtsub1$type) == "P1") {
	    dtsub1$far = 0-dtsub1$far
	}


	subsampleid1 <- paste0(sampleid, unique(dtsub1$type))
	subsampleid2 <- paste0(sampleid, unique(dtsub2$type))
	dtsub1smoothFAR <- smoothFAR(dtsub1)
	dtsub2smoothFAR <- smoothFAR(dtsub2)
	dtsub1seg <- segment(dtsub1smoothFAR, subsampleid1)
	dtsub2seg <- segment(dtsub2smoothFAR, subsampleid2)
	dtsub1segres.raw <- dtsub1seg$output
	dtsub2segres.raw <- dtsub2seg$output

	dtsub1segres <- data.frame(ID=character(), chrom=character(), loc.start=numeric(), loc.end=numeric(), num.mark=numeric(), seg.mean=numeric(), stringsAsFactors=FALSE)
	dtsub2segres <- data.frame(ID=character(), chrom=character(), loc.start=numeric(), loc.end=numeric(), num.mark=numeric(), seg.mean=numeric(), stringsAsFactors=FALSE)

	#breakdown segmentation if the distance between two snps are too far; threshold: 3MB
	for (i in (1:nrow(dtsub1segres.raw))) {
	     dtsub1eval <- dtsub1[(dtsub1$chrom==dtsub1segres.raw$chrom[i] & dtsub1$pos>=dtsub1segres.raw$loc.start[i] & dtsub1$pos<=dtsub1segres.raw$loc.end[i]),]$pos
	     leneval <- which(diff(dtsub1eval)>3000000)
	     lastone <- length(dtsub1eval)

	     if (length(leneval) != 0) {
		 leneval <- c(leneval, lastone)
		 addinline <- dtsub1segres.raw[i,]
		 dtsub1segres.raw$loc.end[i] <- dtsub1eval[leneval[1]]
		 dtsub1segres.raw$num.mark[i] <- leneval[1]
		 dtsub1segres <- rbind(dtsub1segres, dtsub1segres.raw[i,])
		 for (j in (1:(length(leneval)-1))) {
		      addinline$loc.start <- dtsub1eval[leneval[j]+1]
		      addinline$loc.end <- dtsub1eval[leneval[j+1]]
		      addinline$num.mark <- leneval[j+1] - leneval[j]
		      dtsub1segres <- rbind(dtsub1segres, addinline)
		 }
	     } else {
		 dtsub1segres <- rbind(dtsub1segres, dtsub1segres.raw[i,])
	     }
	}

	for (i in (1:nrow(dtsub2segres.raw))) {
	     dtsub2eval <- dtsub2[(dtsub2$chrom==dtsub2segres.raw$chrom[i] & dtsub2$pos>=dtsub2segres.raw$loc.start[i] & dtsub2$pos<=dtsub2segres.raw$loc.end[i]),]$pos
	     leneval <- which(diff(dtsub2eval)>3000000)
	     lastone <- length(dtsub2eval)

	     if (length(leneval) != 0) {
		 leneval <- c(leneval, lastone)
		 addinline <- dtsub2segres.raw[i,]
		 dtsub2segres.raw$loc.end[i] <- dtsub2eval[leneval[1]]
		 dtsub2segres.raw$num.mark[i] <- leneval[1]
		 dtsub2segres <- rbind(dtsub2segres, dtsub2segres.raw[i,])
		 for (j in (1:(length(leneval)-1))) {
		      addinline$loc.start <- dtsub2eval[leneval[j]+1]
		      addinline$loc.end <- dtsub2eval[leneval[j+1]]
		      addinline$num.mark <- leneval[j+1] - leneval[j]
		      dtsub2segres <- rbind(dtsub2segres, addinline)
		 }
	     } else {
		 dtsub2segres <- rbind(dtsub2segres, dtsub2segres.raw[i,])
	     }
	}

	for (i in (1:nrow(dtsub1segres))) {
	    segval = dtsub1segres$seg.mean[i]
            checktype = "PM"
	    checktype = substr(dtsub1segres$ID[i], nchar(dtsub1segres$ID[i])-1, nchar(dtsub1segres$ID[i]))

	    if (checktype == "P1") {
		if (segval <= norm.base.ffupper && segval >= norm.base.fflower) {
		    dtsub1segres$HL[i] <- "PH1"
		} else if (segval >= norm.T2.fflower && segval <= norm.T2.ffupper) {
		    dtsub1segres$HL[i] <- "PH2"
		} else {
		    dtsub1segres$HL[i] <- "inconclusive"
		}
	    } else if (checktype == "M1") {
		#if male fetus, X chromosome evaluation separately
		if (ffsex=="male" && dtsub1segres$chrom[i]=="X") {
		    if (segval >= ffest*0.7 && segval <= ffest*1.3) {
			dtsub1segres$HL[i] <- "MH1"
		    } else if (segval >= -ffest*1.3 && segval <= -ffest*0.7) {
			dtsub1segres$HL[i] <- "MH2"
		    } else {
			dtsub1segres$HL[i] <- "inconclusive"
		    }

		} else {
		    if (segval <= norm.base.ffupper && segval >= norm.base.fflower) {
			dtsub1segres$HL[i] <- "MH2"
		    } else if (segval >= norm.T1.fflower && segval <= norm.T1.ffupper) {
			dtsub1segres$HL[i] <- "MH1"
		    } else {
			dtsub1segres$HL[i] <- "inconclusive"
		    }
		}
	    }
	}


	for (i in (1:nrow(dtsub2segres))) {
	    segval = dtsub2segres$seg.mean[i]
            checktype = "PM"
	    checktype = substr(dtsub2segres$ID[i], nchar(dtsub2segres$ID[i])-1, nchar(dtsub2segres$ID[i]))
           
	    if (checktype == "P2") {
		if (segval <= norm.base.ffupper && segval >= norm.base.fflower) {
		    dtsub2segres$HL[i] <- "PH2"
		} else if (segval >= norm.T1.fflower && segval <= norm.T1.ffupper) {
		    dtsub2segres$HL[i] <- "PH1"
		} else {
		    dtsub2segres$HL[i] <- "inconclusive"
		}
	    } else if (checktype == "M2") {
	        #male chrX
                if (ffsex=="male" && dtsub2segres$chrom[i]=="X") {
                    if (segval >= ffest*0.7 && segval <= ffest*1.3) {
                        dtsub2segres$HL[i] <- "MH1"
                    } else if (segval >= -ffest*1.3 && segval <= -ffest*0.7) {
                        dtsub2segres$HL[i] <- "MH2"
                    } else {
                        dtsub2segres$HL[i] <- "inconclusive"
                    }

                } else {
	    
		    if (segval <= norm.base.ffupper && segval >= norm.base.fflower) {
		        dtsub2segres$HL[i] <- "MH1"
		    } else if (segval >= norm.T2.fflower && segval <= norm.T2.ffupper) {
		        dtsub2segres$HL[i] <- "MH2"
		    } else {
		        dtsub2segres$HL[i] <- "inconclusive"
		    }
	        }
	    }
	}

	dtsub1$segfull <- rep(dtsub1segres$seg.mean, dtsub1segres$num.mark)
	dtsub1$HLfull <- rep(dtsub1segres$HL, dtsub1segres$num.mark)
	dtsub2$segfull <- rep(dtsub2segres$seg.mean, dtsub2segres$num.mark)
	dtsub2$HLfull <- rep(dtsub2segres$HL, dtsub2segres$num.mark)

	dtoutfull <- rbind(dtsub1, dtsub2)
	dtoutfull <- dtoutfull[,c("chrom","pos","siteinfer","far","segfull", "type", "HLfull")]
	colnames(dtoutfull) <- c("chr", "pos", "siteinfer", "FAR", "CBSmean", "type", "CBShomolog")

	dtoutsum <- rbind(dtsub1segres, dtsub2segres)
	sortindex <- order(dtoutsum$chrom, dtoutsum$ID, dtoutsum$loc.start)
	dtoutsum <- dtoutsum[sortindex,]

	fo0 <- gsub(".tsv", ".CBSseg.tsv", infh)
	fo1 <- gsub(".tsv", ".CBSseg.sum.tsv", infh)
	write.table(dtoutfull, fo0, quote=FALSE, row.names=FALSE, sep="\t")
	write.table(dtoutsum, fo1, quote=FALSE, row.names=FALSE, sep="\t")
    }
}
