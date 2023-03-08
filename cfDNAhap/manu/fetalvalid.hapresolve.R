args <- commandArgs(TRUE)
fetalfh <- args[1]
workdir <- args[2]
sample <- args[3]

#fetalpat <- args[1]
#fetalmat <- args[2]
#workdir <- args[3]
#sample <- args[4]


readinfetalfh <- function(fetalfh) {
    fh <- read.delim(fetalfh, header=FALSE, sep="\t")
    colnames(fh) <- c("chromosome", "pos", "rsid", "ref", "alt", "FPH1", "FPH2", "FMH1", "FMH2", "FT1", "FT2", "pathap", "mathap")
    start0 <- fh$pos-1
    ft <- cbind(fh[,1], start0, fh[,c(2,12,13)])
    colnames(ft) <- c("chromosome", "start", "end", "pathap", "mathap")
    #ft$chromosome <- paste0("chr", ft$chromosome)
    return(ft)
}

fetalpatimpute <- function(fetalfh) {
    fh <- readinfetalfh(fetalfh)
    ftpathap <- fh[,c(1,2,3,4)]
    imputepathap <- data.frame(chromosome=character(), start=numeric(), end=numeric(), pathap=character(), stringsAsFactors=FALSE)

    uchrom <- unique(ftpathap$chromosome)
    for (ic in uchrom) {
        imputeval <- "na"
	preimputeval <- "na"
        dt.chr <- subset(ftpathap, chromosome == ic)
	dt.chr$pathap <- as.character(dt.chr$pathap)
	
	for (i in 1:nrow(dt.chr)) {

            if (dt.chr$pathap[i] != "Undetermined" & dt.chr$pathap[i] != "." & dt.chr$pathap[i] != preimputeval) {
	        imputeval <- dt.chr$pathap[i]
		preimputeval <- imputeval
	    }

	    if (dt.chr$pathap[i] != "PH1" & dt.chr$pathap[i] != "PH2") {
	        dt.chr$pathap[i] = imputeval
	    }

	}
	imputepathap <- rbind(imputepathap, dt.chr)
    }
    imputepathap <- subset(imputepathap, pathap == "PH1" | pathap == "PH2")
    fo <- paste0(workdir, "/", sample, ".fetal.genome.impute.pat.hap")

    write.table(imputepathap, fo, row.names=FALSE,quote=FALSE, sep="\t")
    return(imputepathap)
}

fetalpatnoimpute <- function(fetalfh) {
    fh <- readinfetalfh(fetalfh)
    ftpathap <- fh[,c(1,2,3,4)]
    dt <- subset(ftpathap, pathap == "PH1" | pathap == "PH2")
    fo <- paste0(workdir, "/", sample, ".fetal.genome.noimpute.pat.hap")

    write.table(dt, fo, row.names=FALSE,quote=FALSE, sep="\t")
    return(dt)
}

resolve_hapblock <- function(fh) {
    uchrom <- unique(fh$chr)
    uchrom <- uchrom[uchrom!="Y"]
    hapblock <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)

    for (ic in uchrom) {
        dt.chr <- subset(fh, chr==ic)
	seg.end.index <- cumsum(rle(dt.chr$resHL)$lengths)
	seg.val <- rle(dt.chr$resHL)$values
	segL <- length(seg.end.index)
	if (segL > 1) seg.start.index <- c(1, seg.end.index[1:(segL-1)]+1) else seg.start.index <- 1
	seg.start <- dt.chr$pos[seg.start.index]
	seg.end <- dt.chr$pos[seg.end.index]
	hapblock.raw <- data.frame(chr=rep(ic, segL), start=seg.start, end=seg.end, resHL=seg.val)

        hapblock.final <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)
	for (i in (1:nrow(hapblock.raw))) {
	    hapblock.eval <- dt.chr[(dt.chr$pos>=hapblock.raw$start[i] & dt.chr$pos<=hapblock.raw$end[i]),]$pos
	    leneval <- which(diff(hapblock.eval)>3000000)
	    lastone <- length(hapblock.eval)

            if (length(leneval) != 0) {
	        leneval <- c(leneval, lastone)
		addinline <- hapblock.raw[i,]
		hapblock.raw$end[i] <- hapblock.eval[leneval[1]]
		hapblock.final <- rbind(hapblock.final, hapblock.raw[i,])
		for (j in (1:(length(leneval)-1))) {
		    addinline$start <- hapblock.eval[leneval[j]+1]
		    addinline$end <- hapblock.eval[leneval[j+1]]
		    hapblock.final <- rbind(hapblock.final, addinline)
		}
	    } else {
	        hapblock.final <- rbind(hapblock.final, hapblock.raw[i,])
	    }
	}
	hapblock <- rbind(hapblock, hapblock.final)
    }
    return(hapblock)
}


fetalmatimpute <- function(fetalfh) {
    fh <- readinfetalfh(fetalfh)
    ftmathap <- fh[,c(1,2,3,5)]
    imputemathap <- data.frame(chromosome=character(), start=numeric(), end=numeric(), mathap=character(), stringsAsFactors=FALSE)

    uchrom <- unique(ftmathap$chromosome)
    for (ic in uchrom) {
	imputeval <- "na"
	preimputeval <- "na"
	dt.chr <- subset(ftmathap, chromosome == ic)
	dt.chr$mathap <- as.character(dt.chr$mathap)
	for (i in 1:nrow(dt.chr)) {
	    if (dt.chr$mathap[i] != "Undetermined" & dt.chr$mathap[i] != "." & dt.chr$mathap[i] != preimputeval) {
		imputeval <- dt.chr$mathap[i]
		preimputeval <- imputeval
	    }
	    
	    if (dt.chr$mathap[i] != "MH1" & dt.chr$mathap[i] != "MH2") {
		dt.chr$mathap[i] = imputeval
	    }

        }
	imputemathap <- rbind(imputemathap, dt.chr)
    }
    imputemathap <- subset(imputemathap, mathap == "MH1" | mathap == "MH2")
    fo <- paste0(workdir, "/", sample, ".fetal.genome.impute.mat.hap")
    write.table(imputemathap, fo, row.names=FALSE, quote=FALSE, sep="\t")
    return(imputemathap)
}

fetalmatnoimpute <- function(fetalfh) {
    fh <- readinfetalfh(fetalfh)
    ftmathap <- fh[,c(1,2,3,5)]
    dt <- subset(ftmathap, mathap == "MH1" | mathap == "MH2")
    fo <- paste0(workdir, "/", sample, ".fetal.genome.noimpute.mat.hap")

    write.table(dt, fo, row.names=FALSE,quote=FALSE, sep="\t")
    return(dt)
}

resolve_fetalpat <- function(fetalfh, workdir, sample) {
    #fh <- fetalpatimpute(fetalfh)
    fh <- fetalpatnoimpute(fetalfh)
    fh$pathap <- gsub("PH1", 1, fh$pathap)
    fh$pathap <- gsub("PH2", 2, fh$pathap)
    hapblock <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)
    uchrom <- unique(fh$chromosome)
    uchrom <- uchrom[uchrom != "X"]

    for (ic in uchrom) {
        dt.chr <- subset(fh, chromosome==ic)
	seg.end.index <- cumsum(rle(dt.chr$pathap)$lengths)
	seg.val <- rle(dt.chr$pathap)$values
        segL <- length(seg.end.index)
        if (segL > 1) seg.start.index <- c(1, seg.end.index[1:(segL-1)]+1) else seg.start.index <-1
        seg.start <- dt.chr$end[seg.start.index]
        seg.end <- dt.chr$end[seg.end.index]
        hapblock.raw <- data.frame(chr=rep(ic, segL), start=seg.start, end=seg.end, resHL=seg.val)

        hapblock <- rbind(hapblock, hapblock.raw)

    }
    fo <- paste0(workdir, "/", sample, ".pat.resolved.hap")
    write.table(hapblock, fo, quote=FALSE, row.names=FALSE, sep="\t")

    fhb2 <- fh[,c(1,3,4)]
    colnames(fhb2) <- c("chr", "pos", "resHL")
    hapblock2 <- resolve_hapblock(fhb2)
    fo2 <- paste0(workdir, "/", sample, ".pat.resolved.2.hap")
    write.table(hapblock2, fo2, quote=FALSE, row.names=FALSE, sep="\t")
}



resolve_fetalmat <- function(fetalfh, workdir, sample) {
    #fh <- fetalmatimpute(fetalfh)
    fh <- fetalmatnoimpute(fetalfh)
    fh$mathap <- gsub("MH1", 1, fh$mathap)
    fh$mathap <- gsub("MH2", 2, fh$mathap)
    hapblock <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)
    uchrom <- unique(fh$chromosome)
    uchrom <- uchrom[uchrom != "21"]

    for (ic in uchrom) {
        dt.chr <- subset(fh, chromosome==ic)
        seg.end.index <- cumsum(rle(dt.chr$mathap)$lengths)
        seg.val <- rle(dt.chr$mathap)$values
        segL <- length(seg.end.index)
        if (segL > 1) seg.start.index <- c(1, seg.end.index[1:(segL-1)]+1) else seg.start.index <-1
        seg.start <- dt.chr$end[seg.start.index]
        seg.end <- dt.chr$end[seg.end.index]
        hapblock.raw <- data.frame(chr=rep(ic, segL), start=seg.start, end=seg.end, resHL=seg.val)

        hapblock <- rbind(hapblock, hapblock.raw)

    }
    fo <- paste0(workdir, "/", sample, ".mat.resolved.hap")
    write.table(hapblock, fo, quote=FALSE, row.names=FALSE, sep="\t")
    
    fhb2 <- fh[,c(1,3,4)]
    colnames(fhb2) <- c("chr", "pos", "resHL")
    hapblock2 <- resolve_hapblock(fhb2)
    fo2 <- paste0(workdir, "/", sample, ".mat.resolved.2.hap")
    write.table(hapblock2, fo2, quote=FALSE, row.names=FALSE, sep="\t")
}

resolve_fetalpat(fetalfh, workdir, sample)
resolve_fetalmat(fetalfh, workdir, sample)
