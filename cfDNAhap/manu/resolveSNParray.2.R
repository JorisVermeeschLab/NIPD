args <- commandArgs(TRUE)
PGD_Itp2 <- args[1]
cfpat <- args[2]
cfmat <- args[3]

resolve_hapblock <- function(fh) {
    uchrom <- unique(fh$chr)
    uchrom <- uchrom[uchrom!="Y"]
    hapblock <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)

    for (ic in uchrom) {
        dt.chr <- subset(fh, chr==ic)
        seg.end.index <- cumsum(rle(dt.chr$resHL)$lengths)
        seg.val <- rle(dt.chr$resHL)$values
        segL <- length(seg.end.index)
        if (segL > 1) seg.start.index <- c(1, seg.end.index[1:(segL-1)]+1) else seg.start.index <-1
        seg.start <- dt.chr$pos[seg.start.index]
        seg.end <- dt.chr$pos[seg.end.index]
        hapblock.raw <- data.frame(chr=rep(ic, segL), start=seg.start, end=seg.end, resHL=seg.val)

        #hapblock.final <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), stringsAsFactors=FALSE)
        #for (i in (1:nrow(hapblock.raw))) {
        #    hapblock.eval <- dt.chr[(dt.chr$pos>=hapblock.raw$start[i] & dt.chr$pos<=hapblock.raw$end[i]),]$pos
	#    leneval <- which(diff(hapblock.eval)>3000000)
	#    lastone <- length(hapblock.eval)

	#    if (length(leneval) != 0) {
	#        leneval <- c(leneval, lastone)
	#        addinline <- hapblock.raw[i,]
	#        hapblock.raw$end[i] <- hapblock.eval[leneval[1]]
	#        hapblock.final <- rbind(hapblock.final, hapblock.raw[i,])
	#        for (j in (1:(length(leneval)-1))) {
	#            addinline$start <- hapblock.eval[leneval[j]+1]
        #		    addinline$end <- hapblock.eval[leneval[j+1]]
	#	    hapblock.final <- rbind(hapblock.final, addinline)		
	#       }
	#    } else {
        #        hapblock.final <- rbind(hapblock.final, hapblock.raw[i,])
	#    }
        #}
	#hapblock <- rbind(hapblock, hapblock.final)
        hapblock <- rbind(hapblock, hapblock.raw)
    }
    return(hapblock)
}

Get_hapblock_snparr <- function(PGD_Itp2) {
    fh <- read.delim(PGD_Itp2, header=TRUE, sep="\t")
    colnames(fh) <- c("chr", "pos", "pat", "mat")
    fhpat <- fh[,c(1,2,3)]
    colnames(fhpat) <- c("chr", "pos", "resHL")
    fhmat <- fh[,c(1,2,4)]
    colnames(fhmat) <- c("chr", "pos", "resHL")
    fhpat <- subset(fhpat, resHL == 1 | resHL == 2)
    fhmat <- subset(fhmat, resHL == 1 | resHL == 2)

    pat.resolve.hap <- resolve_hapblock(fhpat)
    mat.resolve.hap <- resolve_hapblock(fhmat)

    fo1 <- paste0(PGD_Itp2, ".pat.resolve.2.block")
    fo2 <- paste0(PGD_Itp2, ".mat.resolve.2.block")
    write.table(pat.resolve.hap, fo1, quote=FALSE, row.names=FALSE, sep="\t")
    write.table(mat.resolve.hap, fo2, quote=FALSE, row.names=FALSE, sep="\t")
}

Get_hapblock_cf <- function(cfpat, cfmat) {
	pat <- read.delim(cfpat, header=TRUE, sep="\t")
	mat <- read.delim(cfmat, header=TRUE, sep="\t")
        patsub <- subset(pat, resHL != "inconclusive")
        matsub <- subset(mat, resHL != "inconclusive")
        patsub$resHL <- gsub("PH1", 1, patsub$resHL)
        patsub$resHL <- gsub("PH2", 2, patsub$resHL)
        matsub$resHL <- gsub("MH1", 1, matsub$resHL)
        matsub$resHL <- gsub("MH2", 2, matsub$resHL)

	pat.resolve.hap <- resolve_hapblock(patsub)
	mat.resolve.hap <- resolve_hapblock(matsub)
		
	fo1 <- paste0(cfpat, ".resolve.2.block")
	fo2 <- paste0(cfmat, ".resolve.2.block")
		
	write.table(pat.resolve.hap, fo1, quote=FALSE, row.names=FALSE, sep="\t")
	write.table(mat.resolve.hap, fo2, quote=FALSE, row.names=FALSE, sep="\t")
}

Get_hapblock_snparr(PGD_Itp2)
Get_hapblock_cf(cfpat, cfmat)
