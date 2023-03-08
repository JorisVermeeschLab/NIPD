####segmentation on sibling phasing
script.dir <- dirname(sys.frame(1)$ofile)
pcfsp <- file.path(script.dir, "fastPCF.R")
source(pcfsp)

args = commandArgs(TRUE)
patfh <- args[1]
matfh <- args[2]
sibgender <- args[3]


Run_PCF <- function(x, winsize=20, kmin=10, gamma=50, penalty=T) {
    sdev = getMad(x, winsize)
    res = selectFastPcf(x, kmin, gamma*sdev, penalty)
    yhat = res$yhat
    return(yhat)
}


Paternal_RAF_Fitting <- function(patfh) {
    fh = read.delim(patfh, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
    fh = fh[,c(1,2,9,14,16)]
    colnames(fh) = c("CHROM", "POS", "SIBRAF", "SH1", "TYPE")
    fh = subset(fh, fh$CHROM != "Y")
    fh$SIBRAF = as.numeric(fh$SIBRAF)
    #flip value when sibling is 1|0
    fh$SIBRAF[fh$SH1==1] = 1-fh$SIBRAF[fh$SH1==1]
    fhP1 = subset(fh, TYPE=="P1")
    fhP2 = subset(fh, TYPE=="P2")
    chrlist = as.character(unique(fh$CHROM))
    tdf = data.frame(chr=character(), pos=numeric(), sibraf=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    for (i in 1:length(chrlist)) {
	fhP1chr = subset(fhP1, CHROM == chrlist[i])
	fhP2chr = subset(fhP2, CHROM == chrlist[i])
        xP1 = fhP1chr[,3]
	xP2 = fhP2chr[,3]
	if (length(xP1) != 0 && length(xP2) != 0) {
            yhatP1 = Run_PCF(xP1)
            yhatP2 = Run_PCF(xP2)
            calP1 = cbind(fhP1chr[,c(1,2)], xP1, yhatP1, rep("P1", length(yhatP1)))
            calP2 = cbind(fhP2chr[,c(1,2)], xP2, yhatP2, rep("P2", length(yhatP2)))
	    colnames(calP1) = colnames(calP2) = c("chr","pos","sibraf","yhat","type")
	    calchr = rbind(calP1, calP2)
	    #Y chromosome are all P2 category by definition
	} else if (sibgender == "male" && length(xP2) != 0) {
	    yhatP2 = Run_PCF(xP2)
	    calP2 = cbind(fhP2chr[,c(1,2)], xP2, yhatP2, rep("P2", length(yhatP2)))
            colnames(calP2) = c("chr","pos","sibraf","yhat","type")
	    calchr = calP2
	}
	tdf = rbind(tdf, calchr)
    }
    fo = gsub(".tsv", ".sib.bp.tsv", patfh)
    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
}


Maternal_RAF_Fitting = function(matfh, sibgender) {
    fh = read.delim(matfh, header=FALSE, sep="\t", comment.char="#")
    fh = fh[,c(1,2,9,14,16)]
    colnames(fh) = c("CHROM", "POS", "SIBRAF", "SH1", "TYPE")
    fh$SIBRAF[fh$SH1==1] = 1-fh$SIBRAF[fh$SH1==1]
    fhM1 = subset(fh, TYPE=="M1")
    fhM2 = subset(fh, TYPE=="M2")
    chrlist = as.character(unique(fh$CHROM))
    tdf = data.frame(chr=character(), pos=numeric(), sibraf=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    for (i in 1:length(chrlist)) {
        if (sibgender == "male" & chrlist[i] == "X") {
	    fhchr = subset(fh, CHROM == chrlist[i])
            x = 2*fhchr[, 3]
            yhat = Run_PCF(x)
            calchr = cbind(fhchr[,c(1,2)], x, yhat, rep("M1", length(yhat)))
            colnames(calchr) = c("chr","pos","sibraf","yhat","type")
	} else {
	    fhM1chr = subset(fhM1, CHROM == chrlist[i])
	    fhM2chr = subset(fhM2, CHROM == chrlist[i])
	    xM1 = fhM1chr[,3]
	    xM2 = fhM2chr[,3]
	    yhatM1 = Run_PCF(xM1)
	    yhatM2 = Run_PCF(xM2)
	    calM1 = cbind(fhM1chr[,c(1,2)], xM1, yhatM1, rep("M1", length(yhatM1)))
	    calM2 = cbind(fhM2chr[,c(1,2)], xM2, yhatM2, rep("M2", length(yhatM2)))
	    colnames(calM1) = colnames(calM2) = c("chr","pos","sibraf","yhat","type")
	    calchr = rbind(calM1, calM2)
        }
	tdf = rbind(tdf, calchr)
    }
    fo = gsub(".tsv", ".sib.bp.tsv", matfh)
    write.table(tdf, fo, quote=FALSE, row.names=FALSE, sep="\t")
}

sibpcf_main <- function(patfh, matfh, sibgender) {
    Paternal_RAF_Fitting(patfh)
    Maternal_RAF_Fitting(matfh, sibgender)   
}

sibpcf_main(patfh, matfh, sibgender)
