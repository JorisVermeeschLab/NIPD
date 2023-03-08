script.dir <- dirname(sys.frame(1)$ofile)
CBStoLoad <- file.path(script.dir, "CBS.R")
source(CBStoLoad)

library(tidyr)

ReadReference <- function(refdir, refname, gcc) {
    reference <- paste0(refdir, "/", refname, ".tsv")
    ref <- read.delim(reference, header=TRUE, sep="\t")
    if (gcc) {
	ref <- ref[, c("INDEX", "CHROM", "START", "END", "SET.GCC.MEAN", "SET.GCC.SD")]
    } else {
        ref <- ref[, c("INDEX", "CHROM", "START", "END", "SET.RAW.MEAN", "SET.RAW.SD")]
    }
    colnames(ref) <- c("INDEX", "CHROM", "START", "END", "SET.MEAN", "SET.SD")
    return(ref)
}

GetTarGC <- function(targetgcfh) {
    gcfh <- read.delim(targetgcfh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    gcfh <- gcfh[,c(1,2,3,6,14)]
    colnames(gcfh) <- c("CHROM","START","END","GC","INDEX")
    return(gcfh)
}

GetSampleTarCount <- function(sampledir, samplename, fragsis, targetgcfh) {
    fo <- paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.noref.tsv")
    if (file.exists(fo)) {
	fragout <- read.delim(fo, header=TRUE, sep="\t")

    } else {
        fragcov <- read.delim(fragsis, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	#get Target, avg_coverage (= total.fragment.coverage/target.size)
        fragcov <- fragcov[,c(1,3)]
	fragnew <- fragcov %>% separate(col=Target, into=c("chrom","start"), sep=":") %>% separate(col=start, into=c("start", "end"), sep="-")
	colnames(fragnew) <- c("CHROM", "START", "END", "MEAN")
	fragnew$START <- as.numeric(fragnew$START)-1
	fragnew$CHROM <- as.character(fragnew$CHROM)

        #apply GC correction
	gcfh <- GetTarGC(targetgcfh)
	fragmerge <- merge(fragnew, gcfh, by=c("CHROM", "START", "END"))
	fragloess <- subset(fragmerge, fragmerge$CHROM != "X" & fragmerge$CHROM != "Y" & fragmerge$MEAN != 0)
        #median coverage per GC level
        targcmedian <- tapply(fragloess$MEAN, fragloess$GC, median)
	gc.cor.factor <- predict(loess(targcmedian ~ as.numeric(names(targcmedian))), fragmerge$GC)
	fragmerge$GCC.MEAN <- fragmerge$MEAN * median(fragmerge$MEAN)/gc.cor.factor
	fragmerge$GCC.MEAN[!is.finite(fragmerge$GCC.MEAN)] <- 0

        #Normalization
	autofragmerge <- subset(fragmerge, CHROM != "X" & CHROM != "Y")
        auto.raw.median <- median(as.numeric(autofragmerge$MEAN), na.rm=TRUE)
	auto.gcc.median <- median(as.numeric(autofragmerge$GCC.MEAN), na.rm=TRUE)
	fragmerge$MEAN.RAW.NORM <- fragmerge$MEAN/auto.raw.median
	fragmerge$MEAN.GCC.NORM <- fragmerge$GCC.MEAN/auto.gcc.median
	fragmerge$MEAN.GCC.NORM[fragmerge$MEAN.GCC.NORM < 0] <- 0

        fragout <- fragmerge[,c("INDEX","CHROM", "START", "END", "MEAN", "MEAN.RAW.NORM", "MEAN.GCC.NORM")]
	fragout <- fragout[order(fragout$INDEX),]
	fo <- paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.noref.tsv")
	write.table(fragout, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)
	}
    return(fragout)
}

GetRelCov <- function(sampledir, samplename, fragsis, targetgcfh, refdir, refname, gcc) {
    ref <- ReadReference(refdir, refname, gcc)
    samplefh <- GetSampleTarCount(sampledir, samplename, fragsis, targetgcfh)
    if (gcc) {
        samplefh <- samplefh[,c("INDEX", "CHROM", "START", "END", "MEAN.GCC.NORM")]
	dt <- merge(samplefh, ref, by=c("INDEX", "CHROM", "START", "END"))
	dt$Z <- (dt$MEAN.GCC.NORM - dt$SET.MEAN)/dt$SET.SD
	dt$LOG2R <- log2(dt$MEAN.GCC.NORM/dt$SET.MEAN)
	dt <- dt[,c("INDEX", "CHROM", "START", "END", "MEAN.GCC.NORM", "Z", "LOG2R")]
	dt <- dt[order(dt$INDEX),]

        ###Segmentation on Z-score
        cbs.Z <- dt[,c("CHROM", "START", "Z", "END", "INDEX")]
	colnames(cbs.Z) <- c("chrom", "pos", "far", "END", "INDEX")
	cbs.Z$chrom <- as.character(cbs.Z$chrom)
	
        ina <- (!is.na(cbs.Z$chrom) & is.finite(cbs.Z$pos))
	sortindex <- which(ina)[order(cbs.Z$chrom[ina], cbs.Z$pos[ina])]
	cbs.Z.cal <- cbs.Z[sortindex,]
        cbs.Z.calsm <- smoothFAR(cbs.Z.cal)

	cbs.Z.seg <- segment(cbs.Z.calsm,"Z")
	cbs.Z.out <- cbs.Z.seg$output
	print(cbs.Z.out)

        cbs.Z.out$chrom[cbs.Z.out$chrom=="24"] <- "Y"
	cbs.Z.out$chrom[cbs.Z.out$chrom=="23"] <- "X"
	dt$Seg.Z <- 0

        for (i in 1:nrow(cbs.Z.out)) {
            dt[(dt$CHROM == cbs.Z.out$chrom[i] & dt$START >= cbs.Z.out$loc.start[i] & dt$START <= cbs.Z.out$loc.end[i]),]$Seg.Z <- cbs.Z.out$seg.mean[i]
	}

        ###Segmentation on log2R
	cbs.log2r <- dt[,c("CHROM", "START", "LOG2R", "END", "INDEX")]
	colnames(cbs.log2r) <- c("chrom", "pos", "far", "END", "INDEX")
	cbs.log2r$chrom <- as.character(cbs.log2r$chrom)

        ina <- (!is.na(cbs.log2r$chrom) & is.finite(cbs.log2r$pos))
	sortindex <- which(ina)[order(cbs.log2r$chrom[ina], cbs.log2r$pos[ina])]
        cbs.log2r.cal <- cbs.log2r[sortindex,]
	cbs.log2r.calsm <- smoothFAR(cbs.log2r.cal)

        cbs.log2r.seg <- segment(cbs.log2r.calsm,"log2R")
	cbs.log2r.out <- cbs.log2r.seg$output
        print(cbs.log2r.out)

        cbs.log2r.out$chrom[cbs.log2r.out$chrom=="24"] <- "Y"
	cbs.log2r.out$chrom[cbs.log2r.out$chrom=="23"] <- "X"
        dt$Seg.LOG2R <- 0
	
	for (i in 1:nrow(cbs.log2r.out)) {
            dt[(dt$CHROM == cbs.log2r.out$chrom[i] & dt$START >= cbs.log2r.out$loc.start[i] & dt$START <= cbs.log2r.out$loc.end[i]),]$Seg.LOG2R <- cbs.log2r.out$seg.mean[i]
        }
	
    } else {
        samplefh <- samplefh[,c("INDEX", "CHROM", "START", "END", "MEAN.RAW.NORM")]
	dt <- merge(samplefh, ref, by=c("INDEX", "CHROM", "START", "END"))
	dt$Z <- (dt$MEAN.RAW.NORM - dt$SET.MEAN)/dt$SET.SD
	dt$LOG2R <- log2(dt$MEAN.RAW.NORM/dt$SET.MEAN)
	dt <- dt[,c("INDEX", "CHROM", "START", "END", "MEAN.RAW.NORM", "Z", "LOG2R")]
	dt <- dt[order(dt$INDEX),]

        ###Segmentation on Z-score
	cbs.Z <- dt[,c("CHROM", "START", "Z", "END", "INDEX")]
        colnames(cbs.Z) <- c("chrom", "pos", "far", "END", "INDEX")
	cbs.Z$chrom <- as.character(cbs.Z$chrom)

        ina <- (!is.na(cbs.Z$chrom) & is.finite(cbs.Z$pos))
	sortindex <- which(ina)[order(cbs.Z$chrom[ina], cbs.Z$pos[ina])]
        cbs.Z.cal <- cbs.Z[sortindex,]
	cbs.Z.calsm <- smoothFAR(cbs.Z.cal)

        cbs.Z.seg <- segment(cbs.Z.calsm,"Z")
	cbs.Z.out <- cbs.Z.seg$output

        cbs.Z.out$chrom[cbs.Z.out$chrom=="24"] <- "Y"
	cbs.Z.out$chrom[cbs.Z.out$chrom=="23"] <- "X"
        dt$Seg.Z <- 0

        for (i in 1:nrow(cbs.Z.out)) {
	    dt[(dt$CHROM == cbs.Z.out$chrom[i] & dt$START >= cbs.Z.out$loc.start[i] & dt$START <= cbs.Z.out$loc.end[i]),]$Seg.Z <- cbs.Z.out$seg.mean[i]
        }

        ###Segmentation on log2R
	cbs.log2r <- dt[,c("CHROM", "START", "LOG2R", "END", "INDEX")]
        colnames(cbs.log2r) <- c("chrom", "pos", "far", "END", "INDEX")
        cbs.log2r$chrom <- as.character(cbs.log2r$chrom)

        ina <- (!is.na(cbs.log2r$chrom) & is.finite(cbs.log2r$pos))
	sortindex <- which(ina)[order(cbs.log2r$chrom[ina], cbs.log2r$pos[ina])]
        cbs.log2r.cal <- cbs.log2r[sortindex,]
	cbs.log2r.calsm <- smoothFAR(cbs.log2r.cal)

        cbs.log2r.seg <- segment(cbs.log2r.calsm,"log2R")
	cbs.log2r.out <- cbs.log2r.seg$output

        cbs.log2r.out$chrom[cbs.log2r.out$chrom=="24"] <- "Y"
	cbs.log2r.out$chrom[cbs.log2r.out$chrom=="23"] <- "X"
        dt$Seg.LOG2R <- 0

        for (i in 1:nrow(cbs.log2r.out)) {
	    dt[(dt$CHROM == cbs.log2r.out$chrom[i] & dt$START >= cbs.log2r.out$loc.start[i] & dt$START <= cbs.log2r.out$loc.end[i]),]$Seg.LOG2R <- cbs.log2r.out$seg.mean[i]
	}
    }

    dt <- dt[order(dt$INDEX),]
    if (gcc) {
        fo <- paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.", refname, ".gcc.tsv")
    } else {
        fo <- paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.", refname, ".raw.tsv")
    }
    write.table(dt, fo, sep="\t", row.names=FALSE, quote=FALSE)
}
