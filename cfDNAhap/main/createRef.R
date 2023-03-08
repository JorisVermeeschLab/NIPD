.libPaths(c("/ddn1/vol1/staging/leuven/stg_00002/cgr/Huiwen/sw/R/3.4.2",.libPaths()))

library(tidyr)
library(data.table)
##create reference set for plasma DNA and gDNA samples

GetTarGC <- function(targetgcfh) {
    gcfh <- read.delim(targetgcfh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    gcfh <- gcfh[,c(1,2,3,6,14)]
    colnames(gcfh) <- c("CHROM","START","END","GC","INDEX")
    return(gcfh)
}

GetSampleList <- function(samplelist) {
    splist <- read.delim(samplelist, header=TRUE, sep="\t")
    return(splist)
}

GetSampleTarCount <- function(sampledir, samplename, fragsis, targetgcfh) {
    fo <- paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.noref.tsv")
    if (file.exists(fo)) {
        fragout <- read.delim(fo, header=TRUE, sep="\t")
	
    } else {
        #samplefh <- paste0(sampledir, "/4_cov/", samplename, "/", samplename, ".sorted.merged.md.filtered.fragcov.sample_interval_summary")

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
        fo = paste0(sampledir, "/5_variant/5_cnv/", samplename, ".target.coverage.noref.tsv")
	write.table(fragout, file = fo, sep = "\t", row.names = FALSE, quote = FALSE)
    }
    return(fragout)
}

UpdateMeanSD <- function(n1, sample1.mean, sample1.sd, n2, sample2.mean) {
    if (n2==1) {
        new.mean <- (sample1.mean*n1+sample2.mean)/(n1+1)
	new.sd <- sqrt(((sample2.mean-new.mean)*(sample2.mean-sample1.mean)+(n1-1)*sample1.sd**2)/n1)
	return(c(new.mean, new.sd))
    } else {
        print("Applies to updating one added value!")
    }
}

GetCombinedCov <- function(samplelist, targetgcfh, outdir, setname) {
    splist <- GetSampleList(samplelist)

    mergedt.raw <- data.frame(INDEX=numeric(), CHROM=character(), START=numeric(), END=numeric(), MEAN.RAW.NORM=numeric(), stringsAsFactors=FALSE)
    mergedt.gcc <- data.frame(INDEX=numeric(), CHROM=character(), START=numeric(), END=numeric(), MEAN.GCC.NORM=numeric(), stringsAsFactors=FALSE)

    for (i in 1:dim(splist)[1]) {
        sampledir <- as.character(splist$SAMPLEDIR[i])
	samplename <- as.character(splist$SAMPLENAME[i])
	fragsis <- as.character(splist$FRAGSIS[i])
	fragout <- GetSampleTarCount(sampledir, samplename, fragsis, targetgcfh)
	if (i==1) {
           mergedt.raw <- fragout[,c("INDEX","CHROM","START","END","MEAN.RAW.NORM")]
	   colnames(mergedt.raw)[colnames(mergedt.raw)=="MEAN.RAW.NORM"] <- paste0("m.raw", i)
	   mergedt.gcc <- fragout[,c("INDEX","CHROM","START","END","MEAN.GCC.NORM")]
	   colnames(mergedt.gcc)[colnames(mergedt.gcc)=="MEAN.GCC.NORM"] <- paste0("m.gcc", i)
        } else {
	   selcolraw <- paste0("m.raw", i)
	   selcolgcc <- paste0("m.gcc", i)
	   colnames(fragout)[colnames(fragout)=="MEAN.RAW.NORM"] <- selcolraw
	   colnames(fragout)[colnames(fragout)=="MEAN.GCC.NORM"] <- selcolgcc
 	   mergedt.raw <- merge(mergedt.raw, fragout[, c("INDEX", "CHROM", "START", "END", selcolraw)], by=c("INDEX", "CHROM", "START", "END"))
	   mergedt.gcc <- merge(mergedt.gcc, fragout[, c("INDEX", "CHROM", "START", "END", selcolgcc)], by=c("INDEX", "CHROM", "START", "END"))
	   #print(head(mergedt.gcc, n=10L))
	}
    }
    mergedt.raw$SET.RAW.MEAN <- apply(mergedt.raw[,-c(1,2,3,4)], 1, mean)
    mergedt.raw$SET.RAW.SD <- apply(mergedt.raw[,-c(1,2,3,4)], 1, sd)
    mergedt.gcc$SET.GCC.MEAN <- apply(mergedt.gcc[,-c(1,2,3,4)], 1, mean)
    mergedt.gcc$SET.GCC.SD <- apply(mergedt.gcc[,-c(1,2,3,4)], 1, sd)
    mergedt.raw <- mergedt.raw[,c("INDEX", "CHROM", "START", "END", "SET.RAW.MEAN", "SET.RAW.SD")]
    mergedt.gcc <- mergedt.gcc[,c("INDEX", "CHROM", "START", "END", "SET.GCC.MEAN", "SET.GCC.SD")]

    outdt <- merge(mergedt.raw, mergedt.gcc, by=c("INDEX", "CHROM", "START", "END"))
    outdt <- outdt[order(outdt$INDEX),]
    fo <- paste0(outdir, "/", setname, ".reference.tsv")
    write.table(outdt, fo, sep="\t", row.names=FALSE, quote=FALSE)
}


