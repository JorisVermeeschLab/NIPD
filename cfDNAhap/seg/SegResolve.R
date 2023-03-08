Resolve_Inheritance_CBSseg <- function(infh1, infh2, ffsex) {
	rawdt <- read.delim(infh1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	rawsum <- read.delim(infh2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
	prodt <- data.frame(chr=as.character(), pos=as.numeric(), siteinfer=as.character(), FAR=as.numeric(), CBSmean=as.numeric(), resHL=as.character(), type=as.character(), stringsAsFactors=FALSE)
	prosum <- data.frame(chr=as.character(), start=as.numeric(), end=as.numeric(), resHL=as.character(), support=as.character(),stringsAsFactors=FALSE)

        if (nrow(rawdt) == 0 & nrow(rawsum) == 0) {
            fo1 <- gsub(".tsv", ".resolved.tsv", infh1)
            fo2 <- gsub(".tsv", ".resolved.sum.tsv", infh1)
            write.table(prodt, fo1, quote=FALSE, row.names=FALSE, sep="\t")
            write.table(prosum, fo2, quote=FALSE, row.names=FALSE, sep="\t")
	} else {	
	    uchrom <- unique(rawdt$chr)
	    uchrom <- uchrom[uchrom!="Y"]
	    #if male fetus
	    #if (ffsex == "male") {
		#uchrom <- uchrom[uchrom!="X"]
		#sexdt <- subset(rawdt, chr=="X" | chr=="Y")
                #sexdt <- subset(rawdt, chr=="Y")
		#sexprodt <- sexdt[,c(1,2,3,4,5,7,6)]
		#colnames(sexprodt) = c("chr", "pos", "siteinfer", "FAR", "CBSmean", "resHL", "type")
		#sexsum <- subset(rawsum, chrom=="X" | chrom=="Y")
                #sexsum <- subset(rawsum, chrom=="Y")
		#sexprosum <- sexsum[,c(2,3,4,7)]
		#colnames(sexprosum) = c("chr", "start", "end", "resHL")
	    #}

	    for (ic in uchrom) {
		    dt.chr <- subset(rawdt, chr==ic)
		    sortindex <- order(dt.chr$pos)

		    newdt.chr.T1 <- data.frame(chr=character(), pos=numeric(), siteinfer=character(), FAR=numeric(), CBSmean=numeric(), type=character(), resHL1=character())
		    newdt.chr.T2 <- data.frame(chr=character(), pos=numeric(), siteinfer=character(), FAR=numeric(), CBSmean=numeric(), type=character(), resHL2=character())
                    resHL1 <- rep(NA, nrow(dt.chr[sortindex, c(1,2,3,4,5,6)]))
		    resHL2 <- rep(NA, nrow(dt.chr[sortindex, c(1,2,3,4,5,6)]))
                    newdt.chr.T1 <- cbind(dt.chr[sortindex, c(1,2,3,4,5,6)], resHL1)
		    newdt.chr.T2 <- cbind(dt.chr[sortindex, c(1,2,3,4,5,6)], resHL2)

		    #check missing subtype
		    #missingPM1 <- subset(newdt.chr.T1, type=="P1" | type=="M1")
		    #missingPM2 <- subset(newdt.chr.T2, type=="P2" | type=="M2")
		    #if (nrow(missingPM1)==0) {
		    #    newdt.chr.T1 <- data.frame
		    #}
		    #newdt.chr.T1$resHL1 <- NULL
		    #newdt.chr.T2$resHL2 <- NULL
		    sum.chr <- subset(rawsum, chrom==ic)
		    

		    for (i in (1:nrow(sum.chr))) {
			    if (substr(sum.chr$ID[i], nchar(sum.chr$ID[i])-1, nchar(sum.chr$ID[i])) == "P1" || substr(sum.chr$ID[i], nchar(sum.chr$ID[i])-1, nchar(sum.chr$ID[i])) == "M1") {
				    segstart <- sum.chr$loc.start[i]
				    segend <- sum.chr$loc.end[i]
				    newdt.chr.T1$resHL1[newdt.chr.T1$pos >= segstart & newdt.chr.T1$pos <= segend] <- as.character(sum.chr$HL[i])
			    } else if (substr(sum.chr$ID[i], nchar(sum.chr$ID[i])-1, nchar(sum.chr$ID[i])) == "P2" || substr(sum.chr$ID[i], nchar(sum.chr$ID[i])-1, nchar(sum.chr$ID[i])) == "M2") {
				    segstart <- sum.chr$loc.start[i]
				    segend <- sum.chr$loc.end[i]
				    newdt.chr.T2$resHL2[newdt.chr.T2$pos >= segstart & newdt.chr.T2$pos <= segend] <- as.character(sum.chr$HL[i])
			    }
		    }
		    #check missing subtype
		    #if (!("resHL2" %in% names(newdt.chr.T2))) {
		    #    print("yes")
		    #    newdt.chr.T2$resHL2 <- NULL
		    #	print(newdt.chr.T2)
		    #}

		    newdt.chr <- merge(newdt.chr.T1, newdt.chr.T2, by=intersect(colnames(newdt.chr.T1), colnames(newdt.chr.T2)), all=TRUE)
		    sortindex <- order(newdt.chr$pos)
		    newdt.chr <- newdt.chr[sortindex,]
		    newdt.chr$resHL <- NULL
		    newdt.chr$support <- NULL
		    
		    for (j in (1:nrow(newdt.chr))) {

			    if (is.na(newdt.chr$resHL1[j]) || is.na(newdt.chr$resHL2[j])) {
				    if (!is.na(newdt.chr$resHL1[j])) {
					    newdt.chr$resHL[j] <- newdt.chr$resHL1[j]
					    newdt.chr$support[j] <- "S"
				    } else if (!is.na(newdt.chr$resHL2[j])) {
					    newdt.chr$resHL[j] <- newdt.chr$resHL2[j]
					    newdt.chr$support[j] <- "S"
				    }
			    } else {
				    if (newdt.chr$resHL1[j] == newdt.chr$resHL2[j]) {
					    newdt.chr$resHL[j] <- newdt.chr$resHL1[j]
					    newdt.chr$support[j] <- "B"
				    } else if (newdt.chr$resHL1[j] == "inconclusive" & newdt.chr$resHL2[j] != "inconclusive") {
				            newdt.chr$resHL[j] <- newdt.chr$resHL2[j]
					    newdt.chr$support[j] <- "S"

				    } else if (newdt.chr$resHL1[j] != "inconclusive" & newdt.chr$resHL2[j] == "inconclusive") {
				            newdt.chr$resHL[j] <- newdt.chr$resHL1[j]
					    newdt.chr$support[j] <- "S"
				    } else {
					    newdt.chr$resHL[j] <- "inconclusive"
					    newdt.chr$support[j] <- "N"
				    }
			    }
		    }
		    newdt.chr <- newdt.chr[,c("chr", "pos", "siteinfer", "FAR", "CBSmean", "resHL", "type", "support")]

		    #seg.end.index <- cumsum(rle(newdt.chr$resHL)$lengths)
		    #seg.val <- rle(newdt.chr$resHL)$values
		    multieval <- paste0(newdt.chr$resHL, ":", newdt.chr$support)
		    seg.end.index <- cumsum(rle(multieval)$lengths)
		    #seg.val <- rle(multieval)$values
		    seg.val1 <- newdt.chr$resHL[seg.end.index]
		    seg.val2 <- newdt.chr$support[seg.end.index]
		    segL <- length(seg.end.index)
		    if (segL > 1) seg.start.index <- c(1, seg.end.index[1:(segL-1)]+1) else seg.start.index <- 1
		    seg.start <- newdt.chr$pos[seg.start.index]
		    seg.end <- newdt.chr$pos[seg.end.index]
		    newsum.chr <- data.frame(chr=rep(ic, segL), start=seg.start, end=seg.end, resHL=seg.val1, support=seg.val2)
		    
                    newsum.chr.final <- data.frame(chr=character(), start=numeric(), end=numeric(), resHL=character(), support=character(), stringsAsFactors=FALSE)
                    #evaluate gap
		    for (i in (1:nrow(newsum.chr))) {
		        newsum.eval <- newdt.chr[(newdt.chr$pos>=newsum.chr$start[i] & newdt.chr$pos<=newsum.chr$end[i]),]$pos
			leneval <- which(diff(newsum.eval)>3000000)
			lastone <- length(newsum.eval)
                        
                        if (length(leneval) != 0) {
			    leneval <- c(leneval, lastone)
			    addinline <- newsum.chr[i,]
		            newsum.chr$end[i] <- newsum.eval[leneval[1]]
			    newsum.chr.final <- rbind(newsum.chr.final, newsum.chr[i,])
                            for (j in (1:(length(leneval)-1))) {
			    	addinline$start <- newsum.eval[leneval[j]+1]
				addinline$end <- newsum.eval[leneval[j+1]]
			        newsum.chr.final <- rbind(newsum.chr.final, addinline)
			    } 
			} else {
                            newsum.chr.final <- rbind(newsum.chr.final, newsum.chr[i,])
			}

		    }

                    prodt <- rbind(prodt, newdt.chr)
		    prosum <- rbind(prosum, newsum.chr.final)
	    }
	    #return(prodt)
	    #if (ffsex == "male") {
		#prodt <- rbind(prodt, sexprodt)
		#prosum <- rbind(prosum, sexprosum)
	    #}
	    fo1 <- gsub(".tsv", ".resolved.tsv", infh1)
	    fo2 <- gsub(".tsv", ".resolved.sum.tsv", infh1)
	    write.table(prodt, fo1, quote=FALSE, row.names=FALSE, sep="\t")
	    write.table(prosum, fo2, quote=FALSE, row.names=FALSE, sep="\t")
     }
}

Merge_Raw_Hetboth_Hetmat <- function(infh1, infh2, infh3) {
    hetpatsum <- read.delim(infh1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    #check if there is heteromat type data
    hetmatraw <- try(read.delim(infh2, header=FALSE, comment.char="#", sep="\t", stringsAsFactors=FALSE), silent=TRUE)
    if (inherits(hetmatraw, "try-error")) {
        hetmatraw <- data.frame(CHROM=as.numeric(), POS=as.numeric(), RSID=as.character(), REF=as.character(), ALT=as.character(), QUAL=as.character(), PATRAF=as.character(), MATRAF=as.character(), SIBRAF=as.character(), PH1=as.character(), PH2=as.character(), MH1=as.character(), MH2=as.character(), SH1=as.character(), SH2=as.character(), TYPE=as.character(), PATDP=as.character(), MATDP=as.character(), SIBDP=as.character())
    } else {
        hetmatraw <- hetmatraw
    }

    hetbothraw <- try(read.delim(infh3, header=FALSE, comment.char="#", sep="\t", stringsAsFactors=FALSE), silent=TRUE)
    if (inherits(hetbothraw, "try-error")) {
        hetbothraw <- data.frame(CHROM=as.numeric(), POS=as.numeric(), RSID=as.character(), REF=as.character(), ALT=as.character(), QUAL=as.character(), PATRAF=as.character(), MATRAF=as.character(), SIBRAF=as.character(), PH1=as.character(), PH2=as.character(), MH1=as.character(), MH2=as.character(), SH1=as.character(), SH2=as.character(), TYPE=as.character(), PATDP=as.character(), MATDP=as.character(), SIBDP=as.character())
    } else {
        hetbothraw <- hetbothraw
    }

    colnames(hetmatraw) <- colnames(hetbothraw) <- c("CHROM", "POS", "RSID", "REF", "ALT", "QUAL", "PATRAF", "MATRAF", "SIBRAF", "PH1", "PH2", "MH1", "MH2", "SH1", "SH2", "TYPE", "PATDP", "MATDP", "SIBDP")
    hetmatboth <- updatehetboth <- data.frame(CHROM=as.numeric(), POS=as.numeric(), RSID=as.character(), REF=as.character(), ALT=as.character(), QUAL=as.character(), PATRAF=as.character(), MATRAF=as.character(), SIBRAF=as.character(), PH1=as.character(), PH2=as.character(), MH1=as.character(), MH2=as.character(), SH1=as.character(), SH2=as.character(), TYPE=as.character(), PATDP=as.character(), MATDP=as.character(), SIBDP=as.character())
	
    hetpatref <- subset(hetpatsum, resHL != "inconclusive" & chr != "X")
    uchrom <- unique(hetpatref$chr)
    for (ic in uchrom) {
	    hetpatref.chr <- subset(hetpatref, chr==ic)
	    hetbothraw.chr <- subset(hetbothraw, CHROM==ic)
	    for (i in 1:nrow(hetpatref.chr)) {
		    refstart <- hetpatref.chr$start[i]
		    refend <- hetpatref.chr$end[i]
		    refHL <- hetpatref.chr$resHL[i]

		    if (refHL == "PH1" & nrow(hetbothraw.chr[hetbothraw.chr$POS >= refstart & hetbothraw.chr$POS <= refend, ]) != 0) {
			    hetbothraw.chr[hetbothraw.chr$POS >= refstart & hetbothraw.chr$POS <= refend, ]$TYPE <- "M1"
		    } else if (refHL == "PH2" & nrow(hetbothraw.chr[hetbothraw.chr$POS >= refstart & hetbothraw.chr$POS <= refend, ]) != 0) {
			    hetbothraw.chr[hetbothraw.chr$POS >= refstart & hetbothraw.chr$POS <= refend, ]$TYPE <- "M2"
		    }
	    }
	    hetbothpro.chr <- subset(hetbothraw.chr, TYPE=="M1" | TYPE=="M2")
	    updatehetboth <- rbind(updatehetboth, hetbothpro.chr)
    }
    hetmatboth <- rbind(hetmatraw, updatehetboth)
    hetmatboth$CHROM <- gsub("X", 23, hetmatboth$CHROM)
    sortindex <- order(as.numeric(hetmatboth$CHROM), hetmatboth$POS)
    hetmatboth <- hetmatboth[sortindex,]
    hetmatboth$CHROM <- gsub(23, "X", hetmatboth$CHROM)
    colnames(hetmatboth) <- c("#CHROM", "POS", "RSID", "REF", "ALT", "QUAL", "PATRAF", "MATRAF", "SIBRAF", "PH1", "PH2", "MH1", "MH2", "SH1", "SH2", "TYPE", "PATDP", "MATDP", "SIBDP")

    fo1 <- gsub("heteroboth.tsv", "heteromatboth.tsv", infh3)
    write.table(hetmatboth, fo1, quote=FALSE, row.names=FALSE, sep="\t")
}
