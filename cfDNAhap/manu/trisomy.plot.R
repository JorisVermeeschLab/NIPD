source("./seg/CBS.R")

list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)
library(ggplot2)

args <- commandArgs(TRUE)
cnvfh <- args[1]
patbaf <- args[2]
matbaf <- args[3]
samplename <- args[4]
chrom <- args[5]
patbafsib <- args[6]
matbafsib <- args[7]
figdir <- args[8]

getcnv <- function(cnvfh) {
    fh <- read.delim(cnvfh, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    fh <- fh[,c("CHROM","START","END","MEAN.GCC.NORM","Z","LOG2R","Seg.Z","Seg.LOG2R")]
    colnames(fh) <- c("chromosome","start","end","MEAN.GCC.NORM","Z","LOG2R","Seg.Z","Seg.LOG2R")
    #fh$chromosome <- paste0("chr", fh$chromosome)
    return(fh)
}

cnvplot <- function(cnvfh,samplename, chrom, figdir) {
    fh <- getcnv(cnvfh)
    plotchr <- paste0("chr",seq(1,22))
    plotchr <- c(plotchr, "chrX")

    if (chrom=="all") {
        fo <- paste0(figdir, samplename, ".CNV.manu.pdf")
        pdf(file=fo, width=9, height=3)
        #png(file=fo, width=800, height=300)
        gtrellis_layout(n_track = 1, category=plotchr, remove_chr_prefix=TRUE, track_axis=TRUE, xaxis = FALSE, asist_ticks = FALSE, border=TRUE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(-3, 3), track_ylab = "Log2 Ratio", add_ideogram_track =FALSE, add_name_track = TRUE, axis_label_fontsize = 5, lab_fontsize = 6, name_fontsize = 7, title=NULL, title_fontsize = 10, xlab = NULL)

        add_points_track(fh, fh$LOG2R, pch=1, size=unit(0.05, "mm"), gp=gpar(col="black"))
        add_points_track(fh, fh$Seg.LOG2R, track = current_track(), pch=15, size=unit(0.5, "mm"), gp=gpar(col="red"))
        dev.off()
    } else {
        fh <- subset(fh, chromosome==chrom)
        plotchr <- paste0("chr",chrom)
        fo <- paste0(figdir, samplename, ".CNV.", chrom, ".manu.pdf")
        pdf(file=fo, width=9, height=3)
        gtrellis_layout(n_track = 1, category=plotchr, remove_chr_prefix=TRUE, track_axis=TRUE, xaxis = TRUE, asist_ticks = FALSE, border=TRUE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(-3, 3), track_ylab = "Log2 Ratio", add_ideogram_track =FALSE, add_name_track = TRUE, axis_label_fontsize = 5, lab_fontsize = 6, name_fontsize = 7, title=NULL, title_fontsize = 10, xlab = NULL)
        add_points_track(fh, fh$LOG2R, pch=1, size=unit(0.05, "mm"), gp=gpar(col="black"))
        add_points_track(fh, fh$Seg.LOG2R, track = current_track(), pch=15, size=unit(0.5, "mm"), gp=gpar(col="red"))
        dev.off()
    }
}

cnvplot(cnvfh, samplename, chrom, figdir)


getphaseorder <- function(baf1, baf2) {
    baf1 <- read.delim(baf1, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
    colnames(baf1) <- c("CHROM", "POS","RSID","REF","ALT","QUAL","PATRAF","MATRAF","SIBRAF1","PH1","PH2","MH1","MH2","SH1","SH2","TYPE","PATDP","MATDP","SIBDP")
    print(head(baf1, n=10L))

    baf2 <- read.delim(baf2, header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
    colnames(baf2) <- c("CHROM", "POS","RSID","REF","ALT","QUAL","PATRAF","MATRAF","SIBRAF2","PH1","PH2","MH1","MH2","SH1","SH2","TYPE","PATDP","MATDP","SIBDP")
    baf2 <- baf2[,c("CHROM", "POS", "SIBRAF2")]
    print(head(baf2, n=10L))
    
    dt <- merge(baf1, baf2, by=c("CHROM", "POS"))
    print(head(dt, n=10L))
    dt <- dt[,c("CHROM", "POS","RSID","REF","ALT","QUAL","PATRAF","MATRAF","SIBRAF2","PH1","PH2","MH1","MH2","SH1","SH2","TYPE","PATDP","MATDP","SIBDP")]
    colnames(dt) <- c("CHROM", "POS","RSID","REF","ALT","QUAL","PATRAF","MATRAF","SIBRAF","PH1","PH2","MH1","MH2","SH1","SH2","TYPE","PATDP","MATDP","SIBDP")
    return(dt)
}



getbaf <- function(baf1, baf2, chrom) {
    fh <- getphaseorder(baf1, baf2)
    #fh <- read.delim(bafh, header=FALSE, sep="\t", comment.char = "#",	stringsAsFactors=FALSE)
    #colnames(fh) <- c("CHROM", "POS","RSID","REF","ALT","QUAL","PATRAF","MATRAF","SIBRAF","PH1","PH2","MH1","MH2","SH1","SH2","TYPE","PATDP","MATDP","SIBDP")
    fh$CHROM <- as.character(fh$CHROM)
    chrom <- as.character(chrom)
    fh <- subset(fh, CHROM==chrom)

    raf <- rep(-1, nrow(fh))
    
    for (i in 1:nrow(fh)) {
        if (fh$PH1[i]==1 && fh$PH2[i]==0) {
            raf[i] <- 1-as.numeric(as.character(fh$SIBRAF[i]))

	} else if (fh$TYPE[i]=="M2" && fh$PH1[i]==0) {
	    raf[i] <- as.numeric(as.character(fh$SIBRAF[i]))
	    fh$TYPE[i] <- "M1"
	} else {
	    raf[i] <- as.numeric(as.character(fh$SIBRAF[i]))
	}
    }

    sfh <- fh[,c("CHROM", "POS", "TYPE")]
    colnames(sfh) <- c("chromosome", "end", "type")
    pos1 <- sfh$end-1
    sfh <- cbind(sfh[,1], pos1, raf, sfh[,2:3])
    colnames(sfh) <- c("chrom", "pos", "far", "end", "type")
 
    sfh1 <- subset(sfh, type=="P1" | type=="M1")
    ina1 <- (!is.na(sfh1$chrom) & is.finite(sfh1$pos))
    sortindex1 <- which(ina1)[order(sfh1$chrom[ina1], sfh1$pos[ina1])]
    sfh.cal1 <- sfh1[sortindex1,]
    sfh.calsm1 <- smoothFAR(sfh.cal1)

    sfh.seg1 <- segment(sfh.calsm1,"RAF")
    sfh.out1 <- sfh.seg1$output
    #print(sfh.out1)
    sfh1$Seg.RAF <- 0

    for (i in 1:nrow(sfh.out1)) {
        sfh1[(sfh1$pos >= sfh.out1$loc.start[i] & sfh1$pos <= sfh.out1$loc.end[i]),]$Seg.RAF <- sfh.out1$seg.mean[i]
    }

    sfh2 <- subset(sfh, type=="P2" | type=="M2")
    ina2 <- (!is.na(sfh2$chrom) & is.finite(sfh2$pos))
    sortindex2 <- which(ina2)[order(sfh2$chrom[ina2], sfh2$pos[ina2])]
    sfh.cal2 <- sfh2[sortindex2,]
    print(head(sfh.cal2, n=10L))
    sfh.calsm2 <- smoothFAR(sfh.cal2)

    sfh.seg2 <- segment(sfh.calsm2,"RAF")
    sfh.out2 <- sfh.seg2$output
    #print(sfh.out2)
    sfh2$Seg.RAF <- 0
    
    for (i in 1:nrow(sfh.out2)) {
        sfh2[(sfh2$pos >= sfh.out2$loc.start[i] & sfh2$pos <= sfh.out2$loc.end[i]),]$Seg.RAF <- sfh.out2$seg.mean[i]
    }

    sfh <- rbind(sfh1, sfh2)
    sfh <- sfh[,c(1,2,4,3,6,5)]
    colnames(sfh) <- c("chromosome","start","end","RAF","Seg.RAF","type")
    print(sfh)
    sfh$chromosome <- paste0("chr",sfh$chromosome)
    return(sfh)
}

bafplot <- function(patbaf, patbafsib, matbaf, matbafsib, samplename, chrom, figdir) {
    patbaf <- getbaf(patbafsib,patbaf,chrom)
    matbaf <- getbaf(matbafsib,matbaf,chrom)
    chrom <- paste0("chr", chrom)

    fo1 <- paste0(figdir, samplename, ".", chrom, ".raf.manu.pdf")
    pdf(file=fo1, width=4.5, height=3)
    gtrellis_layout(n_track=2, category=chrom, track_axis=TRUE, xaxis=FALSE, asist_ticks=FALSE, border=TRUE, padding=unit(c(2,2,4,2),"mm"), track_ylim=c(-0.1,1.1,-0.1,1.1), track_ylab=c("Pat RAF","Mat RAF"), add_ideogram_track=FALSE, add_name_track=TRUE, axis_label_fontsize=5, lab_fontsize=6, name_fontsize=7, title=NULL, title_fontsize=10, xlab=NULL, remove_chr_prefix=FALSE)

    add_track(patbaf, panel_fun = function(gr) {
        x <- patbaf$start
	y <- patbaf$RAF
	grid.lines(unit(c(0, 1), "npc"), unit(c(1, 1), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	grid.lines(unit(c(0, 1), "npc"), unit(c(0.5, 0.5), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	#grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patbaf$type=="P1", "#c82027", "#344D90")))
    })
    add_points_track(patbaf, patbaf$RAF, pch=1, track=current_track(), size=unit(0.5, "mm"), gp=gpar(col=ifelse(patbaf$type=="P1", "#b8e186", "#276419")))
    add_points_track(patbaf, patbaf$Seg.RAF, pch=16, track=current_track(), size=unit(1,"mm"), gp=gpar(col="black"))

    add_track(matbaf, panel_fun = function(gr) {
        x <- matbaf$start
	y <- matbaf$RAF
        grid.lines(unit(c(0, 1), "npc"), unit(c(1, 1), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	grid.lines(unit(c(0, 1), "npc"), unit(c(0.5, 0.5), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
        grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
    })
    add_points_track(matbaf, matbaf$RAF, pch=1, track=current_track(), size=unit(0.5, "mm"), gp=gpar(col=ifelse(matbaf$type=="M1", "#c51b7d", "#f1b6da")))
    add_points_track(matbaf, matbaf$Seg.RAF, pch=16, track=current_track(), size=unit(1,"mm"), gp=gpar(col="black"))
    dev.off()
}

#bafplot(patbaf, matbaf, samplename, chrom)
bafplot(patbaf, patbafsib, matbaf, matbafsib, samplename, chrom, figdir)

#patbafh <- getbaf(patbaf,chrom)

#fo1 <- paste0(figdir, samplename, ".", chrom, ".patraf.manu.pdf")
#pdf(file=fo1, width=4.5, height=3)
#ggplot(patbafh, aes(start, RAF))+geom_point(size=0.5, aes(color=type))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major = element_line(linetype = "dashed"),panel.grid.minor = element_blank(),legend.position="none")+scale_color_manual(values = c("#b8e186", "#276419"))+labs(y="Reference Allele Frequency", x="Chromosome 12")
#dev.off()

#matbafh <- getbaf(matbaf,chrom)

#fo2 <- paste0(figdir, samplename, ".", chrom, ".matraf.manu.pdf")
#pdf(file=fo2, width=4.5, height=3)
#ggplot(matbafh, aes(start, RAF))+geom_point(size=0.5, aes(color=type))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.major = element_line(linetype = "dashed"),panel.grid.minor = element_blank(),legend.position="none")+scale_color_manual(values = c("#c51b7d", "#f1b6da"))+labs(y="Reference Allele Frequency", x="Chromosome 12")
#dev.off()
