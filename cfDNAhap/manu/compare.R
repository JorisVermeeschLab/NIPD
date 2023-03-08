list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args <- commandArgs(TRUE)
patfh <- args[1]
matfh <- args[2]
patres <- args[3]
matres <- args[4]
samplename <- args[5]
ffest_scale <- args[6] #input percent - if FF 10%, then input 10
chrom <- args[7] #chr1
snpP1 <- args[8]
snpP2 <- args[9]
snpM1 <- args[10]
snpM2 <- args[11]
snphap <- args[12]
diseaseStart <- args[13]
diseaseEnd <- args[14]
gene <- args[15]
figdir <- args[16]

Get_SNPratio <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 <- fh$pos-1
        fh <- cbind(fh[,1], pos1, fh[,2:6])
        colnames(fh) <- c("chromosome", "start", "end", "siteinf", "FAR", "yhat", "type")
        fh$chromosome <- paste0("chr", fh$chromosome)
	fh$FAR <- 100*as.numeric(fh$FAR)
	fh$yhat <- 100*as.numeric(fh$yhat)
    } else {
        fh = data.frame(chromosome=character(), start=numeric(), end=numeric(), siteinf=character(), FAR=numeric(), yhat=numeric(), type=character(), stringsAsFactors=FALSE)
    }

    return(fh)
}

Get_SegRes <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        colnames(fh) <- c("chromosome", "start", "end", "resHL", "support")
        fh$chromosome <- paste0("chr", fh$chromosome)
        fh <- subset(fh, resHL != "inconclusive")
    } else {
        fh <- data.frame(chromosome=character(), start=numeric(), end=numeric(), resHL=character(), support=character(), stringsAsFactors=FALSE)
    }

    return(fh)
}

Get_SNPcate <- function(fh, cate) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 <- fh$Position - 1
        fh <- cbind(fh[,1], pos1, fh[,2:3])
        colnames(fh) <- c("chromosome", "start", "end", "BAF")
        fh$type <- rep(cate, nrow(fh))
        fh$chromosome <- paste0("chr", fh$chromosome)
    } else {
        fh <- data.frame(chromosome=character(), start=numeric(), end=numeric(), BAF=numeric(), type=character(), stringsAsFactor=FALSE)
    }

    return(fh)
}

Get_SNPhap <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        pos1 <- fh$Position - 1
        fh <- cbind(fh[,1], pos1, fh[,2:4])
        colnames(fh) <- c("chromosome", "start", "end", "pat", "mat")
        fh$chromosome <- paste0("chr", fh$chromosome)
    } else {
        fh <- data.frame(chromosome=character(), start=numeric(), end=numeric(), pat=character(), mat=character(), stringsAsFactor=FALSE)
    }
    return(fh)
}

snparrcomp <- function(patfh, matfh, patres, matres, samplename, ffest_scale, chrom, snpP1, snpP2, snpM1, snpM2, snphap, diseaseStart, diseaseEnd, gene, figdir) {
    ffest_scale <- as.numeric(ffest_scale)
    ffest_scale_fit <- ffest_scale*2

    diseaseStart <- as.numeric(diseaseStart)
    diseaseEnd <- as.numeric(diseaseEnd)
    genepos <- diseaseStart+(diseaseEnd-diseaseStart)/2

    patfh <- Get_SNPratio(patfh)
    matfh <- Get_SNPratio(matfh)

    patres <- Get_SegRes(patres)
    matres <- Get_SegRes(matres)

    patfhsub <- subset(patfh, chromosome==chrom)
    matfhsub <- subset(matfh, chromosome==chrom)
    patressub <- subset(patres, chromosome==chrom)
    matressub <- subset(matres, chromosome==chrom)

    snphap <- Get_SNPhap(snphap)

    snphapsub <- subset(snphap, chromosome==chrom)
    snphapsubpat <- subset(snphapsub, pat==1 | pat==2)
    snphapsubmat <- subset(snphapsub, mat==1 | mat==2)

    chrnumber <- gsub("chr","",chrom)
    chrname <- paste("Chromosome", chrnumber)

    fo <- paste0(figdir, samplename, ".", chrom, ".diseaselocus.manu.pdf")
    pdf(file=fo, width=6, height=4)
    

    gtrellis_layout(n_track = 14, category=chrom, xaxis = FALSE, track_axis = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE), border=FALSE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(0, 1, 0, 1, 0, 1, 0, 0.5, 0, 0.5, 0, 1, 0, 1, 0, 0.5, 0, 1, 0, 1, -ffest_scale_fit, ffest_scale_fit, 0, 1, -ffest_scale_fit, ffest_scale_fit, 0, 1), track_height = unit.c(unit(0.5, "null"), unit(0.5, "null"), unit(0.3, "null"), unit(0.4, "null"), unit(0.3, "null"), unit(0.5, "null"), unit(0.3, "null"), unit(0.3, "null"), unit(0.3, "null"), unit(0.5, "null"), unit(2.5, "null"), unit(0.3, "null"), unit(2.5, "null"), unit(0.3, "null")), track_ylab = c("", "", "", "", "", "", "", "", "", "", "Pat.FAR.seg (%)", "", "Mat.FAR.seg (%)", ""), axis_label_fontsize = 6, lab_fontsize = 7, name_fontsize = 9, title=NULL, title_fontsize = 10, xlab=NULL)
    #gtrellis_layout(n_track = 9, category=chrom, track_axis = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE), border=FALSE, padding = unit(c(2, 2, 4, 2), "mm"), track_ylim = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, -ffest_scale_fit, ffest_scale_fit, 0, 1, -ffest_scale_fit, ffest_scale_fit, 0, 1), track_height = unit.c(unit(0.5, "null"), unit(0.3, "null"), unit(0.3, "null"), unit(0.3, "null"), unit(0.5, "null"), unit(2, "null"), unit(0.3, "null"), unit(2, "null"), unit(0.3, "null")), track_ylab = c("", "", "", "", "", "Pat.FAR(%)", "", "Mat.FAR(%)",""), add_ideogram_track = TRUE, add_name_track = TRUE, axis_label_fontsize = 5, lab_fontsize = 6, name_fontsize = 8, title=samplename, title_fontsize = 10, xlab = NULL)

    #1
    add_track(panel_fun = function(gr) {
        grid.text(c(samplename), gp=gpar(fontsize=9))
    })
    #2
    add_track(panel_fun = function(gr) {
        grid.text(c(chrname), gp=gpar(fontsize=10))
    })
    #3
    add_ideogram_track(track=3)
    #4
    add_track(panel_fun = function(gr) {
        grid.text("25MB", x=unit(25000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("50MB", x=unit(50000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("75MB", x=unit(75000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("100MB", x=unit(100000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("125MB", x=unit(125000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("150MB", x=unit(150000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("175MB", x=unit(175000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("200MB", x=unit(200000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
        grid.text("225MB", x=unit(225000000, "native"), y=unit(0.5, "npc"), gp = gpar(col = "black", fontsize=6))
	grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	#grid.text(c(gene), x=unit(genepos, "native"), y=unit(0.5, "npc"), gp=gpar(col="#efb509", fontsize=7))
																	   		      		      
    })
    #5
    add_track(panel_fun = function(gr) {
        grid.text(c(gene), x=unit(genepos, "native"), y=unit(0.5, "npc"), gp=gpar(col="#efb509", fontsize=8, fontface="italic"))
    })
    #6
    add_track(panel_fun = function(gr) {
        grid.text(c("PGT embryo"), gp=gpar(fontsize=9))
	grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
    })

    #7
    if (nrow(snphapsubpat)!=0) {
        add_rect_track(snphapsubpat, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(snphapsubpat$pat==2, "cornflowerblue", "blue")))
    } else {
        add_track(panel_fun = function(gr) {
       
	grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	    })
    }


    #8
    add_track(panel_fun = function(gr) {     
        grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
    })

    #9
    if (nrow(snphapsubmat)!=0) {
        add_rect_track(snphapsubmat, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(snphapsubmat$mat==1, "red", "pink")))
    } else {
        add_track(panel_fun = function(gr) {
          
        grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    }
    #10
    add_track(panel_fun = function(gr) {
      	grid.text(c("cffDNA"), gp=gpar(fontsize=9))
	grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
    })

    #11
    if (nrow(patfhsub)!=0) {
        add_track(patfhsub, panel_fun = function(gr) {
            x <- patfhsub$start
            y <- patfhsub$yhat
            grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(patfhsub$type=="P1", "#c82027", "#344D90")))
	    grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
            grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    } else {
        add_track(panel_fun = function(gr) {
            grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
	    grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	    grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    }
    #12
    if (nrow(patressub)!=0) {
        add_rect_track(patressub, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(patressub$resHL=="PH2", "cornflowerblue", "blue")))
    } else {
        add_track(panel_fun = function(gr) {
          
        grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    }
    #13
    if (nrow(matfhsub)!=0) {
        add_track(matfhsub, panel_fun = function(gr) {
            x <- matfhsub$start
            y <- matfhsub$yhat
            grid.lines(unit(c(0, 1), "npc"), unit(c(ffest_scale, ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(0, 0), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.lines(unit(c(0, 1), "npc"), unit(c(-ffest_scale, -ffest_scale), "native"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.8))
            grid.points(x, y, pch = 16, size = unit(0.8, "mm"), gp = gpar(col=ifelse(matfhsub$type=="M1", "#c82027", "#344D90")))
	    grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	    grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    } else {
        add_track(panel_fun = function(gr) {
	    grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
	    grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    }
    #14
    if (nrow(matressub)!=0) {
        add_rect_track(matressub, h1=0.05, h2=0.95, gp = gpar(col=NA, fill=ifelse(matressub$resHL=="MH1", "red", "pink")))
    } else {
        add_track(panel_fun = function(gr) {
	    grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
            grid.lines(unit(c(diseaseEnd, diseaseEnd), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
        })
    }
    
    dev.off()
}

snparrcomp(patfh, matfh, patres, matres, samplename, ffest_scale, chrom, snpP1, snpP2, snpM1, snpM2, snphap, diseaseStart, diseaseEnd, gene, figdir)


