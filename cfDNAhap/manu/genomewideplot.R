list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args = commandArgs(TRUE)
patres <- args[1]
matres <- args[2]
scpat <- args[3]
scmat <- args[4]
samplename <- args[5]
indsup <- args[6] #yes or no
figdir <- args[7]


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

Get_SCres <- function(fh) {
    fh <- read.delim(fh, header=TRUE, sep="\t")
    if (nrow(fh) != 0) {
        colnames(fh) <- c("chromosome", "start", "end", "resHL")
        fh$chromosome <- paste0("chr", fh$chromosome)
        fh <- subset(fh, resHL != "inconclusive")
    } else {
        fh <- data.frame(chromosome=character(), start=numeric(), end=numeric(), resHL=numeric(), stringsAsFactors=FALSE)
    }
    return(fh)
}

cffGW <- function(patres, matres, scpat, scmat, samplename, indsup, figdir) {
    patres <- Get_SegRes(patres)
    matres <- Get_SegRes(matres)

    patresB <- subset(patres, support=="B")
    patresS <- subset(patres, support=="S")
    matresB <- subset(matres, support=="B")
    matresS <- subset(matres, support=="S")
    
    scpat <- Get_SCres(scpat)
    scpat <- subset(scpat, chromosome != "chrX")
    scmat <- Get_SCres(scmat)
    plotchr <- paste0("chr", seq(1,22))
    plotchr <- c(plotchr, "chrX")

    fo <- paste0(figdir, samplename, ".", indsup, ".genomewide.manu.pdf")
    pdf(file=fo, width=11.69, height=8.27)
    gtrellis_layout(n_track = 6, ncol = 2, category=plotchr, xaxis=FALSE, track_axis = FALSE, xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE, asist_ticks = FALSE, add_ideogram_track = TRUE, ideogram_track_height =unit(2, "mm"), xlab = NULL, title=samplename, title_fontsize = 12)

    if (nrow(scpat) != 0) {
        add_rect_track(scpat, h1=0.5, h2=0, gp=gpar(col=NA, fill=ifelse(scpat$resHL==2, "cornflowerblue", "blue")))
    } else {
        add_track(panel_fun = function(gr){})
    }
    
    if (nrow(patres) != 0) {
        add_rect_track(patresB, h1=0.5, h2=0, gp = gpar(col=NA, fill=ifelse(patresB$resHL=="PH2", "cornflowerblue", "blue")))
        add_rect_track(patresS, track=current_track(), h1=0.5, h2=0, gp = gpar(col=NA, fill=ifelse(patresS$resHL=="PH2", "cornflowerblue", "blue")))
	if (indsup=="yes") {
	    add_rect_track(patresS, track=current_track(), h1=0.6, h2=0.7, gp=gpar(col=NA, fill="#818181"))
        }
    } else {
        add_track(panel_fun = function(gr){})
    }
    
    add_track(panel_fun = function(gr) {
        grid.text(c(""), gp=gpar(fontsize=10))
        #grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
    })

    if (nrow(scmat) != 0) {
        add_rect_track(scmat, h1=0.5, h2=0, gp=gpar(col=NA, fill=ifelse(scmat$resHL==1, "red", "pink")))
    } else {
        add_track(panel_fun = function(gr){})
    }

    if (nrow(matres) != 0) {
        add_rect_track(matresB, h1=0.5, h2=0, gp = gpar(col=NA, fill=ifelse(matresB$resHL=="MH1", "red", "pink")))
        add_rect_track(matresS, track=current_track(), h1=0.5, h2=0, gp = gpar(col=NA, fill=ifelse(matresS$resHL=="MH1", "red", "pink")))
        if (indsup=="yes") {
            add_rect_track(matresS, track=current_track(), h1=0.6, h2=0.7, gp=gpar(col=NA, fill="#818181"))
        }
    } else {
        add_track(panel_fun = function(gr){})
    }

    #add_track(track=current_track(), panel_fun=function(gr) {
         #grid.lines(unit(c(matresS$start, matresS$end), "native"), unit(c(0, 1), "npc"), gp = gpar(col = "#d9d9d9", lty = 2, lwd=0.5))
    #})

    add_track(panel_fun = function(gr) {
        grid.text(c(""), gp=gpar(fontsize=10))
        #grid.lines(unit(c(diseaseStart, diseaseStart), "native"), unit(c(0, 1), "npc"), gp = gpar(col="#efb509", lty=2, lwd=1.2))
    })
    
    add_track(track = 3, clip = FALSE, panel_fun = function(gr) {
        chr = get_cell_meta_data("name")
        grid.text(chr, x = 0, y = 0, just = c("left", "bottom"), gp=gpar(fontsize=12))
    })
    dev.off()
}


cffGW(patres, matres, scpat, scmat, samplename, indsup, figdir)