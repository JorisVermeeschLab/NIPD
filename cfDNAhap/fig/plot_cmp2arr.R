list.of.packages <- c("gtrellis", "grid", "RColorBrewer", "circlize")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(gtrellis)
library(RColorBrewer)
library(grid)
library(circlize)

args = commandArgs(TRUE)
cfpatfh = args[1]
cfmatfh = args[2]
snppatfh = args[3]
snpmatfh = args[4]
workdir = args[5]
samplename = args[6]

readincf <- function(cffh) {
    cf <- read.delim(cffh, header=FALSE, sep="\t", comment.char="#")
    cf <- cf[1:(nrow(cf)-2),]
    colnames(cf) <- c("chr", "start", "end", "blocksize", "HL")
    cf$chr <- paste0("chr", cf$chr)
    return(cf)
}

readinsnp <- function(snpfh) {
    snp <- read.delim(snpfh, header=FALSE, sep="\t")
    snp <- snp[,1:6]
    colnames(snp) <- c("chr", "start", "end", "length", "cumlen", "HL")
    snp$chr <- paste0("chr", snp$chr)
    return(snp)
}

cmpplot <- function(cfpatfh, cfmatfh, snppatfh, snpmatfh, workdir, samplename) {
    cfpat <- readincf(cfpatfh)
    cfmat <- readincf(cfmatfh)
    snppat <- readinsnp(snppatfh)
    snpmat <- readinsnp(snpmatfh)
    
    plotchr <- paste0("chr", seq(1,22))
    fo <- paste0(workdir, "/", samplename, ".comp2arr.pdf")
    pdf(file=fo, width=11.7, height=7)
    gtrellis_layout(n_track = 4, ncol = 2, category=plotchr, track_axis = FALSE, xpadding = c(0.1, 0), gap = unit(4, "mm"), border = FALSE, asist_ticks = FALSE, add_ideogram_track = TRUE, ideogram_track_height = unit(2, "mm"))
    add_rect_track(snppat, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(snppat$HL==1, "#344D90", "#c82027")))
    add_rect_track(cfpat, h1 = 0.5, h2 = 0, gp = gpar(col=NA, fill = ifelse(cfpat$HL=="PH1", "#344D90", "#c82027")))
    add_rect_track(snpmat, h1=0.5, h2=0, gp = gpar(col=NA, fill = ifelse(snpmat$HL==1, "#efb509", "#3f681c")))
    add_rect_track(cfmat, h1=0.5, h2=0,  gp = gpar(col=NA, fill = ifelse(cfmat$HL=="MH1", "#efb509", "#3f681c")))
    add_track(track = 2, clip = FALSE, panel_fun = function(gr) {
        chr = get_cell_meta_data("name")
        grid.text(chr, x = 0, y = 0, just = c("left", "bottom"))
    })
    dev.off()
}

cmpplot(cfpatfh, cfmatfh, snppatfh, snpmatfh, workdir, samplename)

