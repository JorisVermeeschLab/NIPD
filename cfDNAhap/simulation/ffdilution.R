args <- commandArgs(TRUE)
hetmatfh <- args[1]
hetpatfh <- args[2]
homodifffh <- args[3]
aercfh <- args[4]
outfh <- args[5]
dilution <- args[6]

dilution <- as.numeric(dilution)

readinhetmat <- function(hetmatfh) {
    hetmat <- read.delim(hetmatfh, header=FALSE, comment.char="#", sep="\t")
    hetmatsite <- hetmat[,c(1,2)]
    colnames(hetmatsite) <- c("contig", "position")
    return(hetmatsite)
}

readinhetpat <- function(hetpatfh) {
    hetpat <- read.delim(hetpatfh, header=FALSE, comment.char="#", sep="\t")
    hetpatsite <- hetpat[,c(1,2)]
    colnames(hetpatsite) <- c("contig", "position")
    return(hetpatsite)
}

readinhomodiff <- function(homodifffh) {
    homodiff <- read.delim(homodifffh, header=FALSE, comment.char="#", sep="\t")
    homodiffsite <- homodiff[,c(1,2)]
    colnames(homodiffsite) <- c("contig", "position")
    return(homodiffsite)
}

readinaerc <- function(aercfh) {
    aerc <- read.delim(aercfh, header=TRUE, sep="\t")
    return(aerc)
}

matpoolexpand <- function(dilution, hetmatsite, aerc) {
    aerc.hetmat <- merge(aerc, hetmatsite, by=c("contig", "position"))
    #expand.refCount <- round(aerc.hetmat$totalCount/2)+aerc.hetmat$refCount
    #expand.altCOunt <- round(aerc.hetmat$totalCOunt/2)+aerc.hetmat$altCount
    dtrow <- dim(aerc.hetmat)[1]
    exp.prop <- 1/dilution-1
    totalCount <- aerc.hetmat$totalCount
    refCount <- aerc.hetmat$refCount
    altCount <- aerc.hetmat$altCount
    #direfcnt <- rep(-1, dtrow)
    #dialtcnt <- rep(-1, dtrow)
 
    for (i in 1:dtrow) {
        set.seed(9)
        equalref <- rep(0, round(totalCount[i]/2))
        equalalt <- rep(1, round(totalCount[i]/2))
        equalpool <- c(equalref, equalalt)
        #randomize order in pool
        equalpool <- equalpool[sample(1:length(equalpool))]
        expandsample <- sample(equalpool, round(totalCount[i]*exp.prop), replace=TRUE)
        
        refpool <- rep(0, refCount[i])
        altpool <- rep(1, altCount[i])
        expandpool <- c(refpool, altpool, expandsample)
        #randomize pool
        expandpool <- expandpool[sample(1:length(expandpool))]

        dilutepool <- sample(expandpool, totalCount[i], replace=TRUE)
        aerc.hetmat$refCount[i] <- length(dilutepool[dilutepool==0])
        aerc.hetmat$altCount[i] <- length(dilutepool[dilutepool==1])
    }
    
    return(aerc.hetmat)
}

patpoolexpand <- function(dilution, hetpatsite, aerc) {
    aerc.hetpat <- merge(aerc, hetpatsite, by=c("contig", "position"))

    dtrow <- dim(aerc.hetpat)[1]
    exp.prop <- (1/dilution-1)/2
    totalCount <- aerc.hetpat$totalCount
    refCount <- aerc.hetpat$refCount
    altCount <- aerc.hetpat$altCount
 
    for (i in 1:dtrow) {
        set.seed(9)
        #expand the pool by sampling the majority allele
        if (refCount[i] >= altCount[i]) {
            refpool <- rep(0, round(refCount[i]*(1+exp.prop)))
            altpool <- rep(1, altCount[i])
        } else {
            refpool <- rep(0, refCount[i])
            altpool <- rep(1, round(altCount[i]*(1+exp.prop)))
        }

        expandpool <- c(refpool, altpool)
        expandpool <- expandpool[sample(1:length(expandpool))]

        dilutepool <- sample(expandpool, totalCount[i], replace=TRUE)
        aerc.hetpat$refCount[i] <- length(dilutepool[dilutepool==0])
        aerc.hetpat$altCount[i] <- length(dilutepool[dilutepool==1])
    }
    
    return(aerc.hetpat)
}

ffdilute_main <- function(dilution, hetmatfh, hetpatfh, homodifffh, aercfh, outfh) {
    aerc <- readinaerc(aercfh)
    hetmatsite <- readinhetmat(hetmatfh)
    hetpatsite <- readinhetpat(hetpatfh)
    homodiffsite <- readinhomodiff(homodifffh)
    aerc.hetmat <- matpoolexpand(dilution, hetmatsite, aerc)
    aerc.hetpat <- patpoolexpand(dilution, hetpatsite, aerc)
    aerc.homodiff <- patpoolexpand(dilution, homodiffsite, aerc)
    new.aerc <- rbind(aerc.hetmat, aerc.hetpat, aerc.homodiff)
    sortindex <- order(new.aerc$contig, new.aerc$position)
    new.aerc <- new.aerc[sortindex,]
    write.table(new.aerc, outfh, sep="\t", quote=FALSE, row.names=FALSE)
}

ffdilute_main(dilution, hetmatfh, hetpatfh, homodifffh, aercfh, outfh)
