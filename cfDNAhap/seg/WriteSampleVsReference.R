script.dir <- dirname(sys.frame(1)$ofile)
refcnvsp <- file.path(script.dir, "refcnv.R")
source(refcnvsp)


args <- commandArgs(TRUE)
sampledir <- args[1]
samplename <- args[2]
fragsis <- args[3]
targetgcfh <- args[4]
refdir <- args[5]
refname <- args[6]
gcc <- args[7] #TRUE or FALSE

GetRelCov(sampledir, samplename, fragsis, targetgcfh, refdir, refname, gcc)
