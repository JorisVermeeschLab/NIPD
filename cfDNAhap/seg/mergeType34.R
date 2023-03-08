args <- commandArgs(TRUE)
SegResolve <- args[1]
hetpatsum <- args[2]
hetmatraw <- args[3]
hetbothraw <- args[4]

source(SegResolve)

Merge_Raw_Hetboth_Hetmat(hetpatsum, hetmatraw, hetbothraw)
