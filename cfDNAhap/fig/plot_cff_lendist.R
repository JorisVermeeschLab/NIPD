list.of.packages <- c("ggplot2", "gridExtra", "grid", "lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)

args = commandArgs(TRUE)
fetal_len <- args[1]
shared_len <- args[2]
lenmin <- args[3]
lenmax <- args[4]

readin_lenfh = function(infh1, infh2) {
    fetallen = read.delim(infh1, header = FALSE, sep = "\t")
    fetalrows = dim(fetallen)[1]
    fetallen = cbind(fetallen, rep("Fetal", fetalrows))
    colnames(fetallen) = c("readID", "length", "Type")
    sharedlen = read.delim(infh2, header = FALSE, sep = "\t")
    sharedrows = dim(sharedlen)[1]
    sharedlen = cbind(sharedlen, rep("Shared", sharedrows))
    colnames(sharedlen) = c("readID", "length", "Type")
    combinedlen = rbind(fetallen, sharedlen)
    return(combinedlen)
}

plot_len_dist = function(infh1, infh2, lenmin, lenmax) {
    combinedlen = readin_lenfh(infh1, infh2)
    lenmin = as.numeric(lenmin)
    lenmax = as.numeric(lenmax)
    combinesub = subset(combinedlen, length > lenmin  & length <= lenmax)
    return(combinesub)
}

combinesub=plot_len_dist(fetal_len, shared_len, lenmin, lenmax)
fo = paste0(fetal_len, ".png")
png(file = fo, width = 560, height=320, res=100)
ggplot(combinesub, aes(x=length, color=Type)) + geom_density(adjust=1/2) + scale_color_manual(values = c("#C60000", "#000B29")) +theme_minimal() + labs(x="Fragment Size", y="Density")+scale_x_continuous(breaks=seq(as.numeric(lenmin), as.numeric(lenmax), 50)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

fo1 = paste0(fetal_len, ".zoomin.png")
png(file = fo1, width = 560, height=320, res=100)
ggplot(combinesub, aes(x=length, color=Type)) + geom_density(adjust=1/2) + scale_color_manual(values = c("#C60000", "#000B29")) +theme_minimal() + labs(x="Fragment Size", y="Density")+scale_x_continuous(limits=c(50, 200), breaks=seq(50, 200, 15)) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
