library(ggplot2)

args <- commandArgs(TRUE)
s1accfh <- args[1]
s1gapfh <- args[2]
s2accfh <- args[3]
s2gapfh <- args[4]
ffdilutefh <- args[5]


getfile <- function(fh) {
    infh <- read.delim(fh, header=TRUE, sep="\t")
    return(infh)
}

ffdilute <- getfile(ffdilutefh)
colnames(ffdilute) <- c("mi", "ci", "ff", "acc", "resolution", "Type", "run")
ffdilute$ff <- as.factor(ffdilute$ff)
ffdilute$resolution <- ffdilute$resolution/1000

pdf("Simulation.FF.acc.pdf",width=6, height=4)
ggplot(ffdilute, aes(x=ff, y=acc, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#0000d8", "#d80000", "#fcbba1"), labels=c("Paternal","Maternal","Maternal+Type4"), name="SNP Type")+scale_fill_manual(values=c("#0000d8", "#d80000", "#fcbba1"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()

pdf("Simulation.FF.gap.pdf",width=6, height=4)
ggplot(ffdilute, aes(x=ff, y=resolution, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#0000d8", "#d80000", "#fcbba1"), labels=c("Paternal", "Maternal", "Maternal+Type4"), name="SNP Type")+scale_fill_manual(values=c("#0000d8", "#d80000", "#fcbba1"),name="SNP Type")+ scale_y_continuous(breaks=seq(0, 1600, 250))+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()


###separate
ffdilutePM <- subset(ffdilute, Type!="Type3+Type4")
pdf("Simulation.FF.acc.PM.pdf",width=6, height=4)
ggplot(ffdilutePM, aes(x=ff, y=acc, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#0000d8", "#d80000"), labels=c("Paternal","Maternal"), name="SNP Type")+scale_fill_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()

pdf("Simulation.FF.gap.PM.pdf",width=6, height=4)
ggplot(ffdilutePM, aes(x=ff, y=resolution, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#0000d8", "#d80000"), labels=c("Paternal", "Maternal"), name="SNP Type")+scale_fill_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+ scale_y_continuous(breaks=seq(0,1600, 250))+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()

ffdiluteMM <- subset(ffdilute, Type!="Type2")
pdf("Simulation.FF.acc.MM.pdf",width=6, height=4)
ggplot(ffdiluteMM, aes(x=ff, y=acc, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#d80000", "#fcbba1"), labels=c("Maternal","Maternal+Type4"), name="SNP Type")+scale_fill_manual(values=c("#d80000", "#fcbba1"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()

pdf("Simulation.FF.gap.MM.pdf",width=6, height=4)
ggplot(ffdiluteMM, aes(x=ff, y=resolution, color=factor(Type)))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_color_manual(values=c("#d80000", "#fcbba1"), labels=c("Maternal", "Maternal+Type4"), name="SNP Type")+scale_fill_manual(values=c("#d80000", "#fcbba1"),name="SNP Type")+ scale_y_continuous(breaks=seq(0,1600, 250))+theme_bw()+theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Fetal Fraction(%)") + scale_x_discrete(labels = c("5.50", "6.40", "7.14", "7.86", "8.50", "9.28", "10.00", "11.40", "12.80", "14.20", "17.00", "19.30"))
dev.off()




s1acc <- getfile(s1accfh)
colnames(s1acc) <- c("prop", "totalreads", "Accuracy", "Type", "rep", "thresh")
s1acc <- subset(s1acc, prop != 10 & prop != 15)
s1acc$Accuracy <- 100*s1acc$Accuracy
s1acc$samprop <- paste0(s1acc$prop," (",s1acc$totalreads, ")")
s1acc$samprop <- as.factor(s1acc$samprop)
s1acc$Type <- factor(s1acc$Type, levels=c("Paternal","Maternal"))

s1accthresh <- subset(s1acc, thresh=="30x")
pdf("Family3_085.subsample.30x.acc.pdf",width=6, height=4)
ggplot(s1accthresh, aes(x=samprop, y=Accuracy, color=Type))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_fill_manual(values=c("#0000d8", "#d80000"), labels=c("Paternal", "Maternal"), name="SNP Type")+scale_color_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Subsampling percentage (Total reads in millions)")
dev.off()

s1acc$combine <- paste0(s1acc$Type, "_", s1acc$thresh)
s1acc$combine <- as.factor(s1acc$combine)

pdf("Family3_085.subsample.thresh.acc.pdf",width=6, height=4)
ggplot(s1acc, aes(x=samprop, y=Accuracy, color=combine))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=combine, color=combine))+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Subsampling percentage (Total reads in millions)")
dev.off()


s2acc <- getfile(s2accfh)
colnames(s2acc) <- c("prop", "totalreads", "Accuracy", "Type", "rep", "thresh")
s2acc <- subset(s2acc, prop != 10)
#s2acc$prop <- as.factor(s2acc$prop)
#s2acc$Type <- as.factor(s2acc$Type)
s2acc$Accuracy <- 100*s2acc$Accuracy
s2acc$samprop <- paste0(s2acc$prop," (",s2acc$totalreads, ")")
s2acc$samprop <- as.factor(s2acc$samprop)
s2acc$Type <- factor(s2acc$Type, levels=c("Paternal","Maternal"))

s2accthresh <- subset(s2acc, thresh=="30x")
pdf("Family1_181.subsample.30x.acc.pdf",width=6, height=4)
ggplot(s2accthresh, aes(x=samprop, y=Accuracy, color=Type))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_fill_manual(values=c("#0000d8", "#d80000"),labels=c("Paternal", "Maternal"),name="SNP Type")+scale_color_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Subsampling percentage (Total reads in millions)")
dev.off()

s2acc$combine <- paste0(s2acc$Type, "_", s2acc$thresh)
s2acc$combine <- as.factor(s2acc$combine)

pdf("Family1_181.subsample.thresh.acc.pdf",width=6, height=4)
ggplot(s2acc, aes(x=samprop, y=Accuracy, color=combine))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=combine, color=combine))+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.2))+labs(y = "Accuracy(%)", x="Subsampling percentage (Total reads in millions)")
dev.off()


s1gap <- getfile(s1gapfh)
colnames(s1gap) <- c("prop", "totalreads", "Gap", "Type", "rep", "thresh")
s1gap <- subset(s1gap, prop != 10 & prop !=15)
#s1gap$prop <- as.factor(s1gap$prop)
#s1gap$Type <- as.factor(s1gap$Type)
s1gap$Gap <- s1gap$Gap/1000
s1gap$samprop <- paste0(s1gap$prop," (",s1acc$totalreads, ")")
s1gap$samprop <- as.factor(s1gap$samprop)
s1gap$Type <- factor(s1gap$Type, levels=c("Paternal","Maternal"))

s1gapthresh <- subset(s1gap, thresh=="30x")
pdf("Family3_085.subsample.30x.gap.pdf",width=6, height=4)
ggplot(s1gapthresh, aes(x=samprop, y=Gap, color=Type))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_fill_manual(values=c("#0000d8", "#d80000"),labels=c("Paternal","Maternal"),name="SNP Type")+scale_color_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Subsampling percentage (Total reads in millions)")
dev.off()

s1gap$combine <- paste0(s1gap$Type, "_", s1gap$thresh)
s1gap$combine <- as.factor(s1gap$combine)

pdf("Family3_085.subsample.thresh.gap.pdf",width=6, height=4)
ggplot(s1gap, aes(x=samprop, y=Gap, color=combine))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=combine, color=combine))+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Subsampling percentage (Total reads in millions)")
dev.off()

s2gap <- getfile(s2gapfh)
colnames(s2gap) <- c("prop", "totalreads", "Gap", "Type", "rep", "thresh")
s2gap <- subset(s2gap, prop != 10)
#s2gap$prop <- as.factor(s2gap$prop)
#s2gap$Type <- as.factor(s2gap$Type)
s2gap$Gap <- s2gap$Gap/1000
s2gap$samprop <- paste0(s2gap$prop," (",s2gap$totalreads, ")")
s2gap$samprop <- as.factor(s2gap$samprop)
s2gap$Type <- factor(s2gap$Type, levels=c("Paternal","Maternal"))

s2gapthresh <- subset(s2gap, thresh=="30x")
pdf("Family1_181.subsample.30x.gap.pdf",width=6, height=4)
ggplot(s2gapthresh, aes(x=samprop, y=Gap, color=Type))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=Type, color=Type))+scale_fill_manual(values=c("#0000d8", "#d80000"),labels=c("Paternal","Maternal"), name="SNP Type")+scale_color_manual(values=c("#0000d8", "#d80000"),name="SNP Type")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Subsampling percentage (Total reads in millions)")
dev.off()

s2gap$combine <- paste0(s2gap$Type, "_", s2gap$thresh)
s2gap$combine <- as.factor(s2gap$combine)

pdf("Family1_181.subsample.thresh.gap.pdf",width=6, height=4)
ggplot(s2gap, aes(x=samprop, y=Gap, color=combine))+geom_boxplot(outlier.colour = NULL, outlier.size = 1)+stat_summary(fun.y=mean, geom="line", aes(group=combine, color=combine))+scale_fill_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c"),name="Category")+theme_bw()+theme(axis.text.x = element_text(size=7, angle=45, hjust=1), axis.text.y = element_text(size=10), axis.line = element_line(colour = "black", size = 0.7, linetype = "solid"), legend.position=c(0.85,0.8))+labs(y = "Crossover Resolution(Kb)", x="Subsampling percentage (Total reads in millions)")
dev.off()
