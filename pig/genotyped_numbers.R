library(ggplot2)


genos <- read.table("~/porcine/5_snpeff/num.animals.genotyped.txt", h=T)
nrow(genos)


ggplot(genos) + geom_histogram(aes(x=num.genotyped), binwidth=1) + theme_bw() + xlab("Number of Animals Genotyped") + ylab("Number of Loci") + annotate("text", label=c("14.75%"), x=10, y=2500) + annotate("text", label=c("24.39%"), x=32, y=2500) + annotate("text", label=c("60.85%"), x=57, y=2500) + geom_vline(xintercept=20, color="grey") + geom_vline(xintercept=40, color="grey") + ggtitle("Histogram showing the total number of successfully genotyped animals at variant loci only")
ggplot(genos) + geom_histogram(aes(x=num.na), binwidth=1)
test <- head(genos, 100)


ggplot(test) + geom_bar(aes(x=location, y=num.genotyped), stat="identity")


one.to.twenty <- subset(genos, num.genotyped < 21)
twentyone.to.forty <- subset(genos, num.genotyped > 20 & num.genotyped <= 40)
forty.plus <- subset(genos, num.genotyped > 40)

nrow(one.to.twenty)
nrow(twentyone.to.forty)
nrow(forty.plus)
