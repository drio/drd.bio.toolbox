#!/bin/env Rscript
#
pdf(file='mapping.dist.pdf', height=5, width=9)
par(mfrow=c(1,1))
colors <- c("red", "black", "green")
tools <- c("BFAST", "BWA", "NOVOALIGN")
l_type <- "l"

args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

# Load data 
#d <- scan(in.file, what=list(mapq=0, bfast_ok=0, bfast_mapped=0, bwa_ok=0, bwa_mapped=0, novo_ok=0, novo_mapped=0, total=0), skip=1)
d <- read.table(in.file, header=T)
d <- as.matrix(d)

#barplot(as.matrix(autos_data), main="Mapping Dist", ylab= "Alignments", beside=TRUE, col=rainbow(6))
barplot(c(d$mapq, OK_bwa.data.accuracy.txt), main="Mapping Dist", ylab= "Alignments", beside=TRUE, col=rainbow(2))
