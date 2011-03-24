#!/usr/bin/env Rscript
#
# zcat coverage.txt.gz | head -10000000 | awk '{print $2}' | sort -n | uniq -c | \ 
# scatter.plot.one.var.R title output.pdf
#
# NOTE: if you use pileup straight use filed $8
#
args=(commandArgs(TRUE));
d <- read.table(file="stdin", colClasses=c("numeric","numeric"), col.names=c('x', 'y'), header=F)
pdf(file=args[2], width=7, height=4)
plot(d$y, d$x, xlab="coverage", ylab="frequency", main=args[1], xaxt="n", xlim=c(1,100))
axis(1, seq(0,100,5))
