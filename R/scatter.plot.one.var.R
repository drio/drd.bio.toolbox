#!/usr/bin/env Rscript
#
# zcat coverage.txt.gz | head -10000000 | awk '{print $2}' | sort -n | uniq -c | \ 
# scatter.plot.one.var.R ylabel xlabel title output.png
#
args=(commandArgs(TRUE));
d <- read.table(file="stdin", colClasses=c("numeric","numeric"), col.names=c('x', 'y'), header=F)
#pdf(file=args[1], height=4, width=14)
#bitmap(file=args[4], type="png256", width=7, height=4, units="in", res=300)
pdf(file=args[4], width=7, height=4)
plot(d$y, d$x, log='x', xlab=args[1], ylab=args[2], main=args[3])
