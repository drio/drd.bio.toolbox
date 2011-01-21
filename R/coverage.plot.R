#!/usr/bin/env Rscript
#
# awk '{print $2" "$4}' pileup | R CMD BATCH out.pdf tittle range.max
#
args=(commandArgs(TRUE));

d <- read.table(file="stdin",sep=" ",header=F)
colnames(d) <- c("pos", "cov")
depth <- mean(d[,"cov"])
rangeto <- length(d[,"pos"])        
                                    
#pdf(file=args[1], height=4, width=14)
bitmap(file=args[1], type="png256", width=14, height=4, units="in", res=300)
plot(d$pos, d$cov, xlab="bp position", ylab="depth", type="n")
lines(d$pos, d$cov, col="red")
title(args[2])
