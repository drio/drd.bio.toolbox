#!/usr/bin/env Rscript
#
# $ cat marmo.set | script.R <output.pdf>
# Input
# id, sample-name,n_callable_positions,substitutions,indels,total_events,per_total_difference
# 0,NLL-607,Nomascus_Leucogenys,2370135147,1707841,8092740,9800581,0.413503
# ...
#
args=(commandArgs(TRUE));
pdf(file=args[1], height=6, width=12)
par(mfrow=c(1,2))  
par(mar=c(7,5,2,2))
colors=c("pink", "seagreen1")

d <- read.table(file="stdin", header=T, sep=",")
# First plot: amount of indels/subs per sample
m <- as.matrix(d[c("substitutions", "indels")])
barplot(t(m), names.arg = d$sample.name, beside=TRUE, las=2, col=colors, cex.axis=.6,
        main="SNPs per sample", border="black", ylab="# of events")
legend("top", c("substitutions", "indels"), cex=.7, fill=colors, bty="n")

# Second plot
# Add a couple of columns with the percentage of difference (subs and indels)
nd <- transform(d, per_subs = (100*substitutions)/n_callable_positions, 
                   per_indels = (100*indels)/n_callable_positions)
# We only need 
m <- as.matrix(nd[c("per_subs","per_indels")])
barplot(t(m), beside=TRUE, las=2, ylab="% of difference", names.arg = d$sample.name,
        main="genomic difference", border="black", col=c("pink", "seagreen1"))
legend("top", c("substitutions", "indels"), cex=.7, fill=colors, bty="n")
