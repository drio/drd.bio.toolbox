## 
# $ head out.txt 
# contig window_number #subs #ins #dels
# scaffold_0 0 7 3 4
# scaffold_0 1 0 2 6
# scaffold_0 2 4 4 8
# awk '{print $2" "$3" "$4" "$5}' ident.data | R CMD BATCH "--vanilla --args in.o_file='out.pdf' in.title='XXXXXXXXXXX'" ./this.R
#
args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

pdf(file=in.o_file, height=3.5, width=14)
#bitmap(file=in.o_file, height=5, width=16, pointsize=14)
par(mfrow=c(1,1))
plotchar <- c(0,0,0)
colors   <- c("red", "blue", "green")
types    <- c("Substitutions", "Insertions", "Deletions")

d <- scan(file="stdin", what=list(w_id=0, subs=0, inss=0, dels=0), skip=1)
xrange <- c(0, length(d$w_id))
yrange <- c(0, max(c(max(d$subs), max(d$inss), max(d$dels))))

plot(xrange, yrange, type="n", xlab="windows", ylab="events")
#grid(7, 5, lwd = 2)
#lines(d$w_id, d$subs, type="p", col=colors[1], lwd=2.5, lty=2.8, pch=0)
lines(d$w_id, d$subs, type="l", col=colors[1])
lines(d$w_id, d$inss, type="l", col=colors[2])
lines(d$w_id, d$dels, type="l", col=colors[3])
legend(xrange[1], yrange[2], types, cex=0.8, col=colors, pch=plotchar)
title(in.title)
