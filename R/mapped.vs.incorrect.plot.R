# 
# R CMD BATCH "--vanilla --args in.file='data_file'" < R_code
#
tools    <- c("BFAST", "BWA", "NOVOALIGN_CS")
colors   <- c("red", "blue", "green")
plotchar <- c(0,1,2)
pdf(file='mapped.vs.incorrect.pdf', height=9, width=9)
par(mfrow=c(2,1))

args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

d <- scan(in.file, what=list(mapq=0, bfast_ok=0, bfast_mapped=0, bwa_ok=0, bwa_mapped=0, novo_ok=0, novo_mapped=0, total=0), skip=1)

xrange <- c(0, 5.8e+7)
yrange <- c(0, 2e+6)

plot(xrange, yrange, type="n", xlab="mapped reads", ylab="wrongly mapped reads")
#grid(7, 5, lwd = 2)
lines(d$bfast_mapped, d$bfast_mapped - d$bfast_ok, type="l", col=colors[1], lwd=2.5, lty=2.8, pch=0 )
lines(d$bwa_mapped  , d$bwa_mapped   - d$bwa_ok, type="l", col=colors[2], lwd=2.5, lty=2.8, pch=1 )
lines(d$novo_mapped , d$novo_mapped  - d$novo_ok, type="l", col=colors[3], lwd=2.5, lty=2.8, pch=2 )
legend(xrange[1], yrange[2], tools, cex=0.8, col=colors, pch=plotchar)

title("mapped reads VS incorrectly mapped reads")

xrange <- c(5.2e+7, 5.9e+7)
yrange <- c(0, 1.2e+6)

plot(xrange, yrange, type="n", xlab="mapped reads", ylab="wrongly mapped reads")
#grid(7, 5, lwd = 2)
lines(d$bfast_mapped, d$bfast_mapped - d$bfast_ok, type="l", col=colors[1], lwd=2.5, lty=2.8, pch=0 )
lines(d$bwa_mapped  , d$bwa_mapped   - d$bwa_ok, type="l", col=colors[2], lwd=2.5, lty=2.8, pch=1 )
lines(d$novo_mapped , d$novo_mapped  - d$novo_ok, type="l", col=colors[3], lwd=2.5, lty=2.8, pch=2 )
legend(xrange[1], yrange[2], tools, cex=0.8, col=colors, pch=plotchar)

title("mapped reads VS incorrectly mapped reads(zoom in)")
