# cat i.txt | awk '{print $2" "$8 }' | uniq | awk '{print $2}'| calc_coverage_stats 10 > cov.txt
# i.txt comes from samtools pileup
#
args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

#pdf(file='./cov.pdf', height=5, width=12)
#bitmap(file="./cov.png", height=5, width=14)
#png(file=in.o_file, 8, 3, units="in", res=150)
bitmap(file=in.o_file, type="png256", width=7, height=2, units="in", res=300)
#par(mfrow=c(1,1))
#par(mar=c(11,5,4,2))
data <- scan(file="stdin", what=list(pos=0, me=0, sd=0), skip=0)
plot(data$pos, data$me, main=in.title, xlab="window number", 
ylab="# of reads", col="red", pch=".", cex=.4, ylim=c(0,20), cex.lab=1.2)
