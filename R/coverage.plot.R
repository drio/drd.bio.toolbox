## 
# awk '{print $2" "$4" "$8}' pu.txt | R CMD BATCH "--vanilla --args in.o_file='out.pdf' in.window='50' in.title='XXXXXXXXXXX'" ./cov.R
#
args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

data <- read.table(file="/dev/stdin",sep=" ",header=F)
colnames(data) <- c("pos", "coverage")
#data <- scan("/dev/stdin", what=list(pos=0, coverage=0), skip=0)
depth <- mean(data[,"coverage"])
# depth now has the mean (overall)coverage
# set the bin-size
window <- in.window
rangefrom <- 0
rangeto <- length(data[,"pos"])

data.smoothed<-runmed(data[,"coverage"],k=window)
#pdf(file=in.o_file, height=5, width=14) 
bitmap(file=in.o_file, height=5, width=16, pointsize=14)
plot(x=data[rangefrom:rangeto,"pos"], y=data.smoothed[rangefrom:rangeto], pch=".", cex=1,xlab="bp position",ylab="depth",type="l", col="red")
title(in.title)
dev.off()
print(data[rangefrom:rangeto, "pos"])
print(data.smoothed[rangefrom:rangeto])
