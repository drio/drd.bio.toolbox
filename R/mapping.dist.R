# 
# cat bfast.eval.accuracy.txt | grep -v "0 / 0" | awk '{print $1" "$5" "$6}' > bfast.accuracy.txt
# cat bwa.eval.accuracy.txt | grep -v "0 / 0" | awk '{print $1" "$5" "$6}' > bwa.accuracy.txt
# cat novo.eval.accuracy.txt | grep -v "0 / 0" | awk '{print $1" "$5" "$6}' > novo.accuracy.txt
# R CMD BATCH "--vanilla --args in.file='data_file'" < R_code
#
colors <- c("blue", "green")
what   <- c("mapped", "correctly mapped")
pdf(file='mapping.dist.pdf', height=10, width=12)
par(mfrow=c(3,1))

l <- scan("bfast.accuracy.txt", what=list(mapq=0, ok=0, mapped=0), skip=0)
b <- rbind(l$ok, l$mapped)
barplot(b, beside = TRUE, names.arg = rev(l$mapq), col=colors, main="BFAST MAPQ distribution")

l <- scan("bwa.accuracy.txt", what=list(mapq=0, ok=0, mapped=0), skip=0)
b <- rbind(l$ok, l$mapped)
barplot(b, beside = TRUE, names.arg = rev(l$mapq), col=colors, main="BWA MAPQ distribution")

l <- scan("novo.accuracy.txt", what=list(mapq=0, ok=0, mapped=0), skip=0)
b <- rbind(l$ok, l$mapped)
barplot(b, beside = TRUE, names.arg = rev(l$mapq), col=colors, main="NOVOALIGN_CS MAPQ distribution")
