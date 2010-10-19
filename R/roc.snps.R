#
# R CMD BATCH "--vanilla --args in.bwa='dfile' in.bf='dfile' \ 
# in.novo='dfile' in.o_file='fname' in.title='title'" < R_code
#
args=(commandArgs(TRUE));
for(i in 1:length(args)){
  eval(parse(text=args[[i]]));
}

pdf(file=in.o_file, height=8, width=5)
par(mfrow=c(1,1))
colors <- c("red", "green", "black")
tools <- c("BWA", "BFAST", "NOVOALIGN")
l_type <- "l"

# Load data from each aligner
# 67 1764 3316 3355 94452569
bwa  <- scan(in.bwa, what=list(qual=0, ok=0, called=0, total=0, g_size=0), skip=1)
bf   <- scan(in.bf, what=list(qual=0, ok=0, called=0, total=0, g_size=0), skip=1)
novo <- scan(in.novo, what=list(qual=0, ok=0, called=0, total=0, g_size=0), skip=1)

# Calculate tpr / fpr
bwa_tpr  <- bwa$ok / bwa$total
bwa_fpr  <- log10(1 -  (bwa$g_size - bwa$total) / ((bwa$called - bwa$ok) + (bwa$g_size - bwa$total)))
novo_tpr <- novo$ok / novo$total
novo_fpr  <- log10(1 -  (novo$g_size - novo$total) / ((novo$called - novo$ok) + (novo$g_size - novo$total)))
bf_tpr <- bf$ok / bf$total
bf_fpr  <- log10(1 -  (bf$g_size - bf$total) / ((bf$called - bf$ok) + (bf$g_size - bf$total)))

xrange <- c(-7,-4)
yrange <- c(0,1)

plot(xrange, yrange, type="n", 
  xlab="log10 False Positive Rate (1 - specificity)",
  ylab="True positive rate (sensitivity)")

lines(x=bwa_fpr, y=bwa_tpr, type=l_type, lwd=1.5, col=colors[1])
lines(x=bf_fpr, y=bf_tpr, type=l_type, lwd=1.5, col=colors[2])
lines(x=novo_fpr, y=novo_tpr, type=l_type, lwd=1.5, col=colors[3])

title(in.title)
legend("bottomright", tools, col=colors, lty=1:1:1, cex=1.0)
