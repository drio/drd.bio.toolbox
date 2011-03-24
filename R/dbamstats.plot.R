#!/usr/bin/env Rscript
#
# $ cat data.txt | dbamstats.R output.pdf n_of_samples
#
# raw_pe raw_single_end mapped_paired mapped_unpaired mapped_single_end dups_mapped_paired dups_mapped_unpaired dups_mapped_single_end
# AWG_NLEU.00_000pA 233134668 177049 132578848 70379806 104 4416663 3928916 8
# AWG_NLEU.00_001pA 217104742 378761 118968719 64872054 15992 3344367 3383879 1546
# AWG_NLEU.00_002pA 238495788 314825 124026867 78926206 27852 9230048 6646557 5822
# AWG_NLEU.00_003pA 220715680 348367 113261856 71642118 31892 6719457 6849331 5128
# AWG_NLEU.00_004pA 226150452 373126 104821550 76632515 148010 10171402 7832392 87669
# AWG_NLEU.00_005pA 213642679 361821 101846992 78300015 107538 9604207 8307776 79286
# AWG_NLEU.00_006pA 294270906 27373940 172272102 83707298 17792578 7886288 5690520 2210902
# AWG_NLEU.00_007pA 237126278 435264 113597208 84282547 101110 1874592 4130303 60104
#
args=(commandArgs(TRUE));
pdf(file=args[1], height=7, width=13) 
par(mfrow=c(1,1))  
par(mar=c(11,5,4,2))
f <- read.table(file="stdin")
l <- length(row.names(f))
m <- as.matrix(f)
colors <- rainbow(args[2])
bp <- barplot(m, main="Mapping stats", beside=TRUE, col=colors, horiz=F, las=2, cex.lab=1, cex.axis=1)
legend("topright", c(row.names(f)), cex=1, bty="n", fill=colors);
