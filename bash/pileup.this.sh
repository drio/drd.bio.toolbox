#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

cat <<-EOF
#
# Usage:
# <script> <bam_file> <fast_ref_genome>
#
# 1. Generate pileup
# -v     : print varias only, -c output soap consensus, -f ref genome (fasta)
# -F 1024: Don't include dups in the snp calling.
# -Q20   : Only reads of mapping qual 20
$samtools view -bF 1024 $1 | $samtools pileup -Q20 -vc -f $2 - > $1.pileup

# 2. Filter for high coverage
#                              con_q   SNP_q max_map_q coverage   a1      a2
#  0      1       2       3       4       5       6       7       8       9
# ---------------------------------------------------------------------------------
# chr17   93      t       A       0       6       60      2       .a      "E
# chr17   212     a       N       0       0       0       1       ^#.-1G  !
# chr17   3018    a       G       6       41      60      5       .$ggG.  .9;;"
#
# nr_a# : number of reads supporting allele #
#                                                                               nr_a1  nr_a2     nr_supporting_other_a
#  0      1       2       3       4       5       6       7       8       9      10      11      12      13       14
# -------------------------------------------------------------------------------------------------------------------
# chr17   212     *       -g/-g   40      0       2       1       -g      *       1       0       0       0       0
#
# Filter string:
#  0 1 2 3 4 5 6 7 8 9
# "U Q d D W G g s i X"
# 
# push(@staging, [$score, $flt, $len, @t]);
#
# w SNP within INT bp around a gap to be filtered [10]
# l window size for filtering adjacent gaps       [30]
# W window size for filtering dense SNPs          [10]
# so, max_dist will be 30 under default values.           
#
# Options: -Q INT    minimum RMS mapping quality for SNPs [25]
#          -q INT    minimum RMS mapping quality for gaps [10]
#          -d INT    minimum read depth [3]
#          -D INT    maximum read depth [100]
#          -S INT    minimum SNP quality []
#          -i INT    minimum indel quality []
# 
#          If you have a snp close to a indel (window: -w INT), filter it out if 
#          the indel has low quality.
#          -G INT    min indel score for nearby SNP filtering [25]
#          -w INT    SNP within INT bp around a gap to be filtered [10]
# 
#          -W INT    window size for filtering dense SNPs [10]
#          -N INT    max number of SNPs in a window [2]
# 
#          -l INT    window size for filtering adjacent gaps [30]
# 
#          -p        print filtered variants
$samtools_perl varFilter -p $1.pileup > $1.pileup.varfilter 2> $1.pileup.filtered_out

# 3. Only report Substitutions and Indels based on a quality threshold
# 75 is the quality threshold for indels and 20 for substitutions
awk '(\$3=="*"&&\$6>=50) || (\$3!="*"&&\$6>=20)' $1.pileup.varfilter > $1.pileup.var_filter.qual_threshold
EOF
