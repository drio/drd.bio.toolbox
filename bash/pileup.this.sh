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
# -v: print varias only, -c output soap consensus, -f ref genome (fasta)
$samtools pileup -vc -f $2 $1 > $1.pileup

# 2. Filter for high coverage
#$samtools_perl varFilter -D100 $1.pileup > $1.pileup.filter

# 3. Only report Substitutions and Indels based on a quality threshold
# 75 is the quality threshold for indels and 50 for substitutions
#awk '(\$3=="*"&&\$6>=75)||(\$3!="*"&&\$6>=50)' $1.pileup.filter > $1.pileup.filter.final
EOF
