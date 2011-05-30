#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  echo "Usage:"
  echo "$0 <bam>"
  exit 0
}

if=$1

[ ".$if" == "." ] && usage "I need an input_bam"
[ ! -d "$ipipe_java" ]  && usage "Couldn't find illumina pipe dir: $ipipe_java"

cat <<EOF
# Fix mate information (BWA does not follow the BAM specs)
#
$java -jar -Xmx4g $PICARD/FixMateInformation.jar \
I=$if \
TMP_DIR=$tmp \
VERBOSITY=ERROR \
VALIDATION_STRINGENCY=$PICARD_VALIDATION

# Fix CIGAR
# Class to modify CIGAR and mapping quality of alignments. 
# 1. For mapped reads, if CIGAR extends beyond the chromosome end, 
# clip CIGAR. (not implemented yet).
# 2 For unmapped reads, reset mapping quality to zero and reset 
# CIGAR to *.
#
# $java -jar -Xmx4g $ipipe_java/FixCIGAR.jar \
# I=$if \
# TMP_DIR=$tmp \
# VERBOSITY=ERROR \
# VALIDATION_STRINGENCY=$PICARD_VALIDATION
EOF
