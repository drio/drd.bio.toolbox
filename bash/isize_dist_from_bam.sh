#!/bin/bash
#
set -e
#set -x

usage()
{
  cat <<-EOF
Usage: `basename $0` <bam_file> <insert_size_window> <output_name> <title>

insert_size_window: For a 2kb library, use 4000.
output_name       : output.png
EOF
  exit 1
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

bam=$1
isize_w=$2
o_file=$3
title=$4

[ ! -f $bam ] && error "Bam not found." && usage
[ ".$isize_w" == "." ] && error "insert size window not provided." && usage
[ ".$o_file" == "." ] && error "output file not provided" && usage
[ ".$title" == "." ] && error "title not provided" && usage
`which gnuplot &>/dev/null` || (error "gnuplot not found" && usage)

#(echo "set t png"; echo "set logscale y"; \
cat <<-EOF
(echo "set t png"; \
echo plot \"-\" using 2:1 title \"$title\"; \
samtools view -f3 $bam | \
awk '{if (\$9 < $isize_w && \$9 > 10) print \$9}' | \
sort -n | uniq -c ) | \
gnuplot > $o_file
EOF
