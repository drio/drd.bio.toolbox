#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

r_f=$1
output_seed=$2
lines_per_split=$3

[ ! -f $r1_f ]                 && error "Need path to fastq read file."
[ ".$output_seed" == "." ]     && error "Output seed necessary."
[ ".$lines_per_split" == "." ] && error "Need # of lines per split."

cat <<-EOF
#
# Usage:
# script <fast_file_read_1> <output_seed> <number_of_lines_per_split>
# 
split -a 2 -d -l $lines_per_split $r_f $output_seed.
EOF
