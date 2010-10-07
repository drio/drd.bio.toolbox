#!/bin/bash
#
set -e
#set -x

help()
{
  echo "Usage:" 
  echo "`basename $0` <ref_genome_ff> <output_seed> <lines_per_split>"
  error "$1"
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

r_f=$1
output_seed=$2
lines_per_split=$3

[ ! -f $r1_f ]                 && help "Need path to fastq read file."
[ ".$output_seed" == "." ]     && help "Output seed necessary."
[ ".$lines_per_split" == "." ] && help "Need # of lines per split."

cat <<-EOF
split -a 2 -d -l $lines_per_split $r_f $output_seed.
EOF
