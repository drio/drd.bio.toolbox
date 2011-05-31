#!/bin/bash
#
set -e
#set -x

usage()
{
  cat <<-EOF
Usage: `basename $0` <bam_file>
EOF
  exit 1
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

bam=$1

if [ ".$bam" == "." ] || [ ! -f $bam ] ;then
  usage
fi

title="isize.dist"
samtools view -f3 $bam | \
awk '{if ($9 < $isize_w && $9 > 10) print $9}' | \
sort -n | uniq -c | \
tee ${title}.txt | \
ruby -ane '\
  BEGIN{puts "["}
  puts "  { \"x\":#{$F[1]}, \"y\":#{$F[0]} },"
  END{puts "]"}' \
> ${title}.json
