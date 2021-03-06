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

title="dist_isize"
samtools view -f3 $bam | \
awk '{if ($9 < $isize_w && $9 > 10) print $9}' | \
sort -T$tmp -S$sort_buffer -n | uniq -c | \
awk '{if ($1 > 100) print}' | \
tee ${title}.txt | \
ruby -ane '\
  BEGIN{@h = "["; @body = []}
    @body << "[ #{$F[1]}, #{$F[0]} ]"
  END{puts @h + @body.join(",\n") + "]"}' \
> ${title}.json
