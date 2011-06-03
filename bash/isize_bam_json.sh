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

title="isize_dist"
samtools view -f3 $bam | \
awk '{if ($9 < $isize_w && $9 > 10) print $9}' | \
sort -n | uniq -c | \
tee ${title}.txt | \
ruby -ane '\
  BEGIN{@h = "["; @body = []}
    @body << "[ #{$F[1]}, #{$F[0]} ]"
  END{puts @h + @body.join(",") + "]"}' \
> ${title}.json
