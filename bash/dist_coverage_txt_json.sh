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

title="dist_coverage"
std_pileup $if | \
awk '{print $3}' | \
sort -n | \
uniq -c | \
sort -k1,1rn | \
tee ${title}.txt | \
ruby -ane '\
  BEGIN{puts "["} 
  puts "{ \"x\":#{$F[1]}, \"y\":#{$F[0]} }," 
  END{puts "]"}' \
> ${title}.json
