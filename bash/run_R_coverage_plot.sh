#!/bin/bash
#
#set -e 

usage()
{
  cat<<-EOF
Usage:
  script <output.pdf> <window_size> <title> <pileup_file>
EOF
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname ${BASH_SOURCE[0]}`
r_script="${d_path}/../R/coverage.plot.R"

if [ ".$1" == "." ] || [ ".$2" == "." ] || [ ".$3" == "." ] || [ ! -f $4 ]
then
  usage
  exit 1
fi

awk '{print $2" "$8}' $4 | R CMD BATCH "--vanilla --slave --args in.o_file='$1' in.window='$2' in.title='$3'" $r_script
