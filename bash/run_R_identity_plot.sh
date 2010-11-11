#!/bin/bash
#
#set -e 

usage()
{
  cat<<-EOF
Usage:
  script <output.pdf> <identity_data_file>
EOF
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname ${BASH_SOURCE[0]}`
r_script="${d_path}/../R/identity.plot.R"

if [ ".$1" == "." ] || [ ! -f $2 ]
then
  usage
  exit 1
fi

awk '{print $1" "$4" "$5" "$6}' $2 | R CMD BATCH "--vanilla --slave --args in.o_file='$1'" $r_script
