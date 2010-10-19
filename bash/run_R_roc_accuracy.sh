#!/bin/bash
#
set -e 

usage()
{
  cat<<-EOF
Usage:
  script <data_file> <output name pdf file>
Examples:
$ 
EOF
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname ${BASH_SOURCE[0]}`
r_script="${d_path}/../R/roc.accuracy.R"

[ ! -f $1 ] && usage && error "data file not provided"
[ ".$2" == "." ] && usage && error "Output name not provided"

R CMD BATCH "--vanilla --args in.file='$1' in.o_file='$2'" $r_script
ls -lach $2
