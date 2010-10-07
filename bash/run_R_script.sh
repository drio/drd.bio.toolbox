#!/bin/bash
#
set -e 

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

d_path=`dirname $0`
data_f=$2
r_script="${d_path}/../R/$1.R"

[ ".$1" == "." ] && error "first argument has to be R script to run."
[ ".$2" == "." ] && error "second argument has to be the input data."
[ ! -f $data_f ] && error "data_file not found."
[ ! -f $r_script ] && error "R script not found: $r_name"

#R CMD BATCH "--vanilla --args in.file='$1'" ${d_path}/../R/${r_name}
R CMD BATCH "--vanilla --args in.file='$data_f'" $r_script
