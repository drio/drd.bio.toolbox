#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
[ ".$ff" != "." ] && bn_ff=`basename $ff`

cat <<-EOF
#
# Usage:
# script <fast_file_ref_genome>
# 
ln -s $ff bwa.$bn_ff
$bwa index -a bwtsw -p bwa.$bn_ff.ss $ff
$bwa index -c -a bwtsw -p bwa.$bn_ff.cs $ff
EOF
