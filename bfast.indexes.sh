#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
threads=$2
bn_ff=`basename $ff`

cat <<EOF
# Usage:
# script <path_fasta_ref_genome> <number_of_threads>
# 
EOF

echo "ln -s $ff `$bn_ff`"
echo "$bfast_bin fasta2brg -A0 -f $bn_ff"
echo "$bfast_bin fasta2brg -A1 -f $bn_ff"
echo ""

i=1
for m in ${h_masks[*]}
do
  echo "$bfast_bin index -n$threads -A0 -f $bn_ff -m $m -w 14 -i $i"
  echo "$bfast_bin index -n$threads -A1 -f $bn_ff -m $m -w 14 -i $i"
  i=$[$i+1]
done
cd ..
