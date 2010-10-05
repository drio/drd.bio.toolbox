#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
[ ".$ff" == "." ] && error "Need path to fast file (bf indexes should be there)."

for n in `seq -f"%02g" 0 20`
do
  fq=`ls *reads*.$n 2>/dev/null`
  if [ -f $fq ]
  then
    echo "$n: $fq"
    s_name="bfast.aln.$fq.sh"
    cmd="bfast match -T$tmp -n$threads -A1 -l -f $ff -r $fq > bf.matches.$n.bmf"
    script_this "$cmd" "$s_name"
    dep=`qsub -N "$s_name" -l "nodes=1:ppn=$threads" $s_name`
    echo "  dep aln: $dep"

    cmd="bfast localalign -A1 -n $threads -f $ff -m bf.matches.$n.bmf > bf.local.$n.baf"
    s_name="bfast.local.$fq.sh"
    script_this "$cmd" "$s_name"
    dep=`qsub -W "depend=afterok:$dep" -N $s_name -l "nodes=1:ppn=$threads" $s_name`
    echo "  dep local: $dep"

    cmd="bfast postprocess -A1 -n $threads -f $ff -i bf.local.$n.baf -U -a3 -O1 > bf.$n.sam"
    s_name="bfast.pp.$fq.sh"
    script_this "$cmd" "$s_name"
    dep=`qsub -W "depend=afterok:$dep" -N $s_name -l "nodes=1:ppn=$threads" $s_name`
    echo "  dep post: $dep"
  fi
done
