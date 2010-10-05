#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
[ ".$ff" == "." ] && error "Need a path and seed to locate bwa indexes."

for n in `seq -f"%02g" 0 10`
do
  fq_r1=`ls *r1*.$n 2>/dev/null`
  fq_r2=`ls *r2*.$n 2>/dev/null`
  if [ -f $fq_r1 ] && [ -f $fq_r2 ]
  then
    echo "$n: $fq_r1 $fq_r2"
    s_name_r1="bwa.aln.$fq_r1.sh"
    cmd="bwa aln -c -t$threads $ff $fq_r1 > $fq_r1.sai"
    script_this "$cmd" "$s_name_r1" 
    dep_r1=`qsub -N "$s_name_r1" -l "nodes=1:ppn=$threads" $s_name_r1`
    echo "  dep1: $dep_r1"

    s_name_r2="bwa.aln.$fq_r2.sh"
    cmd="bwa aln -c -t$threads $ff $fq_r2 > $fq_r2.sai"
    script_this "$cmd" "$s_name_r2" 
    dep_r2=`qsub -N "$s_name_r2" -l "nodes=1:ppn=$threads" $s_name_r2`
    echo "  dep2: $dep_r2"

    cmd="bwa sampe $ff $fq_r1.sai $fq_r2.sai $fq_r1 $fq_r2 > $n.sam"
    s_name="sampe.$n.sh"
    script_this "$cmd" "$s_name"
    id=`qsub -W "depend=afterok:$dep_r1:$dep_r2" -N $s_name -l "nodes=1:ppn=1" $s_name`
    echo "  sampe id: $id"
  fi
done
