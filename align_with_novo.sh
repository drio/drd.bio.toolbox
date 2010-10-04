#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

i=$1
[ ".$i" == "." ] && error "Need full path to novo index"

for n in `seq -f"%02g" 0 20`
do
  fq_r1=`ls *read1*.$n 2>/dev/null`
  fq_r2=`ls *read2*.$n 2>/dev/null`
  if [ -f $fq_r1 ] && [ -f $fq_r2 ]
  then
    echo "$n"
    s_name="novo.aln.$n.sh"
    #cmd="$novo_cs/novoalignCS -c $threads -d $i -f $fq_r1 $fq_r2 -o SAM > novo.$n.sam"
    cmd="$novo_cs/novoalignCS -d $i -f $fq_r1 $fq_r2 -o SAM > novo.$n.sam"
    script_this "$cmd" "$s_name"
    qsub -N "$s_name" -l "nodes=1:ppn=1" $s_name
  fi
done
