#!/bin/bash
#
# bam.stats.sh: Generates stats from a bam file
#
set -e
#set -x

usage()
{
  [ ".$1" != "." ] && echo "ERROR: $1"
  cat <<-EOF
Usage:
  `basename $0` <sam/bam file> [n_ends] [min_map_qual]
  n_ends      : 1 or 2 [1]
  min_map_qual: default [0]
EOF
  exit 0
}

ctrl_c()
{
  echo "ctrl+c ... bye ..."
  rm -f ./tr*.$$ ./mr*.$$
}

trap ctrl_c SIGINT

if=$1
[ ".$2"  == "." ]  && n_ends=1 || n_ends=$2
[ ".$3"  == "." ]  && mmq=0 || mmq=$3
[ ".$if" == "." ] && usage "I need a bam/sam file"
if [ "$n_ends" != 1 ] && [ "$n_ends" != 2 ] 
then
  usage "Number of ends: 1 or 2"
fi

echo "n_ends       : $n_ends" >&2
echo "min map qual : $mmq" >&2

if [ `echo ${if#*.} | grep "sam$"` ]
then
  echo "SAM file detected" >&2
  bam=0
  stf="-S"
else
  echo "BAM file detected" >&2
  bam=1
  stf=""
fi

if [ $n_ends == 1 ]
then
  samtools view $stf -q $mmq $if      | wc -l > ./tr1.$$ &
  samtools view $stf -q $mmq -F 4 $if | wc -l > ./mr1.$$ &
else
  samtools view $stf -q $mmq -f 64       $if | wc -l > ./tr1.$$ &
  samtools view $stf -q $mmq -f 128      $if | wc -l > ./tr2.$$ &
  samtools view $stf -q $mmq -f 64 -F 4  $if | wc -l > ./mr1.$$ &
  samtools view $stf -q $mmq -f 128 -F 4 $if | wc -l > ./mr2.$$ &
fi

sleep 1
echo "Waiting for samtools instances to complete ..." >&2
wait

tr1=`cat ./tr1.$$`
mr1=`cat ./mr1.$$`
if [ $n_ends == 2 ]
then
  tr2=`cat ./tr2.$$`
  mr2=`cat ./mr2.$$`
fi

rm -f tr1.$$ mr1.$$ tr2.$$ mr2.$$

if [ $n_ends == 1 ]
then
  echo -ne "Total  # read1:\t $[$tr1] \n" 
  echo -ne "mapped # read1:\t $[$mr1] (`echo "scale=2 ; ($mr1*100)/$tr1 " | bc`%) \n"
else
  echo -ne "Total  # reads:\t $[$tr1+$tr2]\n"
  echo -ne "Total  # read1:\t $[$tr1]\n"
  echo -ne "Total  # read2:\t $[$tr2]\n"
  echo -ne "mapped # read1:\t $[$mr1] (`echo "scale=2 ; ($mr1*100)/$tr1 " | bc`%) \n"
  echo -ne "mapped # read2:\t $[$mr2] (`echo "scale=2 ; ($mr2*100)/$tr2 " | bc`%) \n"
fi
