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

usage: `basename $0` options

Extract various stats from a SAM/BAM file

  -i <sam/bam file> : input SAM/BAM file
  -e n_ends         : number of ends [1]
  -m min_map_qual   : minimum mapping quality to account for [0]

Example(s):

  $ bam.stats.sh -i mybam.sort.dups.bam -e 2 -m 1 > stats.01003310991_3.bf2.n3.merged.sort.dups.1mq.txt
EOF
  exit 0
}

ctrl_c()
{
  echo "ctrl+c ... bye ..."
  rm -f ./tr*.$$ ./mr*.$$ ./d*.$$
}

trap ctrl_c SIGINT

if=""
n_ends=1
mmq=0
while getopts "i:e:m:h" OPTION
do
  case $OPTION in
    i)
      if=$OPTARG
      ;;
    e)
      n_ends=$OPTARG
      ;;
    m)
      mmq=$OPTARG
      ;;
    h)
      usage
      exit 0
      ;;
    ?)
      usage
      exit 2
      ;;
  esac
done
shift $(($OPTIND - 1))

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
  samtools view $stf -q $mmq $if         | wc -l > ./tr1.$$ &
  samtools view $stf -q $mmq -F 4 $if    | wc -l > ./mr1.$$ &
  samtools view $stf -q $mmq -f 1024 $if | wc -l > ./d1.$$ &
else
  samtools view $stf -q $mmq -f 64       $if | wc -l > ./tr1.$$ &
  samtools view $stf -q $mmq -f 128      $if | wc -l > ./tr2.$$ &
  samtools view $stf -q $mmq -f 64 -F 4  $if | wc -l > ./mr1.$$ &
  samtools view $stf -q $mmq -f 128 -F 4 $if | wc -l > ./mr2.$$ &
  samtools view $stf -q $mmq -f 1088     $if | wc -l > ./d1.$$  &
  samtools view $stf -q $mmq -f 1152     $if | wc -l > ./d2.$$  &
fi

sleep 1
echo "Waiting for samtools instances to complete ..." >&2
wait

tr1=`cat ./tr1.$$`
mr1=`cat ./mr1.$$`
d1=`cat ./d1.$$`
if [ $n_ends == 2 ]
then
  tr2=`cat ./tr2.$$`
  mr2=`cat ./mr2.$$`
  d2=`cat ./d2.$$`
fi

rm -f tr*.$$ mr*.$$ d*.$$

if [ $n_ends == 1 ]
then
  echo -ne "# Total  read1:\t $[$tr1] \n"
  echo -ne "# mapped read1:\t $[$mr1] (`echo "scale=2 ; ($mr1*100)/$tr1 " | bc`%) \n"
  echo -ne "# dups   read1:\t $[$d1] (`echo "scale=2 ; ($d1*100)/$tr1 " | bc`%) \n"
else
  echo -ne "# Total  reads:\t $[$tr1+$tr2]\n\n"

  echo -ne "# Total  read1:\t $[$tr1]\n"
  echo -ne "# mapped read1:\t $[$mr1] (`echo "scale=2 ; ($mr1*100)/$tr1 " | bc`%) \n"
  echo -ne "# dups   read1:\t $[$d1] (`echo "scale=2 ; ($d1*100)/$tr1 " | bc`%) \n\n"

  echo -ne "# Total  read2:\t $[$tr2]\n"
  echo -ne "# mapped read2:\t $[$mr2] (`echo "scale=2 ; ($mr2*100)/$tr2 " | bc`%) \n"
  echo -ne "# dups   read1:\t $[$d2] (`echo "scale=2 ; ($d2*100)/$tr1 " | bc`%) \n"
fi
