#!/bin/bash
#
set -e
#set -x

usage()
{
  cat <<-EOF
Usage:
  `basename $0` <path_fasta_genome_ref> <number_of_ends> <rl_1> <rl_2> <des_n_reads>
EOF
  exit 1
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
n_ends=$2
rl_1=$3
rl_2=$4
reads_per_split=$5

[ ! -f $ff ] && usage && error "Path to ref genome file not provided."
[ ".$n_ends" == "." ] && usage && error "Number of ends not provided."
[ ".$rl_1" == "." ] && usage && error "rl_1 not provided."
[ ".$rl_2" == "." ] && usage && error "rl_2 not provided."
[ ".$reads_per_split" == "." ] && \
usage && error "wanted # reads perl split not provided."

genome_size=`cat $ff | \
$ruby -ne 'BEGIN{@t=0}; @t = @t + $_.chomp.size unless $_ =~ /^>/; END{puts @t.to_i}'`
echo "Genome size = $genome_size bp"

[ $n_ends == 1 ] && rl=$rl_1 || rl=`echo "$rl_1+$rl_2" | bc`
for i in 1 2 5 10 15 20 25 30 35 40 45 50
do
  n_rs=`echo "($genome_size)/($rl+($n_ends-1)) * $i" | bc`
  n_reads_per_split=`echo "$n_rs / $reads_per_split" | bc`
  echo -e "${i}x = $n_rs reads. ; $n_reads_per_split splits. ($reads_per_split reads per split)"
done
