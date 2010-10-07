#!/bin/bash
#
set -e
#set -x

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

ff=$1
n_reads=$2
rl_1=$3
rl_2=$4
e_rate=$5

cat <<-EOF
# Usage: 
# script <fast_file> <n_reads> <read_length_1> <read_lenght_2> <error_rate>
#
$dwgsim -c -e $e_rate -N $n_reads -1 $rl_1 -2 $rl_2 $ff -d 2500 -s 250 \
$bwa_read1 $bwa_read2 $bf_reads 2> /dev/null > $true_snps
$novo_bf_2_novo $bf_reads novo
EOF
