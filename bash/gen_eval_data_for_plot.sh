#!/bin/bash
#
#set -e
#set -x

usage()
{
  cat <<-EOF
Usage:
  `basename $0` <pileup_file> <dwgsim_true_snps> <fast_file_ref_genome>
Example:
$ for f in *.final; do  gen_eval_data_for_plot.sh  \$f ../reads/true.snps.txt /data/next_gen_1/drio_scratch/macaque/genomes/macaque_chr17.fa; don
EOF
  exit 1
}

source "`dirname ${BASH_SOURCE[0]}`/common.sh"

pu_f=$1
true_f=$2
ff=$3

[ ! -f "$pu_f" ]   && usage && error "Pileup file not found."
[ ! -f "$true_f" ] && usage && error "SNPs from dwgsim not found."
[ ! -f "$ff" ]    && usage && error "fasta file not found."

genome_size=`cat $ff | \
$ruby -ne 'BEGIN{@t=0}; @t = @t + $_.chomp.size unless $_ =~ /^>/; END{puts @t.to_i}'`
echo "Genome size = $genome_size bp"

dwgsim_pileup_eval.pl $pu_f $true_f > $pu_f."eval".all
cat $pu_f."eval".all | grep "DELS" | sed 's/x//g' | awk "{print \$2\" \"\$6\" \"\$8\" \"\$9\" \"$genome_size}" > $pu_f."eval".dels
cat $pu_f."eval".all | grep "INSS" | sed 's/x//g' | awk "{print \$2\" \"\$6\" \"\$8\" \"\$9\" \"$genome_size}" > $pu_f."eval".inss
cat $pu_f."eval".all | grep "SNPS" | sed 's/x//g' | awk "{print \$2\" \"\$6\" \"\$8\" \"\$9\" \"$genome_size}" > $pu_f."eval".snps
